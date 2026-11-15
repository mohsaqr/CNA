# Orchestration: one .run() for every analysis.
#
# cna_run_analysis() checks the inputs, builds (or restores) the network(s),
# then runs each section. Expensive results are cached in hidden state holders
# whose clearWith lists only what changes them (see tools/generate_jamovi.R),
# so a section recomputes only when one of its own inputs changed. A section
# that fails reports the error in its own heading and the others still run.

#' Generate dynamic instructions HTML for an analysis kind.
#' @param module_type "event_single", "event_cluster", "binary_single",
#'   "binary_cluster" or "delimited".
#' @return HTML string.
#' @noRd
cna_instructions_html <- function(module_type) {
  fields <- switch(module_type,
    event_single  = "<li><b>Event Data</b>: assign <b>Action</b> (required), <b>Actor</b>, <b>Time</b>, <b>Order</b> (optional). Add a <b>Group Variable</b> to build one network per group.</li>",
    event_cluster = "<li><b>Event Data</b>: assign <b>Action</b> and <b>Actor</b> (both required for clustering), <b>Time</b>, <b>Order</b> (optional).</li>",
    binary_single = "<li><b>Binary Data</b>: assign at least two <b>One-Hot Columns</b>, and optionally <b>Actor</b> and <b>Session</b>. Add a <b>Group Variable</b> to build one network per group.</li>",
    binary_cluster = "<li><b>Binary Data</b>: assign at least two <b>One-Hot Columns</b>, and optionally <b>Actor</b> and <b>Session</b>.</li>",
    delimited     = "<li><b>Delimited Data</b>: assign one or more <b>Item Columns</b> whose cells list items separated by the <b>Separator</b> (for example keywords \"network; graph\"), or one item per row grouped by a <b>Transaction ID</b> (the Separator is then not used).</li>")
  mode <- switch(module_type,
    event_single  = , binary_single = "<li><b>Single / Group</b>: one network from all data, or one per group when a Group Variable is set.</li>",
    event_cluster = "<li><b>Clustering</b>: clusters the sessions and builds one network per cluster. Tick <i>Run clustering</i> to start; <i>Choosing k</i> works before that.</li>",
    binary_cluster = "<li><b>Clustering</b>: clusters the one-hot windows and builds one network per cluster. Tick <i>Run clustering</i> to start; <i>Choosing k</i> works before that.</li>",
    delimited     = "<li><b>Co-occurrence</b>: items listed together in a row (or a transaction) co-occur. <i>Similarity</i> normalises the counts (van Eck &amp; Waltman, 2009).</li>")
  paste0(
    '<div style="border: 2px solid #e6f4fe; border-radius: 15px; padding: 15px; ',
    'background-color: #e6f4fe; margin-top: 10px;"><div style="text-align:justify;"><ul>',
    mode, fields,
    "<li>Engine: <b>Nestimate</b> (networks, inference) and <b>cograph</b> (centrality, communities, plots).</li>",
    "</ul></div></div>")
}

# Results that always stay visible (messages and headings).
.cna_message_items <- c("instructions", "errorText", "clusteringNote", "analysisTitle")

#' Hide every output (used while inputs are incomplete).
#' @noRd
.cna_hide_outputs <- function(results, keep = character()) {
  items <- setdiff(results$itemNames, c(.cna_message_items, keep))
  invisible(lapply(items, function(nm) results[[nm]]$setVisible(FALSE)))
}

.cna_show_error <- function(results, message) {
  results$errorText$setContent(message)
  results$errorText$setVisible(TRUE)
}

#' A message for missing inputs, "" for "nothing assigned yet", NULL if ready.
#' @noRd
.cna_input_problem <- function(options, spec) {
  switch(spec$data,
    event = {
      if (is.null(options$buildModel_variables_long_action)) ""
      else if (spec$cluster && is.null(options$buildModel_variables_actor))
        "Clustering requires an Actor variable."
    },
    onehot = {
      n <- length(options$buildModel_variables_onehot)
      if (n == 0L) ""
      else if (n < 2L) "Please select at least 2 one-hot columns."
      else if (spec$cluster && is.null(options$buildModel_variables_actor) &&
               is.null(options$buildModel_variables_session))
        "Clustering one-hot data needs an Actor or a Session: they define the units that are clustered."
    },
    delimited = {
      if (!length(options$buildModel_variables_field)) ""
      else if (is.null(options$buildModel_variables_by) &&
               !nzchar(trimws(options$buildModel_sep %||% "")))
        "Give a Separator for the item columns, or a Transaction ID to group one item per row."
    })
}

#' A friendlier message for build errors (blamed on the Time variable only
#' when one is assigned and the engine's message is about times).
#' @noRd
.cna_build_message <- function(e, options) {
  msg <- conditionMessage(e)
  time_set <- !is.null(options[["buildModel_variables_long_time"]])
  if (time_set && grepl("time|date|posix|format", msg, ignore.case = TRUE))
    paste("Could not read the Time variable:", msg)
  else paste("Could not build the network:", msg)
}

#' The network(s): restored from the state holder, or built and cached.
#' @noRd
.cna_get_model <- function(self, kind) {
  holder <- self$results$buildModelContent
  if (!is.null(holder$state)) return(holder$state)
  model <- tryCatch(cna_build_model(self$data, self$options, kind),
                    error = function(e) {
                      .cna_show_error(self$results, .cna_build_message(e, self$options))
                      NULL
                    })
  if (!is.null(model)) holder$setState(model)
  model
}

#' Restore a cached computation or compute and cache it.
#' @noRd
.cna_cached <- function(holder, compute) {
  if (!is.null(holder$state)) return(holder$state)
  value <- compute()
  holder$setState(value)
  value
}

#' One jamovi cell value from a data.frame cell.
#' @noRd
.cna_cell <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.numeric(x) && !is.finite(x)) return(if (is.na(x)) NaN else x)
  x
}

#' Add the rows of a tidy data.frame to a jamovi table (columns matched by
#' name; extra data.frame columns are ignored). The data.frame is also kept as
#' the table's state.
#' @noRd
.cna_fill <- function(table, df) {
  table$setState(df)
  if (is.null(df) || !nrow(df)) return(invisible(table))
  cols <- intersect(names(df), vapply(table$columns, function(col) col$name, character(1)))
  invisible(lapply(seq_len(nrow(df)), function(i)
    table$addRow(rowKey = i, values = lapply(as.list(df[i, cols, drop = FALSE]), .cna_cell))))
}

#' Add numeric columns `keys` (titled `titles`) to a table.
#' @noRd
.cna_add_columns <- function(table, keys, titles = keys, type = "number") {
  invisible(Map(function(k, t) table$addColumn(name = k, title = t, type = type),
                keys, titles))
}

.cna_no_edges_note <- paste(
  "The network has no edges between different nodes (only self co-occurrence).",
  "With one-hot data, a wider window lets codes in neighbouring rows co-occur.")

.cna_any <- function(options, names) any(vapply(names, function(n) isTRUE(options[[n]]), logical(1)))

#' Grid image size for n panels.
#' @noRd
.cna_size_for <- function(image, n, panel_w = 420L, panel_h = 380L,
                          min_w = 600L, min_h = 450L) {
  if (is.null(image) || n <= 1L) return(invisible())
  grid <- .cna_grid(n)
  image$setSize(max(min_w, grid[[2L]] * panel_w), max(min_h, grid[[1L]] * panel_h))
}

#' Does the data hold several units to resample (actors, or sessions for
#' one-hot data)? Without them the whole data set is one unit: bootstrap
#' intervals have zero width, stability coefficients are 0 and every
#' permutation p-value is 1 (verified with Nestimate 1.0.7), so those
#' analyses are not run.
#' @noRd
.cna_has_units <- function(options, spec) {
  switch(spec$data,
    event = !is.null(options$buildModel_variables_actor),
    onehot = !is.null(options$buildModel_variables_actor) ||
      !is.null(options$buildModel_variables_session),
    FALSE)
}

.cna_units_message <- function(spec) {
  need <- if (identical(spec$data, "onehot")) "an Actor or a Session" else "an Actor"
  paste0("This analysis resamples units of the data, so it needs ", need,
         ". Without one the data form a single unit: bootstrap intervals would ",
         "have zero width, stability coefficients would be 0 and every ",
         "permutation p-value would be 1.")
}

#' Why a resampling section is not run here, or NULL.
#'
#' * no units (Actor for events, Actor or Session for one-hot): the data are
#'   one unit, so intervals would have zero width, CS = 0 and every
#'   permutation p = 1 (Nestimate >= 1.0.7 refuses this itself);
#' * a scaled model: resampling is defined on raw counts only. Scaled models
#'   are never resampled (a design decision, 2026-11-06), although Nestimate
#'   >= 1.0.7 could scale each resample like the network.
#' @param section "bootstrap", "permutation" or "stability".
#' @noRd
.cna_resampling_problem <- function(ctx, section) {
  o <- ctx$options
  spec <- ctx$spec
  if (!.cna_has_units(o, spec)) return(.cna_units_message(spec))
  if (!is.null(.cna_scaling(o)))
    return(paste("Bootstrap, permutation and stability are run on raw counts.",
                 "Set Scaling to 'No scaling' for these analyses."))
  NULL
}

#' Stop a resampling section that cannot be valid: explain in each table's
#' note and hide its plots. Returns TRUE when the section may run.
#' @noRd
.cna_resampling_ok <- function(ctx, section, tables, plots) {
  problem <- .cna_resampling_problem(ctx, section)
  if (is.null(problem)) return(TRUE)
  invisible(lapply(tables, function(nm) ctx$results[[nm]]$setNote("units", problem)))
  invisible(lapply(plots, function(nm) ctx$results[[nm]]$setVisible(FALSE)))
  FALSE
}

# ---------------------------------------------------------------------------
# Sections. Each takes the run context; ctx$model is the network(s).
# ---------------------------------------------------------------------------

.cna_section_model <- function(ctx) {
  o <- ctx$options; r <- ctx$results
  n_missing <- attr(ctx$model, "n_missing_group")
  notes <- c(attr(ctx$model, "cna_notes"),
             if (isTRUE(n_missing > 0)) sprintf(
               "%d row%s with a missing Group value were left out.", n_missing,
               if (n_missing == 1L) "" else "s"))
  if (length(notes)) r$buildModelTitle$setContent(paste(notes, collapse = "\n"))
  if (isTRUE(o$buildModel_show_matrix)) {
    m <- cna_matrix_df(ctx$model)
    states <- setdiff(names(m), c("group", "state"))
    keys <- paste0("m", seq_along(states))
    names(m)[match(states, names(m))] <- keys
    .cna_add_columns(r$buildModel_matrixTable, keys, states)
    .cna_fill(r$buildModel_matrixTable, m)
  }
  if (isTRUE(o$buildModel_show_summary))
    .cna_fill(r$buildModel_summaryTable, cna_summary_df(ctx$model))
  # only analyses that declare the mixed-network options may read them
  # (jmvcore errors on an undeclared option), so check the kind first
  if (isTRUE(ctx$spec$mixed) && isTRUE(o$buildModel_show_mixed)) {
    pairs <- .cna_cached(r$buildModel_mixedState, function()
      cna_build_mixed(ctx$data, o, ctx$kind, ctx$model))
    n <- length(pairs)
    if (identical(o$buildModel_mixed_type, "overlay")) .cna_size_for(r$buildModel_mixed, n)
    else r$buildModel_mixed$setSize(900L, max(450L, 420L * n))
  }
  n <- ctx$n_nets
  .cna_size_for(r$buildModel_plot, n)
  .cna_size_for(r$buildModel_histo, n)
}

.cna_section_centrality <- function(ctx) {
  o <- ctx$options; tbl <- ctx$results$centralityTable
  if (!.cna_any(o, c("centrality_show_table", "centrality_show_plot"))) return()
  measures <- cna_selected_measures(o)
  cent <- if (!is.null(tbl$state)) tbl$state else
    cna_centrality_df(ctx$model, measures, loops = isTRUE(o$centrality_loops),
                      normalized = isTRUE(o$centrality_normalize),
                      invert = isTRUE(o$centrality_invert))
  .cna_add_columns(tbl, measures, .cna_measure_titles[measures])
  .cna_fill(tbl, cent)
  img <- ctx$results$centrality_plot
  if (ctx$n_nets > 1L && !identical(o$centrality_plot_style, "distribution"))
    img$setSize(650L, max(450L, 330L * ctx$n_nets))
}

.cna_section_stability <- function(ctx) {
  o <- ctx$options; r <- ctx$results
  if (!.cna_any(o, c("stability_show_table", "stability_show_plot"))) return()
  if (!.cna_resampling_ok(ctx, "stability", c("stabilityTable", "reliabilityTable"),
                          "stability_plot")) return()
  st <- .cna_cached(r$stabilityState, function() cna_stability_compute(
    ctx$model, cna_selected_measures(o), iter = o$stability_iter,
    threshold = o$stability_threshold, certainty = o$stability_certainty,
    method = o$stability_method, loops = isTRUE(o$centrality_loops),
    invert = isTRUE(o$centrality_invert),
    seed = ctx$seed))
  if (isTRUE(o$stability_show_table)) {
    .cna_fill(r$stabilityTable, cna_stability_df(st, ctx$multi))
    .cna_fill(r$reliabilityTable, cna_reliability_df(st))
    failures <- cna_stability_failures(st)
    if ("cs" %in% names(failures)) r$stabilityTable$setNote("cs", failures[["cs"]])
    if ("cd" %in% names(failures)) r$stabilityTable$setNote("cd", failures[["cd"]])
    if ("nr" %in% names(failures)) r$reliabilityTable$setNote("nr", failures[["nr"]])
    missing <- cna_stability_missing(st, ctx$multi)
    if (length(missing))
      r$stabilityTable$setNote("constant", sprintf(paste(
        "No coefficient for %s: the measure does not vary across nodes in at",
        "least one network, so its stability is undefined."),
        paste(missing, collapse = ", ")))
  }
  if (ctx$n_nets > 1L) r$stability_plot$setSize(700L, max(450L, 300L * ctx$n_nets))
}

.cna_section_community <- function(ctx) {
  o <- ctx$options; r <- ctx$results
  if (!.cna_any(o, c("community_show_table", "community_show_plot"))) return()
  coms <- .cna_cached(r$communityState, function() cna_community_compute(
    ctx$model, o$community_methods, gamma = o$community_gamma, seed = ctx$seed))
  if (isTRUE(o$community_show_table)) .cna_fill(r$communityTable, cna_community_df(coms))
  .cna_size_for(r$community_plot, ctx$n_nets)
}

.cna_section_edge_betweenness <- function(ctx) {
  o <- ctx$options; r <- ctx$results
  if (!.cna_any(o, c("edgebetweenness_show_table", "edgebetweenness_show_plot"))) return()
  if (isTRUE(o$edgebetweenness_show_table)) {
    df <- cna_edge_betweenness_df(cna_edge_betweenness_compute(
      ctx$model, isTRUE(o$edgebetweenness_invert)))
    .cna_fill(r$edgeBetweennessTable, df)
    if (is.null(df)) r$edgeBetweennessTable$setNote("none", .cna_no_edges_note)
  }
  .cna_size_for(r$edgeBetweenness_plot, ctx$n_nets)
}

.cna_section_cliques <- function(ctx) {
  o <- ctx$options; r <- ctx$results
  if (!.cna_any(o, c("cliques_show_table", "cliques_show_plot", "cliques_show_complex"))) return()
  complexes <- .cna_cached(r$cliquesState, function()
    cna_cliques_compute(ctx$model, o$cliques_size, o$cliques_threshold))
  cliques <- cna_cliques_df(complexes, ctx$model, o$cliques_size)
  if (isTRUE(o$cliques_show_table)) {
    .cna_fill(r$cliquesTable, cliques)
    if (is.null(cliques))
      r$cliquesTable$setNote("none", sprintf(
        "No cliques of %d nodes with every edge weight at least %s.",
        o$cliques_size, format(o$cliques_threshold)))
  }
  if (isTRUE(o$cliques_show_complex)) {
    cx <- cna_complex_df(complexes)
    .cna_fill(r$cliquesComplexTable, cx$rows)
    invisible(Map(function(note, k) r$cliquesComplexTable$setNote(paste0("n", k), note),
                  cx$notes, seq_along(cx$notes)))
  }
  if (!is.null(cliques)) {
    per_group <- table(cliques$group)
    n_shown <- sum(pmin(per_group, o$cliques_max_plot))
    .cna_size_for(r$cliques_plot, n_shown, panel_w = 300L, panel_h = 280L)
  }
}

.cna_section_backbone <- function(ctx) {
  o <- ctx$options; r <- ctx$results
  if (!.cna_any(o, c("backbone_show_table", "backbone_show_plot"))) return()
  bb <- .cna_cached(r$backboneState, function() cna_backbone_compute(
    ctx$model, method = o$backbone_method, level = o$backbone_level,
    threshold = o$backbone_threshold, lowest = o$backbone_lowest))
  if (isTRUE(o$backbone_show_table)) {
    df <- cna_backbone_df(bb)
    .cna_fill(r$backboneTable, df)
    if (is.null(df)) r$backboneTable$setNote("none", .cna_no_edges_note)
  }
  .cna_size_for(r$backbone_plot, ctx$n_nets)
}

.cna_section_links <- function(ctx) {
  o <- ctx$options; tbl <- ctx$results$linksTable
  if (!isTRUE(o$links_show_table)) return()
  links <- cna_links_df(ctx$model, method = o$links_method,
                        exclude_existing = isTRUE(o$links_exclude_existing),
                        threshold = o$links_threshold, top_n = o$links_top_n)
  .cna_fill(tbl, links)
  if (is.null(links)) {
    one_node <- all(vapply(.cna_nets(ctx$model), function(net) nrow(net$weights) < 2L,
                           logical(1)))
    tbl$setNote("none", if (one_node)
      "A network with a single node has no links to predict."
    else paste(
      "No absent edges to predict: every pair of nodes already co-occurs above",
      "the threshold. Raise the threshold, or untick 'Only absent edges' to",
      "score the existing ones."))
  }
}

.cna_section_rules <- function(ctx) {
  o <- ctx$options; r <- ctx$results
  if (!.cna_any(o, c("rules_show_table", "rules_show_plot"))) return()
  rules <- .cna_cached(r$rulesState, function() cna_rules_compute(
    ctx$model, min_support = o$rules_min_support, min_confidence = o$rules_min_confidence,
    min_lift = o$rules_min_lift, max_length = o$rules_max_length))
  if (isTRUE(o$rules_show_table)) {
    df <- cna_rules_df(rules, o$rules_max_rows)
    .cna_fill(r$rulesTable, df)
    if (is.null(df)) r$rulesTable$setNote("none",
      "No rules meet the support, confidence and lift limits; lower them to find more.")
  }
  if (ctx$n_nets > 1L) r$rules_plot$setSize(700L, max(450L, 330L * ctx$n_nets))
}

.cna_section_bootstrap <- function(ctx) {
  o <- ctx$options; r <- ctx$results
  if (!.cna_any(o, c("bootstrap_show_table", "bootstrap_show_plot"))) return()
  if (!.cna_resampling_ok(ctx, "bootstrap", "bootstrapTable", "bootstrap_plot")) return()
  threshold <- if (identical(o$bootstrap_method, "threshold") && isTRUE(o$bootstrap_set_threshold))
    o$bootstrap_threshold
  bs <- .cna_cached(r$bootstrapState, function() cna_bootstrap_compute(
    ctx$model, iter = o$bootstrap_iteration, level = o$bootstrap_level,
    inference = o$bootstrap_method,
    range = c(o$bootstrap_range_low, o$bootstrap_range_up), threshold = threshold,
    ci_method = o$bootstrap_ci_method, seed = ctx$seed))
  if (isTRUE(o$bootstrap_show_table))
    .cna_fill(r$bootstrapTable, cna_bootstrap_df(
      bs, ctx$multi, significant_only = isTRUE(o$bootstrap_table_significant_only),
      show_all = isTRUE(o$bootstrap_table_show_all), max_rows = o$bootstrap_table_max_rows))
  .cna_size_for(r$bootstrap_plot, ctx$n_nets)
}

.cna_section_permutation <- function(ctx) {
  o <- ctx$options; r <- ctx$results
  if (!.cna_any(o, c("permutation_show_table", "permutation_show_plot"))) return()
  if (!.cna_resampling_ok(ctx, "permutation", "permutationTable", "permutation_plot")) return()
  pt <- .cna_cached(r$permutationState, function() cna_permutation_compute(
    ctx$model, iter = o$permutation_iter, alpha = o$permutation_level,
    paired = isTRUE(o$permutation_paired), adjust = o$permutation_adjust,
    seed = ctx$seed))
  if (isTRUE(o$permutation_show_table)) {
    .cna_fill(r$permutationTable, cna_permutation_df(
      pt, show_all = isTRUE(o$permutation_table_show_all),
      max_rows = o$permutation_table_max_rows))
    if (!identical(o$permutation_adjust, "none"))
      r$permutationTable$setNote("adjust", sprintf(
        "p-values adjusted (%s) over the unique edges of each comparison.",
        o$permutation_adjust))
  }
  .cna_size_for(r$permutation_plot, length(pt), min_h = 600L)
}

#' Resolve the two networks to compare (group levels or cluster numbers).
#' @noRd
.cna_compare_pair <- function(model, options, cluster) {
  groups <- names(model)
  if (cluster) {
    i <- options$compare_cluster_i; j <- options$compare_cluster_j
    if (i > length(groups) || j > length(groups))
      .cna_stop(sprintf("There are %d clusters; choose numbers from 1 to %d.",
                        length(groups), length(groups)))
  } else {
    gi <- options$compare_group_i %||% groups[[1L]]
    gj <- options$compare_group_j %||% groups[[min(2L, length(groups))]]
    if (!nzchar(gi)) gi <- groups[[1L]]
    if (!nzchar(gj)) gj <- groups[[min(2L, length(groups))]]
    i <- match(gi, groups); j <- match(gj, groups)
  }
  if (is.na(i) || is.na(j) || i == j) .cna_stop("Please select two different networks.")
  list(i = i, j = j, group_i = groups[[i]], group_j = groups[[j]])
}

.cna_section_compare <- function(ctx) {
  o <- ctx$options; r <- ctx$results
  shows <- c("compare_show_summary", "compare_show_network", "compare_show_edges",
             "compare_show_plot", "compare_show_diff_plot")
  if (!.cna_any(o, shows)) return()
  state <- .cna_cached(r$compareState, function() {
    pair <- .cna_compare_pair(ctx$model, o, ctx$spec$cluster)
    c(pair, list(result = cna_compare_compute(ctx$model, pair$i, pair$j,
                                              scaling = o$compare_scaling)))
  })
  r$compareTitle$setContent(paste("Comparing:", state$group_i, "vs", state$group_j))
  if (isTRUE(o$compare_show_summary))
    .cna_fill(r$compareSummaryTable, cna_compare_summary_df(state$result))
  if (isTRUE(o$compare_show_network)) {
    r$compareNetworkTable$getColumn("group_i")$setTitle(state$group_i)
    r$compareNetworkTable$getColumn("group_j")$setTitle(state$group_j)
    .cna_fill(r$compareNetworkTable, cna_compare_network_df(state$result))
  }
  if (isTRUE(o$compare_show_edges)) {
    r$compareEdgeTable$getColumn("weight_i")$setTitle(state$group_i)
    r$compareEdgeTable$getColumn("weight_j")$setTitle(state$group_j)
    .cna_fill(r$compareEdgeTable, cna_compare_edges_df(state$result, o$compare_edges_max_rows))
  }
  if (identical(o$compare_plot_type, "all")) r$compare_plot$setSize(900L, 700L)
}

.cna_section_cluster_diag <- function(ctx) {
  o <- ctx$options
  if (isTRUE(o$cluster_show_diagnostics))
    .cna_fill(ctx$results$clusterDiagTable, cna_cluster_diagnostics_df(ctx$model))
}

#' Choosing k runs before (and without) the clustering itself.
#' @noRd
.cna_section_cluster_choice <- function(self, kind) {
  o <- self$options; r <- self$results
  if (!.cna_any(o, c("cluster_show_choice", "cluster_show_choice_plot"))) return()
  cc <- .cna_cached(r$clusterChoiceState, function() {
    seqs <- if (.cna_kinds[[kind]]$data == "event") cna_event_sequences(self$data, o)
            else cna_onehot_units(self$data, o)
    cna_cluster_choice_compute(seqs, o$cluster_choice_kmin, o$cluster_choice_kmax,
                               dissimilarity = o$clustering_dissimilarity,
                               method = o$clustering_method, seed = .cna_seed(o))
  })
  if (isTRUE(o$cluster_show_choice))
    .cna_fill(r$clusterChoiceTable, cna_cluster_choice_df(cc))
}

# Section registry: heading result, runner, and whether it needs several
# networks. Sections whose heading is absent from an analysis are skipped.
.cna_sections <- list(
  list(title = "buildModelTitle",      run = .cna_section_model),
  list(title = "clusterTitle",         run = .cna_section_cluster_diag),
  list(title = "centralityTitle",      run = .cna_section_centrality),
  list(title = "stabilityTitle",       run = .cna_section_stability),
  list(title = "communityTitle",       run = .cna_section_community),
  list(title = "edgeBetweennessTitle", run = .cna_section_edge_betweenness),
  list(title = "cliquesTitle",         run = .cna_section_cliques),
  list(title = "backboneTitle",        run = .cna_section_backbone),
  list(title = "linksTitle",           run = .cna_section_links),
  list(title = "rulesTitle",           run = .cna_section_rules),
  list(title = "bootstrapTitle",       run = .cna_section_bootstrap),
  list(title = "permutationTitle",     run = .cna_section_permutation, multi = TRUE,
       items = c("permutationTitle", "permutationState", "permutationTable", "permutation_plot")),
  list(title = "compareTitle",         run = .cna_section_compare, multi = TRUE,
       items = c("compareTitle", "compareState", "compareSummaryTable", "compareNetworkTable",
                 "compareEdgeTable", "compare_plot", "compare_diff_plot")))

#' Run one section, reporting a failure in its heading.
#' @noRd
.cna_run_section <- function(ctx, section) {
  heading <- ctx$results[[section$title]]
  if (is.null(heading)) return(invisible())
  if (isTRUE(section$multi) && !ctx$multi) {
    invisible(lapply(section$items, function(nm) ctx$results[[nm]]$setVisible(FALSE)))
    return(invisible())
  }
  # headings keep their text across runs (Preformatted content is restored),
  # so clear any earlier message before the section reports afresh
  heading$setContent("")
  tryCatch(section$run(ctx), error = function(e) {
    heading$setContent(paste("This analysis could not be completed:", conditionMessage(e)))
  })
  invisible()
}

#' The single .run() of every CNA analysis.
#' @param self The jamovi analysis object.
#' @param kind One of names(.cna_kinds).
#' @return NULL, invisibly; the results are written into `self$results`.
#' @noRd
cna_run_analysis <- function(self, kind) {
  options <- self$options
  results <- self$results
  spec <- .cna_kinds[[kind]]
  results$instructions$setContent(cna_instructions_html(kind))

  problem <- .cna_input_problem(options, spec)
  if (!is.null(problem)) {
    .cna_hide_outputs(results)
    if (nzchar(problem)) .cna_show_error(results, problem)
    return(invisible())
  }

  if (spec$cluster) {
    choice_items <- c("clusterTitle", "clusterChoiceState", "clusterChoiceTable",
                      "clusterChoice_plot")
    results$clusterTitle$setContent("")
    tryCatch(.cna_section_cluster_choice(self, kind), error = function(e)
      results$clusterTitle$setContent(paste("Choosing k could not be completed:",
                                            conditionMessage(e))))
    if (!isTRUE(options$clustering_run)) {
      results$clusteringNote$setContent(paste0(
        '<div style="border: 2px solid #e6f4fe; border-radius: 15px; padding: 15px; ',
        'background-color: #e6f4fe; margin-top: 10px;"><p style="text-align:center; ',
        'font-weight: bold; color: #1a5276; font-size: 14px;">Tick "Run clustering" ',
        'to cluster and build the networks. "Choosing k" (under Cluster Diagnostics) ',
        'helps pick the number of clusters first.</p></div>'))
      results$clusteringNote$setVisible(TRUE)
      .cna_hide_outputs(results, keep = choice_items)
      return(invisible())
    }
  }

  model <- .cna_get_model(self, kind)
  if (is.null(model)) {
    .cna_hide_outputs(results)
    return(invisible())
  }

  multi <- inherits(model, "netobject_group")
  ctx <- list(model = model, options = options, results = results, data = self$data,
              kind = kind, spec = spec, multi = multi,
              n_nets = if (multi) length(model) else 1L, seed = .cna_seed(options))
  invisible(lapply(.cna_sections, function(s) .cna_run_section(ctx, s)))
}
