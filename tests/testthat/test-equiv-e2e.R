# End-to-end equivalence through the real jamovi runners.
#
# Each of the five analyses runs through its generated jmvcore runner under
# seeded random option combinations drawn from its .a.yaml (inside the min/max
# limits). Every table it fills must equal a reference built here by calling
# Nestimate and cograph DIRECTLY (never CNA's cna_* engine), so option wiring,
# model building and table assembly are checked together across the option
# space. On top of equality each run checks the table contracts: row caps,
# documented sort orders, filters, the group column ("" and hidden for one
# network, the group or cluster name otherwise), column types, hidden tables.
#
# Rules the references follow:
# * stochastic parts get the run's seed (seed_set is always on): engine `seed=`
#   arguments, and set.seed() around cluster_choice(), as the module does;
# * a section the engine cannot run (a cluster of one unit, clusters with
#   different states) must report the engine's own message in its heading;
#   resampling without units, or on a scaled model, must be refused in a note;
# * rows are compared in the documented order (a total order: ties are broken
#   by the names), except the edge differences, compared per unordered pair.
#
# Section 7 holds minimal repros of the defects the random runs found; they
# assert the intended behaviour and guard against regressions. A test titled
# "(defect)" describes behaviour the module does not have yet.

# ---------------------------------------------------------------------------
# Fixtures (synthetic, deterministic). Equal-sized groups so paired
# permutation tests have matched units.
# ---------------------------------------------------------------------------

e2e_events <- function(n_actor = 15L, per = 12L, groups = c("G1", "G2", "G3"), seed = 101L) {
  withr::with_seed(seed, {
    actors <- sprintf("s%02d", seq_len(n_actor))
    grp <- rep(groups, length.out = n_actor)
    base <- as.POSIXct("2021-03-01 09:00:00", tz = "UTC")
    do.call(rbind, Map(function(a, g) {
      tm <- base + sort(sample.int(3000L, per))
      data.frame(Actor = a, Action = sample(LETTERS[1:5], per, replace = TRUE),
                 Time = tm, Order = seq_len(per), Group = g, stringsAsFactors = FALSE)
    }, actors, grp))
  }) |> (\(d) { rownames(d) <- NULL; d })()
}

E2E_CODES <- paste0("code", 1:5)

e2e_onehot <- function(n_actor = 18L, rows = 4L, groups = c("G1", "G2", "G3"), seed = 202L) {
  withr::with_seed(seed, {
    actors <- sprintf("a%02d", seq_len(n_actor))
    grp <- rep(groups, length.out = n_actor)
    d <- do.call(rbind, Map(function(a, g) {
      n <- 2L * rows
      m <- matrix(stats::rbinom(n * length(E2E_CODES), 1L, 0.45), ncol = length(E2E_CODES),
                  dimnames = list(NULL, E2E_CODES))
      data.frame(m, Actor = a, Session = rep(c("s1", "s2"), each = rows), Group = g,
                 stringsAsFactors = FALSE)
    }, actors, grp))
    rownames(d) <- NULL
    d
  })
}

e2e_delimited <- function(n = 40L, seed = 303L) {
  items <- c("network", "graph", "matrix", "algebra", "model", "cluster")
  tags <- c("graph", "theory", "model", "data")
  withr::with_seed(seed, data.frame(
    doc = sprintf("d%02d", seq_len(n)),
    keywords = vapply(seq_len(n), \(i) paste(sample(items, sample(2:4, 1L)), collapse = "; "),
                      character(1)),
    tags = vapply(seq_len(n), \(i) paste(sample(tags, sample(1:2, 1L)), collapse = "; "),
                  character(1)),
    stringsAsFactors = FALSE))
}

# One item per row with a transaction ID (the "by" layout of the same data);
# `two_fields` adds a second item column with padded, blank and missing items
# and two rows without an ID.
e2e_delimited_long <- function(wide, two_fields = FALSE) {
  items <- strsplit(wide$keywords, "; ", fixed = TRUE)
  long <- data.frame(doc = rep(wide$doc, lengths(items)), item = unlist(items),
                     stringsAsFactors = FALSE)
  if (!two_fields) return(long)
  long$item2 <- rep(c(" theory", "", NA, "data ", "graph"), length.out = nrow(long))
  long$doc[c(3L, 17L)] <- NA
  long
}

# ---------------------------------------------------------------------------
# Option space
# ---------------------------------------------------------------------------

# Every option's default from the analysis's .a.yaml (Variables: NULL).
e2e_yaml_defaults <- function(analysis) {
  opts <- Filter(\(o) !identical(o$type, "Data"), read_def(analysis, "a")$options)
  vals <- lapply(opts, \(o) {
    v <- o$default
    if (identical(o$type, "Integer") && !is.null(v)) v <- as.integer(v)
    v
  })
  stats::setNames(vals, vapply(opts, `[[`, "", "name"))
}

E2E_MEASURE_OPTS <- c("Strength", "Degree", "Betweenness", "Closeness", "Harmonic",
                      "Eigenvector", "PageRank", "Coreness", "Clustering", "Diffusion",
                      "BetweennessRSP", "Eccentricity", "Constraint", "EffectiveSize")

# Table switches: each table of a drawn combination is shown with p = 0.85, so
# the hidden-table contract is exercised too.
E2E_TABLE_FLAGS <- c("buildModel_show_matrix", "buildModel_show_summary", "centrality_show_table",
                     "edgebetweenness_show_table", "community_show_table", "cliques_show_table",
                     "cliques_show_complex", "backbone_show_table", "links_show_table",
                     "rules_show_table", "bootstrap_show_table", "stability_show_table",
                     "permutation_show_table", "compare_show_summary", "compare_show_network",
                     "compare_show_edges", "cluster_show_diagnostics", "cluster_show_choice")

#' One seeded random draw over the option space of `analysis` (respecting the
#' .a.yaml min/max). Returns list(data, opts) with every declared option set.
e2e_draw <- function(analysis, i) withr::with_seed(9100L + i * 37L + 1000L * match(analysis, ANALYSES), {
  pick <- function(x) x[[sample.int(length(x), 1L)]]
  coin <- function(p = 0.5) stats::runif(1L) < p
  o <- list(seed_set = TRUE, seed_value = sample.int(9999L, 1L),
            buildModel_scaling = pick(c("noScaling", "noScaling", "noScaling", "minmax", "max", "rank")),
            centrality_loops = coin(), centrality_normalize = coin(), centrality_invert = coin(0.7),
            edgebetweenness_invert = coin(0.7),
            community_methods = pick(c("louvain", "leiden", "walktrap", "fast_greedy", "infomap",
                                       "label_propagation", "edge_betweenness",
                                       "leading_eigenvector", "spinglass", "optimal")),
            community_gamma = pick(c(0.5, 1, 1.5, 2)),
            cliques_size = pick(2:4), cliques_threshold = pick(c(0, 0, 0.3, 10)),
            backbone_method = pick(c("disparity", "threshold", "lowest")),
            backbone_level = pick(c(0.05, 0.3, 0.5, 0.9)),
            backbone_threshold = pick(c(0, 0.2, 10)), backbone_lowest = pick(c(0.05, 0.2, 0.5)),
            links_method = pick(c("resource_allocation", "adamic_adar", "common_neighbors",
                                  "jaccard", "preferential_attachment", "katz")),
            links_exclude_existing = coin(), links_threshold = pick(c(0, 0.3, 15)),
            links_top_n = pick(c(1L, 3L, 5L, 20L)),
            rules_min_support = pick(c(0.05, 0.1, 0.2, 0.4)),
            rules_min_confidence = pick(c(0.1, 0.3, 0.5, 0.8)),
            rules_min_lift = pick(c(0, 0.8, 1, 1.1)), rules_max_length = pick(2:4),
            rules_max_rows = pick(c(1L, 5L, 15L, 50L)),
            bootstrap_iteration = pick(c(20L, 30L)), bootstrap_level = pick(c(0.05, 0.1)),
            bootstrap_method = pick(c("stability", "threshold")),
            bootstrap_range_low = pick(c(0.6, 0.75, 0.9)), bootstrap_range_up = pick(c(1.1, 1.25, 1.5)),
            bootstrap_set_threshold = coin(), bootstrap_threshold = pick(c(0, 5, 15)),
            bootstrap_ci_method = pick(c("percentile", "basic")),
            bootstrap_table_max_rows = pick(c(1L, 7L, 20L)),
            bootstrap_table_show_all = coin(0.3), bootstrap_table_significant_only = coin(0.3),
            stability_iter = 50L, stability_threshold = pick(c(0.5, 0.7, 0.8)),
            stability_certainty = pick(c(0.9, 0.95)),
            stability_method = pick(c("pearson", "spearman", "kendall")),
            permutation_iter = pick(c(20L, 30L)), permutation_paired = coin(0.3),
            permutation_level = pick(c(0.05, 0.1)),
            permutation_adjust = pick(c("none", "holm", "hochberg", "hommel", "bonferroni",
                                        "BH", "BY", "fdr")),
            permutation_table_max_rows = pick(c(3L, 10L, 50L)),
            permutation_table_show_all = coin(0.3),
            compare_scaling = pick(c("none", "minmax", "max", "rank")),
            compare_edges_max_rows = pick(c(2L, 8L, 50L)),
            clustering_run = TRUE, clustering_k = pick(2:3),
            clustering_dissimilarity = pick(c("hamming", "osa", "lv", "dl", "lcs", "qgram",
                                              "cosine", "jaccard", "jw")),
            clustering_method = pick(c("pam", "ward.D", "ward.D2", "complete", "average",
                                       "single", "mcquitty", "median", "centroid")),
            cluster_choice_kmin = 2L, cluster_choice_kmax = pick(2:4))
  o[E2E_TABLE_FLAGS] <- lapply(E2E_TABLE_FLAGS, \(f) coin(0.85))
  measures <- E2E_MEASURE_OPTS[sample.int(length(E2E_MEASURE_OPTS), sample(1:4, 1L))]
  o[paste0("centrality_", E2E_MEASURE_OPTS)] <- as.list(E2E_MEASURE_OPTS %in% measures)
  pair <- sample.int(o$clustering_k, 2L)
  o$compare_cluster_i <- pair[[1L]]
  o$compare_cluster_j <- pair[[2L]]

  grouped <- coin(0.6)
  levels3 <- c("G1", "G2", "G3")
  gpair <- sample(levels3, 2L)
  data <- switch(analysis,
    CNAEventSingle = , CNAEventCluster = {
      clocks <- pick(c("time", "time", "order", "both"))
      o$buildModel_variables_long_action <- "Action"
      o$buildModel_variables_actor <- if (identical(analysis, "CNAEventCluster") || coin(0.85)) "Actor"
      o$buildModel_variables_long_time <- if (clocks %in% c("time", "both")) "Time"
      o$buildModel_variables_long_order <- if (clocks %in% c("order", "both")) "Order"
      o$buildModel_threshold <- pick(c(120L, 600L, 3000L))
      o$buildModel_variables_group <- if (grouped) "Group"
      e2e_events()
    },
    CNABinarySingle = , CNABinaryCluster = {
      units <- if (identical(analysis, "CNABinaryCluster")) pick(c("both", "actor"))
               else pick(c("both", "both", "actor", "session", "none"))
      o$buildModel_variables_onehot <- if (coin(0.7)) E2E_CODES else
        sort(E2E_CODES[sample.int(5L, sample(3:4, 1L))])
      o$buildModel_variables_actor <- if (units %in% c("both", "actor")) "Actor"
      o$buildModel_variables_session <- if (units %in% c("both", "session")) "Session"
      o$buildModel_window <- pick(1:3)
      o$buildModel_window_type <- pick(c("tumbling", "sliding"))
      o$buildModel_variables_group <- if (grouped) "Group"
      e2e_onehot()
    },
    CNADelimited = {
      layout <- pick(c("sep", "sep", "two_fields", "by", "by_two"))
      o$buildModel_similarity <- pick(c("none", "none", "jaccard", "cosine", "inclusion",
                                        "association", "dice", "equivalence"))
      o$buildModel_min_occur <- pick(1:3)
      o$buildModel_diagonal <- coin()
      o$buildModel_edge_threshold <- pick(c(0, 0, 0.2, 2))
      o$buildModel_top_n <- pick(c(0L, 0L, 4L, 8L))
      o$buildModel_sep <- ";"
      wide <- e2e_delimited()
      if (layout %in% c("by", "by_two")) {
        o$buildModel_variables_field <- if (identical(layout, "by_two")) c("item", "item2") else "item"
        o$buildModel_variables_by <- "doc"
        e2e_delimited_long(wide, two_fields = identical(layout, "by_two"))
      } else {
        o$buildModel_variables_field <- if (identical(layout, "two_fields")) c("keywords", "tags") else "keywords"
        o$buildModel_variables_by <- NULL
        wide
      }
    })
  o$compare_group_i <- if (isTRUE(grouped) && !is.null(o$buildModel_variables_group)) gpair[[1L]]
  o$compare_group_j <- if (isTRUE(grouped) && !is.null(o$buildModel_variables_group)) gpair[[2L]]
  declared <- names(e2e_yaml_defaults(analysis))
  full <- e2e_yaml_defaults(analysis)
  for_analysis <- o[intersect(names(o), declared)]
  full[names(for_analysis)] <- for_analysis
  # Variable options left unassigned are NULL (the runner has no default for them)
  vars <- grep("_variables_|^compare_group_", declared, value = TRUE)
  full[setdiff(vars, names(for_analysis))] <- list(NULL)
  list(data = data, opts = full)
})

# Diagnosed third-party warnings, muffled by message and nothing else:
# Nestimate's advice for one-hot data without windows or sessions, the ggplot2
# deprecation jmvcore triggers when a results object is created, igraph's
# note on weighted edge-betweenness communities, and Nestimate's notes that a
# measure or the edge weights do not vary (no stability to assess).
E2E_KNOWN_WARNINGS <- c("produces sparse networks", "argument of `element_line()` is deprecated",
                        "Membership vector will be selected based on the highest modularity",
                        "All centrality measures have zero variance",
                        "Original edge vector has zero variance")
e2e_known <- function(expr) withCallingHandlers(expr,
  warning = function(w)
    if (any(vapply(E2E_KNOWN_WARNINGS, \(k) grepl(k, conditionMessage(w), fixed = TRUE), logical(1))))
      invokeRestart("muffleWarning"),
  # Nestimate's note for a single-sequence network (the runner refuses
  # resampling itself, with its own message)
  message = function(m)
    if (grepl("network with one long sequence", conditionMessage(m), fixed = TRUE))
      invokeRestart("muffleMessage"))

#' A full option list for a hand-made run: .a.yaml defaults, every table shown,
#' small iteration counts, a seed, then `...` (names must be declared options).
e2e_opts <- function(analysis, ..., tables = TRUE) {
  full <- e2e_yaml_defaults(analysis)
  vars <- grep("_variables_|^compare_group_", names(full), value = TRUE)
  full[vars] <- list(NULL)
  flags <- intersect(E2E_TABLE_FLAGS, names(full))
  full[flags] <- list(isTRUE(tables))
  small <- list(seed_set = TRUE, seed_value = 11L, bootstrap_iteration = 20L,
                permutation_iter = 20L, stability_iter = 50L, clustering_run = TRUE,
                cluster_choice_kmax = 3L)
  small <- small[intersect(names(small), names(full))]
  full[names(small)] <- small
  over <- list(...)
  stopifnot("every override is a declared option" = all(names(over) %in% names(full)))
  full[names(over)] <- over
  full
}

#' Every numeric option inside its .a.yaml min/max.
e2e_in_range <- function(analysis, opts) {
  defs <- Filter(\(d) d$type %in% c("Integer", "Number"), read_def(analysis, "a")$options)
  all(vapply(defs, \(d) {
    v <- opts[[d$name]]
    (is.null(d$min) || v >= d$min) && (is.null(d$max) || v <= d$max)
  }, logical(1)))
}

#' Run an analysis through its generated jmvcore runner with every option set.
e2e_run <- function(analysis, data, opts)
  e2e_known(do.call(cna(analysis), c(list(data = data), opts)))

# ---------------------------------------------------------------------------
# Reference: networks built with direct Nestimate calls.
# ---------------------------------------------------------------------------

e2e_scaling <- function(o) if (identical(o$buildModel_scaling, "noScaling")) NULL else o$buildModel_scaling
e2e_mode <- function(o) if (identical(o$buildModel_window_type, "sliding")) "overlapping" else "non-overlapping"
e2e_compact <- function(x) Filter(Negate(is.null), x)

# Plain named list of the networks of a model ("" for a single network).
e2e_nets <- function(model) {
  if (inherits(model, "netobject_group"))
    return(stats::setNames(lapply(seq_along(model), \(k) model[[k]]), names(model)))
  stats::setNames(list(model), "")
}

# Event Actor and Group reach Nestimate as text (the module documents this, and
# it fixes the unit order resampling sees whatever the column's storage type).
e2e_event_args <- function(data, o) {
  cols <- unlist(o[c("buildModel_variables_actor", "buildModel_variables_long_time",
                     "buildModel_variables_long_action", "buildModel_variables_long_order")])
  data <- e2e_as_text(data, o$buildModel_variables_actor)
  e2e_compact(list(data = data[cols], action = o$buildModel_variables_long_action,
                   actor = o$buildModel_variables_actor, time = o$buildModel_variables_long_time,
                   order = o$buildModel_variables_long_order,
                   time_threshold = o$buildModel_threshold))
}

e2e_as_text <- function(data, cols) {
  data[cols] <- lapply(data[cols], as.character)
  data
}

#' Rows with a Group value, stably sorted by group, the Group as text. The
#' group order (which Nestimate takes from first appearance) is the factor's
#' levels, else the values in numeric order when all are numbers, else in
#' radix (C-locale) order; never the row order.
e2e_by_group <- function(data, group) {
  raw <- data[[group]]
  txt <- as.character(raw)
  data <- data[!is.na(txt) & nzchar(trimws(txt)), , drop = FALSE]
  raw <- data[[group]]
  txt <- as.character(raw)
  values <- unique(txt)
  as_number <- utils::type.convert(values, as.is = TRUE)
  levels <- if (is.factor(raw)) intersect(levels(raw), values)
            else if (is.numeric(as_number)) values[order(as_number)]
            else values[order(values, method = "radix")]
  out <- data[order(match(txt, levels), method = "radix"), , drop = FALSE]
  out[[group]] <- as.character(out[[group]])
  out
}

# One-hot Actor and Session reach Nestimate as text, rows sorted by group.
e2e_onehot_args <- function(data, o, with_group = TRUE) {
  group <- if (with_group) o$buildModel_variables_group
  data <- e2e_as_text(data, c(o$buildModel_variables_actor, o$buildModel_variables_session))
  if (!is.null(group)) data <- e2e_by_group(data, group)
  cols <- c(o$buildModel_variables_onehot, o$buildModel_variables_actor,
            o$buildModel_variables_session, group)
  e2e_compact(list(data = data[cols], codes = o$buildModel_variables_onehot, format = "onehot",
                   actor = o$buildModel_variables_actor, session = o$buildModel_variables_session,
                   window_size = o$buildModel_window, mode = e2e_mode(o), group = group))
}

#' The clustering of the analysis and its unit data (event sequences, or the
#' one-hot network whose actor x session units are clustered).
e2e_cluster_input <- function(analysis, data, o) {
  if (identical(analysis, "CNAEventCluster")) {
    args <- e2e_event_args(data, o)
    return(do.call(Nestimate::prepare, args)$sequence_data)
  }
  do.call(Nestimate::build_network, c(e2e_onehot_args(data, o, with_group = FALSE),
                                      method = "co_occurrence"))
}

e2e_clustering <- function(input, o)
  Nestimate::build_clusters(input, k = o$clustering_k,
                            dissimilarity = o$clustering_dissimilarity,
                            method = o$clustering_method, seed = o$seed_value)

#' Cluster networks as groups: rows labelled "Cluster k" and sorted by k (so
#' the networks come as Cluster 1, 2, ...), every cluster sharing the states,
#' the scaling applied.
e2e_cluster_groups <- function(args, labels, scaling) {
  k <- as.integer(sub("^Cluster ", "", labels))
  args$data$.cluster <- labels
  args$data <- args$data[order(k, method = "radix"), , drop = FALSE]
  do.call(Nestimate::build_network, c(args, e2e_compact(list(
    group = ".cluster", method = "co_occurrence", scaling = scaling))))
}

#' The reference model: a netobject or a netobject_group.
e2e_ref_model <- function(analysis, data, o) {
  scaling <- e2e_scaling(o)
  switch(analysis,
    CNAEventSingle = {
      group <- o$buildModel_variables_group
      if (!is.null(group)) data <- e2e_by_group(data, group)
      args <- e2e_event_args(data, o)
      if (!is.null(group)) {
        args$data[[group]] <- data[[group]]
        args$group <- group
      }
      do.call(Nestimate::build_network, c(args, e2e_compact(list(method = "co_occurrence",
                                                                 scaling = scaling))))
    },
    CNABinarySingle = do.call(Nestimate::build_network,
                              c(e2e_onehot_args(data, o), e2e_compact(list(
                                method = "co_occurrence", scaling = scaling)))),
    CNAEventCluster = {
      seqs <- e2e_cluster_input(analysis, data, o)
      cl <- e2e_clustering(seqs, o)
      model <- e2e_cluster_groups(list(data = seqs), paste("Cluster", cl$assignments), scaling)
      structure(model, e2e_clustering = cl)
    },
    CNABinaryCluster = {
      # each row takes its actor x session unit's cluster (assignments follow
      # the units' first-seen order), so the windows stay the transactions and
      # units, as in Binary Single. Sliding windows: a unit with fewer rows than
      # the window has no complete window, is not clustered and is left out.
      cl <- e2e_clustering(e2e_cluster_input(analysis, data, o), o)
      args <- e2e_onehot_args(data, o, with_group = FALSE)
      unit_cols <- c(o$buildModel_variables_actor, o$buildModel_variables_session)
      key <- do.call(paste, c(unname(as.list(args$data[unit_cols])), sep = "\r"))
      units <- unique(key)
      if (identical(e2e_mode(o), "overlapping")) {
        rows <- tabulate(match(key, units), length(units))
        units <- units[rows >= o$buildModel_window]
        args$data <- args$data[key %in% units, , drop = FALSE]
        key <- key[key %in% units]
      }
      model <- e2e_cluster_groups(args, paste("Cluster", cl$assignments[match(key, units)]), scaling)
      structure(model, e2e_clustering = cl)
    },
    CNADelimited = {
      by <- o$buildModel_variables_by
      fields <- o$buildModel_variables_field
      top_n <- o$buildModel_top_n
      df <- data[c(fields, by)]
      if (!is.null(by)) {
        # a Transaction ID: the items of every item column pooled into one
        # column (named like the first), trimmed, blank / NA items and rows
        # without an ID dropped
        items <- trimws(unlist(lapply(data[fields], as.character), use.names = FALSE))
        ids <- rep(data[[by]], length(fields))
        keep <- !is.na(ids) & !is.na(items) & nzchar(items)
        df <- stats::setNames(data.frame(ids[keep], items[keep], stringsAsFactors = FALSE),
                              c(by, fields[[1L]]))
        fields <- fields[[1L]]
      }
      Nestimate::cooccurrence(
        df, field = fields, by = by, sep = if (is.null(by)) o$buildModel_sep,
        similarity = o$buildModel_similarity, threshold = o$buildModel_edge_threshold,
        min_occur = o$buildModel_min_occur, diagonal = o$buildModel_diagonal,
        top_n = if (top_n > 0L) top_n)
    })
}

# ---------------------------------------------------------------------------
# Reference tables: direct Nestimate / cograph calls on the reference model,
# assembled in the documented row order and capped as the options ask.
# ---------------------------------------------------------------------------

# Checkbox -> cograph measure (what each centrality checkbox promises).
E2E_COGRAPH <- c(Strength = "strength", Degree = "degree", Betweenness = "betweenness",
                 Closeness = "closeness", Harmonic = "harmonic", Eigenvector = "eigenvector",
                 PageRank = "pagerank", Coreness = "coreness", Clustering = "transitivity",
                 Diffusion = "diffusion", BetweennessRSP = "rsp_betweenness",
                 Eccentricity = "eccentricity", Constraint = "constraint",
                 EffectiveSize = "effective_size")
E2E_TITLES <- c(Strength = "Strength", Degree = "Degree", Betweenness = "Betweenness",
                Closeness = "Closeness", Harmonic = "Harmonic closeness",
                Eigenvector = "Eigenvector", PageRank = "PageRank", Coreness = "Coreness",
                Clustering = "Clustering", Diffusion = "Diffusion",
                BetweennessRSP = "RSP betweenness", Eccentricity = "Eccentricity",
                Constraint = "Constraint", EffectiveSize = "Effective size")

e2e_measures <- function(o) E2E_MEASURE_OPTS[vapply(E2E_MEASURE_OPTS, \(m)
  isTRUE(o[[paste0("centrality_", m)]]), logical(1))]

# cograph's column for a measure ("_all" for undirected ones).
e2e_col <- function(ct, m) intersect(c(m, paste0(m, "_all")), names(ct))[[1L]]

# Row-bind per-network frames with a leading group column (NULL when empty).
e2e_bind <- function(nets, fun) {
  parts <- Filter(\(d) !is.null(d) && nrow(d) > 0L,
                  unname(Map(\(net, g) {
                    d <- fun(net, g)
                    if (is.null(d) || !nrow(d)) return(NULL)
                    data.frame(group = g, d, check.names = FALSE, stringsAsFactors = FALSE)
                  }, nets, names(nets))))
  if (!length(parts)) return(NULL)
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

# Unique undirected non-zero off-diagonal edges in matrix order.
e2e_edges <- function(w) {
  w <- unclass(w)
  idx <- which(upper.tri(w) & w != 0, arr.ind = TRUE)
  data.frame(from = rownames(w)[idx[, 1L]], to = colnames(w)[idx[, 2L]],
             weight = as.numeric(w[idx]), stringsAsFactors = FALSE)
}

e2e_sort <- function(df, ...) {
  if (is.null(df)) return(NULL)
  out <- df[do.call(order, list(...)), , drop = FALSE]
  rownames(out) <- NULL
  out
}

e2e_quiet <- function(expr) withCallingHandlers(expr, message = function(m)
  if (grepl("Weights inverted", conditionMessage(m), fixed = TRUE) ||
      grepl("centralities computed excluding loops", conditionMessage(m), fixed = TRUE))
    invokeRestart("muffleMessage"))

# Each reference table: a data.frame (expected rows, documented order) or NULL
# (no rows: the table must explain itself in a note).
#' Centralities of a network (netobject or weight matrix) as the checkboxes
#' promise them: a named list, one numeric vector per measure, in node order.
#' cograph computes each; two measures have their own rule:
#' * Eccentricity with inverted weights is a path length, so it is computed
#'   on the matrix of inverse weights (1 / w, 0 where there is no edge) with
#'   the weights taken as they are;
#' * Effective size never counts self-loops.
e2e_cent <- function(x, ms, loops, normalized, invert) {
  w <- if (is.matrix(x)) x else unclass(x$weights)
  calc <- function(input, m, lp, inv) {
    args <- list(input, measures = E2E_COGRAPH[[m]], loops = lp, normalized = normalized,
                 invert_weights = inv)
    if (is.matrix(input)) args$directed <- FALSE
    ct <- e2e_quiet(do.call(cograph::centrality, args))
    stats::setNames(as.numeric(ct[[e2e_col(ct, E2E_COGRAPH[[m]])]]), as.character(ct$node))
  }
  # Diffusion is tna's, range-normalised to 0-1 (Nestimate's default); the
  # stability reference takes it from Nestimate's built-in, so a matrix never
  # asks for it here
  diffusion <- function() {
    d <- e2e_quiet(Nestimate::net_centrality(x, measures = "Diffusion", loops = loops,
                                             normalize_diffusion = TRUE))
    stats::setNames(as.numeric(d$Diffusion), as.character(d$state))
  }
  # Clustering: Zhang & Horvath on w / max(w) = cograph's "onnela" transitivity
  # (tna's formula on w + t(w)) divided by 2 x the largest off-diagonal weight
  clustering <- function() {
    off <- w
    diag(off) <- 0
    ct <- e2e_quiet(do.call(cograph::centrality, c(list(x, measures = "transitivity",
      transitivity_type = "onnela"), if (is.matrix(x)) list(directed = FALSE))))
    v <- as.numeric(ct[[e2e_col(ct, "transitivity")]])
    stats::setNames(if (max(off) > 0) v / (2 * max(off)) else rep(NaN, length(v)),
                    as.character(ct$node))
  }
  values <- lapply(ms, \(m) {
    if (identical(m, "Diffusion")) diffusion()
    else if (identical(m, "Clustering")) clustering()
    else if (identical(m, "Eccentricity") && isTRUE(invert)) calc(ifelse(w > 0, 1 / w, 0), m, loops, FALSE)
    else if (identical(m, "EffectiveSize")) calc(x, m, FALSE, invert)
    else calc(x, m, loops, invert)
  })
  lapply(stats::setNames(values, ms), \(v) unname(v[rownames(w)]))
}

e2e_ref_centrality <- function(nets, o) {
  ms <- e2e_measures(o)
  e2e_bind(nets, \(net, g) {
    ct <- e2e_cent(net, ms, loops = o$centrality_loops, normalized = o$centrality_normalize,
                   invert = o$centrality_invert)
    data.frame(state = rownames(net$weights), ct, check.names = FALSE, stringsAsFactors = FALSE)
  })
}

e2e_ref_matrix <- function(nets) {
  states <- unique(unlist(lapply(nets, \(n) rownames(n$weights))))
  e2e_bind(nets, \(net, g) {
    w <- unclass(net$weights)
    full <- matrix(NA_real_, nrow(w), length(states), dimnames = list(rownames(w), states))
    full[, colnames(w)] <- w
    data.frame(state = rownames(w), full, check.names = FALSE, stringsAsFactors = FALSE)
  })
}

e2e_ref_summary <- function(nets) e2e_bind(nets, \(net, g) {
  s <- summary(net)
  data.frame(metric = as.character(s$metric), value = as.numeric(s$value),
             stringsAsFactors = FALSE)
})

e2e_ref_edge_betweenness <- function(nets, o) e2e_bind(nets, \(net, g) {
  eb <- Nestimate::net_edge_betweenness(net, invert = o$edgebetweenness_invert)
  e <- e2e_edges(net$weights)
  if (!nrow(e)) return(NULL)
  e$betweenness <- as.numeric(unclass(eb$weights)[cbind(e$from, e$to)])
  e2e_sort(e, -e$betweenness, e$from, e$to)
})

e2e_ref_community <- function(nets, o) e2e_bind(nets, \(net, g) {
  cm <- if (identical(o$community_methods, "spinglass"))
    cograph::communities(net, method = "spinglass", gamma = o$community_gamma, seed = o$seed_value)
  else cograph::communities(net, method = o$community_methods, resolution = o$community_gamma,
                            seed = o$seed_value)
  data.frame(state = as.character(cm$node), community = as.integer(cm$community),
             stringsAsFactors = FALSE)
})

e2e_ref_cliques <- function(nets, o) e2e_bind(nets, \(net, g) {
  sc <- Nestimate::build_simplicial(net, type = "clique", threshold = o$cliques_threshold,
                                    max_dim = o$cliques_size - 1L)
  members <- Filter(\(s) length(s) == o$cliques_size, sc$simplices)
  if (!length(members)) return(NULL)
  w <- unclass(net$weights)
  pair_w <- lapply(members, \(s) { sub <- w[sc$nodes[s], sc$nodes[s]]; sub[upper.tri(sub)] })
  d <- data.frame(members = vapply(members, \(s) paste(sc$nodes[s], collapse = ", "), ""),
                  mean_weight = vapply(pair_w, mean, 0), min_weight = vapply(pair_w, min, 0),
                  stringsAsFactors = FALSE)
  d <- e2e_sort(d, -d$mean_weight, d$members)
  data.frame(clique = seq_len(nrow(d)), d)
})

e2e_ref_complex <- function(nets, o) e2e_bind(nets, \(net, g) {
  size <- max(o$cliques_size, nrow(net$weights))
  sc <- Nestimate::build_simplicial(net, type = "clique", threshold = o$cliques_threshold,
                                    max_dim = size - 1L)
  betti <- withCallingHandlers(Nestimate::betti_numbers(sc),
                               warning = function(w) invokeRestart("muffleWarning"))
  f <- as.integer(sc$f_vector)
  data.frame(metric = c(sprintf("Cliques of %d node%s", seq_along(f),
                                ifelse(seq_along(f) == 1L, "", "s")),
                        "Euler characteristic", sprintf("Betti number b%d", seq_along(betti) - 1L),
                        "Maximum q (Q-analysis)"),
             value = c(f, Nestimate::euler_characteristic(sc), as.numeric(betti),
                       as.numeric(Nestimate::q_analysis(sc)$max_q)),
             stringsAsFactors = FALSE)
})

e2e_ref_backbone <- function(nets, o) e2e_bind(nets, \(net, g) {
  p <- unclass(Nestimate::net_prune(net, method = o$backbone_method, threshold = o$backbone_threshold,
                                    lowest = o$backbone_lowest, level = o$backbone_level)$weights)
  e <- e2e_edges(net$weights)
  if (!nrow(e)) return(NULL)
  # undirected: an edge survives when either direction does
  kept <- p[cbind(e$from, e$to)] != 0 | p[cbind(e$to, e$from)] != 0
  e$status <- ifelse(kept, "kept", "removed")
  e2e_sort(e, e$status, -e$weight, e$from, e$to)
})

e2e_ref_links <- function(nets, o) e2e_bind(nets, \(net, g) {
  w <- unclass(net$weights)
  if (nrow(w) < 2L) return(NULL)                 # one node: no pair to score
  w[abs(w) <= o$links_threshold] <- 0
  diag(w) <- 0                                   # a node is not its own neighbour
  # Jaccard compares neighbour sets: binary adjacency; the others use weights
  p <- Nestimate::predict_links(w, methods = o$links_method,
                                weighted = !identical(o$links_method, "jaccard"),
                                top_n = o$links_top_n,
                                exclude_existing = o$links_exclude_existing)$predictions
  if (is.null(p) || !nrow(p)) return(NULL)
  data.frame(from = as.character(p$from), to = as.character(p$to), score = as.numeric(p$score),
             rank = as.integer(p$rank), existing = ifelse(p$existing, "yes", "no"),
             stringsAsFactors = FALSE)
})

e2e_ref_rules <- function(nets, o) e2e_bind(nets, \(net, g) {
  r <- Nestimate::association_rules(net, min_support = o$rules_min_support,
                                    min_confidence = o$rules_min_confidence,
                                    min_lift = o$rules_min_lift, max_length = o$rules_max_length)$rules
  if (is.null(r) || !nrow(r)) return(NULL)
  r <- e2e_sort(r, -r$lift, -r$support, r$antecedent, r$consequent)
  r <- utils::head(r, o$rules_max_rows)
  data.frame(antecedent = as.character(r$antecedent), consequent = as.character(r$consequent),
             support = r$support, confidence = r$confidence, lift = r$lift,
             conviction = r$conviction, count = as.integer(r$count), stringsAsFactors = FALSE)
})

e2e_ref_bootstrap <- function(model, nets, o) {
  threshold <- if (identical(o$bootstrap_method, "threshold") && isTRUE(o$bootstrap_set_threshold))
    o$bootstrap_threshold
  bs <- Nestimate::bootstrap_network(model, iter = o$bootstrap_iteration, ci_level = o$bootstrap_level,
                                     inference = o$bootstrap_method,
                                     consistency_range = c(o$bootstrap_range_low, o$bootstrap_range_up),
                                     edge_threshold = threshold, ci_method = o$bootstrap_ci_method,
                                     seed = o$seed_value)
  parts <- if (inherits(model, "netobject_group")) stats::setNames(lapply(names(model), \(g) bs[[g]]), names(model))
           else list(bs) |> stats::setNames("")
  d <- e2e_bind(parts, \(b, g) {
    s <- b$summary
    data.frame(from = as.character(s$from), to = as.character(s$to), weight = s$weight,
               # threshold inference has no consistency range: the columns stay empty
               p_value = s$p_value, cr_lower = s$cr_lower %||% NA_real_,
               cr_upper = s$cr_upper %||% NA_real_,
               ci_lower = s$ci_lower, ci_upper = s$ci_upper,
               significant = ifelse(s$sig %in% TRUE, "Yes", "No"), stringsAsFactors = FALSE)
  })
  d <- e2e_sort(d, d$significant != "Yes", d$p_value, -abs(d$weight), d$group, d$from, d$to)
  if (isTRUE(o$bootstrap_table_significant_only)) d <- d[d$significant == "Yes", , drop = FALSE]
  if (!isTRUE(o$bootstrap_table_show_all)) d <- utils::head(d, o$bootstrap_table_max_rows)
  if (!nrow(d)) NULL else { rownames(d) <- NULL; d }
}

e2e_ref_permutation <- function(model, o) {
  pt <- Nestimate::permutation(model, iter = o$permutation_iter, alpha = o$permutation_level,
                               paired = o$permutation_paired, adjust = "none", seed = o$seed_value)
  cmp <- stats::setNames(lapply(names(pt), \(nm) pt[[nm]]), names(pt))
  d <- e2e_bind(cmp, \(x, nm) {
    s <- x$summary
    # one test per unique undirected edge: the correction runs over exactly those
    data.frame(edge_name = paste(s$from, "–", s$to), diff_true = s$diff,
               effect_size = s$effect_size,
               p_value = stats::p.adjust(s$p_value, method = o$permutation_adjust),
               stringsAsFactors = FALSE)
  })
  names(d)[names(d) == "group"] <- "group_comparison"
  d <- e2e_sort(d, d$p_value, -abs(d$diff_true), d$group_comparison, d$edge_name)
  if (!isTRUE(o$permutation_table_show_all)) d <- utils::head(d, o$permutation_table_max_rows)
  d
}

e2e_compare_pair <- function(analysis, model, o) {
  if (analysis %in% c("CNAEventCluster", "CNABinaryCluster"))
    return(c(o$compare_cluster_i, o$compare_cluster_j))
  match(c(o$compare_group_i, o$compare_group_j), names(model))
}

e2e_ref_compare <- function(analysis, model, o) {
  ij <- e2e_compare_pair(analysis, model, o)
  cmp <- Nestimate::compare_model(model, i = ij[[1L]], j = ij[[2L]], scaling = o$compare_scaling)
  e <- cmp$edge_metrics
  # unique unordered pairs, diagonal included
  key <- paste(pmin(e$source, e$target), pmax(e$source, e$target), sep = "|")
  e <- e[!duplicated(key), , drop = FALSE]
  list(summary = data.frame(category = as.character(cmp$summary_metrics$category),
                            metric = as.character(cmp$summary_metrics$metric),
                            value = as.numeric(cmp$summary_metrics$value), stringsAsFactors = FALSE),
       network = data.frame(metric = as.character(cmp$network_metrics$metric),
                            group_i = as.numeric(cmp$network_metrics$x),
                            group_j = as.numeric(cmp$network_metrics$y), stringsAsFactors = FALSE),
       edges = data.frame(key = paste(pmin(e$source, e$target), pmax(e$source, e$target), sep = "|"),
                          weight_i = e$weight_x, weight_j = e$weight_y,
                          difference = e$raw_difference,
                          relative_difference = e$relative_difference, stringsAsFactors = FALSE),
       title = paste("Comparing:", names(model)[ij[[1L]]], "vs", names(model)[ij[[2L]]]))
}

e2e_attempt <- function(expr) tryCatch(expr, error = function(e) e)

e2e_ref_stability <- function(model, o) {
  ms <- e2e_measures(o)
  # Diffusion through Nestimate's built-in (tna's, range-normalised, as the
  # table); every other measure through centrality_fn on each resampled
  # matrix, computed as the table computes it (normalisation does not apply)
  own <- setdiff(ms, "Diffusion")
  fn <- if (length(own)) function(mat) stats::setNames(
    e2e_cent(mat, own, loops = o$centrality_loops, normalized = FALSE, invert = o$centrality_invert),
    paste0("cg_", own))
  cs <- e2e_attempt(Nestimate::centrality_stability(
    model, measures = ifelse(ms == "Diffusion", "Diffusion", paste0("cg_", ms)),
    iter = o$stability_iter, threshold = o$stability_threshold,
    certainty = o$stability_certainty, method = o$stability_method, centrality_fn = fn,
    loops = o$centrality_loops, normalize_diffusion = TRUE, seed = o$seed_value))
  cd <- e2e_attempt(Nestimate::casedrop_reliability(
    model, iter = o$stability_iter, threshold = o$stability_threshold,
    certainty = o$stability_certainty, method = o$stability_method, seed = o$seed_value))
  nr <- e2e_attempt(Nestimate::network_reliability(model, iter = o$stability_iter, seed = o$seed_value))
  grouped <- inherits(model, "netobject_group")
  nets <- e2e_nets(model)
  parts <- function(x) if (grouped) stats::setNames(lapply(names(model), \(g) x[[g]]), names(model))
                       else stats::setNames(list(x), "")
  # a measure constant across the nodes of the full network has no stability:
  # it is left out of that network's rows
  constant <- lapply(nets, \(net) {
    v <- e2e_cent(net, ms, loops = o$centrality_loops, normalized = FALSE, invert = o$centrality_invert)
    ms[vapply(v, \(x) { x <- x[is.finite(x)]; length(x) < 2L || isTRUE(all.equal(min(x), max(x))) },
              logical(1))]
  })
  returned <- if (!inherits(cs, "error")) Map(\(x, g) {
    m <- sub("^cg_", "", names(x$cs))
    keep <- !m %in% constant[[g]]
    stats::setNames(as.numeric(x$cs)[keep], unname(E2E_TITLES[m[keep]]))
  }, parts(cs), names(nets))
  cs_rows <- if (!is.null(returned)) e2e_bind(returned, \(v, g)
    if (length(v)) data.frame(measure = names(v), cs = unname(v), stringsAsFactors = FALSE))
  cd_rows <- if (!inherits(cd, "error")) e2e_bind(parts(cd), \(x, g)
    data.frame(measure = "Edge weights", cs = as.numeric(x$cs), stringsAsFactors = FALSE))
  rows <- rbind(cs_rows, cd_rows)
  if (!is.null(rows)) rows <- e2e_sort(rows, match(rows$group, names(nets)), rows$measure == "Edge weights")
  rel <- if (!inherits(nr, "error") && nrow(nr$summary)) {
    s <- nr$summary
    # one network: group "" like every other table (Nestimate names it by method)
    data.frame(group = if (grouped) as.character(s$model) else "",
               metric = as.character(s$metric), mean = as.numeric(s$mean),
               sd = as.numeric(s$sd), stringsAsFactors = FALSE)
  }
  list(stability = rows, reliability = rel,
       failed = c(cs = inherits(cs, "error"), cd = inherits(cd, "error"), nr = inherits(nr, "error")),
       # titles missing from at least one network's rows (the "constant" note)
       missing = if (is.null(returned)) character() else
         setdiff(unname(E2E_TITLES[ms]), Reduce(intersect, lapply(returned, names))))
}

e2e_ref_cluster_diag <- function(cl) {
  d <- as.data.frame(Nestimate::cluster_diagnostics(cl))
  data.frame(cluster = paste("Cluster", d$cluster), size = as.integer(d$size),
             pct = as.numeric(d$pct), mean_within_dist = as.numeric(d$mean_within_dist),
             silhouette = as.numeric(d$sil_mean), stringsAsFactors = FALSE)
}

e2e_ref_cluster_choice <- function(analysis, data, o) {
  cc <- withr::with_seed(o$seed_value, Nestimate::cluster_choice(
    e2e_cluster_input(analysis, data, o), k = seq.int(o$cluster_choice_kmin, o$cluster_choice_kmax),
    dissimilarity = o$clustering_dissimilarity, method = o$clustering_method))
  data.frame(k = as.integer(cc$k), silhouette = as.numeric(cc$silhouette),
             mean_within_dist = as.numeric(cc$mean_within_dist),
             min_size = as.integer(cc$min_size), max_size = as.integer(cc$max_size),
             size_ratio = as.numeric(cc$size_ratio), stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Comparing a run with its reference
# ---------------------------------------------------------------------------

# Tables that say in a note why they have no rows (others may be empty by a filter).
E2E_EMPTY_NOTE <- c("edgeBetweennessTable", "cliquesTable", "backboneTable", "linksTable",
                    "rulesTable")

e2e_notes <- function(res, name) vapply(res[[name]]$notes, \(n) n$note, character(1))

# A table as the runner filled it, with NaN read as NA and plain row names.
e2e_table <- function(res, name) {
  d <- table_df(res, name)
  d[] <- lapply(d, \(x) { if (is.double(x)) x[is.nan(x)] <- NA_real_; x })
  rownames(d) <- NULL
  d
}

e2e_norm <- function(d) {
  d[] <- lapply(d, \(x) { if (is.factor(x)) x <- as.character(x); if (is.double(x)) x[is.nan(x)] <- NA_real_; x })
  rownames(d) <- NULL
  d
}

# Column types promised by the .r.yaml (dynamic columns are numbers).
e2e_types_ok <- function(analysis, name, d) {
  item <- Filter(\(x) identical(x$name, name), read_def(analysis, "r")$items)[[1L]]
  declared <- stats::setNames(vapply(item$columns, \(c) c$type %||% "number", ""),
                              vapply(item$columns, `[[`, "", "name"))
  types <- declared[names(d)]
  types[is.na(types)] <- "number"
  all(vapply(seq_along(d), \(k) switch(types[[k]], text = is.character(d[[k]]),
                                       integer = is.integer(d[[k]]),
                                       number = is.numeric(d[[k]]), FALSE), logical(1)))
}

#' Assert one table against its expected outcome.
#' @param expect list(kind = "hidden" | "rows" | "empty" | "refused" | "error", ...)
e2e_expect <- function(res, analysis, name, expect, label) {
  lab <- paste(label, name)
  tbl <- res[[name]]
  switch(expect$kind,
    hidden = expect_false(isTRUE(tbl$visible), label = paste(lab, "hidden")),
    error = expect_match(res[[expect$heading]]$content,
                         paste("could not be completed:", conditionMessage(expect$error)),
                         fixed = TRUE, label = paste(lab, "heading reports the engine error")),
    refused = {
      expect_equal(tbl$rowCount, 0, label = paste(lab, "refused: no rows"))
      expect_match(paste(e2e_notes(res, name), collapse = " "), expect$pattern,
                   label = paste(lab, "refusal note"))
    },
    empty = {
      expect_equal(tbl$rowCount, 0, label = paste(lab, "no rows"))
      if (name %in% E2E_EMPTY_NOTE)
        expect_gt(length(tbl$notes), 0, label = paste(lab, "explains why it is empty"))
    },
    rows = {
      expect_true(isTRUE(tbl$visible), label = paste(lab, "visible"))
      got <- e2e_table(res, name)
      want <- e2e_norm(expect$df)
      if (!nrow(got)) {
        fail(sprintf("%s: expected %d rows, the table is empty", lab, nrow(want)))
        return(invisible())
      }
      if (!is.null(expect$rename)) names(got) <- expect$rename(names(got))
      # one network: the group column is hidden (asDF leaves it out)
      if ("group" %in% names(want) && !"group" %in% names(got) && all(want$group == "")) want$group <- NULL
      expect_true(e2e_types_ok(analysis, name, table_df(res, name)), label = paste(lab, "column types"))
      expect_setequal(names(got), names(want))
      common <- intersect(names(want), names(got))
      expect_equal(got[common], want[common], tolerance = 1e-8, label = lab)
    })
}

e2e_outcome <- function(on, compute, heading) {
  if (!isTRUE(on)) return(list(kind = "hidden"))
  val <- e2e_attempt(compute())
  if (inherits(val, "error")) return(list(kind = "error", error = val, heading = heading))
  if (is.null(val) || !nrow(val)) list(kind = "empty") else list(kind = "rows", df = val)
}

e2e_units <- function(analysis, o) switch(analysis,
  CNAEventSingle = , CNAEventCluster = !is.null(o$buildModel_variables_actor),
  CNABinarySingle = , CNABinaryCluster = !is.null(o$buildModel_variables_actor) ||
    !is.null(o$buildModel_variables_session),
  FALSE)

# Why resampling is refused for these options (regex for the note), or NULL.
e2e_refusal <- function(analysis, o) {
  if (!e2e_units(analysis, o)) return("resamples units of the data, so it needs")
  if (!is.null(e2e_scaling(o))) return("run on raw counts")
  NULL
}

#' Expected outcome of every table of a run.
e2e_expectations <- function(analysis, data, o, res) {
  model <- e2e_ref_model(analysis, data, o)
  nets <- e2e_nets(model)
  multi <- inherits(model, "netobject_group")
  refusal <- e2e_refusal(analysis, o)
  matrix_titles <- function(nms) vapply(nms, \(n) if (grepl("^m[0-9]+$", n))
    res$buildModel_matrixTable$getColumn(n)$title else n, "")
  out <- list(
    buildModel_matrixTable = c(e2e_outcome(o$buildModel_show_matrix, \() e2e_ref_matrix(nets),
                                           "buildModelTitle"), list(rename = matrix_titles)),
    buildModel_summaryTable = e2e_outcome(o$buildModel_show_summary, \() e2e_ref_summary(nets), "buildModelTitle"),
    centralityTable = e2e_outcome(o$centrality_show_table, \() e2e_ref_centrality(nets, o), "centralityTitle"),
    edgeBetweennessTable = e2e_outcome(o$edgebetweenness_show_table, \() e2e_ref_edge_betweenness(nets, o),
                                       "edgeBetweennessTitle"),
    communityTable = e2e_outcome(o$community_show_table, \() e2e_ref_community(nets, o), "communityTitle"),
    cliquesTable = e2e_outcome(o$cliques_show_table, \() e2e_ref_cliques(nets, o), "cliquesTitle"),
    cliquesComplexTable = e2e_outcome(o$cliques_show_complex, \() e2e_ref_complex(nets, o), "cliquesTitle"),
    backboneTable = e2e_outcome(o$backbone_show_table, \() e2e_ref_backbone(nets, o), "backboneTitle"),
    linksTable = e2e_outcome(o$links_show_table, \() e2e_ref_links(nets, o), "linksTitle"),
    rulesTable = e2e_outcome(o$rules_show_table, \() e2e_ref_rules(nets, o), "rulesTitle"))

  resample <- function(flag, compute, heading, table) {
    if (!isTRUE(o[[flag]])) return(list(kind = "hidden"))
    if (!is.null(refusal)) return(list(kind = "refused", pattern = refusal))
    e2e_outcome(TRUE, compute, heading)
  }
  out$bootstrapTable <- resample("bootstrap_show_table", \() e2e_ref_bootstrap(model, nets, o),
                                 "bootstrapTitle")
  out$permutationTable <- if (!multi) list(kind = "hidden") else
    resample("permutation_show_table", \() e2e_ref_permutation(model, o), "permutationTitle")
  if (isTRUE(o$stability_show_table) && is.null(refusal)) {
    st <- e2e_ref_stability(model, o)
    out$stabilityTable <- if (is.null(st$stability)) list(kind = "empty") else list(kind = "rows", df = st$stability)
    out$reliabilityTable <- if (is.null(st$reliability)) list(kind = "empty") else list(kind = "rows", df = st$reliability)
    out$.stability <- st
  } else {
    out$stabilityTable <- out$reliabilityTable <-
      if (isTRUE(o$stability_show_table)) list(kind = "refused", pattern = refusal) else list(kind = "hidden")
  }
  if (multi && any(vapply(c("compare_show_summary", "compare_show_network", "compare_show_edges"),
                          \(f) isTRUE(o[[f]]), logical(1)))) {
    cmp <- e2e_attempt(e2e_ref_compare(analysis, model, o))
    if (inherits(cmp, "error")) out$.compare_error <- list(kind = "error", error = cmp, heading = "compareTitle")
    pick <- function(flag, part) {
      if (!isTRUE(o[[flag]])) return(list(kind = "hidden"))
      if (inherits(cmp, "error")) return(list(kind = "error", error = cmp, heading = "compareTitle"))
      list(kind = "rows", df = cmp[[part]])
    }
    out$compareSummaryTable <- pick("compare_show_summary", "summary")
    out$compareNetworkTable <- pick("compare_show_network", "network")
    out$.compare <- cmp
  } else {
    out$compareSummaryTable <- out$compareNetworkTable <- list(kind = "hidden")
  }
  if (analysis %in% c("CNAEventCluster", "CNABinaryCluster")) {
    out$clusterDiagTable <- e2e_outcome(o$cluster_show_diagnostics,
                                        \() e2e_ref_cluster_diag(attr(model, "e2e_clustering")),
                                        "clusterTitle")
    out$clusterChoiceTable <- e2e_outcome(o$cluster_show_choice,
                                          \() e2e_ref_cluster_choice(analysis, data, o), "clusterTitle")
  }
  list(tables = out, model = model, multi = multi)
}

#' Check every table of a run: reference equality plus the table contracts.
e2e_check <- function(analysis, data, o, label) {
  res <- e2e_run(analysis, data, o)
  ex <- e2e_known(e2e_expectations(analysis, data, o, res))
  tables <- ex$tables
  errors <- Filter(\(e) identical(e$kind, "error"), tables)
  expected_fail <- unique(vapply(errors, `[[`, "", "heading"))
  failures <- section_failures(res)
  expect_setequal(names(failures), expected_fail)
  expect_false(isTRUE(res$errorText$visible), label = paste(label, "no input error"))
  present <- intersect(names(tables), res$itemNames)
  invisible(lapply(present, \(nm) e2e_expect(res, analysis, nm, tables[[nm]], label)))
  e2e_expect_compare_edges(res, o, tables$.compare, label)
  e2e_expect_stability_notes(res, tables$.stability, label)
  e2e_contracts(res, analysis, o, ex, label)
  if (identical(analysis, "CNABinaryCluster") && is.null(e2e_scaling(o))) {
    # the row-to-cluster mapping of the reference: each cluster's network is
    # the one Nestimate builds from the clustering itself
    direct <- e2e_nets(Nestimate::build_network(attr(ex$model, "e2e_clustering"),
                                                method = "co_occurrence"))
    # (Nestimate keeps only the states a cluster uses; the rest must be 0 here)
    invisible(Map(\(net, g) {
      w <- unclass(net$weights)
      d <- unclass(direct[[g]]$weights)
      used <- rownames(d)
      expect_equal(w[used, used], d[used, used], label = paste(label, g, "== build_network(<net_clustering>)"))
      expect_true(all(w[setdiff(rownames(w), used), ] == 0), label = paste(label, g, "unused states are 0"))
    }, e2e_nets(ex$model), names(e2e_nets(ex$model))))
  }
  invisible(list(res = res, ex = ex))
}

# Edge differences: the top rows by |difference| (ties broken by name, so
# compared as values per unordered pair, and as the multiset of |difference|).
e2e_expect_compare_edges <- function(res, o, cmp, label) {
  if (!"compareEdgeTable" %in% res$itemNames) return(invisible())
  if (is.null(cmp) || !isTRUE(o$compare_show_edges)) {
    expect_false(isTRUE(res$compareEdgeTable$visible) && res$compareEdgeTable$rowCount > 0,
                 label = paste(label, "compareEdgeTable hidden"))
    return(invisible())
  }
  if (inherits(cmp, "error")) return(invisible())
  got <- e2e_table(res, "compareEdgeTable")
  want <- cmp$edges
  expect_equal(nrow(got), min(o$compare_edges_max_rows, nrow(want)),
               label = paste(label, "compareEdgeTable rows (max_rows)"))
  key <- paste(pmin(got$from, got$to), pmax(got$from, got$to), sep = "|")
  ref <- want[match(key, want$key), c("weight_i", "weight_j", "difference", "relative_difference")]
  rownames(ref) <- NULL
  expect_equal(got[names(ref)], e2e_norm(ref), tolerance = 1e-8,
               label = paste(label, "compareEdgeTable values"))
  expect_equal(sort(abs(got$difference), decreasing = TRUE),
               utils::head(sort(abs(want$difference), decreasing = TRUE), nrow(got)),
               label = paste(label, "compareEdgeTable keeps the largest differences"))
  expect_equal(res$compareTitle$content, cmp$title, label = paste(label, "compare title"))
}

e2e_expect_stability_notes <- function(res, st, label) {
  if (is.null(st)) return(invisible())
  keys_cs <- vapply(res$stabilityTable$notes, \(n) n$key, "")
  keys_nr <- vapply(res$reliabilityTable$notes, \(n) n$key, "")
  expect_equal(c(cs = "cs" %in% keys_cs, cd = "cd" %in% keys_cs, nr = "nr" %in% keys_nr),
               st$failed, label = paste(label, "stability failure notes"))
  expect_equal("constant" %in% keys_cs, length(st$missing) > 0 && !st$failed[["cs"]],
               label = paste(label, "constant-measure note"))
}

#' Contracts every table must keep regardless of the reference: the group
#' column, the row caps, the documented sort orders.
e2e_contracts <- function(res, analysis, o, ex, label) {
  groups <- names(e2e_nets(ex$model))
  # the rows as filled (the state keeps the group column even when it is hidden)
  tbl <- function(nm) {
    st <- if (nm %in% res$itemNames) res[[nm]]$state
    if (is.data.frame(st) && isTRUE(res[[nm]]$visible)) st else data.frame()
  }
  lab <- function(...) paste(label, ...)
  with_group <- c("buildModel_matrixTable", "buildModel_summaryTable", "centralityTable",
                  "edgeBetweennessTable", "communityTable", "cliquesTable", "cliquesComplexTable",
                  "backboneTable", "linksTable", "rulesTable", "bootstrapTable", "stabilityTable",
                  "reliabilityTable")
  invisible(lapply(intersect(with_group, res$itemNames), \(nm) {
    state <- res[[nm]]$state
    if (!is.data.frame(state) || !nrow(state)) return()
    if (ex$multi) {
      expect_true(all(state$group %in% groups), label = lab(nm, "group column names the networks"))
      expect_true(isTRUE(res[[nm]]$getColumn("group")$visible), label = lab(nm, "group column shown"))
    } else {
      expect_true(all(state$group == ""), label = lab(nm, "group is '' for one network"))
      expect_false(isTRUE(res[[nm]]$getColumn("group")$visible), label = lab(nm, "group column hidden"))
    }
  }))
  sorted_within <- function(d, by) all(vapply(split(by, d$group), \(v) !is.unsorted(v), logical(1)))
  bs <- tbl("bootstrapTable")
  if (nrow(bs)) {
    if (!isTRUE(o$bootstrap_table_show_all))
      expect_lte(nrow(bs), o$bootstrap_table_max_rows, label = lab("bootstrap max_rows"))
    if (isTRUE(o$bootstrap_table_significant_only))
      expect_true(all(bs$significant == "Yes"), label = lab("bootstrap significant_only"))
    expect_false(is.unsorted(bs$significant != "Yes"), label = lab("bootstrap: significant first"))
    expect_true(all(vapply(split(bs$p_value, bs$significant), \(p) !is.unsorted(p), logical(1))),
                label = lab("bootstrap: then by p-value"))
  }
  pm <- tbl("permutationTable")
  if (nrow(pm)) {
    if (!isTRUE(o$permutation_table_show_all))
      expect_lte(nrow(pm), o$permutation_table_max_rows, label = lab("permutation max_rows"))
    expect_false(is.unsorted(pm$p_value), label = lab("permutation: smallest p first"))
  }
  ru <- tbl("rulesTable")
  if (nrow(ru)) {
    expect_true(all(table(ru$group) <= o$rules_max_rows), label = lab("rules max_rows per network"))
    expect_true(sorted_within(ru, -ru$lift), label = lab("rules: highest lift first"))
  }
  lk <- tbl("linksTable")
  if (nrow(lk)) {
    expect_true(all(table(lk$group) <= o$links_top_n), label = lab("links top_n per network"))
    expect_true(sorted_within(lk, lk$rank), label = lab("links: by rank"))
    if (isTRUE(o$links_exclude_existing)) expect_true(all(lk$existing == "no"), label = lab("links: only absent"))
  }
  eb <- tbl("edgeBetweennessTable")
  if (nrow(eb)) expect_true(sorted_within(eb, -eb$betweenness), label = lab("edge betweenness: highest first"))
  bb <- tbl("backboneTable")
  if (nrow(bb)) expect_true(sorted_within(bb, bb$status == "removed"), label = lab("backbone: kept first"))
  cq <- tbl("cliquesTable")
  if (nrow(cq)) {
    expect_true(sorted_within(cq, -cq$mean_weight), label = lab("cliques: strongest first"))
    expect_true(all(lengths(strsplit(cq$members, ", ", fixed = TRUE)) == o$cliques_size),
                label = lab("cliques: of the chosen size"))
  }
  ce <- tbl("compareEdgeTable")
  if (nrow(ce)) expect_false(is.unsorted(-abs(ce$difference)), label = lab("compare edges: largest first"))
}

# ---------------------------------------------------------------------------
# Helpers for the hand-made runs below
# ---------------------------------------------------------------------------

e2e_read <- function(name)
  utils::read.csv(file.path(repo_root(), "data", name), stringsAsFactors = FALSE)

E2E_SRL_CODES <- c("Plan", "Monitor", "Evaluate", "Reflect", "Read", "Summarize",
                   "Practice", "Review")

# Every visible table of a run, as filled.
e2e_snapshot <- function(res) {
  tables <- Filter(\(nm) inherits(res[[nm]], "Table") && isTRUE(res[[nm]]$visible), res$itemNames)
  stats::setNames(lapply(tables, \(nm) e2e_table(res, nm)), tables)
}

e2e_event_opts <- function(analysis = "CNAEventSingle", ...)
  e2e_opts(analysis, buildModel_variables_long_action = "Action",
           buildModel_variables_actor = "Actor", buildModel_variables_long_time = "Time",
           buildModel_threshold = 600L, ...)

# ---------------------------------------------------------------------------
# 1. Random option combinations, every analysis
# ---------------------------------------------------------------------------

E2E_DRAWS <- 8L

invisible(lapply(ANALYSES, \(analysis) test_that(
  sprintf("%s: %d random option combinations == direct Nestimate/cograph", analysis, E2E_DRAWS), {
    skip_on_cran()
    invisible(lapply(seq_len(E2E_DRAWS), \(i) {
      d <- e2e_draw(analysis, i)
      label <- sprintf("%s #%d", analysis, i)
      expect_true(e2e_in_range(analysis, d$opts), label = paste(label, "options inside the .a.yaml limits"))
      e2e_check(analysis, d$data, d$opts, label)
    }))
  })))

# ---------------------------------------------------------------------------
# 2. Table contracts at their limits
# ---------------------------------------------------------------------------

test_that("row caps, show-all and significant-only hold at their limits", {
  dat <- e2e_events()
  grouped <- \(...) e2e_event_opts(buildModel_variables_group = "Group", compare_group_i = "G3",
                                   compare_group_j = "G1", bootstrap_method = "stability", ...)
  n_edges <- 5 * 6 / 2                       # 5 actions, unique pairs incl. the diagonal
  capped <- e2e_check("CNAEventSingle", dat, grouped(
    bootstrap_table_max_rows = 1L, permutation_table_max_rows = 1L, rules_max_rows = 1L,
    compare_edges_max_rows = 1L, links_top_n = 1L, links_exclude_existing = FALSE,
    rules_min_support = 0.05), "limits capped")$res
  expect_equal(capped$bootstrapTable$rowCount, 1)
  expect_equal(capped$permutationTable$rowCount, 1)
  expect_equal(capped$compareEdgeTable$rowCount, 1)
  expect_equal(as.vector(table(table_df(capped, "linksTable")$group)), c(1L, 1L, 1L))
  expect_equal(as.vector(table(table_df(capped, "rulesTable")$group)), c(1L, 1L, 1L))

  everything <- e2e_check("CNAEventSingle", dat, grouped(
    bootstrap_table_show_all = TRUE, permutation_table_show_all = TRUE,
    bootstrap_table_max_rows = 1L, permutation_table_max_rows = 1L), "limits show_all")$res
  expect_equal(everything$bootstrapTable$rowCount, 3 * n_edges)      # 3 groups
  expect_equal(everything$permutationTable$rowCount, 3 * n_edges)    # 3 comparisons
  expect_setequal(unique(table_df(everything, "permutationTable")$group_comparison),
                  c("G1 vs G2", "G1 vs G3", "G2 vs G3"))

  significant <- e2e_check("CNAEventSingle", dat, grouped(
    bootstrap_table_show_all = TRUE, bootstrap_table_significant_only = TRUE), "limits significant")$res
  sig <- table_df(significant, "bootstrapTable")
  expect_true(all(sig$significant == "Yes"))
  expect_equal(nrow(sig), sum(table_df(everything, "bootstrapTable")$significant == "Yes"))
})

# ---------------------------------------------------------------------------
# 3. Bundled datasets through the analyses they fit
# ---------------------------------------------------------------------------

test_that("bundled event datasets: every table == direct Nestimate/cograph", {
  skip_on_cran()
  srl <- e2e_read("SRL_events.csv")
  srl_opts <- \(analysis, ...) e2e_opts(analysis, buildModel_variables_long_action = "Action",
                                        buildModel_variables_actor = "Actor",
                                        buildModel_variables_long_time = "Time",
                                        buildModel_variables_long_order = "Order",
                                        buildModel_threshold = 900L, ...)
  e2e_check("CNAEventSingle", srl, srl_opts("CNAEventSingle", buildModel_variables_group = "Group",
                                            compare_group_i = "Strategic", compare_group_j = "Surface"),
            "SRL_events (groups)")
  e2e_check("CNAEventCluster", srl, srl_opts("CNAEventCluster", clustering_k = 2L), "SRL_events (clusters)")
  reg <- regulation_data(n_actors = 80, seed = 3)
  e2e_check("CNAEventSingle", reg, e2e_opts("CNAEventSingle", buildModel_variables_long_action = "Action",
                                            buildModel_variables_actor = "Actor",
                                            buildModel_variables_long_time = "Time_Posix",
                                            buildModel_variables_group = "Group", compare_group_i = "B",
                                            compare_group_j = "C", buildModel_threshold = 900L),
            "Regulation_long (ISO 8601 times)")
})

test_that("bundled one-hot and delimited datasets: every table == direct Nestimate/cograph", {
  skip_on_cran()
  srl <- e2e_read("SRL_onehot.csv")
  e2e_check("CNABinarySingle", srl, e2e_opts("CNABinarySingle", buildModel_variables_onehot = E2E_SRL_CODES,
                                             buildModel_variables_actor = "Actor",
                                             buildModel_variables_group = "Group", buildModel_window = 2L,
                                             compare_group_i = "Strategic", compare_group_j = "Surface"),
            "SRL_onehot (groups)")
  e2e_check("CNABinaryCluster", srl, e2e_opts("CNABinaryCluster", buildModel_variables_onehot = E2E_SRL_CODES,
                                              buildModel_variables_actor = "Actor", buildModel_window = 2L,
                                              buildModel_window_type = "sliding"),
            "SRL_onehot (clusters, sliding)")
  la <- e2e_read("LearningActivities.csv")
  e2e_check("CNABinaryCluster", la, e2e_opts("CNABinaryCluster",
                                             buildModel_variables_onehot = c("Reading", "Video", "Forum",
                                                                             "Quiz", "Coding", "Review"),
                                             buildModel_variables_actor = "Student", clustering_k = 3L),
            "LearningActivities (clusters, numeric Student)")
  hai <- e2e_read("HumanAI_turns.csv")
  e2e_check("CNADelimited", hai, e2e_opts("CNADelimited", buildModel_variables_field = "Codes",
                                          buildModel_similarity = "association"),
            "HumanAI_turns (delimited codes)")
  e2e_check("CNADelimited", e2e_read("Regulation_sessions.csv"),
            e2e_opts("CNADelimited", buildModel_variables_field = "Actions", buildModel_min_occur = 5L),
            "Regulation_sessions (delimited actions)")
})

# ---------------------------------------------------------------------------
# 4. Storage types
# ---------------------------------------------------------------------------

test_that("events: Actor/Group as numbers, factors or text and Time as POSIXct, ISO text or seconds give the same tables", {
  reg <- regulation_data(n_actors = 30, seed = 4)
  o <- e2e_opts("CNAEventSingle", buildModel_variables_long_action = "Action",
                buildModel_variables_actor = "Actor", buildModel_variables_long_time = "Time_Posix",
                buildModel_variables_group = "Group", compare_group_i = "A", compare_group_j = "B",
                buildModel_threshold = 900L)
  base <- e2e_snapshot(e2e_run("CNAEventSingle", reg, o))        # integer Actor, ISO text time
  expect_gt(length(base), 10)
  variants <- list(
    actor_text = transform(reg, Actor = as.character(Actor)),
    actor_factor = transform(reg, Actor = factor(Actor)),
    group_factor = transform(reg, Group = factor(Group)),
    time_posixct = transform(reg, Time_Posix = as.POSIXct(Time_Posix, format = "%Y-%m-%dT%H:%M:%OSZ",
                                                          tz = "UTC")),
    time_seconds = transform(reg, Time_Posix = Time_Unix))
  invisible(Map(\(d, nm) expect_equal(e2e_snapshot(e2e_run("CNAEventSingle", d, o)), base, label = nm),
                variants, names(variants)))
})

test_that("one-hot: Actor/Group as factors and 0/1 codes as nominal labels give the same tables", {
  srl <- e2e_read("SRL_onehot.csv")
  o <- e2e_opts("CNABinarySingle", buildModel_variables_onehot = E2E_SRL_CODES,
                buildModel_variables_actor = "Actor", buildModel_variables_group = "Group",
                buildModel_window = 2L, compare_group_i = "Strategic", compare_group_j = "Surface")
  base <- e2e_snapshot(e2e_run("CNABinarySingle", srl, o))
  as_labels <- srl
  as_labels[E2E_SRL_CODES] <- lapply(srl[E2E_SRL_CODES], \(x) factor(x, levels = c(0, 1)))
  variants <- list(actor_factor = transform(srl, Actor = factor(Actor)),
                   group_factor = transform(srl, Group = factor(Group)), codes_labels = as_labels)
  invisible(Map(\(d, nm) expect_equal(e2e_snapshot(e2e_run("CNABinarySingle", d, o)), base, label = nm),
                variants, names(variants)))
})

test_that("one-hot Actor stored as integers or as the same labels in text gives the same tables", {
  # Nestimate orders units by their sorted values (integers 1, 2, ..., 10, text
  # "1", "10", "2"), so the Actor must reach it as text, as for events, or the
  # same seed draws different resamples (regression)
  srl <- e2e_read("SRL_onehot.csv")
  srl$Actor <- match(srl$Actor, unique(srl$Actor))
  o <- e2e_opts("CNABinarySingle", buildModel_variables_onehot = E2E_SRL_CODES,
                buildModel_variables_actor = "Actor", buildModel_window = 2L)
  as_integer <- e2e_snapshot(e2e_run("CNABinarySingle", srl, o))
  as_text <- e2e_snapshot(e2e_run("CNABinarySingle", transform(srl, Actor = as.character(Actor)), o))
  invisible(lapply(names(as_text), \(nm) expect_equal(as_integer[[nm]], as_text[[nm]], label = nm)))
})

# ---------------------------------------------------------------------------
# 5. Invariant: row order
# ---------------------------------------------------------------------------

# The rule: row order is never part of the definition for events with a Time
# (the order within an actor is the time order) or for delimited data (a row is
# a transaction), so any shuffle keeps every table, resampling included. For
# one-hot data with windows wider than one row, the order of rows within an
# actor x session IS the definition (it decides the windows), so only shuffles
# that interleave the units while keeping each unit's own row order apply.
# Groups are kept out here; their order is tested separately below.
test_that("shuffling rows leaves every table unchanged where row order is not the definition", {
  ev <- e2e_events()
  ev_opts <- e2e_event_opts()
  shuffled <- withr::with_seed(21L, ev[sample.int(nrow(ev)), ])
  expect_equal(e2e_snapshot(e2e_run("CNAEventSingle", shuffled, ev_opts)),
               e2e_snapshot(e2e_run("CNAEventSingle", ev, ev_opts)))

  oh <- e2e_onehot()
  oh_opts <- e2e_opts("CNABinarySingle", buildModel_variables_onehot = E2E_CODES,
                      buildModel_variables_actor = "Actor", buildModel_variables_session = "Session",
                      buildModel_window = 2L)
  unit <- paste(oh$Actor, oh$Session)
  position <- stats::ave(seq_along(unit), unit, FUN = seq_along)
  draw <- withr::with_seed(22L, stats::runif(length(unique(unit))))[match(unit, unique(unit))]
  interleaved <- oh[order(position, draw), ]
  expect_false(identical(interleaved$Actor, oh$Actor))
  expect_equal(e2e_snapshot(e2e_run("CNABinarySingle", interleaved, oh_opts)),
               e2e_snapshot(e2e_run("CNABinarySingle", oh, oh_opts)))

  dl <- e2e_delimited()
  dl_opts <- e2e_opts("CNADelimited", buildModel_variables_field = c("keywords", "tags"))
  expect_equal(e2e_snapshot(e2e_run("CNADelimited", withr::with_seed(23L, dl[sample.int(nrow(dl)), ]), dl_opts)),
               e2e_snapshot(e2e_run("CNADelimited", dl, dl_opts)))
})

test_that("group order and comparison direction do not depend on row order", {
  # Nestimate orders groups by first appearance; the module sorts the rows by
  # group first, so G2's rows coming first must not reorder the grouped tables
  # nor turn "G1 vs G2" into "G2 vs G1" with every difference negated (regression)
  dat <- e2e_events(groups = c("G1", "G2"), n_actor = 10L)
  o <- e2e_event_opts(buildModel_variables_group = "Group", compare_group_i = "G1",
                      compare_group_j = "G2", bootstrap_method = "stability")
  g2_first <- dat[order(dat$Group != "G2", seq_len(nrow(dat))), ]
  reordered <- e2e_run("CNAEventSingle", g2_first, o)
  expect_equal(e2e_snapshot(reordered), e2e_snapshot(e2e_run("CNAEventSingle", dat, o)))
  expect_equal(unique(table_df(reordered, "permutationTable")$group_comparison), "G1 vs G2")
  # a factor's levels set the order
  by_levels <- e2e_run("CNAEventSingle", transform(dat, Group = factor(Group, levels = c("G2", "G1"))), o)
  expect_equal(unique(table_df(by_levels, "centralityTable")$group), c("G2", "G1"))
})

# ---------------------------------------------------------------------------
# 6. Error path
# ---------------------------------------------------------------------------

test_that("a column assigned both as a code and as the Actor is a classed error the runner explains", {
  oh <- e2e_onehot()
  o <- e2e_opts("CNABinarySingle", buildModel_variables_onehot = c(E2E_CODES, "Actor"),
                buildModel_variables_actor = "Actor", tables = FALSE)
  err <- tryCatch(cna("cna_build_model")(oh, o, "binary_single"), cna_bad_input = \(e) e)
  expect_s3_class(err, "cna_bad_input")
  expect_error(cna("cna_build_model")(oh, o, "binary_single"), class = "cna_bad_input")
  res <- e2e_run("CNABinarySingle", oh, c(o[setdiff(names(o), "centrality_show_table")],
                                          list(centrality_show_table = TRUE)))
  expect_true(res$errorText$visible)
  expect_identical(res$errorText$content, paste("Could not build the network:", conditionMessage(err)))
  expect_false(res$centralityTable$visible)
})

# ---------------------------------------------------------------------------
# 7. Regressions for the defects the random runs found
# ---------------------------------------------------------------------------

test_that("threshold inference fills the bootstrap table; its consistency-range columns stay empty", {
  # Nestimate's threshold-inference summary has no cr_lower / cr_upper, and
  # the table assembly fails on them ("arguments imply differing number of rows")
  dat <- e2e_events()
  o <- e2e_event_opts(tables = FALSE, bootstrap_show_table = TRUE, bootstrap_method = "threshold",
                      bootstrap_table_show_all = TRUE)
  res <- e2e_run("CNAEventSingle", dat, o)
  expect_identical(res$bootstrapTitle$content, "")
  net <- e2e_ref_model("CNAEventSingle", dat, o)
  ref <- Nestimate::bootstrap_network(net, iter = 20L, ci_level = 0.05, inference = "threshold",
                                      seed = 11L)$summary
  got <- table_df(res, "bootstrapTable")
  expect_equal(nrow(got), nrow(ref))
  if (nrow(got)) {
    expect_true(all(is.na(got$cr_lower)) && all(is.na(got$cr_upper)))
    expect_equal(got$p_value[match(paste(ref$from, ref$to), paste(got$from, got$to))], ref$p_value)
  }
})

test_that("the reliability table of one network has group '' like every other table", {
  res <- e2e_run("CNAEventSingle", e2e_events(), e2e_event_opts(tables = FALSE, stability_show_table = TRUE))
  expect_gt(res$reliabilityTable$rowCount, 0)
  expect_equal(unique(res$stabilityTable$state$group), "")
  expect_equal(unique(res$reliabilityTable$state$group), "")      # Nestimate names it "co_occurrence"
})

test_that("event clusters lacking some actions share every state and can be compared and permuted", {
  # two clean clusters, A-B-C and C-D-E: built as groups, both networks carry
  # all five states (zero rows for the unused ones), so compare_model() and
  # permutation() accept them (regression: built from the clustering, each had
  # only its own states and both sections failed)
  dat <- do.call(rbind, lapply(seq_len(12L), \(k) data.frame(
    Actor = sprintf("s%02d", k),
    Action = rep(if (k <= 6L) c("A", "B", "C") else c("C", "D", "E"), length.out = 9L),
    stringsAsFactors = FALSE)))
  o <- e2e_opts("CNAEventCluster", buildModel_variables_long_action = "Action",
                buildModel_variables_actor = "Actor", tables = FALSE, compare_show_summary = TRUE,
                compare_show_edges = TRUE, permutation_show_table = TRUE, buildModel_show_matrix = TRUE)
  out <- e2e_check("CNAEventCluster", dat, o, "clusters with different actions")
  used <- lapply(e2e_nets(out$ex$model), \(net) {
    w <- unclass(net$weights)
    sort(rownames(w)[rowSums(w) > 0])
  })
  expect_equal(unname(used), list(c("A", "B", "C"), c("C", "D", "E")))
  invisible(lapply(e2e_nets(out$ex$model), \(net) expect_setequal(rownames(net$weights), LETTERS[1:5])))
  expect_identical(out$res$compareTitle$content, "Comparing: Cluster 1 vs Cluster 2")
  expect_gt(out$res$compareSummaryTable$rowCount, 0)
  expect_gt(out$res$permutationTable$rowCount, 0)
})

test_that("delimited with a Transaction ID pools every item column; rows without an ID are left out and counted", {
  long <- e2e_delimited_long(e2e_delimited(), two_fields = TRUE)
  o <- e2e_opts("CNADelimited", buildModel_variables_field = c("item", "item2"),
                buildModel_variables_by = "doc")
  out <- e2e_check("CNADelimited", long, o, "delimited by two item columns")
  expect_match(out$res$buildModelTitle$content, "2 rows without a Transaction ID were left out.",
               fixed = TRUE)
  # items of the second column, trimmed, are states; blanks and NAs are not
  expect_true(all(c("theory", "data") %in% rownames(out$ex$model$weights)))
  expect_false(any(c("", " theory", "data ") %in% rownames(out$ex$model$weights)))
})

test_that("the stability table keeps the networks' order when a network has no stable measure (fixed)", {
  # Degree is constant in the complete networks of G1 and G3 (no CS row for
  # them) but not in G2, whose actors use only A and B. Fixed 2026-11-06:
  # cna_stability_df() ordered the groups by first appearance among the bound
  # rows (CS rows first), so the table read G2, G1, G3 while every other table
  # reads G1, G2, G3.
  dat <- e2e_events()
  dat$Action[dat$Group == "G2"] <- rep(c("A", "B"), length.out = sum(dat$Group == "G2"))
  o <- e2e_event_opts(tables = FALSE, stability_show_table = TRUE, buildModel_variables_group = "Group",
                      centrality_Strength = FALSE, centrality_Betweenness = FALSE,
                      centrality_Degree = TRUE)
  res <- e2e_run("CNAEventSingle", dat, o)
  expect_equal(unique(res$stabilityTable$state$group), c("G1", "G2", "G3"))
})

test_that("a one-node network has no links to predict: explained in a note, not a raw error", {
  # Nestimate::predict_links() fails on a 1 x 1 matrix; such networks arise from
  # small clusters and from data that use a single code
  dat <- data.frame(Actor = rep(c("a", "b", "c"), each = 4), Action = "A", stringsAsFactors = FALSE)
  res <- e2e_run("CNAEventSingle", dat, e2e_opts("CNAEventSingle", buildModel_variables_long_action = "Action",
                                                 buildModel_variables_actor = "Actor", tables = FALSE,
                                                 links_show_table = TRUE))
  expect_identical(res$linksTitle$content, "")
  expect_equal(res$linksTable$rowCount, 0)
  expect_match(paste(e2e_notes(res, "linksTable"), collapse = " "), "single node has no links")
})
