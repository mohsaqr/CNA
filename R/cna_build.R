# Network builders: jamovi data + options -> Nestimate network(s).
#
# Every builder returns a `netobject` (one network) or a `netobject_group`
# (one network per group or cluster). Groups use Nestimate's native `group=`
# argument and clusters use build_clusters() + build_network(), so downstream
# code never loops over groups by hand: a netobject_group flows through every
# analysis and cograph::splot() draws it as a grid.

# Analysis kinds and what they are built from.
# `mixed`: whether the transition + co-occurrence view is offered.
.cna_kinds <- list(
  event_single   = list(data = "event",     cluster = FALSE, mixed = TRUE),
  event_cluster  = list(data = "event",     cluster = TRUE,  mixed = TRUE),
  binary_single  = list(data = "onehot",    cluster = FALSE, mixed = TRUE),
  binary_cluster = list(data = "onehot",    cluster = TRUE,  mixed = TRUE),
  delimited      = list(data = "delimited", cluster = FALSE, mixed = FALSE))

#' Scaling option -> Nestimate `scaling` (NULL = raw counts).
#' @noRd
.cna_scaling <- function(options) {
  scaling <- options$buildModel_scaling
  if (is.null(scaling) || identical(scaling, "noScaling")) NULL else scaling
}

#' The reproducibility seed, or NULL when "Set seed" is off.
#' @noRd
.cna_seed <- function(options) {
  if (isTRUE(options$seed_set)) as.integer(options$seed_value) else NULL
}

#' Evaluate `expr` with a fixed seed, restoring the caller's RNG state after.
#' Used for engine functions that take no `seed` argument.
#' @noRd
.cna_with_seed <- function(seed, expr) {
  if (is.null(seed)) return(expr)
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  on.exit(
    if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else rm(".Random.seed", envir = globalenv()),
    add = TRUE)
  set.seed(seed)
  expr
}

#' Event (long-format) columns as Nestimate expects them.
#'
#' Action, actor and group become character. Time is passed through untouched
#' (factors become character) so Nestimate::prepare() parses it: converting with
#' as.POSIXct() here would drop the time of day from ISO 8601 stamps such as
#' "2025-01-01T08:27:07Z" and merge each day into one session. A numeric Order
#' column is kept numeric so events sort 1, 2, ..., 10 rather than "1", "10", "2".
#' Rows with a missing Group value are dropped and counted (attribute
#' `n_missing_group`).
#' @noRd
.cna_event_frame <- function(data, options, with_group = TRUE) {
  vars <- list(action = options$buildModel_variables_long_action,
               actor  = options$buildModel_variables_actor,
               time   = options$buildModel_variables_long_time,
               order  = options$buildModel_variables_long_order,
               group  = if (with_group) options$buildModel_variables_group)
  cols <- unlist(vars[c("actor", "time", "action", "order", "group")], use.names = FALSE)
  df <- as.data.frame(data)[cols]

  as_text <- function(x) if (is.factor(x)) as.character(x) else x
  df[[vars$action]] <- as.character(df[[vars$action]])
  if (!is.null(vars$actor)) df[[vars$actor]] <- as.character(df[[vars$actor]])
  if (!is.null(vars$time))  df[[vars$time]]  <- as_text(df[[vars$time]])
  if (!is.null(vars$order)) df[[vars$order]] <- .cna_numeric_if_possible(df[[vars$order]])
  df <- .cna_drop_missing_group(df, vars$group)
  attr(df, "cna_vars") <- vars
  df
}

#' A column whose every value reads as a number becomes numeric; anything else
#' is returned unchanged.
#' @noRd
.cna_numeric_if_possible <- function(x) {
  if (is.numeric(x)) return(x)
  txt <- trimws(as.character(x))
  number <- "^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$"
  if (all(is.na(txt) | grepl(number, txt))) as.numeric(txt) else x
}

#' Drop rows whose Group value is missing (counted in `n_missing_group`), and
#' put the rows in group order. Nestimate orders groups by first appearance,
#' so otherwise the group order, and the direction of every comparison
#' ("G1 vs G2"), would follow the row order. Groups follow the factor levels
#' (jamovi's level order), otherwise numeric or alphabetical order. The sort is
#' stable: each group keeps its own row order (one-hot windows depend on it).
#' @noRd
.cna_drop_missing_group <- function(df, group) {
  if (is.null(group)) return(df)
  raw <- df[[group]]
  g <- as.character(raw)
  missing <- is.na(g) | !nzchar(trimws(g))
  values <- unique(g[!missing])
  levels <- if (is.factor(raw)) intersect(levels(raw), values)
            else values[order(.cna_numeric_if_possible(values), method = "radix")]
  keep <- which(!missing)
  keep <- keep[order(match(g[keep], levels), method = "radix")]
  out <- df[keep, , drop = FALSE]
  out[[group]] <- g[keep]
  attr(out, "n_missing_group") <- sum(missing)
  out
}

#' Event co-occurrence network: single, or one per group via `group=`.
#' @return A `netobject` or `netobject_group`.
#' @noRd
cna_build_event_nest <- function(data, options) {
  df <- .cna_event_frame(data, options)
  vars <- attr(df, "cna_vars")
  args <- Filter(Negate(is.null), list(
    data = df, method = "co_occurrence",
    action = vars$action, actor = vars$actor, time = vars$time,
    order = vars$order, time_threshold = options$buildModel_threshold,
    scaling = .cna_scaling(options), group = vars$group))
  net <- do.call(Nestimate::build_network, args)
  attr(net, "n_missing_group") <- attr(df, "n_missing_group")
  net
}

#' One 0/1 indicator column from a jamovi column: numbers 0/1, logicals, or
#' labels "0"/"1"/"TRUE"/"FALSE" (jamovi often stores 0/1 as a nominal factor).
#' Missing values stay NA (Nestimate counts them as inactive). Anything else is
#' an error naming the column.
#' @noRd
.cna_as_binary <- function(x, name) {
  if (is.logical(x)) return(as.integer(x))
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) {
    key <- toupper(trimws(x))
    map <- c("0" = 0L, "1" = 1L, "FALSE" = 0L, "TRUE" = 1L)
    bad <- unique(x[!is.na(x) & !key %in% names(map)])
    if (length(bad))
      .cna_stop(sprintf("One-hot column '%s' must hold 0/1 values; it contains %s.",
                        name, paste(sQuote(utils::head(bad, 3), FALSE), collapse = ", ")),
                "cna_not_binary")
    return(unname(map[key]))
  }
  bad <- unique(x[!is.na(x) & !x %in% c(0, 1)])
  if (length(bad))
    .cna_stop(sprintf("One-hot column '%s' must hold 0/1 values; it contains %s.",
                      name, paste(utils::head(bad, 3), collapse = ", ")),
              "cna_not_binary")
  as.integer(x)
}

#' The one-hot data as Nestimate should receive it, plus notes for the user.
#'
#' Validates the code columns (0/1; jamovi labels and logicals accepted) and
#' drops rows whose Actor, Session or Group is missing, saying how many
#' (Nestimate >= 1.0.7 refuses such rows). Only the assigned columns are passed,
#' each with its role declared, so no column takes a role from its name. Rows
#' keep their order: one-hot windows and transitions follow the row order within
#' each actor/session.
#' @return list(data, codes, actor, session, group, notes).
#' @noRd
.cna_onehot_frame <- function(data, options, with_group = TRUE) {
  codes <- options$buildModel_variables_onehot
  roles <- list(actor = options$buildModel_variables_actor,
                session = options$buildModel_variables_session,
                group = if (with_group) options$buildModel_variables_group)
  overlap <- intersect(unlist(roles), codes)
  if (length(overlap))
    .cna_stop(sprintf("'%s' is assigned both as a one-hot column and as Actor, Session or Group.",
                      overlap[[1L]]), "cna_bad_input")

  df <- as.data.frame(data)[c(codes, unlist(roles))]
  df[codes] <- Map(.cna_as_binary, df[codes], codes)
  n_na_cells <- sum(is.na(df[codes]))

  units <- unlist(roles[c("actor", "session")])
  # as for events: integer IDs and the same IDs as text are the same units,
  # in the same order (Nestimate orders units by their sorted values)
  df[units] <- lapply(df[units], as.character)
  missing_unit <- if (length(units)) Reduce(`|`, lapply(df[units], is.na)) else logical(nrow(df))
  df <- df[!missing_unit, , drop = FALSE]
  df <- .cna_drop_missing_group(df, roles$group)

  notes <- c(
    if (!is.null(roles$session) && is.null(roles$actor))
      paste("Rows with the same Session value form one sequence across actors; if",
            "session numbers restart for each actor, assign the Actor as well."),
    if (n_na_cells > 0L) sprintf("%d missing one-hot value%s counted as 0 (inactive).",
                                 n_na_cells, if (n_na_cells == 1L) " was" else "s were"),
    if (any(missing_unit)) sprintf("%d row%s with a missing Actor or Session %s left out.",
                                   sum(missing_unit), if (sum(missing_unit) == 1L) "" else "s",
                                   if (sum(missing_unit) == 1L) "was" else "were"),
    .cna_sliding_note(df, units, options))
  list(data = df, codes = codes, actor = roles$actor, session = roles$session,
       group = roles$group, notes = notes, n_missing_group = attr(df, "n_missing_group"))
}

#' In sliding (overlapping) mode Nestimate counts only complete windows (as
#' tna's import_onehot() does), so an actor/session with fewer rows than the
#' window contributes nothing at all (tumbling mode keeps the last, partial
#' window). Say how many.
#' @noRd
.cna_sliding_note <- function(df, units, options) {
  window <- .cna_window(options)
  if (window <= 1L || !identical(.cna_window_mode(options), "overlapping")) return(NULL)
  rows <- if (length(units)) as.vector(table(interaction(df[units], drop = TRUE))) else nrow(df)
  short <- sum(rows < window)
  if (!short) return(NULL)
  sprintf(paste("%d of %d %s fewer rows than the sliding window (%d) and contribute",
                "nothing; tumbling windows or a smaller window keep them."),
          short, length(rows), if (length(units)) "actors/sessions have" else "data set has",
          window)
}

#' Arguments for Nestimate::build_network() on one-hot data (any method).
#' @noRd
.cna_onehot_args <- function(frame, options) {
  Filter(Negate(is.null), list(
    data = frame$data, codes = frame$codes, format = "onehot",
    actor = frame$actor, session = frame$session,
    window_size = .cna_window(options), mode = .cna_window_mode(options),
    scaling = .cna_scaling(options), group = frame$group))
}

#' One-hot (binary) co-occurrence network: single or per group.
#' @noRd
cna_build_onehot_nest <- function(data, options) {
  frame <- .cna_onehot_frame(data, options)
  args <- .cna_onehot_args(frame, options)
  net <- do.call(Nestimate::build_network, c(args, method = "co_occurrence"))
  attr(net, "n_missing_group") <- frame$n_missing_group
  attr(net, "cna_notes") <- frame$notes
  attr(net, "cna_onehot_args") <- args   # the mixed view and rules reuse them
  net
}

.cna_window <- function(options) {
  window <- options$buildModel_window
  if (is.null(window) || window < 1) 1L else as.integer(window)
}

.cna_window_mode <- function(options) {
  if (identical(options$buildModel_window_type, "sliding")) "overlapping"
  else "non-overlapping"
}

#' Event sequences (wide), as clustering and cluster_choice() take them.
#' @noRd
cna_event_sequences <- function(data, options) {
  df <- .cna_event_frame(data, options, with_group = FALSE)
  vars <- attr(df, "cna_vars")
  args <- Filter(Negate(is.null), list(
    data = df, actor = vars$actor, action = vars$action, time = vars$time,
    order = vars$order, time_threshold = options$buildModel_threshold))
  do.call(Nestimate::prepare, args)$sequence_data
}

#' The one-hot network that clustering and Choosing k take: unscaled and
#' ungrouped. build_clusters() and cluster_choice() cluster one windowed
#' sequence per actor x session of it, with the analysis's windows (tumbling or
#' sliding), and return the assignments in first-seen actor x session order.
#' @noRd
cna_onehot_units <- function(data, options,
                             frame = .cna_onehot_frame(data, options, with_group = FALSE)) {
  args <- .cna_onehot_args(frame, options)
  args$scaling <- NULL
  do.call(Nestimate::build_network, c(args, method = "co_occurrence"))
}

#' Cluster sequences and build one co-occurrence network per cluster.
#'
#' The networks are built as groups (one group per cluster) from the clustered
#' sequences, so that, as with Group, every cluster's network has every action
#' (zero rows for the actions a cluster never uses) and scaling applies.
#' build_network(<net_clustering>) instead gives each cluster only its own
#' actions, which compare_model() and permutation() refuse, and ignores
#' `scaling` (Nestimate 1.0.7). On the actions a cluster uses, the weights are
#' the same.
#' @return A `netobject_group` keyed "Cluster 1", "Cluster 2", ... carrying the
#'   `clustering` attribute that cluster_diagnostics() reads.
#' @noRd
cna_cluster_networks <- function(seqs, options) {
  clustering <- .cna_clusters(seqs, options)
  seqs$.cna_cluster <- paste("Cluster", clustering$assignments)
  args <- Filter(Negate(is.null), list(
    data = seqs, method = "co_occurrence", group = ".cna_cluster",
    scaling = .cna_scaling(options)))
  nets <- .cna_order_clusters(do.call(Nestimate::build_network, args))
  attr(nets, "clustering") <- clustering
  nets
}

cna_build_event_cluster_nest <- function(data, options) {
  if (!isTRUE(options$clustering_run)) return(NULL)
  cna_cluster_networks(cna_event_sequences(data, options), options)
}

#' Cluster with the analysis's distance, method, k and seed.
#' @param x Event sequences (wide), or the one-hot network of cna_onehot_units().
#' @noRd
.cna_clusters <- function(x, options) {
  Nestimate::build_clusters(
    x, k = options$clustering_k, dissimilarity = options$clustering_dissimilarity,
    method = options$clustering_method, seed = .cna_seed(options))
}

#' Cluster the actor x session units of a one-hot network, then build each
#' cluster's network from its own rows like a group.
#'
#' Unscaled, the result equals build_network(<net_clustering>) (verified on
#' Nestimate 1.0.7, a1754ab); building the clusters as groups also applies the
#' scaling, which build_network(<net_clustering>) ignores, and keeps the rows
#' each network rests on for the mixed view and resampling. Every row gets its
#' unit's cluster through the first-seen actor x session order of the
#' assignments.
#' @noRd
cna_build_onehot_cluster_nest <- function(data, options) {
  if (!isTRUE(options$clustering_run)) return(NULL)
  frame <- .cna_onehot_frame(data, options, with_group = FALSE)
  clustering <- .cna_clusters(cna_onehot_units(data, options, frame), options)
  units <- c(frame$actor, frame$session)
  key <- do.call(paste, c(unname(as.list(frame$data[units])), sep = "\r"))
  unit_order <- unique(key)
  # Nestimate clusters the units that have a complete window. With sliding
  # windows a unit with fewer rows than the window has none (the model note
  # counts them): it gets no cluster, and its rows add nothing to any network.
  if (identical(.cna_window_mode(options), "overlapping")) {
    n_rows <- as.vector(table(key)[unit_order])
    unit_order <- unit_order[n_rows >= .cna_window(options)]
    frame$data <- frame$data[key %in% unit_order, , drop = FALSE]
    key <- key[key %in% unit_order]
  }
  stopifnot("one cluster assignment per actor/session" =
              length(unit_order) == length(clustering$assignments))
  frame$data$.cna_cluster <- paste("Cluster", clustering$assignments[match(key, unit_order)])
  frame$group <- ".cna_cluster"
  args <- .cna_onehot_args(frame, options)
  nets <- .cna_order_clusters(do.call(Nestimate::build_network,
                                      c(args, method = "co_occurrence")))
  attr(nets, "clustering") <- clustering
  attr(nets, "cna_onehot_args") <- args
  attr(nets, "cna_notes") <- frame$notes
  nets
}

#' Order a cluster netobject_group as Cluster 1, 2, ..., 10 (not 1, 10, 2).
#' @noRd
.cna_order_clusters <- function(nets) {
  if (!all(grepl("^Cluster [0-9]+$", names(nets)))) return(nets)
  ord <- order(as.integer(sub("^Cluster ", "", names(nets))))
  keep <- attributes(nets)
  out <- unclass(nets)[ord]
  attributes(out) <- c(list(names = names(nets)[ord]), keep[setdiff(names(keep), "names")])
  out
}

#' Delimited / bipartite co-occurrence via Nestimate::cooccurrence().
#'
#' Item columns holding several items separated by `sep` (keywords, codes), or
#' one item per row grouped into transactions by `by`. Similarity normalises the
#' raw counts (van Eck & Waltman, 2009).
#' @return A `netobject` without source data (no resampling analyses).
#' @noRd
cna_build_delimited_nest <- function(data, options) {
  fields <- options$buildModel_variables_field
  by <- options$buildModel_variables_by
  # cooccurrence() reads field + sep as delimited cells and then ignores `by`,
  # so with a Transaction ID the separator is not passed (one item per row).
  sep <- if (is.null(by)) options$buildModel_sep
  if (is.null(sep) || !nzchar(trimws(sep))) sep <- NULL
  df <- as.data.frame(data)[c(fields, by)]
  df[fields] <- lapply(df[fields], as.character)
  notes <- NULL
  if (!is.null(by)) {
    df <- .cna_transaction_items(df, fields, by)
    notes <- attr(df, "cna_notes")
    fields <- fields[[1L]]
  }
  top_n <- options$buildModel_top_n
  net <- Nestimate::cooccurrence(
    df, field = fields, by = by, sep = sep,
    similarity = options$buildModel_similarity,
    threshold = options$buildModel_edge_threshold,
    min_occur = as.integer(options$buildModel_min_occur),
    diagonal = isTRUE(options$buildModel_diagonal),
    top_n = if (is.null(top_n) || top_n <= 0) NULL else as.integer(top_n))
  attr(net, "cna_notes") <- notes
  net
}

#' Transaction-ID input: one item per cell, pooled over the Item columns into
#' one column (named like the first), trimmed, with blank and missing items
#' and rows without a Transaction ID left out; delimited mode treats its items
#' the same way. Nestimate 1.0.7's cooccurrence(by =) takes the values
#' verbatim (" a" and "a" are two items) and fails on several Item columns.
#' @return Data frame with columns `by` and `fields[[1]]`; a "cna_notes"
#'   attribute counts the rows without a Transaction ID.
#' @noRd
.cna_transaction_items <- function(df, fields, by) {
  items <- trimws(unlist(lapply(df[fields], as.character), use.names = FALSE))
  long <- data.frame(rep(df[[by]], length(fields)), items, stringsAsFactors = FALSE)
  names(long) <- c(by, fields[[1L]])
  no_id <- sum(is.na(df[[by]]))
  long <- long[!is.na(long[[by]]) & !is.na(items) & nzchar(items), , drop = FALSE]
  rownames(long) <- NULL
  attr(long, "cna_notes") <- if (no_id)
    sprintf("%d row%s without a Transaction ID %s left out.", no_id,
            if (no_id == 1L) "" else "s", if (no_id == 1L) "was" else "were")
  long
}

#' Transition + co-occurrence pairs for the mixed-network view.
#'
#' The co-occurrence side is the analysis's own network. The transition side is
#' a frequency network built from the same data with the same arguments, so
#' both share units, windows, scaling and groups:
#'   events  - build_network(<each network's sequences>, method = "frequency");
#'   one-hot - build_network(<same one-hot arguments>, method = "frequency"),
#'             i.e. Nestimate's wtna transitions within each actor x session.
#' @return Named list, one element per network ("" for a single network), of
#'   list(transition, cooccurrence) netobjects.
#' @noRd
cna_build_mixed <- function(data, options, kind, model) {
  nets <- .cna_nets(model)
  if (identical(.cna_kinds[[kind]]$data, "event")) {
    scaling <- .cna_scaling(options)
    return(lapply(nets, function(net) list(
      transition = Nestimate::build_network(net$data, method = "frequency",
                                            scaling = scaling),
      cooccurrence = net)))
  }
  transitions <- .cna_nets(.cna_order_clusters(do.call(
    Nestimate::build_network, c(attr(model, "cna_onehot_args"), method = "frequency"))))
  stopifnot("transition and co-occurrence networks pair up" =
              setequal(names(transitions), names(nets)))
  Map(function(tr, co) list(transition = tr, cooccurrence = co),
      transitions[match(names(nets), names(transitions))], nets)
}

#' Build the analysis's network(s) from data and options.
#' @noRd
cna_build_model <- function(data, options, kind) {
  switch(kind,
    event_single   = cna_build_event_nest(data, options),
    event_cluster  = cna_build_event_cluster_nest(data, options),
    binary_single  = cna_build_onehot_nest(data, options),
    binary_cluster = cna_build_onehot_cluster_nest(data, options),
    delimited      = cna_build_delimited_nest(data, options),
    stop(errorCondition(sprintf("unknown analysis kind '%s'", kind),
                        class = "cna_bad_input", call = NULL)))
}
