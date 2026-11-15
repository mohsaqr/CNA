# Analyses as pure functions: network(s) + parameters -> results.
#
# `*_compute()` functions call Nestimate / cograph and return the engine's own
# objects (what the plots draw). `*_df()` functions turn those into a tidy
# data.frame, one row per observation, with a `group` column ("" for a single
# network, the group or cluster name otherwise). Nothing here touches jamovi,
# so every function can be tested against a direct engine call.

# ---------------------------------------------------------------------------
# Shared helpers
# ---------------------------------------------------------------------------

#' The networks of a model as a plain named list ("" for a single network).
#' @noRd
.cna_nets <- function(model) {
  if (inherits(model, "netobject_group")) {
    nets <- unclass(model)
    attributes(nets) <- list(names = names(model))
    return(nets)
  }
  stats::setNames(list(model), "")
}

#' The parts of a grouped engine result as a plain named list.
#' @noRd
.cna_parts <- function(x, multi) {
  if (!multi) return(stats::setNames(list(x), ""))
  parts <- unclass(x)
  attributes(parts) <- list(names = names(x))
  parts
}

#' Row-bind per-network data.frames, dropping empty ones.
#' @noRd
.cna_bind <- function(parts) {
  parts <- Filter(function(p) !is.null(p) && nrow(p) > 0L, unname(parts))
  if (!length(parts)) return(NULL)
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' Apply `fun(net, group)` to each network and bind the results with a
#' leading `group` column.
#' @noRd
.cna_per_net <- function(model, fun) {
  nets <- .cna_nets(model)
  .cna_bind(Map(function(net, g) {
    df <- fun(net, g)
    if (is.null(df) || !nrow(df)) return(NULL)
    data.frame(group = rep(g, nrow(df)), df, check.names = FALSE,
               stringsAsFactors = FALSE)
  }, nets, names(nets)))
}

#' Unique undirected edges of a weight matrix: from, to, weight (upper
#' triangle, off-diagonal, non-zero).
#' @noRd
.cna_edges <- function(w) {
  idx <- which(upper.tri(w) & w != 0, arr.ind = TRUE)
  data.frame(from = rownames(w)[idx[, 1L]], to = colnames(w)[idx[, 2L]],
             weight = as.numeric(w[idx]), stringsAsFactors = FALSE)
}

#' Muffle two informational notes the engines print on every call (they flood
#' the console inside resampling loops): cograph's "Weights inverted ..." and
#' Nestimate's "centralities computed excluding loops ...". Every other
#' message passes through.
#' @noRd
.cna_quiet_notes <- function(expr) {
  notes <- c("Weights inverted", "centralities computed excluding loops")
  withCallingHandlers(expr, message = function(m) {
    if (any(vapply(notes, grepl, logical(1), x = conditionMessage(m), fixed = TRUE)))
      invokeRestart("muffleMessage")
  })
}

#' Evaluate `expr`, collecting its warnings instead of printing them.
#' @return list(value, warnings).
#' @noRd
.cna_collect_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = unique(warnings))
}

.cna_stop <- function(message, class = "cna_bad_input")
  stop(errorCondition(message, class = class, call = NULL))

# ---------------------------------------------------------------------------
# Network description
# ---------------------------------------------------------------------------

#' Weight matrix as a data.frame: group, state, then one column per state
#' (the union of states across groups; NA where a group lacks a state).
#' @noRd
cna_matrix_df <- function(model) {
  nets <- .cna_nets(model)
  states <- unique(unlist(lapply(nets, function(n) rownames(n$weights))))
  .cna_per_net(model, function(net, g) {
    w <- net$weights
    full <- matrix(NA_real_, nrow(w), length(states),
                   dimnames = list(rownames(w), states))
    full[, colnames(w)] <- unclass(w)
    data.frame(state = rownames(w), as.data.frame(full, optional = TRUE),
               check.names = FALSE, stringsAsFactors = FALSE)
  })
}

#' Network-level metrics (Nestimate's summary()): group, metric, value.
#' @noRd
cna_summary_df <- function(model) {
  .cna_per_net(model, function(net, g) {
    s <- summary(net)
    data.frame(metric = as.character(s$metric), value = as.numeric(s$value),
               stringsAsFactors = FALSE)
  })
}

# ---------------------------------------------------------------------------
# Centrality
# ---------------------------------------------------------------------------

# Centrality option suffix -> cograph::centrality() measure.
.cna_measure_map <- c(
  Strength = "strength", Degree = "degree", Betweenness = "betweenness",
  Closeness = "closeness", Harmonic = "harmonic", Eigenvector = "eigenvector",
  PageRank = "pagerank", Coreness = "coreness", Clustering = "transitivity",
  Diffusion = "diffusion", BetweennessRSP = "rsp_betweenness",
  Eccentricity = "eccentricity", Constraint = "constraint",
  EffectiveSize = "effective_size")

.cna_measure_titles <- c(
  Strength = "Strength", Degree = "Degree", Betweenness = "Betweenness",
  Closeness = "Closeness", Harmonic = "Harmonic closeness",
  Eigenvector = "Eigenvector", PageRank = "PageRank", Coreness = "Coreness",
  Clustering = "Clustering", Diffusion = "Diffusion",
  BetweennessRSP = "RSP betweenness", Eccentricity = "Eccentricity",
  Constraint = "Constraint", EffectiveSize = "Effective size")

#' Measures ticked in the options, in display order.
#' @noRd
cna_selected_measures <- function(options) {
  measures <- names(.cna_measure_map)
  measures[vapply(measures, function(m) isTRUE(options[[paste0("centrality_", m)]]),
                  logical(1))]
}

#' The cograph column that holds `measure` (undirected measures carry "_all").
#' @noRd
.cna_cograph_column <- function(ct, measure) {
  hit <- intersect(c(measure, paste0(measure, "_all")), names(ct))
  if (!length(hit))
    .cna_stop(sprintf("cograph returned no column for measure '%s'", measure),
              "cna_engine_contract")
  hit[[1L]]
}

#' Centralities of one network with cograph::centrality().
#'
#' Diffusion is tna's diffusion centrality (row sums of W + W^2 + ... + W^n),
#' range-normalised to 0-1 as Nestimate does by default
#' (`net_centrality(normalize_diffusion = TRUE)`): on raw counts the power
#' series runs into the billions. cograph's "diffusion" is a different,
#' binary measure (Kandhway & Kuri), constant on dense co-occurrence networks.
#' A matrix input (a stability resample) never asks for it: stability uses
#' Nestimate's built-in Diffusion, the same definition.
#'
#' Clustering is the weighted coefficient of Zhang & Horvath (2005), the one
#' tna uses, on the weights divided by their maximum: the formula is defined
#' for weights in 0-1, and on raw counts it grows with the counts (tna's value
#' equals 2 x max(w) x this one, tna summing both directions of the matrix).
#' cograph's default "transitivity" is unweighted and is 1 for every node of a
#' complete network, which most co-occurrence networks are.
#'
#' Two cograph measures are computed apart from the others (cograph 2.7.12):
#' * cograph computes eccentricity on the raw weights even with
#'   `invert_weights = TRUE`, unlike its other path measures. With `invert`,
#'   Eccentricity is computed on the lengths 1/w with inversion off, so every
#'   path measure uses the same lengths.
#' * Burt's effective size sums over the alters of ego other than ego itself,
#'   so a self-loop cannot change it; cograph counts self-loops as redundant
#'   ties when `loops = TRUE`. EffectiveSize is computed without loops.
#' @param x A netobject, or a weight matrix (the stability resamples).
#' @return Named list: `node` (character) and one numeric vector per measure
#'   (named by the option suffix), in node order.
#' @references Burt, R. S. (1992). Structural holes. Harvard University Press.
#'
#'   Banerjee, A., Chandrasekhar, A. G., Duflo, E., & Jackson, M. O. (2013).
#'   The diffusion of microfinance. Science, 341, 1236498. (Diffusion)
#'
#'   Zhang, B., & Horvath, S. (2005). A general framework for weighted gene
#'   co-expression network analysis. Statistical Applications in Genetics and
#'   Molecular Biology, 4(1), 17. (Clustering)
#' @noRd
.cna_centralities <- function(x, measures, loops = FALSE, normalized = FALSE,
                              invert = TRUE) {
  run <- function(input, ms, loops, inv) {
    cg <- unname(.cna_measure_map[ms])
    args <- list(input, measures = cg, loops = loops, normalized = normalized,
                 invert_weights = inv)
    if (is.matrix(input)) args$directed <- FALSE
    ct <- .cna_quiet_notes(do.call(cograph::centrality, args))
    c(list(node = as.character(ct$node)),
      stats::setNames(lapply(cg, function(m) as.numeric(ct[[.cna_cograph_column(ct, m)]])), ms))
  }
  diffusion <- intersect(measures, "Diffusion")
  if (length(diffusion) && is.matrix(x))
    .cna_stop("Diffusion is computed on a network, not a matrix", "cna_engine_contract")
  if (is.matrix(x) && is.null(rownames(x)))   # node names to align the parts
    dimnames(x) <- rep(list(as.character(seq_len(nrow(x)))), 2L)
  clustering <- intersect(measures, "Clustering")
  by_length <- if (invert) intersect(measures, "Eccentricity") else character()
  loop_free <- if (loops) intersect(measures, "EffectiveSize") else character()
  rest <- setdiff(measures, c(diffusion, clustering, by_length, loop_free))
  w <- if (is.matrix(x)) x else unclass(x$weights)
  tna_diffusion <- function() {
    d <- .cna_quiet_notes(Nestimate::net_centrality(
      x, measures = "Diffusion", loops = loops, normalize_diffusion = TRUE))
    list(node = as.character(d$state), Diffusion = as.numeric(d$Diffusion))
  }
  parts <- c(
    if (length(rest)) list(run(x, rest, loops, invert)),
    if (length(diffusion)) list(tna_diffusion()),
    if (length(clustering)) list(list(node = rownames(w),
                                      Clustering = .cna_weighted_clustering(w))),
    if (length(loop_free)) list(run(x, loop_free, FALSE, invert)),
    if (length(by_length)) list(run(ifelse(w > 0, 1 / w, 0), by_length, loops, FALSE)))
  node <- parts[[1L]]$node
  values <- unlist(lapply(parts, function(p)
    lapply(p[setdiff(names(p), "node")], function(v) v[match(node, p$node)])),
    recursive = FALSE)
  c(list(node = node), values[measures])
}

#' Weighted clustering coefficient (Zhang & Horvath, 2005) on w / max(w),
#' diagonal excluded:
#'   C_i = sum_jk w_ij w_jk w_ki / ((sum_j w_ij)^2 - sum_j w_ij^2),
#' in 0-1. NaN where it is undefined (fewer than two neighbours, no edges).
#' @noRd
.cna_weighted_clustering <- function(w) {
  w <- unname(unclass(w))
  diag(w) <- 0
  top <- max(w)
  if (!isTRUE(top > 0)) return(rep(NaN, nrow(w)))
  w <- w / top
  triangles <- rowSums((w %*% w) * w)
  strength <- rowSums(w)
  triangles / (strength^2 - rowSums(w^2))
}

#' Centralities via cograph::centrality(): group, state, one column per
#' selected measure (named by the option suffix, e.g. "Strength").
#' @param invert Invert weights for path-based measures, so a strong
#'   co-occurrence is a short distance (Nestimate's default).
#' @noRd
cna_centrality_df <- function(model, measures, loops = FALSE,
                              normalized = FALSE, invert = TRUE) {
  if (!length(measures))
    .cna_stop("Select at least one centrality measure.", "cna_no_measures")
  .cna_per_net(model, function(net, g) {
    ct <- .cna_centralities(net, measures, loops = loops, normalized = normalized,
                            invert = invert)
    out <- data.frame(state = ct$node, stringsAsFactors = FALSE)
    out[measures] <- ct[measures]
    out
  })
}

#' Case-dropping stability of centralities (CS-coefficients), edge-weight
#' case-drop reliability and split-half reliability (Nestimate).
#'
#' Every measure is computed on each resample by the same definition as the
#' centrality table, so the stability is that of the values the table shows:
#' Diffusion by Nestimate's built-in (tna's, range-normalised, as the table),
#' every other measure by the table's own function (.cna_centralities).
#' Nestimate uses a built-in whenever a requested name matches one (its
#' Clustering is tna's, on raw counts), so those are passed under a "cg_"
#' prefix. Normalising does not change a correlation, so it is not applied.
#' @return list(cs, cd, nr, labels, requested, multi): the engine objects
#'   (`cs` with display titles as measure names), a map from those titles
#'   back to option names, the requested options, and whether the model has
#'   several networks.
#' @references Epskamp, S., Borsboom, D., & Fried, E. I. (2018). Estimating
#'   psychological networks and their accuracy. Behavior Research Methods, 50,
#'   195-212.
#' @noRd
cna_stability_compute <- function(model, measures, iter, threshold, certainty,
                                  method = "pearson", loops = FALSE,
                                  invert = TRUE, seed = NULL) {
  if (!length(measures))
    .cna_stop("Select at least one centrality measure.", "cna_no_measures")
  # Diffusion: Nestimate's built-in, the table's definition (tna, range-
  # normalised); every other measure: the table's own function under "cg_".
  own <- setdiff(measures, "Diffusion")
  nest_measures <- ifelse(measures == "Diffusion", "Diffusion", paste0("cg_", measures))
  constant <- lapply(.cna_nets(model), .cna_constant_measures, measures = measures,
                     loops = loops, invert = invert)
  centrality_fn <- if (length(own)) function(mat)
    stats::setNames(.cna_centralities(mat, own, loops = loops, invert = invert)[own],
                    paste0("cg_", own))
  # The three procedures are independent: one that the engine cannot run on
  # this network (see .cna_attempt) does not hide the results of the others.
  cs <- .cna_tidy_stability(.cna_attempt(Nestimate::centrality_stability(
    model, measures = nest_measures, iter = iter, threshold = threshold,
    certainty = certainty, method = method, centrality_fn = centrality_fn,
    loops = loops, normalize_diffusion = TRUE, seed = seed)), constant)
  cd <- .cna_attempt(Nestimate::casedrop_reliability(
    model, iter = iter, threshold = threshold, certainty = certainty,
    method = method, seed = seed))
  nr <- .cna_attempt(Nestimate::network_reliability(model, iter = iter, seed = seed))
  list(cs = cs, cd = cd, nr = nr,
       labels = stats::setNames(measures, unname(.cna_measure_titles[measures])),
       requested = measures, multi = inherits(model, "netobject_group"))
}

#' Measures (option names) constant across the nodes of `net`: a constant
#' cannot be correlated, so its stability is undefined.
#' @noRd
.cna_constant_measures <- function(net, measures, loops, invert) {
  ct <- .cna_centralities(net, measures, loops = loops, invert = invert)
  measures[vapply(ct[measures], function(v) {
    v <- v[is.finite(v)]
    length(v) < 2L || isTRUE(all.equal(min(v), max(v)))
  }, logical(1))]
}

#' Tidy a stability result: measures get their display titles (they are
#' computed under "cg_<option>" names, see cna_stability_compute()), and a
#' measure constant in a network is removed from that network's result.
#' Nestimate drops such a measure when others vary, but reports CS = 0 for
#' every measure when all are constant; either way there is no coefficient.
#' The table, the missing-measure note and Nestimate's plot read these fields.
#' @param constant Named list (by network) of constant option names.
#' @noRd
.cna_tidy_stability <- function(cs, constant) {
  if (.cna_failed(cs)) return(cs)
  if (inherits(cs, "net_stability_group")) {
    out <- Map(.cna_tidy_stability, unclass(cs), lapply(names(cs), function(g) constant[g]))
    attributes(out) <- attributes(cs)
    return(out)
  }
  title <- function(nm) unname(.cna_measure_titles[sub("^cg_", "", nm)])
  drop <- unname(.cna_measure_titles[unlist(constant, use.names = FALSE)])
  keep <- function(x) x[!names(x) %in% drop]
  names(cs$cs) <- title(names(cs$cs))
  names(cs$correlations) <- title(names(cs$correlations))
  cs$measures <- setdiff(title(cs$measures), drop)
  cs$cs <- keep(cs$cs)
  cs$correlations <- keep(cs$correlations)
  cs
}

#' The engine's result, or the error it raised (kept so it can be reported).
#' @noRd
.cna_attempt <- function(expr) tryCatch(expr, error = function(e) e)

.cna_failed <- function(x) inherits(x, "error")

#' Why a stability procedure is missing, phrased for the table note.
#' @param hint Optional sentence appended to each message.
#' @return Named character vector (cs, cd, nr) of the procedures that failed.
#' @noRd
cna_stability_failures <- function(st, hint = NULL) {
  what <- c(cs = "Centrality stability", cd = "Edge-weight case-drop reliability",
            nr = "Split-half reliability")
  failed <- Filter(function(k) .cna_failed(st[[k]]), names(what))
  stats::setNames(vapply(failed, function(k)
    trimws(paste(sprintf("%s could not be computed by Nestimate (%s).",
                         what[[k]], conditionMessage(st[[k]])), hint %||% "")),
    character(1)), failed)
}



#' Requested measures the stability analysis returned no coefficient for.
#' Nestimate leaves out a measure that is constant across nodes (for example
#' Degree in a complete network): a constant cannot be correlated, so its
#' stability is undefined. The table says so instead of dropping it silently.
#' @noRd
cna_stability_missing <- function(st, multi) {
  if (.cna_failed(st$cs)) return(character())
  returned <- lapply(.cna_parts(st$cs, multi), function(cs)
    .cna_measure_label(names(cs$cs), st$labels))
  setdiff(unname(.cna_measure_titles[st$requested]), Reduce(intersect, returned))
}

#' CS-coefficients: group, measure, cs (plus an "Edge weights" row per group).
#' @noRd
cna_stability_df <- function(st, multi) {
  cs_rows <- if (!.cna_failed(st$cs)) {
    parts <- .cna_parts(st$cs, multi)
    .cna_bind(Map(function(cs, g) {
      if (!length(cs$cs)) return(NULL)   # every measure constant in this network
      data.frame(group = g, measure = .cna_measure_label(names(cs$cs), st$labels),
                 cs = as.numeric(cs$cs), stringsAsFactors = FALSE)
    }, parts, names(parts)))
  }
  cd_rows <- if (!.cna_failed(st$cd)) {
    parts <- .cna_parts(st$cd, multi)
    .cna_bind(Map(function(cd, g)
      data.frame(group = g, measure = "Edge weights", cs = as.numeric(cd$cs),
                 stringsAsFactors = FALSE),
      parts, names(parts)))
  }
  rows <- .cna_bind(list(cs_rows, cd_rows))
  if (is.null(rows)) return(NULL)
  # networks in the model's order (a network whose measures are all constant
  # has only its edge-weight row), edge weights last within each network
  groups <- unique(unlist(lapply(list(st$cs, st$cd), function(x)
    if (!.cna_failed(x)) names(.cna_parts(x, multi)))))
  rows <- rows[order(match(rows$group, groups), rows$measure == "Edge weights"), ]
  rownames(rows) <- NULL
  rows
}

#' Display title of a stability measure (Nestimate name -> option -> title).
#' @noRd
.cna_measure_label <- function(nest_names, labels) {
  option <- unname(labels[nest_names])
  option[is.na(option)] <- nest_names[is.na(option)]
  title <- unname(.cna_measure_titles[option])
  title[is.na(title)] <- option[is.na(title)]
  title
}

#' Split-half reliability: group, metric, mean, sd ("" group for one network;
#' Nestimate labels a single network's rows with its method name).
#' @noRd
cna_reliability_df <- function(st) {
  if (.cna_failed(st$nr)) return(NULL)
  s <- st$nr$summary
  if (is.null(s) || !nrow(s)) return(NULL)
  data.frame(group = if (isTRUE(st$multi) && "model" %in% names(s)) as.character(s$model) else "",
             metric = as.character(s$metric), mean = as.numeric(s$mean),
             sd = as.numeric(s$sd), stringsAsFactors = FALSE)
}

# ---------------------------------------------------------------------------
# Communities
# ---------------------------------------------------------------------------

.cna_community_methods <- c("louvain", "leiden", "walktrap", "fast_greedy",
                            "infomap", "label_propagation", "edge_betweenness",
                            "leading_eigenvector", "spinglass", "optimal")

#' Community detection with cograph::communities(), one result per network.
#' @param gamma Resolution: `resolution` for Louvain and Leiden, `gamma` for
#'   Spinglass (cograph passes `resolution` only to the first two).
#' @return Named list of `cograph_communities`.
#' @noRd
cna_community_compute <- function(model, method, gamma = 1, seed = NULL) {
  if (!method %in% .cna_community_methods)
    .cna_stop(sprintf("Unknown community method '%s'.", method))
  lapply(.cna_nets(model), function(net) {
    if (identical(method, "optimal") && nrow(net$weights) > 50L)
      .cna_stop(paste("Optimal community detection grows exponentially with",
                      "the number of nodes; use it for at most 50 nodes."),
                "cna_too_large")
    if (identical(method, "spinglass"))
      cograph::communities(net, method = method, gamma = gamma, seed = seed)
    else
      cograph::communities(net, method = method, resolution = gamma, seed = seed)
  })
}

#' group, state, community.
#' @noRd
cna_community_df <- function(coms) {
  .cna_bind(Map(function(cm, g)
    data.frame(group = g, state = as.character(cm$node),
               community = as.integer(cm$community), stringsAsFactors = FALSE),
    coms, names(coms)))
}

# ---------------------------------------------------------------------------
# Edge betweenness
# ---------------------------------------------------------------------------

#' Edge betweenness (Nestimate::net_edge_betweenness), one per network.
#' @noRd
cna_edge_betweenness_compute <- function(model, invert = TRUE) {
  lapply(.cna_nets(model), function(net)
    Nestimate::net_edge_betweenness(net, invert = invert))
}

#' group, from, to, weight (co-occurrence), betweenness; highest first.
#' @noRd
cna_edge_betweenness_df <- function(ebs) {
  .cna_bind(Map(function(eb, g) {
    w <- eb$source_weights
    edges <- .cna_edges(w)
    if (!nrow(edges)) return(NULL)
    edges$betweenness <- as.numeric(eb$weights[cbind(edges$from, edges$to)])
    edges <- edges[order(-edges$betweenness, edges$from, edges$to), ]
    data.frame(group = g, edges, stringsAsFactors = FALSE)
  }, ebs, names(ebs)))
}

# ---------------------------------------------------------------------------
# Cliques and the clique complex
# ---------------------------------------------------------------------------

.cna_clique_limit <- 1e6

#' Upper bound on the number of cliques of up to `size` nodes in the graph of
#' edges with weight >= threshold: a j-clique contains j vertices whose
#' degree is at least j - 1, so there are at most sum_v C(deg_v, j - 1) / j.
#' @noRd
.cna_clique_bound <- function(w, threshold, size) {
  a <- abs(unclass(w))
  diag(a) <- 0
  a <- pmax(a, t(a))
  deg <- rowSums(a > 0 & a >= threshold)
  sum(vapply(seq_len(size), function(j) sum(choose(deg, j - 1L)) / j, numeric(1)))
}

#' Clique complexes (Nestimate::build_simplicial, edges with weight >=
#' threshold), one per network.
#'
#' `listing` holds every clique of up to `size` nodes (for the clique table and
#' plot). `full` is the clique complex for the summary, built up to the largest
#' clique size the enumeration bound allows: usually the whole complex, so its
#' Betti numbers are final; when it had to stop short, Nestimate's warning about
#' the truncation is reported with the summary.
#' @return Named list (one per network) of list(listing, full).
#' @noRd
cna_cliques_compute <- function(model, size, threshold) {
  lapply(.cna_nets(model), function(net) {
    bound <- .cna_clique_bound(net$weights, threshold, size)
    if (bound > .cna_clique_limit)
      .cna_stop(sprintf(paste(
        "This network can hold up to %s cliques of %d or fewer nodes at",
        "threshold %s; raise the threshold or lower the size."),
        format(round(bound), big.mark = ","), size, format(threshold)),
        "cna_too_many_cliques")
    n <- nrow(net$weights)
    feasible <- Filter(function(k) .cna_clique_bound(net$weights, threshold, k) <= .cna_clique_limit,
                       seq.int(size, max(size, n)))
    full_size <- max(feasible)
    list(listing = Nestimate::build_simplicial(net, type = "clique", threshold = threshold,
                                               max_dim = as.integer(size) - 1L),
         full = Nestimate::build_simplicial(net, type = "clique", threshold = threshold,
                                            max_dim = as.integer(full_size) - 1L))
  })
}

#' Cliques of exactly `size` nodes: group, clique, members, mean_weight,
#' min_weight; strongest (highest mean weight) first.
#' @noRd
cna_cliques_df <- function(complexes, model, size) {
  nets <- .cna_nets(model)
  .cna_bind(Map(function(cx, net, g) {
    sc <- cx$listing
    members <- Filter(function(s) length(s) == size, sc$simplices)
    if (!length(members)) return(NULL)
    w <- unclass(net$weights)
    labels <- sc$nodes
    edge_stats <- vapply(members, function(s) {
      sub <- w[labels[s], labels[s], drop = FALSE]
      e <- sub[upper.tri(sub)]
      c(mean = mean(e), min = min(e))
    }, numeric(2))
    df <- data.frame(group = g,
                     members = vapply(members, function(s) paste(labels[s], collapse = ", "),
                                      character(1)),
                     mean_weight = edge_stats["mean", ], min_weight = edge_stats["min", ],
                     stringsAsFactors = FALSE)
    df <- df[order(-df$mean_weight, df$members), ]
    df$clique <- seq_len(nrow(df))
    df[c("group", "clique", "members", "mean_weight", "min_weight")]
  }, complexes, nets, names(nets)))
}

#' Clique-complex summary: group, metric, value (clique counts by size, Euler
#' characteristic, Betti numbers, maximum q of the Q-analysis), from the
#' `full` complex of cna_cliques_compute().
#' @return list(rows, notes); notes carry Nestimate's warnings (e.g. a
#'   complex too large to build completely).
#' @references Atkin, R. H. (1974). Mathematical structure in human affairs.
#' @noRd
cna_complex_df <- function(complexes) {
  parts <- Map(function(cx, g) {
    sc <- cx$full
    betti <- .cna_collect_warnings(Nestimate::betti_numbers(sc))
    q <- Nestimate::q_analysis(sc)
    f <- as.integer(sc$f_vector)
    rows <- data.frame(
      group = g,
      metric = c(sprintf("Cliques of %d node%s", seq_along(f), ifelse(seq_along(f) == 1L, "", "s")),
                 "Euler characteristic",
                 sprintf("Betti number b%d", seq_along(betti$value) - 1L),
                 "Maximum q (Q-analysis)"),
      value = c(f, Nestimate::euler_characteristic(sc), as.numeric(betti$value),
                as.numeric(q$max_q)),
      stringsAsFactors = FALSE)
    list(rows = rows, notes = betti$warnings)
  }, complexes, names(complexes))
  list(rows = .cna_bind(lapply(parts, `[[`, "rows")),
       notes = unique(unlist(lapply(parts, `[[`, "notes"))))
}

# ---------------------------------------------------------------------------
# Backbone extraction
# ---------------------------------------------------------------------------

#' Backbone via Nestimate::net_prune(). The network is undirected, so an edge
#' is kept when either of its two directions survives (for the disparity
#' filter: when it is significant for at least one endpoint, Serrano, Boguna &
#' Vespignani, 2009). net_prune() never disconnects the network.
#' @return Named list of list(weights, keep, backbone).
#' @references Serrano, M. A., Boguna, M., & Vespignani, A. (2009). Extracting
#'   the multiscale backbone of complex weighted networks. PNAS, 106,
#'   6483-6488.
#' @noRd
cna_backbone_compute <- function(model, method = "disparity", level = 0.5,
                                 threshold = 0, lowest = 0.05) {
  lapply(.cna_nets(model), function(net) {
    pruned <- Nestimate::net_prune(net, method = method, threshold = threshold,
                                   lowest = lowest, level = level)
    w <- unclass(net$weights)
    p <- unclass(pruned$weights)
    keep <- ((p != 0) | (t(p) != 0)) & (w != 0)
    list(weights = w, keep = keep, backbone = w * keep)
  })
}

#' group, from, to, weight, status ("kept" / "removed"); kept edges first.
#' @noRd
cna_backbone_df <- function(backbones) {
  .cna_bind(Map(function(b, g) {
    edges <- .cna_edges(b$weights)
    if (!nrow(edges)) return(NULL)
    kept <- b$keep[cbind(edges$from, edges$to)]
    edges$status <- ifelse(kept, "kept", "removed")
    edges <- edges[order(edges$status, -edges$weight, edges$from, edges$to), ]
    data.frame(group = g, edges, stringsAsFactors = FALSE)
  }, backbones, names(backbones)))
}

# ---------------------------------------------------------------------------
# Link prediction
# ---------------------------------------------------------------------------

#' Link prediction (Nestimate::predict_links): group, from, to, score, rank,
#' existing. Weights at or below `threshold` count as absent edges. The
#' diagonal (a code co-occurring with itself) is not a link: predict_links
#' would count it in the degrees and neighbourhoods and shift every score.
#' @references Liben-Nowell, D., & Kleinberg, J. (2007). The link-prediction
#'   problem for social networks. JASIST, 58, 1019-1031.
#' @noRd
cna_links_df <- function(model, method = "resource_allocation",
                         exclude_existing = TRUE, threshold = 0, top_n = 20L) {
  .cna_per_net(model, function(net, g) {
    w <- unclass(net$weights)
    if (nrow(w) < 2L) return(NULL)   # one node: no pair to score
    w[abs(w) <= threshold] <- 0
    diag(w) <- 0
    # Jaccard compares neighbour sets; Nestimate 1.0.7's weighted Jaccard
    # mixes a weighted intersection with a count-based union (scores of 0
    # for pairs that share neighbours), so it is given the binary adjacency.
    pl <- Nestimate::predict_links(w, methods = method,
                                   weighted = !identical(method, "jaccard"),
                                   top_n = as.integer(top_n),
                                   exclude_existing = exclude_existing)
    p <- pl$predictions
    if (is.null(p) || !nrow(p)) return(NULL)
    data.frame(from = as.character(p$from), to = as.character(p$to),
               score = as.numeric(p$score), rank = as.integer(p$rank),
               existing = ifelse(p$existing, "yes", "no"),
               stringsAsFactors = FALSE)
  })
}

# ---------------------------------------------------------------------------
# Association rules
# ---------------------------------------------------------------------------

#' Association rules (Nestimate::association_rules, Apriori) per network.
#' @references Agrawal, R., & Srikant, R. (1994). Fast algorithms for mining
#'   association rules. VLDB.
#' @noRd
cna_rules_compute <- function(model, min_support = 0.1, min_confidence = 0.5,
                              min_lift = 1, max_length = 5L) {
  rules <- function(x) Nestimate::association_rules(
    x, min_support = min_support, min_confidence = min_confidence,
    min_lift = min_lift, max_length = as.integer(max_length))
  # On a one-hot network the codes are the items and each window of each
  # actor x session is a transaction (Nestimate >= 1.0.7; the tests check it
  # against an independent transaction count).
  lapply(.cna_nets(model), rules)
}

#' group, antecedent, consequent, support, confidence, lift, conviction,
#' count; highest lift first, at most `max_rows` per network.
#' @noRd
cna_rules_df <- function(rules, max_rows = 50L) {
  .cna_bind(Map(function(r, g) {
    df <- r$rules
    if (is.null(df) || !nrow(df)) return(NULL)
    df <- df[order(-df$lift, -df$support, df$antecedent, df$consequent), ]
    df <- utils::head(df, max_rows)
    data.frame(group = g, antecedent = as.character(df$antecedent),
               consequent = as.character(df$consequent),
               support = df$support, confidence = df$confidence, lift = df$lift,
               conviction = df$conviction, count = as.integer(df$count),
               stringsAsFactors = FALSE)
  }, rules, names(rules)))
}

# ---------------------------------------------------------------------------
# Bootstrap and permutation
# ---------------------------------------------------------------------------

#' Edge bootstrap (Nestimate::bootstrap_network). `threshold = NULL` lets
#' Nestimate choose the threshold for threshold inference.
#' @noRd
cna_bootstrap_compute <- function(model, iter, level, inference = "stability",
                                  range = c(0.75, 1.25), threshold = NULL,
                                  ci_method = "percentile", seed = NULL) {
  Nestimate::bootstrap_network(
    model, iter = iter, ci_level = level, inference = inference,
    consistency_range = range, edge_threshold = threshold,
    ci_method = ci_method, seed = seed)
}

#' group, from, to, weight, p_value, cr/ci bounds, significant. Significant
#' edges first, then by p-value; filtered and capped as the table options ask.
#' Threshold inference has no consistency range (Nestimate leaves the cr_*
#' columns out), so those cells are NA there.
#' @noRd
cna_bootstrap_df <- function(bs, multi, significant_only = FALSE,
                             show_all = FALSE, max_rows = 20L) {
  parts <- .cna_parts(bs, multi)
  summ <- .cna_bind(Map(function(b, g) {
    s <- b$summary
    if (is.null(s) || !nrow(s)) return(NULL)
    data.frame(group = g, from = as.character(s$from), to = as.character(s$to),
               weight = s$weight, p_value = s$p_value,
               cr_lower = s$cr_lower %||% NA_real_, cr_upper = s$cr_upper %||% NA_real_,
               ci_lower = s$ci_lower, ci_upper = s$ci_upper,
               significant = ifelse(s$sig %in% TRUE, "Yes", "No"),
               stringsAsFactors = FALSE)
  }, parts, names(parts)))
  if (is.null(summ)) return(NULL)
  summ <- summ[order(summ$significant != "Yes", summ$p_value, -abs(summ$weight),
                     summ$group, summ$from, summ$to), ]
  if (significant_only) summ <- summ[summ$significant == "Yes", , drop = FALSE]
  if (!show_all) summ <- utils::head(summ, max_rows)
  summ
}

#' Edge permutation test between every pair of groups (Nestimate::permutation).
#'
#' Nestimate adjusts p-values over all n x n cells of the difference matrix,
#' which for an undirected network counts every edge twice (Nestimate 1.0.7).
#' The test is run unadjusted and the correction applied here over the unique
#' edges each comparison reports; the adjusted values are written back into
#' `p_values` and `diff_sig` so the plot's stars agree with the table.
#' @noRd
cna_permutation_compute <- function(model, iter, alpha = 0.05, paired = FALSE,
                                    adjust = "none", seed = NULL) {
  pt <- Nestimate::permutation(model, iter = iter, alpha = alpha,
                               paired = paired, adjust = "none", seed = seed)
  cna_adjust_permutation(pt, adjust, alpha)
}

#' @noRd
cna_adjust_permutation <- function(pt, adjust, alpha) {
  if (identical(adjust, "none")) return(pt)
  adjusted <- lapply(pt, function(cmp) {
    s <- cmp$summary
    s$p_value <- stats::p.adjust(s$p_value, method = adjust)
    s$sig <- s$p_value < alpha
    p <- cmp$p_values
    p[cbind(s$from, s$to)] <- s$p_value
    if (!isTRUE(cmp$x$directed)) p[cbind(s$to, s$from)] <- s$p_value
    cmp$p_values <- p
    cmp$diff_sig <- cmp$diff * (p < alpha)
    cmp$summary <- s
    cmp$adjust <- adjust
    cmp
  })
  attributes(adjusted) <- attributes(pt)
  adjusted
}

#' group_comparison, edge_name, diff_true, effect_size, p_value; smallest
#' p-value first.
#' @noRd
cna_permutation_df <- function(pt, show_all = FALSE, max_rows = 50L) {
  rows <- .cna_bind(Map(function(cmp, nm) {
    s <- cmp$summary
    if (is.null(s) || !nrow(s)) return(NULL)
    data.frame(group_comparison = nm,
               edge_name = paste(s$from, "–", s$to),
               diff_true = s$diff, effect_size = s$effect_size,
               p_value = s$p_value, stringsAsFactors = FALSE)
  }, pt, names(pt)))
  if (is.null(rows)) return(NULL)
  rows <- rows[order(rows$p_value, -abs(rows$diff_true), rows$group_comparison,
                     rows$edge_name), ]
  if (!show_all) rows <- utils::head(rows, max_rows)
  rows
}

# ---------------------------------------------------------------------------
# Comparing two networks
# ---------------------------------------------------------------------------

#' Compare two networks of a group or cluster model (Nestimate::compare_model).
#' @noRd
cna_compare_compute <- function(model, i, j, scaling = "none") {
  Nestimate::compare_model(model, i = i, j = j, scaling = scaling)
}

#' category, metric, value.
#' @noRd
cna_compare_summary_df <- function(cmp) {
  s <- cmp$summary_metrics
  data.frame(category = as.character(s$category), metric = as.character(s$metric),
             value = as.numeric(s$value), stringsAsFactors = FALSE)
}

#' metric, group_i, group_j.
#' @noRd
cna_compare_network_df <- function(cmp) {
  s <- cmp$network_metrics
  data.frame(metric = as.character(s$metric), group_i = as.numeric(s$x),
             group_j = as.numeric(s$y), stringsAsFactors = FALSE)
}

#' Edge-level differences over the unique undirected pairs (diagonal
#' included): from, to, weight_i, weight_j, difference, relative_difference;
#' largest absolute difference first.
#' @noRd
cna_compare_edges_df <- function(cmp, max_rows = 50L) {
  e <- cmp$edge_metrics
  states <- unique(c(e$source, e$target))
  pos_from <- match(e$source, states)
  pos_to <- match(e$target, states)
  e <- e[pos_from <= pos_to, , drop = FALSE]
  df <- data.frame(from = as.character(e$source), to = as.character(e$target),
                   weight_i = e$weight_x, weight_j = e$weight_y,
                   difference = e$raw_difference,
                   relative_difference = e$relative_difference,
                   stringsAsFactors = FALSE)
  df <- df[order(-abs(df$difference), df$from, df$to), ]
  utils::head(df, max_rows)
}

# ---------------------------------------------------------------------------
# Clustering diagnostics
# ---------------------------------------------------------------------------

#' cluster, size, pct, mean_within_dist, silhouette.
#' @noRd
cna_cluster_diagnostics_df <- function(model) {
  d <- as.data.frame(Nestimate::cluster_diagnostics(attr(model, "clustering")))
  data.frame(cluster = paste("Cluster", d$cluster), size = as.integer(d$size),
             pct = as.numeric(d$pct), mean_within_dist = as.numeric(d$mean_within_dist),
             silhouette = as.numeric(d$sil_mean), stringsAsFactors = FALSE)
}

#' Fit k = kmin..kmax clusterings and compare them (Nestimate::cluster_choice).
#' @references Rousseeuw, P. J. (1987). Silhouettes: a graphical aid to the
#'   interpretation and validation of cluster analysis. J. Comput. Appl. Math.,
#'   20, 53-65.
#' @noRd
cna_cluster_choice_compute <- function(seqs, kmin, kmax, dissimilarity = "hamming",
                                       method = "pam", seed = NULL) {
  if (kmin > kmax) .cna_stop("The smallest k must not exceed the largest k.")
  .cna_with_seed(seed, Nestimate::cluster_choice(
    seqs, k = seq.int(kmin, kmax), dissimilarity = dissimilarity, method = method))
}

#' k, silhouette, mean_within_dist, min_size, max_size, size_ratio.
#' @noRd
cna_cluster_choice_df <- function(cc) {
  data.frame(k = as.integer(cc$k), silhouette = as.numeric(cc$silhouette),
             mean_within_dist = as.numeric(cc$mean_within_dist),
             min_size = as.integer(cc$min_size), max_size = as.integer(cc$max_size),
             size_ratio = as.numeric(cc$size_ratio), stringsAsFactors = FALSE)
}
