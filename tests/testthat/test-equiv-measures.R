# Network MEASURES on a built network, checked against igraph called directly
# or against an independent base-R implementation of the published formula.
# Nothing here compares CNA with cograph / Nestimate themselves: the reference
# side never calls either package (they only BUILD the networks under test).
#
# Definitions CNA reports (read from cograph 2.7.12 / Nestimate 1.0.7):
# * Centrality (cograph::centrality, undirected, weighted = TRUE):
#   - Strength / Degree: a self-loop counts twice (igraph convention).
#   - Betweenness, Closeness, Harmonic: shortest paths on 1/w with "invert",
#     on w without it. Closeness = 1 / sum of distances to the nodes reached
#     (igraph). Normalised: closeness x (nodes reached), harmonic / (n - 1).
#   - Every other measure is normalised by dividing by its maximum.
#   - Eigenvector: principal eigenvector, max = 1 (Bonacich 1972). PageRank:
#     damping 0.85, dangling nodes teleport uniformly (Brin & Page 1998).
#     With loops these two put w_ii (not 2 w_ii, as igraph does) on the
#     diagonal.
#   - Coreness and Effective size are unweighted.
#   - Clustering (2026-11-06): weighted, Zhang & Horvath (2005) on w / max(w)
#     without the diagonal, in 0-1 (NaN below two neighbours); computed by CNA
#     (cograph's default transitivity is unweighted, 1 on complete networks).
#   - Diffusion (2026-11-06): tna's, the row sums of W + W^2 + ... + W^n
#     (Banerjee et al. 2013), range-normalised to 0-1 (NaN when constant, as
#     tna's ranger()); it comes from Nestimate, not cograph.
#   - Constraint: Burt (1992) with p_ij = w_ij / s_i.
#   - RSP betweenness: Kivimaki et al. (2016), beta = 0.01, cost 1/w,
#     component-local; the invert switch does not reach it (cograph documents
#     this).
# * Communities: igraph's algorithm on the graph INCLUDING the diagonal as
#   self-loops, with the co-occurrence counts as weights. Leiden uses the CPM
#   objective (cograph's default), and spinglass uses the largest component.
# * Backbone (Nestimate::net_prune): the disparity filter counts the diagonal
#   in s_i and k_i (tna's convention), and an edge stays when it is
#   significant for either endpoint. Threshold and lowest never disconnect the
#   network. The cut-off for "lowest" is a quantile over every positive cell
#   (both triangles plus the diagonal), as Nestimate documents.
# * Link prediction on the off-diagonal matrix (weighted, except Jaccard,
#   which compares neighbour sets): CN = sum_z w_xz w_zy;
#   RA = sum_z w_xz w_zy / k_z; AA = sum_z w_xz w_zy / log k_z; PA = k_x k_y;
#   Katz = ((I - bA)^-1 - I) on the 0/1 adjacency, b = 0.9 / rho(A).
#
# Blocks marked DEFECT assert the correct value and fail until the upstream
# package is fixed (cograph or Nestimate; docs/EQUIVALENCE_AUDIT.md). Blocks
# marked "(fixed)" pin a defect CNA has fixed.

# ---------------------------------------------------------------------------
# Networks
# ---------------------------------------------------------------------------

# Evaluate `expr` under a seed without leaving the global RNG changed.
with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
            rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed)
  expr
}

# A symmetric matrix of non-negative integer counts as a co-occurrence netobject.
mat_net <- function(m) Nestimate::build_network(m, method = "co_occurrence")

# Same network with the diagonal removed (integer-count networks only).
strip_diag <- function(net) {
  m <- unclass(net$weights)
  diag(m) <- 0
  mat_net(m)
}

named_matrix <- function(nodes, edges, diagonal = NULL) {
  m <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
  m[cbind(edges$from, edges$to)] <- edges$w
  m[cbind(edges$to, edges$from)] <- edges$w
  if (!is.null(diagonal)) diag(m) <- diagonal
  m
}

# Random event networks: 4-8 codes, 4-10 actors, short to long windows (sparse
# with isolates and several components through to complete).
mc_event_net <- function(seed, scaling = "noScaling") with_seed(seed, {
  k <- sample(4:8, 1L)
  n_actor <- sample(4:10, 1L)
  per <- sample(5:12, 1L)
  thr <- sample(c(120, 400, 3000), 1L)
  cna("cna_build_event_nest")(
    make_event_data(n_actor = n_actor, actions = LETTERS[seq_len(k)], per = per, seed = seed),
    opts_event(buildModel_threshold = thr, buildModel_scaling = scaling))
})

reg_net <- function(n_actors, seed) cna("cna_build_event_nest")(
  regulation_data(n_actors = n_actors, seed = seed),
  opts_event(buildModel_variables_long_time = "Time_Posix", buildModel_threshold = 900))

hand_nets <- list(
  # triangle + pair + an isolate that has a self-loop; diagonal elsewhere too
  two_components = mat_net(named_matrix(
    letters[1:6], data.frame(from = c("a", "b", "a", "d"), to = c("b", "c", "c", "e"),
                             w = c(4, 1, 2, 3)), diagonal = c(2, 0, 1, 5, 0, 7))),
  path3 = mat_net(named_matrix(c("a", "b", "c"),
                               data.frame(from = c("a", "b"), to = c("b", "c"), w = c(10, 1)))),
  pair2 = mat_net(named_matrix(c("a", "b"), data.frame(from = "a", to = "b", w = 5),
                               diagonal = c(3, 0))),
  star = mat_net(named_matrix(c("s", "l1", "l2", "l3", "l4"),
                              data.frame(from = c("s", "s", "s", "s", "l1"),
                                         to = c("l1", "l2", "l3", "l4", "l2"),
                                         w = c(1, 2, 3, 4, 2)))),
  # all weights tied: a 5-cycle with one chord (several equal shortest paths)
  tied_cycle = mat_net(named_matrix(LETTERS[1:5],
                                    data.frame(from = c("A", "B", "C", "D", "E", "A"),
                                               to = c("B", "C", "D", "E", "A", "C"),
                                               w = rep(2, 6)), diagonal = rep(2, 5))),
  # planted communities: two 4-cliques joined by one weak tie
  two_cliques = mat_net(named_matrix(
    c(paste0("p", 1:4), paste0("q", 1:4)),
    data.frame(from = c("p1", "p1", "p1", "p2", "p2", "p3", "q1", "q1", "q1", "q2", "q2", "q3", "p1"),
               to   = c("p2", "p3", "p4", "p3", "p4", "p4", "q2", "q3", "q4", "q3", "q4", "q4", "q1"),
               w    = c(9, 11, 8, 10, 12, 9, 10, 8, 11, 9, 12, 10, 1)),
    diagonal = c(3, 0, 2, 1, 0, 4, 1, 2))),
  # three triangles on a ring of weak ties
  three_groups = mat_net(named_matrix(
    paste0("g", 1:9),
    data.frame(from = c("g1", "g1", "g2", "g4", "g4", "g5", "g7", "g7", "g8", "g3", "g6", "g9"),
               to   = c("g2", "g3", "g3", "g5", "g6", "g6", "g8", "g9", "g9", "g4", "g7", "g1"),
               w    = c(7, 6, 8, 9, 7, 6, 8, 6, 7, 1, 2, 1)))))

# Sparse event networks: many codes, short windows (isolates, several components).
sparse_net <- function(seed) with_seed(seed, cna("cna_build_event_nest")(
  make_event_data(n_actor = 5, actions = LETTERS[1:9], per = 6, seed = seed),
  opts_event(buildModel_threshold = 90)))

# The helper fixtures call set.seed(); with_seed() restores the RNG afterwards.
pool <- with_seed(20261106, c(
  stats::setNames(lapply(seq_len(25), mc_event_net), sprintf("event_seed%02d", seq_len(25))),
  stats::setNames(lapply(100 + seq_len(5), sparse_net), sprintf("sparse_%d", 100 + seq_len(5))),
  list(reg_20 = reg_net(20, 1), reg_60 = reg_net(60, 2), reg_150 = reg_net(150, 3),
       delim_counts = cna("cna_build_delimited_nest")(make_delimited_data(), opts_delimited()),
       delim_cosine = cna("cna_build_delimited_nest")(
         make_delimited_data(n = 30, seed = 8),
         opts_delimited(buildModel_similarity = "cosine", buildModel_diagonal = FALSE)),
       fixture = fixture_event_net()),
  hand_nets))

weights_of <- function(net) unclass(net$weights)
n_edges_off <- function(net) {
  w <- weights_of(net)
  sum(w[upper.tri(w)] != 0)
}
has_diag <- function(net) any(diag(weights_of(net)) != 0)

# ---------------------------------------------------------------------------
# Independent references (igraph and base R only)
# ---------------------------------------------------------------------------

# Undirected igraph graph straight from the weight matrix; edges in canonical
# order (upper triangle with the diagonal, column-major).
ref_graph <- function(w, loops = FALSE) {
  w <- unclass(w)
  if (!loops) diag(w) <- 0
  idx <- which(upper.tri(w, diag = TRUE) & w != 0, arr.ind = TRUE)
  g <- igraph::make_empty_graph(nrow(w), directed = FALSE)
  g <- igraph::set_vertex_attr(g, "name", value = rownames(w))
  igraph::add_edges(g, as.vector(t(idx)), weight = w[idx])
}

lengths_of <- function(g, invert) {
  w <- igraph::E(g)$weight
  if (invert) 1 / w else w
}

max_scale <- function(v, on) {
  if (!on || !any(!is.na(v))) return(v)
  m <- max(v, na.rm = TRUE)
  if (m > 0) v / m else v
}

# PageRank as the solution of x = d M x + (1 - d) / n with dangling rows
# teleporting uniformly (Brin & Page 1998).
ref_pagerank <- function(w, damping = 0.85) {
  n <- nrow(w)
  s <- rowSums(w)
  dangling <- s == 0
  p <- w / ifelse(dangling, 1, s)
  m <- t(p) + outer(rep(1 / n, n), as.numeric(dangling))
  as.numeric(solve(diag(n) - damping * m, rep((1 - damping) / n, n)))
}

# Principal (Perron) eigenvector of a non-negative symmetric matrix, max = 1.
ref_eigenvector <- function(w) {
  if (all(w == 0)) return(rep(1, nrow(w)))
  e <- eigen(unname(w), symmetric = TRUE)
  v <- abs(e$vectors[, which.max(e$values)])
  v / max(v)
}

# Zhang & Horvath (2005) weighted clustering on w / max(w), diagonal excluded,
# written as the explicit double sum over neighbour pairs j != k.
ref_zhang_horvath <- function(w) {
  w <- unname(unclass(w))
  diag(w) <- 0
  if (!any(w > 0)) return(rep(NaN, nrow(w)))
  w <- w / max(w)
  vapply(seq_len(nrow(w)), function(i) {
    pairs <- outer(w[i, ], w[i, ])
    diag(pairs) <- 0
    sum(pairs * w) / sum(pairs)
  }, numeric(1))
}

# tna's diffusion centrality: row sums of W + W^2 + ... + W^n (Banerjee et al.
# 2013), range-normalised to 0-1; a constant vector has no range (NaN).
ref_diffusion <- function(w) {
  w <- unname(unclass(w))
  powers <- Reduce(function(p, i) p %*% w, seq_len(ncol(w)), diag(ncol(w)), accumulate = TRUE)
  r <- rowSums(Reduce(`+`, powers[-1L]))
  span <- max(r) - min(r)
  if (!isTRUE(span > 1e-10 * max(abs(r)))) return(rep(NaN, length(r)))
  (r - min(r)) / span
}

# Burt (1992): c_i = sum_{j in N(i)} (p_ij + sum_{q != i, j} p_iq p_qj)^2.
# Undefined (NaN) for an isolate; with loops an isolate that has a self-loop
# scores 0, as in igraph::constraint().
ref_constraint <- function(w, loops = FALSE) {
  a <- unclass(w)
  self <- diag(a) != 0
  diag(a) <- 0
  p <- a / ifelse(rowSums(a) > 0, rowSums(a), 1)
  vapply(seq_len(nrow(a)), function(i) {
    nb <- which(a[i, ] != 0)
    if (!length(nb)) return(if (loops && self[i]) 0 else NaN)
    sum(vapply(nb, function(j) {
      q <- setdiff(seq_len(nrow(a)), c(i, j))
      (p[i, j] + sum(p[i, q] * p[q, j]))^2
    }, numeric(1)))
  }, numeric(1))
}

# Burt (1992) effective size for binary ties: n - 2t / n, with n the contacts
# of ego (ego excluded) and t the ties among them. An isolate scores 0.
ref_effective_size <- function(w) {
  a <- (unclass(w) != 0) * 1
  diag(a) <- 0
  vapply(seq_len(nrow(a)), function(i) {
    nb <- which(a[i, ] != 0)
    k <- length(nb)
    if (!k) return(0)
    k - sum(a[nb, nb]) / k
  }, numeric(1))
}

# Simple RSP betweenness, Kivimaki et al. (2016):
# bet = diag(Z (Z^div - n D_Z^-1)^T Z), Z = (I - W)^-1, W = P_ref * exp(-beta C),
# C = 1/w, computed per connected component (pairs without a path add nothing).
ref_rsp <- function(w, beta = 0.01) {
  a <- unclass(w)
  diag(a) <- 0
  comp <- igraph::components(ref_graph(a))$membership
  parts <- lapply(split(seq_len(nrow(a)), comp), function(idx) {
    nc <- length(idx)
    if (nc < 2L) return(stats::setNames(0, idx))
    ac <- a[idx, idx, drop = FALSE]
    wm <- (ac / rowSums(ac)) * exp(-beta * ifelse(ac > 0, 1 / ac, 0))
    z <- solve(diag(nc) - wm)
    x <- 1 / z - nc * diag(1 / diag(z), nc)
    stats::setNames(diag(z %*% t(x) %*% z), idx)
  })
  v <- unlist(unname(parts))
  out <- numeric(nrow(a))
  out[as.integer(names(v))] <- v
  out
}

# Every centrality CNA offers except Eccentricity (tested on its own below).
REF_MEASURES <- c("Strength", "Degree", "Coreness", "Betweenness", "BetweennessRSP",
                  "Closeness", "Harmonic", "Eigenvector", "PageRank", "Diffusion",
                  "Clustering", "Constraint", "EffectiveSize")

ref_centrality <- function(w, loops = FALSE, normalized = FALSE, invert = TRUE) {
  w <- unclass(w)
  wl <- w
  if (!loops) diag(wl) <- 0
  g <- ref_graph(w, loops)          # the graph the measure sees
  gs <- ref_graph(w, FALSE)         # loop-free, for paths
  len <- lengths_of(gs, invert)
  sc <- function(v) max_scale(v, normalized)
  data.frame(
    state = rownames(w),
    Strength = sc(igraph::strength(g)),
    Degree = sc(igraph::degree(g)),
    Coreness = sc(igraph::coreness(g)),
    Betweenness = sc(igraph::betweenness(gs, weights = len)),
    BetweennessRSP = sc(ref_rsp(w)),
    Closeness = igraph::closeness(gs, weights = len, normalized = normalized),
    Harmonic = igraph::harmonic_centrality(gs, weights = len, normalized = normalized),
    Eigenvector = sc(ref_eigenvector(wl)),
    PageRank = sc(ref_pagerank(wl)),
    Diffusion = ref_diffusion(wl),   # already 0-1; "normalize" leaves it
    Clustering = ref_zhang_horvath(w),   # already 0-1; "normalize" leaves it
    Constraint = sc(ref_constraint(w, loops)),
    EffectiveSize = sc(ref_effective_size(w)),
    row.names = NULL, stringsAsFactors = FALSE)
}

# CNA's centrality table vs the reference, measure by measure, aligned by name.
expect_centralities <- function(net, label, measures = REF_MEASURES, ...) {
  out <- cna("cna_centrality_df")(net, measures, ...)
  ref <- ref_centrality(net$weights, ...)
  expect_setequal(out$state, ref$state)
  ord <- match(ref$state, out$state)
  invisible(lapply(measures, function(m)
    expect_equal(unname(out[[m]][ord]), unname(ref[[m]]), info = paste(label, m))))
}

edge_key <- function(a, b) paste(pmin(a, b), pmax(a, b), sep = "|")

# ---------------------------------------------------------------------------
# 1. Centralities
# ---------------------------------------------------------------------------

test_that("13 centralities == igraph / base-R formulas (defaults: no loops, inverted)", {
  invisible(Map(expect_centralities, pool, names(pool)))
})

test_that("normalized = TRUE: closeness/harmonic as igraph, the rest divided by the max", {
  invisible(Map(function(net, nm) expect_centralities(net, nm, normalized = TRUE),
                pool, names(pool)))
})

test_that("invert = FALSE: path measures take the raw weights as lengths", {
  invisible(Map(function(net, nm) expect_centralities(net, nm, invert = FALSE),
                pool, names(pool)))
})

test_that("loops = TRUE: degree/strength count a loop twice, spectral measures once", {
  invisible(Map(function(net, nm) expect_centralities(net, nm, measures = REF_MEASURES,
                                                     loops = TRUE),
                pool, names(pool)))
})

test_that("second reference: igraph's eigenvector, PageRank and constraint (loop-free)", {
  invisible(Map(function(net, nm) {
    out <- cna("cna_centrality_df")(net, c("Eigenvector", "PageRank", "Constraint"))
    g <- ref_graph(net$weights)
    ord <- match(igraph::V(g)$name, out$state)
    expect_equal(out$PageRank[ord], unname(igraph::page_rank(g)$vector), info = nm)
    expect_equal(out$Constraint[ord], unname(igraph::constraint(g)), info = nm)
    if (igraph::is_connected(g))
      expect_equal(out$Eigenvector[ord], unname(igraph::eigen_centrality(g)$vector), info = nm)
  }, pool, names(pool)))
})

test_that("Eccentricity, invert = FALSE == igraph eccentricity on raw lengths", {
  invisible(Map(function(net, nm) {
    out <- cna("cna_centrality_df")(net, "Eccentricity", invert = FALSE)
    g <- ref_graph(net$weights)
    ref <- igraph::eccentricity(g, weights = lengths_of(g, FALSE))
    expect_equal(out$Eccentricity[match(names(ref), out$state)], unname(ref), info = nm)
  }, pool, names(pool)))
})

# Fixed in CNA 2026-11-06; the block pins the fix. Was: (CNA option mapping) the UI promises "Invert weights for path-based
# measures" and lists eccentricity as a path measure, but cograph documents
# that eccentricity ignores invert_weights, so with the default Invert = TRUE
# a strong co-occurrence counts as a LONG distance. Minimal repro: path
# a-b (10), b-c (1). Inverted lengths 0.1 and 1 give eccentricities
# (1.1, 1, 1.1); CNA reports (11, 10, 11).
test_that("Eccentricity with invert = TRUE uses 1/w lengths (fixed)", {
  net <- hand_nets$path3
  out <- cna("cna_centrality_df")(net, "Eccentricity", invert = TRUE)
  expect_equal(out$Eccentricity, c(1.1, 1, 1.1))
  invisible(Map(function(net, nm) {
    out <- cna("cna_centrality_df")(net, "Eccentricity", invert = TRUE)
    g <- ref_graph(net$weights)
    ref <- igraph::eccentricity(g, weights = lengths_of(g, TRUE))
    expect_equal(out$Eccentricity[match(names(ref), out$state)], unname(ref), info = nm)
  }, pool[c("fixture", "reg_60", "delim_counts", "two_components")],
  c("fixture", "reg_60", "delim_counts", "two_components")))
})

# DEFECT (cograph .cg_betweenness): ties between shortest paths are detected
# with an ABSOLUTE tolerance of 1e-15, so mathematically equal path lengths
# that differ by a few ulps at larger magnitude count as different. Max-scaled
# weights make the inverted lengths large. Minimal repro (counts, then
# scaling = "max"): A-B 3, B-C 6, A-C 2, C-D 997, D-E 1. With inverted lengths,
# A-B-C (997/3 + 997/6) and A-C (997/2) tie exactly, so B lies on half of the
# geodesics from A to C, D and E: betweenness(B) = 1.5 (igraph agrees). CNA
# reports 0.25.
test_that("DEFECT: betweenness counts floating-point ties among scaled weights", {
  m <- named_matrix(LETTERS[1:5], data.frame(from = c("A", "B", "A", "C", "D"),
                                             to = c("B", "C", "C", "D", "E"),
                                             w = c(3, 6, 2, 997, 1)))
  net <- Nestimate::build_network(m, method = "co_occurrence", scaling = "max")
  out <- cna("cna_centrality_df")(net, "Betweenness")
  expect_equal(out$Betweenness, c(0, 1.5, 4, 3, 0))
  g <- ref_graph(net$weights)
  expect_equal(unname(igraph::betweenness(g, weights = lengths_of(g, TRUE))), c(0, 1.5, 4, 3, 0))
})

test_that("Monte Carlo: betweenness on max-scaled event networks == igraph", {
  nets <- lapply(seq_len(25), mc_event_net, scaling = "max")
  invisible(Map(function(net, s) {
    out <- cna("cna_centrality_df")(net, "Betweenness")
    g <- ref_graph(net$weights)
    ref <- igraph::betweenness(g, weights = lengths_of(g, TRUE))
    expect_equal(out$Betweenness[match(names(ref), out$state)], unname(ref),
                 info = sprintf("max-scaled seed %d", s))
  }, nets, seq_along(nets)))
})

# Fixed in CNA 2026-11-06; the block pins the fix. Was: (cograph calculate_effective_size; CNA now computes it
# without loops) with loops = TRUE the neighbour
# lists contain the node itself (twice for an undirected loop), so an alter
# with a self-loop counts as "redundant with itself". Burt's sums exclude ego
# and j, so a self-loop cannot change effective size. Minimal repro, the
# two_components net: b has contacts {a, c}, one tie among them, so
# ES(b) = 2 - 2/2 = 1. CNA gives 1 with loops = FALSE and 0 with loops = TRUE.
test_that("effective size ignores self-loops (Burt) when loops = TRUE (fixed)", {
  invisible(Map(function(net, nm) {
    out <- cna("cna_centrality_df")(net, "EffectiveSize", loops = TRUE)
    ref <- ref_effective_size(net$weights)
    expect_equal(out$EffectiveSize, ref, info = nm)
  }, pool[c("two_components", "fixture", "reg_60")], c("two_components", "fixture", "reg_60")))
})

test_that("RSP betweenness: the invert switch does not change it (cost is always 1/w)", {
  net <- pool$reg_60
  f <- cna("cna_centrality_df")
  expect_equal(f(net, "BetweennessRSP", invert = TRUE)$BetweennessRSP,
               f(net, "BetweennessRSP", invert = FALSE)$BetweennessRSP)
})

# ---------------------------------------------------------------------------
# 2. Edge betweenness
# ---------------------------------------------------------------------------

expect_edge_betweenness <- function(net, label, invert) {
  out <- cna("cna_edge_betweenness_df")(cna("cna_edge_betweenness_compute")(net, invert = invert))
  if (!n_edges_off(net)) return(expect_null(out, label = label))
  g <- ref_graph(net$weights)
  ends <- igraph::ends(g, igraph::E(g))
  ref <- stats::setNames(igraph::edge_betweenness(g, weights = lengths_of(g, invert),
                                                  directed = FALSE),
                         edge_key(ends[, 1], ends[, 2]))
  expect_setequal(edge_key(out$from, out$to), names(ref))
  expect_equal(out$betweenness, unname(ref[edge_key(out$from, out$to)]), info = label)
}

test_that("edge betweenness == igraph across the pool, inverted and raw lengths", {
  invisible(Map(function(net, nm) {
    expect_edge_betweenness(net, nm, invert = TRUE)
    expect_edge_betweenness(net, nm, invert = FALSE)
  }, pool, names(pool)))
})

test_that("edge betweenness == igraph on max-scaled networks; pool covers disconnected ones", {
  nets <- lapply(seq_len(25), mc_event_net, scaling = "max")
  invisible(Map(function(net, s) expect_edge_betweenness(net, sprintf("scaled %d", s), TRUE),
                nets, seq_along(nets)))
  disconnected <- Filter(function(net) !igraph::is_connected(ref_graph(net$weights)), pool)
  expect_gte(length(disconnected), 8)
})

# ---------------------------------------------------------------------------
# 3. Communities
# ---------------------------------------------------------------------------

same_partition <- function(a, b) isTRUE(all(outer(a, a, `==`) == outer(b, b, `==`)))

# igraph's edge-betweenness clustering says which cut of its dendrogram it
# returns; this informational warning is raised for CNA and the reference.
muffle_eb_note <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (grepl("Membership vector will be selected", conditionMessage(w), fixed = TRUE))
    invokeRestart("muffleWarning")
})

# The igraph call each UI method maps to, with cograph's documented defaults.
ref_partition <- function(g, method, gamma, seed) with_seed(seed, switch(
  method,
  louvain = igraph::cluster_louvain(g, resolution = gamma),
  leiden = igraph::cluster_leiden(g, resolution = gamma, objective_function = "CPM",
                                  beta = 0.01, n_iterations = 2),
  walktrap = igraph::cluster_walktrap(g, steps = 4),
  fast_greedy = igraph::cluster_fast_greedy(g),
  infomap = igraph::cluster_infomap(g, nb.trials = 10),
  label_propagation = igraph::cluster_label_prop(g, mode = "out"),
  edge_betweenness = muffle_eb_note(igraph::cluster_edge_betweenness(g, directed = FALSE)),
  leading_eigenvector = igraph::cluster_leading_eigen(g),
  spinglass = igraph::cluster_spinglass(g, spins = 25, gamma = gamma),
  optimal = igraph::cluster_optimal(g)))

expect_communities <- function(net, label, method, gamma = 1, seed = 1L) {
  g <- ref_graph(net$weights, loops = TRUE)
  out <- muffle_eb_note(cna("cna_community_df")(
    cna("cna_community_compute")(net, method, gamma = gamma, seed = seed)))
  mem <- igraph::membership(ref_partition(g, method, gamma, seed))
  info <- sprintf("%s / %s / gamma %s / seed %d", label, method, gamma, seed)
  expect_setequal(out$state, names(mem))
  cna_mem <- out$community[match(names(mem), out$state)]
  expect_true(same_partition(cna_mem, as.integer(mem)), info = info)
  expect_equal(igraph::compare(cna_mem, as.integer(mem), method = "nmi"), 1, info = info)
  expect_equal(igraph::modularity(g, cna_mem, weights = igraph::E(g)$weight),
               igraph::modularity(g, as.integer(mem), weights = igraph::E(g)$weight), info = info)
}

connected_pool <- Filter(function(net) igraph::is_connected(ref_graph(net$weights)), pool)

test_that("all 10 community methods == igraph::cluster_* called directly (seeded)", {
  methods <- cna(".cna_community_methods")
  nets <- connected_pool[intersect(c("fixture", "reg_20", "reg_60", "reg_150", "delim_counts",
                                     "delim_cosine", "star", "tied_cycle", "path3", "pair2",
                                     "two_cliques", "three_groups"),
                                   names(connected_pool))]
  expect_gte(length(nets), 6)
  invisible(lapply(methods, function(m) Map(function(net, nm) {
    expect_communities(net, nm, m, seed = 3L)
  }, nets, names(nets))))
})

test_that("Monte Carlo: every method on random event networks, two seeds", {
  skip_on_cran()
  methods <- setdiff(cna(".cna_community_methods"), "spinglass")
  nets <- pool[sprintf("event_seed%02d", seq_len(12))]
  invisible(lapply(methods, function(m) Map(function(net, nm) {
    expect_communities(net, nm, m, seed = 1L)
    expect_communities(net, nm, m, seed = 42L)
  }, nets, names(nets))))
})

test_that("resolution reaches Louvain/Leiden and gamma reaches spinglass (== igraph)", {
  nets <- connected_pool[c("reg_60", "reg_150", "fixture", "two_cliques", "three_groups")]
  invisible(lapply(c(0.5, 1, 2, 5, 40), function(gm) Map(function(net, nm) {
    expect_communities(net, nm, "louvain", gamma = gm, seed = 5L)
    expect_communities(net, nm, "leiden", gamma = gm, seed = 5L)
  }, nets, names(nets))))
  invisible(lapply(c(0.5, 2), function(gm)
    expect_communities(nets$reg_60, "reg_60", "spinglass", gamma = gm, seed = 5L)))
})

test_that("Leiden (CPM at gamma = 1) puts all codes in one community on count networks", {
  net <- pool$reg_60
  out <- cna("cna_community_df")(cna("cna_community_compute")(net, "leiden", gamma = 1, seed = 1L))
  expect_length(unique(out$community), 1L)
})

test_that("spinglass on a disconnected network: largest component only, with a warning", {
  net <- hand_nets$two_components
  expect_warning(out <- cna("cna_community_df")(
    cna("cna_community_compute")(net, "spinglass", seed = 2L)), "largest component")
  expect_setequal(out$state, c("a", "b", "c"))
  g <- ref_graph(net$weights, loops = TRUE)
  comp <- igraph::components(g)
  sub <- igraph::induced_subgraph(g, which(comp$membership == which.max(comp$csize)))
  mem <- igraph::membership(ref_partition(sub, "spinglass", 1, 2L))
  expect_true(same_partition(out$community[match(names(mem), out$state)], as.integer(mem)))
})

# ---------------------------------------------------------------------------
# 4. Cliques and the clique complex
# ---------------------------------------------------------------------------

thresholded_graph <- function(w, thr) {
  a <- unclass(w)
  diag(a) <- 0
  igraph::graph_from_adjacency_matrix((a > 0 & a >= thr) * 1, mode = "undirected")
}

clique_keys <- function(sets) sort(vapply(sets, function(v) paste(sort(v), collapse = ", "), ""))

test_that("clique enumeration == igraph::cliques for sizes 2-5 across thresholds", {
  nets <- pool[c(sprintf("event_seed%02d", seq_len(10)), "reg_60", "delim_counts", "fixture",
                 "two_components", "tied_cycle")]
  invisible(Map(function(net, nm) {
    w <- weights_of(net)
    pos <- w[upper.tri(w) & w > 0]
    thrs <- unique(c(0, if (length(pos)) stats::quantile(pos, c(0.5, 0.75), names = FALSE)))
    lapply(thrs, function(thr) lapply(2:5, function(k) {
      out <- cna("cna_cliques_df")(cna("cna_cliques_compute")(net, size = k, threshold = thr), net, k)
      ref <- igraph::cliques(thresholded_graph(w, thr), min = k, max = k)
      ref_keys <- clique_keys(lapply(ref, names))
      info <- sprintf("%s thr %.3f size %d", nm, thr, k)
      if (!length(ref_keys)) return(expect_null(out, info = info))
      expect_equal(clique_keys(strsplit(out$members, ", ", fixed = TRUE)), ref_keys, info = info)
      stats <- t(vapply(strsplit(out$members, ", ", fixed = TRUE), function(v) {
        e <- w[v, v][upper.tri(w[v, v])]
        c(mean(e), min(e))
      }, numeric(2)))
      expect_equal(unname(cbind(out$mean_weight, out$min_weight)), unname(stats), info = info)
    }))
  }, nets, names(nets)))
})

test_that("clique complex: f-vector, b0 and max q == igraph clique counts / components", {
  nets <- pool[c(sprintf("event_seed%02d", seq_len(10)), "reg_60", "two_components", "star")]
  invisible(Map(function(net, nm) {
    w <- weights_of(net)
    g <- thresholded_graph(w, 0)
    sm <- cna("cna_complex_df")(cna("cna_cliques_compute")(net, size = 3L, threshold = 0))$rows
    f <- sm$value[grepl("^Cliques of", sm$metric)]
    ref_f <- vapply(seq_len(igraph::clique_num(g)), function(k)
      as.numeric(length(igraph::cliques(g, min = k, max = k))), numeric(1))
    expect_equal(f, ref_f, info = nm)
    expect_equal(sm$value[sm$metric == "Betti number b0"], igraph::components(g)$no, info = nm)
    expect_equal(sm$value[sm$metric == "Maximum q (Q-analysis)"], igraph::clique_num(g) - 1,
                 info = nm)
  }, nets, names(nets)))
})

# ---------------------------------------------------------------------------
# 5. Backbone
# ---------------------------------------------------------------------------

# Serrano, Boguna & Vespignani (2009): alpha_ij = (1 - w_ij / s_i)^(k_i - 1);
# an undirected edge stays when alpha < level from either endpoint.
# `count_diagonal`: whether the diagonal enters s_i and k_i.
serrano_keep <- function(w, level, count_diagonal) {
  a <- unclass(w)
  if (!count_diagonal) diag(a) <- 0
  s <- rowSums(a)
  k <- rowSums(a > 0)
  alpha <- (1 - a / ifelse(s > 0, s, 1))^(k - 1)
  off <- unclass(w)
  diag(off) <- 0
  keep <- (alpha < level | t(alpha) < level) & off > 0
  dimnames(keep) <- dimnames(off)
  keep
}

backbone_keep <- function(net, ...) {
  k <- cna("cna_backbone_compute")(net, ...)[[1]]$keep
  diag(k) <- FALSE
  k
}

LEVELS <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9)

test_that("disparity filter == Serrano et al. (2009) on loop-free networks, all levels", {
  nets <- c(list(delim_cosine = pool$delim_cosine,
                 delim_nodiag = with_seed(5, cna("cna_build_delimited_nest")(
                   make_delimited_data(n = 60, seed = 5), opts_delimited(buildModel_diagonal = FALSE))),
                 star = hand_nets$star, path3 = hand_nets$path3),
            lapply(pool[c(sprintf("event_seed%02d", seq_len(10)), "reg_60", "fixture")], strip_diag))
  invisible(Map(function(net, nm) lapply(LEVELS, function(lv)
    expect_equal(backbone_keep(net, method = "disparity", level = lv),
                 serrano_keep(net$weights, lv, count_diagonal = FALSE),
                 info = sprintf("%s level %.2f", nm, lv))), nets, names(nets)))
})

test_that("disparity filter with self-co-occurrences: the diagonal counts in s_i and k_i", {
  nets <- Filter(has_diag, pool)
  expect_gt(length(nets), 10)
  invisible(Map(function(net, nm) lapply(LEVELS, function(lv)
    expect_equal(backbone_keep(net, method = "disparity", level = lv),
                 serrano_keep(net$weights, lv, count_diagonal = TRUE),
                 info = sprintf("%s level %.2f", nm, lv))), nets, names(nets)))
  # the convention matters: on the fixture at level 0.5 it differs from the
  # loop-free formula
  fx <- pool$fixture
  expect_false(identical(serrano_keep(fx$weights, 0.5, TRUE), serrano_keep(fx$weights, 0.5, FALSE)))
})

# Threshold/lowest backbones of a connected network: every edge above the
# cut-off stays, the backbone stays connected, and a kept edge at or below the
# cut-off is a bridge (removing it would disconnect).
expect_cut_backbone <- function(net, keep, cut, label) {
  w <- weights_of(net)
  off <- w
  diag(off) <- 0
  expect_true(all(keep[off > cut]), info = label)
  expect_false(any(keep[off == 0]), info = label)
  bb <- ref_graph(off * keep)
  expect_true(igraph::is_connected(bb), info = label)
  low <- which(upper.tri(off) & keep & off <= cut, arr.ind = TRUE)
  if (nrow(low)) {
    ends <- igraph::ends(bb, igraph::bridges(bb))
    expect_true(all(edge_key(rownames(off)[low[, 1]], colnames(off)[low[, 2]]) %in%
                      edge_key(ends[, 1], ends[, 2])), info = label)
  }
  # when the edges above the cut-off already connect the network, it is exact
  if (igraph::is_connected(ref_graph(off * (off > cut))))
    expect_equal(unname(keep), unname((off > cut) & off > 0), info = label)
}

test_that("threshold backbone on connected networks: weight <= threshold removed unless a bridge", {
  invisible(Map(function(net, nm) {
    w <- weights_of(net)
    pos <- w[upper.tri(w) & w > 0]
    lapply(stats::quantile(pos, c(0.1, 0.3, 0.5, 0.8), names = FALSE), function(thr)
      expect_cut_backbone(net, backbone_keep(net, method = "threshold", threshold = thr), thr,
                          sprintf("%s threshold %.3f", nm, thr)))
  }, connected_pool, names(connected_pool)))
})

test_that("lowest backbone: cut-off = quantile over all positive cells (diagonal included)", {
  invisible(Map(function(net, nm) lapply(c(0.05, 0.2, 0.5), function(p) {
    w <- weights_of(net)
    cut <- stats::quantile(w[w > 0], probs = p, names = FALSE)
    expect_cut_backbone(net, backbone_keep(net, method = "lowest", lowest = p), cut,
                        sprintf("%s lowest %.2f", nm, p))
  }), connected_pool, names(connected_pool)))
})

# DEFECT (Nestimate .prune_default, inherited from tna): the connectivity guard
# removes an edge only if the WHOLE network is weakly connected afterwards, so
# a network that is already disconnected (one isolated code is enough) loses
# no edge at all. Minimal repro: a-b 5, a-c 4, b-c 1, c-d 6, d-e 1 plus an
# isolated f, threshold 2: b-c (1) is not a bridge and must go (d-e stays,
# it is one). CNA keeps every edge.
test_that("DEFECT: threshold backbone also prunes networks that have an isolated node", {
  m <- named_matrix(letters[1:6], data.frame(from = c("a", "a", "b", "c", "d"),
                                             to = c("b", "c", "c", "d", "e"),
                                             w = c(5, 4, 1, 6, 1)))
  df <- cna("cna_backbone_df")(cna("cna_backbone_compute")(mat_net(m), method = "threshold",
                                                           threshold = 2))
  status <- stats::setNames(df$status, edge_key(df$from, df$to))
  expect_equal(unname(status[edge_key("b", "c")]), "removed")
  expect_equal(unname(status[edge_key("d", "e")]), "kept")
  # same network without the isolate: b-c is removed (passes today)
  df2 <- cna("cna_backbone_df")(cna("cna_backbone_compute")(mat_net(m[1:5, 1:5]),
                                                            method = "threshold", threshold = 2))
  expect_equal(df2$status[edge_key(df2$from, df2$to) == edge_key("b", "c")], "removed")
})

test_that("backbone table: one row per off-diagonal edge, weights unchanged", {
  invisible(Map(function(net, nm) {
    df <- cna("cna_backbone_df")(cna("cna_backbone_compute")(net, method = "disparity", level = 0.3))
    if (!n_edges_off(net)) return(expect_null(df))
    w <- weights_of(net)
    expect_equal(nrow(df), n_edges_off(net), info = nm)
    expect_equal(df$weight, unname(w[cbind(df$from, df$to)]), info = nm)
  }, pool, names(pool)))
})

# ---------------------------------------------------------------------------
# 6. Link prediction
# ---------------------------------------------------------------------------

LINK_METHODS <- c("resource_allocation", "adamic_adar", "common_neighbors", "jaccard",
                  "preferential_attachment", "katz")

# Scores of every node pair by the published formula, on a loop-free matrix.
ref_link_scores <- function(w, method) {
  w <- unclass(w)
  a <- (w != 0) * 1
  k <- rowSums(a)
  s <- switch(method,
    common_neighbors = w %*% w,
    resource_allocation = w %*% diag(ifelse(k > 0, 1 / k, 0), nrow(w)) %*% w,
    adamic_adar = w %*% diag(ifelse(k > 1, 1 / log(k), 0), nrow(w)) %*% w,
    jaccard = {
      inter <- a %*% a
      uni <- outer(k, k, "+") - inter
      ifelse(uni > 0, inter / uni, 0)
    },
    preferential_attachment = outer(k, k),
    katz = {
      rho <- max(abs(eigen(a, symmetric = TRUE, only.values = TRUE)$values))
      if (rho > 0) solve(diag(nrow(a)) - (0.9 / rho) * a) - diag(nrow(a)) else 0 * a
    })
  dimnames(s) <- dimnames(w)
  s
}

absent_pairs <- function(w) {
  idx <- which(upper.tri(w) & w == 0, arr.ind = TRUE)
  data.frame(from = rownames(w)[idx[, 1]], to = colnames(w)[idx[, 2]], stringsAsFactors = FALSE)
}

expect_links <- function(net, label, method, threshold = 0) {
  out <- cna("cna_links_df")(net, method = method, threshold = threshold, top_n = 1000L)
  w <- weights_of(net)
  w[w <= threshold] <- 0
  pairs <- absent_pairs(w)
  info <- sprintf("%s / %s / threshold %s", label, method, threshold)
  if (!nrow(pairs)) return(expect_null(out, info = info))
  ref <- ref_link_scores(w, method)
  expect_setequal(edge_key(out$from, out$to), edge_key(pairs$from, pairs$to))
  expect_true(all(out$existing == "no"), info = info)
  expect_equal(out$score, unname(ref[cbind(out$from, out$to)]), info = info)
  expect_equal(out$rank, seq_len(nrow(out)), info = info)
  expect_false(is.unsorted(rev(out$score)), info = info)
}

loopfree_pool <- c(lapply(pool[c(sprintf("event_seed%02d", seq_len(25)), "reg_20", "reg_60",
                                 "delim_counts", "fixture", "two_components", "star",
                                 "tied_cycle")], strip_diag),
                   list(delim_cosine = pool$delim_cosine, path3 = hand_nets$path3))

test_that("link scores == formulas on loop-free networks (5 methods, two thresholds)", {
  methods <- setdiff(LINK_METHODS, "jaccard")
  invisible(lapply(methods, function(m) Map(function(net, nm) {
    expect_links(net, nm, m)
    w <- weights_of(net)
    thr <- stats::quantile(w[upper.tri(w) & w > 0], 0.4, names = FALSE)
    expect_links(net, nm, m, threshold = thr)
  }, loopfree_pool, names(loopfree_pool))))
})

test_that("on 0/1 networks the scores are the set definitions (Liben-Nowell & Kleinberg 2007)", {
  binary <- lapply(loopfree_pool[sprintf("event_seed%02d", seq_len(12))], function(net) {
    m <- (weights_of(net) != 0) * 1
    mat_net(m)
  })
  invisible(Map(function(net, nm) {
    w <- weights_of(net)
    nb <- lapply(seq_len(nrow(w)), function(i) which(w[i, ] != 0))
    set_score <- function(method, x, y) {
      common <- intersect(nb[[x]], nb[[y]])
      deg <- lengths(nb)[common]
      switch(method,
        common_neighbors = length(common),
        resource_allocation = sum(1 / deg),
        adamic_adar = sum(ifelse(deg > 1, 1 / log(deg), 0)),
        preferential_attachment = length(nb[[x]]) * length(nb[[y]]),
        jaccard = if (length(union(nb[[x]], nb[[y]])))
          length(common) / length(union(nb[[x]], nb[[y]])) else 0)
    }
    lapply(c("common_neighbors", "resource_allocation", "adamic_adar",
             "preferential_attachment"), function(m) {
      out <- cna("cna_links_df")(net, method = m, top_n = 1000L)
      if (is.null(out)) return(expect_equal(nrow(absent_pairs(w)), 0L))
      x <- match(out$from, rownames(w))
      y <- match(out$to, rownames(w))
      ref <- vapply(seq_along(x), function(i) set_score(m, x[i], y[i]), numeric(1))
      expect_equal(out$score, ref, info = paste(nm, m))
    })
  }, binary, names(binary)))
})

# Fixed in CNA 2026-11-06; the block pins the fix. Was: (Nestimate .lp_jaccard with weighted = TRUE, which CNA
# used; CNA now passes the binary adjacency for Jaccard):
# the intersection is a sum of weight PRODUCTS but the union is a count of
# binary neighbours, so on co-occurrence counts the denominator goes negative
# and the score is set to 0 for every pair. The Jaccard index
# is |N(x) & N(y)| / |N(x) | N(y)|, which Nestimate documents as "ratio of
# shared neighbors to total neighbors". Minimal repro: a-b 2, b-c 3, c-d 4,
# a-c 1. Pairs a-d and b-d each share c out of two neighbours: Jaccard 1/2.
# CNA reports 0.
test_that("Jaccard link score is the neighbour-set ratio, within [0, 1] (fixed)", {
  m <- named_matrix(c("a", "b", "c", "d"), data.frame(from = c("a", "b", "c", "a"),
                                                     to = c("b", "c", "d", "c"),
                                                     w = c(2, 3, 4, 1)))
  out <- cna("cna_links_df")(mat_net(m), method = "jaccard")
  expect_equal(out$score[match(edge_key(c("a", "b"), "d"), edge_key(out$from, out$to))],
               c(0.5, 0.5))
  invisible(Map(function(net, nm) expect_links(net, nm, "jaccard"),
                loopfree_pool[c("event_seed03", "reg_20", "delim_counts", "star")],
                c("event_seed03", "reg_20", "delim_counts", "star")))
})

test_that("only absent links are scored; exclude_existing = FALSE scores every pair", {
  invisible(Map(function(net, nm) {
    w <- weights_of(net)
    out <- cna("cna_links_df")(net, method = "common_neighbors", exclude_existing = FALSE,
                               top_n = 1000L)
    n <- nrow(w)
    expect_equal(nrow(out), n * (n - 1) / 2, info = nm)
    expect_equal(out$existing == "yes", unname(w[cbind(out$from, out$to)] != 0), info = nm)
    absent <- cna("cna_links_df")(net, method = "common_neighbors", top_n = 1000L)
    expect_equal(NROW(absent), nrow(absent_pairs(w)), info = nm)
  }, pool[c("event_seed01", "event_seed05", "reg_20", "two_components", "star")],
  c("event_seed01", "event_seed05", "reg_20", "two_components", "star")))
})

# Fixed in CNA 2026-11-06; the block pins the fix. Was: (CNA passed the matrix with its diagonal to predict_links, and
# Nestimate counts a self-loop in degrees and in Katz walks): a code's
# co-occurrence with itself changes the scores of pairs of OTHER codes, which
# no link-prediction definition allows (N(x) never contains x). Minimal repro:
# the a-b-c-d network above with self-co-occurrences 5 on b and 4 on c.
# Preferential attachment of b-d becomes 3 instead of 2. Resource allocation
# of b-d drops from 4 to 3 (c, the shared neighbour, gains a degree), and both
# Katz scores change. Common neighbours of absent pairs are unaffected.
test_that("self-co-occurrence does not change other pairs' link scores (fixed)", {
  m <- named_matrix(c("a", "b", "c", "d"), data.frame(from = c("a", "b", "c", "a"),
                                                     to = c("b", "c", "d", "c"),
                                                     w = c(2, 3, 4, 1)))
  m_loop <- m
  diag(m_loop) <- c(0, 5, 4, 0)
  invisible(lapply(setdiff(LINK_METHODS, "jaccard"), function(meth) {
    plain <- cna("cna_links_df")(mat_net(m), method = meth)
    looped <- cna("cna_links_df")(mat_net(m_loop), method = meth)
    key <- edge_key(plain$from, plain$to)
    expect_equal(looped$score[match(key, edge_key(looped$from, looped$to))], plain$score,
                 info = meth)
  }))
})

# ---------------------------------------------------------------------------
# 7. Network summary and weight matrix
# ---------------------------------------------------------------------------

summary_value <- function(s, metric) s$value[s$metric == metric]

test_that("summary metrics == independent formulas (Nestimate's conventions named)", {
  invisible(Map(function(net, nm) {
    s <- cna("cna_summary_df")(net)
    w <- weights_of(net)
    off <- w
    diag(off) <- 0
    gs <- ref_graph(w)
    n <- nrow(w)
    expect_equal(summary_value(s, "Node Count"), n, info = nm)
    # edge count: unique non-zero pairs INCLUDING self-co-occurrences
    expect_equal(summary_value(s, "Edge Count"), sum(w[upper.tri(w, diag = TRUE)] != 0), info = nm)
    # mean distance: co-occurrence counts as lengths (not inverted), over
    # the connected ordered pairs
    md <- if (igraph::ecount(gs)) igraph::mean_distance(gs, weights = igraph::E(gs)$weight,
                                                        unconnected = TRUE) else NaN
    expect_equal(summary_value(s, "Mean Distance"), md, info = nm)
    # strength and degree count the diagonal once
    expect_equal(summary_value(s, "Mean Out-Strength"), mean(rowSums(w)), info = nm)
    expect_equal(summary_value(s, "SD Out-Strength"), stats::sd(rowSums(w)), info = nm)
    expect_equal(summary_value(s, "Mean In-Strength"), mean(colSums(w)), info = nm)
    expect_equal(summary_value(s, "Mean Out-Degree"), mean(rowSums(w != 0)), info = nm)
    expect_equal(summary_value(s, "SD Out-Degree"), stats::sd(rowSums(w != 0)), info = nm)
    # Freeman (1979) degree centralization on the loop-free graph; with
    # n <= 2 the maximum is 0 and Nestimate reports 0 (igraph: NaN)
    centr <- if (n > 2) igraph::centr_degree(gs, loops = FALSE, normalized = TRUE)$centralization else 0
    expect_equal(summary_value(s, "Centralization (Out-Degree)"), centr, info = nm)
    expect_equal(summary_value(s, "Reciprocity"), 1, info = nm)
  }, pool, names(pool)))
})

test_that("density of a loop-free network == edges / (n (n - 1) / 2) (Wasserman & Faust)", {
  nets <- Filter(Negate(has_diag), pool)
  expect_gt(length(nets), 2)
  invisible(Map(function(net, nm)
    expect_equal(summary_value(cna("cna_summary_df")(net), "Network Density"),
                 igraph::edge_density(ref_graph(net$weights)), info = nm), nets, names(nets)))
})

# DEFECT (Nestimate .summary_metrics_from_weights): the numerator counts
# self-co-occurrences as edges while the denominator n(n-1)/2 counts only
# pairs of distinct nodes, so density is inflated (and then capped at 1).
# Minimal repro: path a-b-c-d-e (3 edges) with diagonal
# (1, 2, 0, 3, 0). Density = 3/10 = 0.3 (0.4 even if loops counted as
# possible edges); CNA reports 0.6.
test_that("DEFECT: density does not count self-co-occurrences as edges", {
  m <- named_matrix(letters[1:5], data.frame(from = c("a", "b", "c"), to = c("b", "c", "d"),
                                             w = c(2, 3, 4)), diagonal = c(1, 2, 0, 3, 0))
  expect_equal(summary_value(cna("cna_summary_df")(mat_net(m)), "Network Density"), 0.3)
  nets <- Filter(has_diag, pool)[c("event_seed01", "event_seed02", "reg_20", "two_components")]
  invisible(Map(function(net, nm)
    expect_equal(summary_value(cna("cna_summary_df")(net), "Network Density"),
                 igraph::edge_density(ref_graph(net$weights)), info = nm), nets, names(nets)))
})

test_that("matrix table: cells == weights, union of states, NA where a group lacks a state", {
  g1 <- mat_net(named_matrix(c("A", "B", "C"), data.frame(from = c("A", "B"), to = c("B", "C"),
                                                          w = c(2, 5)), diagonal = c(1, 0, 4)))
  g2 <- mat_net(named_matrix(c("A", "B", "D"), data.frame(from = "A", to = "D", w = 7)))
  model <- structure(list(G1 = g1, G2 = g2), class = c("netobject_group", "list"))
  m <- cna("cna_matrix_df")(model)
  expect_equal(names(m), c("group", "state", "A", "B", "C", "D"))
  expect_equal(nrow(m), 6L)
  invisible(Map(function(net, g) {
    sub <- subset(m, group == g)
    w <- weights_of(net)
    cols <- colnames(w)
    expect_equal(unname(as.matrix(sub[cols])), unname(w[sub$state, cols]), info = g)
    expect_true(all(is.na(as.matrix(sub[setdiff(c("A", "B", "C", "D"), cols)]))), info = g)
  }, list(g1, g2), c("G1", "G2")))
  # a real group model: every cell equals that group's weight
  gm <- with_seed(42, fixture_group_net())
  mm <- cna("cna_matrix_df")(gm)
  invisible(lapply(names(gm), function(g) {
    sub <- subset(mm, group == g)
    w <- weights_of(gm[[g]])
    expect_equal(unname(as.matrix(sub[colnames(w)])), unname(w[sub$state, colnames(w)]))
  }))
})

# ---------------------------------------------------------------------------
# 8. Invariants and error paths
# ---------------------------------------------------------------------------

test_that("relabelling the nodes permutes every measure accordingly", {
  all_measures <- names(cna(".cna_measure_map"))
  invisible(Map(function(net, nm) with_seed(7, {
    w <- weights_of(net)
    perm <- sample(nrow(w))
    pnet <- mat_net(w[perm, perm])
    a <- cna("cna_centrality_df")(net, all_measures)
    b <- cna("cna_centrality_df")(pnet, all_measures)
    expect_equal(unname(as.matrix(b[match(a$state, b$state), all_measures])),
                 unname(as.matrix(a[all_measures])), info = nm)
    eb_a <- cna("cna_edge_betweenness_df")(cna("cna_edge_betweenness_compute")(net))
    eb_b <- cna("cna_edge_betweenness_df")(cna("cna_edge_betweenness_compute")(pnet))
    expect_equal(eb_b$betweenness[match(edge_key(eb_a$from, eb_a$to), edge_key(eb_b$from, eb_b$to))],
                 eb_a$betweenness, info = nm)
    expect_equal(backbone_keep(pnet, method = "disparity", level = 0.3)[rownames(w), rownames(w)],
                 backbone_keep(net, method = "disparity", level = 0.3), info = nm)
  }), pool[c("event_seed04", "event_seed09", "reg_60", "delim_counts", "two_components")],
  c("event_seed04", "event_seed09", "reg_60", "delim_counts", "two_components")))
})

test_that("scaling weights by c = 3: strength, closeness, harmonic scale; the rest do not", {
  # (Diffusion is not scale-free: W^k grows as c^k)
  scale_free <- c("Degree", "Coreness", "Betweenness", "Eigenvector", "PageRank",
                  "Clustering", "Constraint", "EffectiveSize")
  invisible(Map(function(net, nm) {
    w <- weights_of(net)
    big <- mat_net(3 * w)
    a <- cna("cna_centrality_df")(net, c("Strength", "Closeness", "Harmonic", scale_free))
    b <- cna("cna_centrality_df")(big, c("Strength", "Closeness", "Harmonic", scale_free))
    expect_equal(b$Strength, 3 * a$Strength, info = nm)
    expect_equal(b$Closeness, 3 * a$Closeness, info = nm)    # inverted lengths shrink by 3
    expect_equal(b$Harmonic, 3 * a$Harmonic, info = nm)
    invisible(lapply(scale_free, function(m) expect_equal(b[[m]], a[[m]], info = paste(nm, m))))
    expect_equal(backbone_keep(big, method = "disparity", level = 0.3),
                 backbone_keep(net, method = "disparity", level = 0.3), info = nm)
    la <- cna("cna_links_df")(net, method = "resource_allocation", top_n = 1000L)
    lb <- cna("cna_links_df")(big, method = "resource_allocation", top_n = 1000L)
    if (!is.null(la)) expect_equal(lb$score, 9 * la$score, info = nm)  # product of two weights
  }, pool[c("event_seed02", "event_seed07", "reg_20", "two_components", "star")],
  c("event_seed02", "event_seed07", "reg_20", "two_components", "star")))
})

test_that("error paths are classed: optimal on > 50 nodes, missing cograph column", {
  nodes <- sprintf("v%02d", seq_len(51))
  ring <- named_matrix(nodes, data.frame(from = nodes, to = c(nodes[-1], nodes[1]), w = 1))
  expect_error(cna("cna_community_compute")(mat_net(ring), "optimal"), class = "cna_too_large")
  expect_error(cna(".cna_cograph_column")(data.frame(node = "a", strength_all = 1), "pagerank"),
               class = "cna_engine_contract")
})
