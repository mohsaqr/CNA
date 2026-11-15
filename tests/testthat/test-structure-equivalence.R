# Structural analyses added in 1.1.0: edge betweenness, cliques and the clique
# complex, backbone extraction, link prediction and association rules. Each is
# checked against the engine call it wraps and, where one exists, against an
# independent computation (igraph) or a mathematical invariant.

# --- edge betweenness -----------------------------------------------------------

test_that("edge betweenness == Nestimate::net_edge_betweenness", {
  net <- fixture_event_net()
  out <- cna("cna_edge_betweenness_df")(cna("cna_edge_betweenness_compute")(net, invert = TRUE))
  ref <- Nestimate::net_edge_betweenness(net, invert = TRUE)$weights
  expect_equal(out$betweenness, unname(ref[cbind(out$from, out$to)]))
  expect_equal(out$weight, unname(net$weights[cbind(out$from, out$to)]))
  # one row per unique undirected edge, highest betweenness first
  w <- net$weights
  expect_equal(nrow(out), sum(w[upper.tri(w)] != 0))
  expect_false(is.unsorted(rev(out$betweenness)))
})

test_that("edge betweenness matches igraph on inverted weights (independent)", {
  net <- fixture_event_net()
  out <- cna("cna_edge_betweenness_df")(cna("cna_edge_betweenness_compute")(net, invert = TRUE))
  w <- unclass(net$weights)
  diag(w) <- 0
  g <- igraph::graph_from_adjacency_matrix(w, mode = "undirected", weighted = TRUE)
  eb <- igraph::edge_betweenness(g, weights = 1 / igraph::E(g)$weight, directed = FALSE)
  ends <- igraph::ends(g, igraph::E(g))
  key <- function(a, b) paste(pmin(a, b), pmax(a, b))
  ref <- stats::setNames(eb, key(ends[, 1], ends[, 2]))
  expect_equal(out$betweenness, unname(ref[key(out$from, out$to)]))
})

# --- cliques ----------------------------------------------------------------------

test_that("cliques == igraph::cliques on the thresholded graph (independent)", {
  net <- fixture_event_net()
  w <- unclass(net$weights)
  thr <- stats::median(w[upper.tri(w)])
  cx <- cna("cna_cliques_compute")(net, size = 3L, threshold = thr)
  out <- cna("cna_cliques_df")(cx, net, 3L)
  a <- w
  diag(a) <- 0
  g <- igraph::graph_from_adjacency_matrix((a >= thr & a > 0) * 1, mode = "undirected")
  ref <- igraph::cliques(g, min = 3, max = 3)
  ref_keys <- sort(vapply(ref, function(v) paste(sort(names(v)), collapse = ", "), ""))
  out_keys <- sort(vapply(strsplit(out$members, ", "), function(v) paste(sort(v), collapse = ", "), ""))
  expect_equal(out_keys, ref_keys)
  # every listed edge reaches the threshold; min/mean recomputed independently
  stats3 <- t(vapply(strsplit(out$members, ", "), function(v) {
    e <- w[v, v][upper.tri(w[v, v])]
    c(mean(e), min(e))
  }, numeric(2)))
  expect_equal(unname(cbind(out$mean_weight, out$min_weight)), unname(stats3))
  expect_true(all(out$min_weight >= thr))
})

test_that("clique complex summary satisfies Euler = alternating f-vector sum", {
  net <- fixture_event_net()
  w <- unclass(net$weights)
  cx <- cna("cna_cliques_compute")(net, size = 3L, threshold = stats::median(w[upper.tri(w)]))
  sm <- cna("cna_complex_df")(cx)$rows
  f <- sm$value[grepl("^Cliques of", sm$metric)]
  euler <- sm$value[sm$metric == "Euler characteristic"]
  expect_equal(euler, sum(f * (-1)^(seq_along(f) - 1)))
  # and Euler = alternating sum of Betti numbers (Euler-Poincare)
  betti <- sm$value[grepl("^Betti", sm$metric)]
  expect_equal(euler, sum(betti * (-1)^(seq_along(betti) - 1)))
})

test_that("an infeasible clique enumeration is refused with a classed error", {
  n <- 60
  w <- matrix(1, n, n, dimnames = list(paste0("v", 1:n), paste0("v", 1:n)))
  net <- list(weights = w)
  expect_error(cna("cna_cliques_compute")(net, size = 10L, threshold = 0),
               class = "cna_too_many_cliques")
})

# --- backbone -----------------------------------------------------------------

test_that("backbone keeps an edge when either direction survives net_prune", {
  net <- fixture_event_net()
  bb <- cna("cna_backbone_compute")(net, method = "disparity", level = 0.5)[[1]]
  pr <- Nestimate::net_prune(net, method = "disparity", level = 0.5)$weights
  expect_equal(bb$keep, ((pr != 0) | (t(pr) != 0)) & (net$weights != 0))
  expect_true(isSymmetric(unname(bb$backbone)))
  df <- cna("cna_backbone_df")(stats::setNames(list(bb), ""))
  w <- net$weights
  expect_equal(nrow(df), sum(w[upper.tri(w)] != 0))
  expect_true(all(df$status %in% c("kept", "removed")))
  # the disparity filter at level 0.5 removes some, but not all, edges here
  expect_true(any(df$status == "kept"))
})

test_that("threshold backbone removes only weights at or below the threshold", {
  net <- fixture_event_net()
  w <- net$weights
  thr <- stats::quantile(w[upper.tri(w)], 0.5)
  df <- cna("cna_backbone_df")(cna("cna_backbone_compute")(net, method = "threshold",
                                                           threshold = thr))
  expect_true(all(df$weight[df$status == "removed"] <= thr))
  expect_true(all(df$status[df$weight > thr] == "kept"))
})

# --- link prediction ------------------------------------------------------------

test_that("link prediction == Nestimate::predict_links on the thresholded matrix", {
  net <- fixture_event_net()
  w <- unclass(net$weights)
  thr <- stats::quantile(w[upper.tri(w)], 0.6)
  out <- cna("cna_links_df")(net, method = "adamic_adar", exclude_existing = TRUE,
                             threshold = thr, top_n = 10L)
  w2 <- w
  w2[abs(w2) <= thr] <- 0
  ref <- Nestimate::predict_links(w2, methods = "adamic_adar", top_n = 10L,
                                  exclude_existing = TRUE)$predictions
  expect_equal(out$from, as.character(ref$from))
  expect_equal(out$score, ref$score)
  expect_true(all(out$existing == "no"))
  # predicted pairs really are absent after thresholding
  expect_true(all(w2[cbind(out$from, out$to)] == 0))
})

test_that("a complete network has no absent links to predict", {
  dat <- regulation_data(n_actors = 300)
  net <- Nestimate::build_network(dat, method = "co_occurrence", actor = "Actor",
                                  action = "Action", time = "Time_Posix")
  expect_null(cna("cna_links_df")(net, exclude_existing = TRUE))
  expect_gt(nrow(cna("cna_links_df")(net, exclude_existing = FALSE, top_n = 5L)), 0)
})

# --- association rules ----------------------------------------------------------

test_that("association rules == Nestimate::association_rules", {
  dat <- regulation_data(n_actors = 300)
  net <- Nestimate::build_network(dat, method = "co_occurrence", actor = "Actor",
                                  action = "Action", time = "Time_Posix")
  rules <- cna("cna_rules_compute")(net, min_support = 0.05, min_confidence = 0.4)
  out <- cna("cna_rules_df")(rules, max_rows = 1000L)
  ref <- Nestimate::association_rules(net, min_support = 0.05, min_confidence = 0.4)$rules
  expect_equal(nrow(out), nrow(ref))
  key <- function(d) paste(d$antecedent, "=>", d$consequent)
  expect_setequal(key(out), key(ref))
  expect_equal(out$lift[match(key(ref), key(out))], ref$lift)
  # lift = confidence / support(consequent) >= the requested minimum
  expect_true(all(out$lift >= 1 - 1e-12))
  expect_true(all(out$support >= 0.05 - 1e-12 & out$confidence >= 0.4 - 1e-12))
  expect_false(is.unsorted(rev(out$lift)))
})

test_that("rules respect max_rows per network", {
  dat <- regulation_data(n_actors = 300)
  net <- Nestimate::build_network(dat, method = "co_occurrence", actor = "Actor",
                                  action = "Action", time = "Time_Posix", group = "Group")
  out <- cna("cna_rules_df")(cna("cna_rules_compute")(net, min_support = 0.05), max_rows = 4L)
  expect_true(all(table(out$group) <= 4))
})

# --- network description ----------------------------------------------------------

test_that("matrix table reproduces the weight matrix", {
  net <- fixture_group_net()
  m <- cna("cna_matrix_df")(net)
  sub <- m[m$group == "G1", ]
  states <- setdiff(names(m), c("group", "state"))
  expect_equal(unname(as.matrix(sub[states])), unname(unclass(net$G1$weights)[sub$state, states]))
})

test_that("network summary == Nestimate summary(), per group", {
  net <- fixture_group_net()
  s <- cna("cna_summary_df")(net)
  ref <- summary(net$G2)
  expect_equal(s$value[s$group == "G2"], as.numeric(ref$value))
  # undirected edge count: unique non-zero pairs incl. the diagonal (Nestimate >= 1.0.7)
  testthat::skip_if(utils::packageVersion("Nestimate") < "1.0.7")
  w <- net$G2$weights
  expect_equal(s$value[s$group == "G2" & s$metric == "Edge Count"],
               sum(w[upper.tri(w, diag = TRUE)] != 0))
})
