# Communities: the table must equal cograph::communities() (seeded), the
# resolution must reach Louvain, and every method in the UI must be accepted.

same_partition <- function(a, b) {
  # equal up to relabelling: the two labellings define the same pairs
  isTRUE(all(outer(a, a, `==`) == outer(b, b, `==`)))
}

test_that("community table == cograph::communities (seeded)", {
  dat <- make_event_data()
  res <- run_event_single(dat, community_show_table = TRUE, community_methods = "louvain",
                          seed_set = TRUE, seed_value = 9L)
  tbl <- table_df(res, "communityTable")
  net <- cna("cna_build_event_nest")(dat, opts_event())
  ref <- cograph::communities(net, method = "louvain", resolution = 1, seed = 9L)
  ord <- match(ref$node, tbl$state)
  expect_true(same_partition(tbl$community[ord], ref$community))
})

test_that("walktrap partition matches igraph directly (independent check)", {
  net <- fixture_event_net()
  coms <- cna("cna_community_compute")(net, "walktrap")[[1]]
  w <- unclass(net$weights)
  diag(w) <- 0
  g <- igraph::graph_from_adjacency_matrix(w, mode = "undirected", weighted = TRUE)
  ref <- igraph::membership(igraph::cluster_walktrap(g))
  expect_true(same_partition(coms$community[match(names(ref), coms$node)], as.integer(ref)))
})

test_that("resolution (gamma) reaches Louvain", {
  net <- regulation_data(n_actors = 300) |>
    (\(d) Nestimate::build_network(d, method = "co_occurrence", actor = "Actor",
                                   action = "Action", time = "Time_Posix"))()
  f <- cna("cna_community_compute")
  low <- f(net, "louvain", gamma = 0.5, seed = 1L)[[1]]
  high <- f(net, "louvain", gamma = 3, seed = 1L)[[1]]
  expect_gte(length(unique(high$community)), length(unique(low$community)))
  expect_false(same_partition(low$community, high$community))
})

test_that("gamma reaches Spinglass (review fix: cograph passes resolution only to Louvain/Leiden)", {
  net <- regulation_data(n_actors = 300) |>
    (\(d) Nestimate::build_network(d, method = "co_occurrence", actor = "Actor",
                                   action = "Action", time = "Time_Posix"))()
  f <- cna("cna_community_compute")
  low <- f(net, "spinglass", gamma = 0.5, seed = 1L)[[1]]
  high <- f(net, "spinglass", gamma = 3, seed = 1L)[[1]]
  ref <- cograph::communities(net, method = "spinglass", gamma = 3, seed = 1L)
  expect_true(same_partition(high$community, ref$community))
  expect_gt(length(unique(high$community)), length(unique(low$community)))
})

test_that("every community method offered in the UI runs", {
  net <- fixture_event_net()
  methods <- cna(".cna_community_methods")
  offered <- vapply(option_def("CNAEventSingle", "community_methods")$options, `[[`, "", "name")
  expect_setequal(offered, methods)
  for (m in methods) {
    out <- cna("cna_community_df")(cna("cna_community_compute")(net, m, seed = 1L))
    expect_equal(nrow(out), nrow(net$weights), info = m)
    expect_false(anyNA(out$community), info = m)
  }
})

test_that("an unknown method is a classed error", {
  expect_error(cna("cna_community_compute")(fixture_event_net(), "spectral"),
               class = "cna_bad_input")
})

test_that("group communities: one partition per group", {
  out <- cna("cna_community_df")(cna("cna_community_compute")(fixture_group_net(), "louvain",
                                                              seed = 2L))
  expect_setequal(unique(out$group), c("G1", "G2"))
})
