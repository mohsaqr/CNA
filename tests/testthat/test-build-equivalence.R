# Builders: every cna_build_* must give exactly what a direct Nestimate call
# with the same arguments gives (CNA is a thin wrapper; the weight matrix is
# the ground truth), plus regressions for the two build bugs fixed in 1.1.0.

test_that("event network == direct Nestimate::build_network", {
  dat <- make_event_data()
  net <- cna("cna_build_event_nest")(dat, opts_event())
  ref <- Nestimate::build_network(dat[c("Time", "Actor", "Action")], method = "co_occurrence",
                                  action = "Action", actor = "Actor", time = "Time",
                                  time_threshold = 3000)
  expect_s3_class(net, "netobject")
  expect_equal(net$weights, ref$weights)
})

test_that("ISO 8601 times keep their time of day (regression, 1.1.0)", {
  # Two sessions on the same day, an hour apart: with the old as.POSIXct()
  # conversion "2025-01-01T08:00:00Z" became midnight and the sessions merged.
  dat <- data.frame(
    Actor = "a",
    Time = c("2025-01-01T08:00:00Z", "2025-01-01T08:01:00Z",
             "2025-01-01T09:30:00Z", "2025-01-01T09:31:00Z"),
    Action = c("A", "B", "C", "D"), stringsAsFactors = FALSE)
  net <- cna("cna_build_event_nest")(dat, opts_event(buildModel_threshold = 900))
  w <- net$weights
  expect_equal(w["A", "B"], 1)              # same session
  expect_equal(w["C", "D"], 1)
  expect_equal(w["A", "C"], 0)              # different sessions never co-occur
  expect_equal(w["B", "D"], 0)
})

test_that("a numeric Order column sorts numerically", {
  f <- cna(".cna_numeric_if_possible")
  expect_identical(f(c("1", "10", "2")), c(1, 10, 2))
  expect_identical(f(c("a", "1")), c("a", "1"))
  expect_identical(f(c(3L, 1L)), c(3L, 1L))
})

test_that("group mode is a netobject_group keyed by group; missing groups are counted", {
  dat <- make_event_data()
  dat$Group[1:3] <- NA
  net <- cna("cna_build_event_nest")(dat, opts_event(buildModel_variables_group = "Group"))
  expect_s3_class(net, "netobject_group")
  expect_setequal(names(net), c("G1", "G2"))
  expect_equal(attr(net, "n_missing_group"), 3L)
})

test_that("scaling is wired through (minmax bounds the weights by 1)", {
  dat <- make_event_data()
  raw <- cna("cna_build_event_nest")(dat, opts_event())
  scaled <- cna("cna_build_event_nest")(dat, opts_event(buildModel_scaling = "minmax"))
  expect_false(isTRUE(all.equal(raw$weights, scaled$weights)))
  expect_lte(max(scaled$weights), 1 + 1e-9)
})

test_that("one-hot network == direct build_network", {
  dat <- make_onehot_data()
  codes <- paste0("code", 1:4)
  net <- cna("cna_build_onehot_nest")(dat, opts_onehot())
  ref <- Nestimate::build_network(dat[c(codes, "Actor", "Session")], method = "co_occurrence",
                                  codes = codes, format = "onehot", actor = "Actor",
                                  session = "Session", window_size = 1L,
                                  mode = "non-overlapping")
  expect_equal(net$weights, ref$weights)
})

test_that("delimited network == Nestimate::cooccurrence with the same arguments", {
  dat <- make_delimited_data()
  for (sim in c("none", "jaccard", "cosine")) {
    net <- cna("cna_build_delimited_nest")(dat, opts_delimited(buildModel_similarity = sim))
    ref <- Nestimate::cooccurrence(dat["keywords"], field = "keywords", sep = ";",
                                   similarity = sim)
    expect_equal(net$weights, ref$weights, info = sim)
  }
  # raw counts: an item's diagonal is the number of rows that list it
  net <- cna("cna_build_delimited_nest")(dat, opts_delimited())
  n_network <- sum(grepl("\\bnetwork\\b", dat$keywords))
  expect_equal(unname(net$weights["network", "network"]), n_network)
})

test_that("a Transaction ID groups one item per row, even with a Separator set (review fix)", {
  # cooccurrence() ignores `by` when `sep` is given, and the Separator option
  # defaults to ";": the builder must not pass it alongside a Transaction ID.
  long <- data.frame(paper = c(1, 1, 1, 2, 2, 3, 3),
                     keyword = c("network", "graph", "matrix", "graph", "algebra",
                                 "network", "algebra"), stringsAsFactors = FALSE)
  net <- cna("cna_build_delimited_nest")(long, opts_delimited(
    buildModel_variables_field = "keyword", buildModel_variables_by = "paper",
    buildModel_sep = ";"))
  ref <- Nestimate::cooccurrence(long, field = "keyword", by = "paper")
  expect_equal(net$weights, ref$weights)
  # independent: network and graph share paper 1 only; graph and algebra paper 2
  expect_equal(unname(net$weights["network", "graph"]), 1)
  expect_equal(unname(net$weights["graph", "algebra"]), 1)
  expect_equal(unname(net$weights["matrix", "algebra"]), 0)
})

test_that("delimited top_n keeps only that many edges", {
  dat <- make_delimited_data()
  net <- cna("cna_build_delimited_nest")(dat, opts_delimited(buildModel_top_n = 3L))
  w <- net$weights
  expect_equal(sum(w[upper.tri(w)] != 0), 3)
})

test_that("clusters == Nestimate::cluster_network when unscaled", {
  seqs <- cna("cna_event_sequences")(make_event_data(n_actor = 12, per = 12), opts_event())
  o <- list(clustering_k = 2L, clustering_dissimilarity = "hamming",
            clustering_method = "pam", buildModel_scaling = "noScaling",
            seed_set = TRUE, seed_value = 11L)
  net <- cna("cna_cluster_networks")(seqs, o)
  ref <- Nestimate::cluster_network(seqs, k = 2, method = "co_occurrence",
                                    cluster_by = "pam", dissimilarity = "hamming", seed = 11)
  expect_s3_class(net, "netobject_group")
  expect_equal(names(net), names(ref))
  invisible(Map(function(a, b) expect_equal(unclass(a$weights), unclass(b$weights)), net, ref))
})

test_that("cluster scaling is honoured (regression, 1.1.0)", {
  # build_network(<net_clustering>) ignores `scaling`; CNA rebuilds each
  # cluster's network from its own sequences with the scaling applied.
  seqs <- cna("cna_event_sequences")(make_event_data(n_actor = 12, per = 12), opts_event())
  o <- list(clustering_k = 2L, clustering_dissimilarity = "hamming",
            clustering_method = "pam", seed_set = TRUE, seed_value = 11L)
  raw <- cna("cna_cluster_networks")(seqs, c(o, buildModel_scaling = "noScaling"))
  scaled <- cna("cna_cluster_networks")(seqs, c(o, buildModel_scaling = "minmax"))
  expect_equal(names(scaled), names(raw))
  expect_true(inherits(attr(scaled, "clustering"), "net_clustering"))
  invisible(Map(function(r, s) {
    expect_lte(max(s$weights), 1 + 1e-9)
    # minmax of the cluster's own raw weights, over the off-diagonal and diagonal alike
    expect_equal(cor(as.vector(r$weights), as.vector(s$weights)), 1, tolerance = 1e-9)
  }, raw, scaled))
})

test_that("unknown analysis kind raises a classed error", {
  expect_error(cna("cna_build_model")(data.frame(), list(), "nope"), class = "cna_bad_input")
})
