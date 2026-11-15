# Mixed network (transition + co-occurrence). Events build the transition
# network from the same sequences as each co-occurrence network; one-hot data
# use wtna(method = "both"). Both displays must draw readable, labelled panels.

test_that("event mixed: transitions == build_network(method = 'frequency') on the same data", {
  dat <- make_event_data()
  o <- opts_event(buildModel_variables_group = "Group")
  model <- cna("cna_build_event_nest")(dat, o)
  pairs <- cna("cna_build_mixed")(dat, o, "event_single", model)
  expect_setequal(names(pairs), c("G1", "G2"))
  ref <- Nestimate::build_network(dat[c("Time", "Actor", "Action", "Group")],
                                  method = "frequency", action = "Action", actor = "Actor",
                                  time = "Time", time_threshold = 3000, group = "Group")
  for (g in c("G1", "G2")) {
    expect_true(pairs[[g]]$transition$directed)
    expect_equal(unclass(pairs[[g]]$transition$weights), unclass(ref[[g]]$weights))
    expect_identical(pairs[[g]]$cooccurrence, model[[g]])
  }
})

test_that("event mixed: a single network is keyed '' and clusters get one pair each", {
  dat <- make_event_data(n_actor = 12, per = 12)
  single <- cna("cna_build_mixed")(dat, opts_event(), "event_single",
                                   cna("cna_build_event_nest")(dat, opts_event()))
  expect_identical(names(single), "")
  seqs <- cna("cna_event_sequences")(dat, opts_event())
  o <- list(clustering_k = 2L, clustering_dissimilarity = "hamming", clustering_method = "pam",
            buildModel_scaling = "noScaling", seed_set = TRUE, seed_value = 4L)
  clusters <- cna("cna_cluster_networks")(seqs, o)
  pairs <- cna("cna_build_mixed")(dat, o, "event_cluster", clusters)
  expect_equal(names(pairs), names(clusters))
  invisible(Map(function(p, net) {
    ref <- Nestimate::build_network(net$data, method = "frequency")
    expect_equal(unclass(p$transition$weights), unclass(ref$weights))
    # every transition is between states of that cluster's co-occurrence network
    expect_true(all(rownames(p$transition$weights) %in% rownames(net$weights)))
  }, pairs, clusters))
})

test_that("one-hot mixed: the analysis's network beside build_network(method = 'frequency')", {
  dat <- make_onehot_data()
  codes <- paste0("code", 1:4)
  o <- opts_onehot(buildModel_window = 2L, buildModel_window_type = "sliding")
  model <- cna("cna_build_onehot_nest")(dat, o)
  pairs <- cna("cna_build_mixed")(dat, o, "binary_single", model)
  expect_identical(names(pairs), "")
  expect_identical(pairs[[1]]$cooccurrence, model)
  ref <- Nestimate::build_network(dat[c(codes, "Actor", "Session")], method = "frequency",
                                  codes = codes, format = "onehot", actor = "Actor",
                                  session = "Session", window_size = 2L, mode = "overlapping")
  expect_equal(pairs[[1]]$transition$weights, ref$weights)
  # and that is wtna's transition count within each actor x session
  w <- Nestimate::wtna(dat[c(codes, "Actor", "Session")], method = "transition", codes = codes,
                       window_size = 2L, mode = "overlapping", actor = c("Actor", "Session"))
  expect_equal(unclass(pairs[[1]]$transition$weights), unclass(w$weights))
})

test_that("one-hot transitions never cross from one actor into the next (fix)", {
  # within each actor only A -> A and B -> B; the boundary a1 (B) -> a2 (A)
  # would be the only B -> A transition if rows were treated as one sequence
  dat <- data.frame(A = c(1, 1, 0, 0, 1, 1), B = c(0, 0, 1, 1, 0, 0),
                    Actor = c("a1", "a1", "a1", "a1", "a2", "a2"))
  o <- list(buildModel_variables_onehot = c("A", "B"), buildModel_variables_actor = "Actor",
            buildModel_variables_session = NULL, buildModel_variables_group = NULL,
            buildModel_window = 1L, buildModel_window_type = "tumbling")
  model <- cna("cna_build_onehot_nest")(dat, o)
  tr <- cna("cna_build_mixed")(dat, o, "binary_single", model)[[1]]$transition$weights
  expect_equal(unname(tr["B", "A"]), 0)
  expect_equal(unname(tr["A", "B"]), 1)
})

test_that("the mixed view is offered where it can be built, and nowhere else", {
  with_mixed <- c("CNAEventSingle", "CNAEventCluster", "CNABinarySingle", "CNABinaryCluster")
  for (a in ANALYSES)
    expect_identical("buildModel_show_mixed" %in% option_names(a), a %in% with_mixed,
                     label = a)
})

test_that("both displays draw labelled panels through the analysis", {
  dat <- make_event_data()
  for (type in c("group", "overlay")) {
    res <- run_event_single(dat, buildModel_variables_group = "Group",
                            compare_group_i = "G1", compare_group_j = "G2",
                            buildModel_show_mixed = TRUE, buildModel_mixed_type = type)
    img <- res$buildModel_mixed
    expect_true(img$visible)
    texts <- svg_texts(cna("cna_render")(img$analysis, "buildModel_mixed"))
    if (type == "group") {
      expect_true(all(c("G1 Transition", "G1 Co-occurrence", "G2 Transition",
                        "G2 Co-occurrence") %in% texts))
      expect_equal(img$height, 2L * 420L)       # one row per group
    } else {
      expect_true(all(c("G1", "G2") %in% texts))
    }
  }
})

test_that("one-hot transition panels are drawn directed, with arrows (cograph styling fix)", {
  testthat::skip_if_not_installed("svglite")
  dat <- make_onehot_data(n = 80)
  res <- run_cna("CNABinarySingle", dat, buildModel_show_mixed = TRUE,
                 buildModel_mixed_type = "group")
  pairs <- res$buildModel_mixedState$state
  expect_true(pairs[[1]]$transition$directed)
  f <- tempfile(fileext = ".svg")
  svglite::svglite(f, 9, 5)
  cna("cna_render")(res$buildModel_mixed$analysis, "buildModel_mixed")
  grDevices::dev.off()
  svg <- readLines(f, warn = FALSE)
  unlink(f)
  # arrowheads are drawn as polygons; the co-occurrence panel has none
  tr <- unclass(pairs[[1]]$transition$weights)
  n_directed_edges <- sum(tr[row(tr) != col(tr)] != 0)
  expect_gte(sum(grepl("<polygon", svg)), n_directed_edges)
})
