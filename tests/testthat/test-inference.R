# Resampling and comparison: bootstrap, permutation (with the corrected
# multiple-testing family), stability, network comparison,
# clustering diagnostics, and seed reproducibility.

# --- bootstrap ------------------------------------------------------------------

test_that("bootstrap table == bootstrap_network summary (seeded)", {
  net <- fixture_event_net()
  bs <- cna("cna_bootstrap_compute")(net, iter = 60L, level = 0.05, seed = 4L)
  ref <- Nestimate::bootstrap_network(net, iter = 60L, ci_level = 0.05, seed = 4L)
  out <- cna("cna_bootstrap_df")(bs, multi = FALSE, show_all = TRUE)
  key <- paste(out$from, out$to)
  rk <- paste(ref$summary$from, ref$summary$to)
  expect_setequal(key, rk)
  expect_equal(out$p_value, ref$summary$p_value[match(key, rk)])
  expect_equal(out$weight, unname(net$weights[cbind(out$from, out$to)]))
  # significant edges listed first, and the filter never adds rows
  sig <- cna("cna_bootstrap_df")(bs, multi = FALSE, significant_only = TRUE, show_all = TRUE)
  expect_true(all(sig$significant == "Yes"))
  expect_lte(nrow(sig), nrow(out))
  expect_false(is.unsorted(out$significant != "Yes"))
})

test_that("bootstrap CI method reaches the engine", {
  net <- fixture_event_net()
  a <- cna("cna_bootstrap_compute")(net, iter = 60L, level = 0.05, ci_method = "percentile", seed = 4L)
  b <- cna("cna_bootstrap_compute")(net, iter = 60L, level = 0.05, ci_method = "basic", seed = 4L)
  expect_false(isTRUE(all.equal(a$summary$ci_lower, b$summary$ci_lower)))
})

# --- permutation ------------------------------------------------------------------

test_that("permutation p-values are corrected over the unique edges (regression, 1.1.0)", {
  net <- fixture_group_net()
  raw <- cna("cna_permutation_compute")(net, iter = 99L, adjust = "none", seed = 5L)
  adj <- cna("cna_permutation_compute")(net, iter = 99L, adjust = "BH", seed = 5L)
  s_raw <- raw[[1]]$summary
  s_adj <- adj[[1]]$summary
  n <- nrow(net$G1$weights)
  # undirected: one test per unordered pair, diagonal included
  expect_equal(nrow(s_raw), n * (n + 1) / 2)
  # independent check: BH over exactly those tests
  expect_equal(s_adj$p_value, stats::p.adjust(s_raw$p_value, "BH"))
  # the plot reads p_values / diff_sig: they must agree with the table
  p <- adj[[1]]$p_values
  expect_equal(unname(p[cbind(s_adj$from, s_adj$to)]), s_adj$p_value)
  expect_true(isSymmetric(unname(p)))
  expect_equal(adj[[1]]$diff_sig, adj[[1]]$diff * (p < adj[[1]]$alpha))
})

test_that("permutation table: one row per edge per comparison, smallest p first", {
  net <- fixture_group_net()
  pt <- cna("cna_permutation_compute")(net, iter = 50L, seed = 5L)
  df <- cna("cna_permutation_df")(pt, show_all = TRUE)
  n <- nrow(net$G1$weights)
  expect_equal(nrow(df), n * (n + 1) / 2)
  expect_false(is.unsorted(df$p_value))
  expect_equal(nrow(cna("cna_permutation_df")(pt, max_rows = 5L)), 5L)
})

# --- stability -----------------------------------------------------------------

test_that("stability covers built-in and cograph-only measures", {
  net <- fixture_event_net(n_actor = 16)
  st <- cna("cna_stability_compute")(net, c("Strength", "Eigenvector", "PageRank"),
                                     iter = 50L, threshold = 0.7, certainty = 0.95, seed = 3L)
  df <- cna("cna_stability_df")(st, multi = FALSE)
  expect_equal(df$measure, c("Strength", "Eigenvector", "PageRank", "Edge weights"))
  expect_true(all(df$cs >= 0 & df$cs <= 0.9))
  expect_length(cna("cna_stability_failures")(st), 0)
  expect_length(cna("cna_stability_missing")(st, multi = FALSE), 0)
  rel <- cna("cna_reliability_df")(st)
  expect_true(all(c("metric", "mean", "sd") %in% names(rel)))
})

test_that("a measure constant across nodes is reported, not dropped silently", {
  net <- fixture_event_net(n_actor = 16)
  w <- net$weights
  testthat::skip_if(any(w[upper.tri(w)] == 0), "fixture is not a complete network")
  # complete network: every node has the same degree
  st <- cna("cna_stability_compute")(net, c("Strength", "Degree"), iter = 50L,
                                     threshold = 0.7, certainty = 0.95, seed = 3L)
  expect_equal(cna("cna_stability_missing")(st, multi = FALSE), "Degree")
  expect_false("Degree" %in% cna("cna_stability_df")(st, multi = FALSE)$measure)
})

test_that("stability is reproducible with a seed", {
  net <- fixture_event_net(n_actor = 16)
  f <- function() cna("cna_stability_df")(cna("cna_stability_compute")(
    net, c("Strength", "Betweenness"), iter = 50L, threshold = 0.7, certainty = 0.95,
    seed = 8L), multi = FALSE)
  expect_equal(f(), f())
})

test_that("a stability procedure the engine cannot run is reported, not hidden", {
  # one procedure failing must not hide the others: inject an engine failure
  net <- fixture_event_net(n_actor = 16)
  st <- cna("cna_stability_compute")(net, "Strength", iter = 50L, threshold = 0.7,
                                     certainty = 0.95, seed = 1L)
  st$nr <- simpleError("engine failed")
  failures <- cna("cna_stability_failures")(st, hint = "HINT")
  expect_identical(names(failures), "nr")
  expect_match(failures[["nr"]], "^Split-half reliability could not be computed by Nestimate \\(engine failed\\)\\. HINT$")
  expect_gt(nrow(cna("cna_stability_df")(st, multi = FALSE)), 0)
})

# --- comparing networks -----------------------------------------------------------

test_that("edge differences equal the weight differences (independent)", {
  net <- fixture_group_net()
  cmp <- cna("cna_compare_compute")(net, 1, 2)
  df <- cna("cna_compare_edges_df")(cmp, max_rows = 1000L)
  n <- nrow(net$G1$weights)
  expect_equal(nrow(df), n * (n + 1) / 2)
  expect_equal(df$weight_i, unname(net$G1$weights[cbind(df$from, df$to)]))
  expect_equal(df$weight_j, unname(net$G2$weights[cbind(df$from, df$to)]))
  expect_equal(df$difference, df$weight_i - df$weight_j)
  expect_false(is.unsorted(rev(abs(df$difference))))
})

test_that("compare pair resolution: groups by level, clusters by number", {
  net <- fixture_group_net()
  f <- cna(".cna_compare_pair")
  p <- f(net, list(compare_group_i = "G2", compare_group_j = "G1"), cluster = FALSE)
  expect_equal(c(p$i, p$j), c(2L, 1L))
  expect_error(f(net, list(compare_group_i = "G1", compare_group_j = "G1"), cluster = FALSE),
               class = "cna_bad_input")
  expect_error(f(net, list(compare_cluster_i = 1L, compare_cluster_j = 5L), cluster = TRUE),
               class = "cna_bad_input")
})

# --- clustering -----------------------------------------------------------------

test_that("cluster diagnostics and choosing k match Nestimate", {
  seqs <- cna("cna_event_sequences")(make_event_data(n_actor = 14, per = 12), opts_event())
  o <- list(clustering_k = 2L, clustering_dissimilarity = "hamming",
            clustering_method = "pam", buildModel_scaling = "noScaling",
            seed_set = TRUE, seed_value = 2L)
  net <- cna("cna_cluster_networks")(seqs, o)
  diag_df <- cna("cna_cluster_diagnostics_df")(net)
  expect_equal(sum(diag_df$size), nrow(seqs))
  expect_equal(sum(diag_df$pct), 100, tolerance = 1e-6)
  cc <- cna("cna_cluster_choice_compute")(seqs, 2L, 4L, seed = 2L)
  df <- cna("cna_cluster_choice_df")(cc)
  expect_equal(df$k, 2:4)
  ref <- Nestimate::cluster_choice(seqs, k = 2:4, dissimilarity = "hamming", method = "pam")
  expect_equal(df$silhouette, as.numeric(ref$silhouette))
  expect_error(cna("cna_cluster_choice_compute")(seqs, 5L, 3L), class = "cna_bad_input")
})

# --- resampling needs units ----------------------------------------------------

test_that("without units to resample, the resampling analyses explain and do not run", {
  # one actor-less sequence: Nestimate would give zero-width intervals, CS = 0
  # and p = 1 for every edge (checked on Nestimate 1.0.7)
  dat <- make_event_data()
  res <- run_event_single(dat, buildModel_variables_actor = NULL,
                          buildModel_variables_group = "Group",
                          compare_group_i = "G1", compare_group_j = "G2",
                          bootstrap_show_table = TRUE, bootstrap_show_plot = TRUE,
                          permutation_show_table = TRUE, stability_show_table = TRUE,
                          bootstrap_iteration = 20L, permutation_iter = 20L)
  for (tbl in c("bootstrapTable", "permutationTable", "stabilityTable")) {
    expect_equal(res[[tbl]]$rowCount, 0, label = tbl)
    expect_match(paste(unlist(lapply(res[[tbl]]$notes, `[[`, "note")), collapse = " "),
                 "needs an Actor", label = tbl)
  }
  expect_false(res$bootstrap_plot$visible)
  expect_null(res$bootstrapState$state)
  # the descriptive sections still run
  res2 <- run_event_single(dat, buildModel_variables_actor = NULL, centrality_show_table = TRUE)
  expect_gt(res2$centralityTable$rowCount, 0)
})

test_that("a permutation without significant edges draws a message, not a blank image", {
  dat <- make_event_data(n_actor = 12)
  res <- run_event_single(dat, buildModel_variables_group = "Group",
                          compare_group_i = "G1", compare_group_j = "G2",
                          permutation_show_plot = TRUE, permutation_iter = 10L,
                          permutation_level = 0.0001)
  pt <- res$permutationState$state
  expect_true(all(vapply(pt, function(cmp) all(cmp$diff_sig == 0), logical(1))))
  f <- tempfile(fileext = ".png")
  grDevices::png(f)
  drew <- cna("cna_render")(res$permutation_plot$analysis, "permutation_plot")
  grDevices::dev.off()
  expect_true(drew)
  expect_true(file.exists(f))
  texts <- svg_texts(cna("cna_render")(res$permutation_plot$analysis, "permutation_plot"))
  expect_true(any(grepl("No edge differs significantly", texts)))
})

# --- seeds ----------------------------------------------------------------------

test_that("Set seed makes the analysis reproducible; a different seed differs", {
  dat <- make_event_data()
  run <- function(seed) table_df(run_event_single(
    dat, bootstrap_show_table = TRUE, bootstrap_iteration = 40L,
    bootstrap_table_show_all = TRUE, seed_set = TRUE, seed_value = seed), "bootstrapTable")
  a <- run(10L)
  expect_equal(a, run(10L))
  expect_false(isTRUE(all.equal(a$ci_lower, run(11L)$ci_lower)))
})

test_that(".cna_with_seed restores the caller's RNG state", {
  set.seed(99)
  before <- .Random.seed
  cna(".cna_with_seed")(1L, stats::runif(3))
  expect_identical(.Random.seed, before)
})

test_that("scaled models do not resample (Nestimate rescales resamples differently)", {
  dat <- make_event_data()
  res <- run_event_single(dat, buildModel_scaling = "minmax", bootstrap_show_table = TRUE,
                          bootstrap_iteration = 20L, stability_show_table = TRUE)
  for (tbl in c("bootstrapTable", "stabilityTable")) {
    expect_equal(res[[tbl]]$rowCount, 0, label = tbl)
    expect_match(paste(unlist(lapply(res[[tbl]]$notes, `[[`, "note")), collapse = " "),
                 "Set Scaling to 'No scaling'", label = tbl)
  }
  expect_null(res$bootstrapState$state)
})

test_that("one-hot: stability, bootstrap and permutation run (by actor, windowed)", {
  dat <- make_onehot_data(n = 96)
  res <- run_cna("CNABinarySingle", dat, buildModel_variables_group = "Group",
                 compare_group_i = "G1", compare_group_j = "G2", buildModel_window = 2L,
                 stability_show_table = TRUE, bootstrap_show_table = TRUE,
                 permutation_show_table = TRUE, bootstrap_iteration = 20L,
                 permutation_iter = 20L, stability_iter = 50L, seed_set = TRUE)
  expect_length(section_failures(res), 0)
  expect_gt(res$stabilityTable$rowCount, 0)
  expect_gt(res$reliabilityTable$rowCount, 0)
  expect_gt(res$bootstrapTable$rowCount, 0)
  expect_gt(res$permutationTable$rowCount, 0)
})

test_that("one-hot stability == Nestimate::centrality_stability on the same network (seeded)", {
  dat <- make_onehot_data(n = 96)
  net <- cna("cna_build_onehot_nest")(dat, opts_onehot(buildModel_window = 2L))
  st <- cna("cna_stability_compute")(net, c("Strength", "Betweenness"), iter = 50L,
                                     threshold = 0.7, certainty = 0.95, seed = 5L)
  expect_length(cna("cna_stability_failures")(st), 0)
  ref <- Nestimate::centrality_stability(net, measures = c("InStrength", "Betweenness"),
                                         iter = 50L, threshold = 0.7, certainty = 0.95,
                                         seed = 5L)
  expect_equal(unname(st$cs$cs), unname(ref$cs))
})
