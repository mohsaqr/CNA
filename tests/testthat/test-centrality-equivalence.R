# Centrality: the table must equal cograph::centrality() with the same
# switches, every toggle must move the numbers, and the option <-> engine
# measure mapping must cover every checkbox the UI offers.

test_that("centrality == cograph::centrality for every offered measure", {
  # Three exceptions, each tested below: Eccentricity (cograph ignores
  # invert_weights for it, so CNA gives it the lengths 1/w), Diffusion (tna's,
  # from Nestimate, range-normalised to 0-1) and Clustering (weighted, Zhang &
  # Horvath on w / max(w)).
  net <- fixture_event_net()
  map <- cna(".cna_measure_map")
  map <- map[!names(map) %in% c("Eccentricity", "Diffusion", "Clustering")]
  out <- cna("cna_centrality_df")(net, names(map), loops = FALSE, normalized = FALSE,
                                  invert = TRUE)
  ref <- suppressMessages(cograph::centrality(net, measures = unname(map), loops = FALSE,
                                              invert_weights = TRUE))
  expect_setequal(out$state, ref$node)
  ord <- match(ref$node, out$state)
  invisible(Map(function(opt, m) {
    col <- intersect(c(m, paste0(m, "_all")), names(ref))[[1]]
    expect_equal(out[[opt]][ord], ref[[col]], info = opt)
  }, names(map), map))
})

test_that("Diffusion == Nestimate::net_centrality(normalize_diffusion = TRUE), 0-1, loops wired", {
  net <- fixture_group_net()
  invisible(lapply(c(FALSE, TRUE), function(loops) {
    out <- cna("cna_centrality_df")(net, c("Strength", "Diffusion"), loops = loops)
    invisible(lapply(names(net), function(g) {
      ref <- suppressMessages(Nestimate::net_centrality(net[[g]], measures = "Diffusion",
                                                         loops = loops, normalize_diffusion = TRUE))
      got <- subset(out, group == g)
      expect_equal(got$Diffusion[match(as.character(ref$state), got$state)], ref$Diffusion,
                   info = sprintf("%s, loops = %s", g, loops))
      expect_equal(range(got$Diffusion), c(0, 1))
    }))
  }))
  expect_false(isTRUE(all.equal(cna("cna_centrality_df")(net, "Diffusion", loops = TRUE)$Diffusion,
                                cna("cna_centrality_df")(net, "Diffusion", loops = FALSE)$Diffusion)))
  expect_error(cna(".cna_centralities")(unclass(net[[1]]$weights), "Diffusion"),
               class = "cna_engine_contract")
})

test_that("Clustering is weighted: tna's coefficient / (2 x max weight), in 0-1, not constant", {
  # tna (and Nestimate's built-in, and cograph's "onnela" type) apply Zhang &
  # Horvath to w + t(w) on raw counts; on w / max(w) the same formula gives
  # tna's value / (2 max w). The unweighted local transitivity is 1 everywhere
  # on this complete network.
  net <- fixture_group_net()
  out <- cna("cna_centrality_df")(net, c("Strength", "Clustering"))
  invisible(lapply(names(net), function(g) {
    w <- unclass(net[[g]]$weights)
    diag(w) <- 0
    tna_like <- suppressMessages(Nestimate::net_centrality(net[[g]], measures = "Clustering"))
    got <- subset(out, group == g)
    expect_equal(got$Clustering[match(as.character(tna_like$state), got$state)],
                 tna_like$Clustering / (2 * max(w)), info = g)
    expect_true(all(got$Clustering > 0 & got$Clustering <= 1), info = g)
    expect_gt(stats::sd(got$Clustering), 0)
  }))
  # loops and normalise leave it as it is; scaling the weights does not change it
  ref <- cna("cna_centrality_df")(net[[1]], "Clustering")$Clustering
  expect_equal(cna("cna_centrality_df")(net[[1]], "Clustering", loops = TRUE,
                                        normalized = TRUE)$Clustering, ref)
  big <- Nestimate::build_network(3 * unclass(net[[1]]$weights), method = "co_occurrence")
  expect_equal(cna("cna_centrality_df")(big, "Clustering")$Clustering, ref)
})

test_that("inverted Eccentricity == cograph's eccentricity on the lengths 1/w (cograph ignores invert)", {
  net <- fixture_group_net()
  out <- cna("cna_centrality_df")(net, "Eccentricity", invert = TRUE)
  invisible(lapply(names(net), function(g) {
    w <- unclass(net[[g]]$weights)
    ref <- suppressMessages(cograph::centrality(ifelse(w > 0, 1 / w, 0), measures = "eccentricity",
                                                loops = FALSE, invert_weights = FALSE,
                                                directed = FALSE))
    got <- subset(out, group == g)
    expect_equal(got$Eccentricity[match(ref$node, got$state)], ref$eccentricity_all %||% ref$eccentricity,
                 info = g)
  }))
  # without inversion the raw weights are the lengths, as cograph computes it
  raw <- cna("cna_centrality_df")(net[[1]], "Eccentricity", invert = FALSE)
  ref <- suppressMessages(cograph::centrality(net[[1]], measures = "eccentricity", loops = FALSE,
                                              invert_weights = FALSE))
  expect_equal(raw$Eccentricity[match(ref$node, raw$state)], ref$eccentricity_all %||% ref$eccentricity)
})

test_that("loops, normalize and invert are wired and move the numbers", {
  net <- fixture_event_net()
  f <- cna("cna_centrality_df")
  base <- f(net, c("Strength", "Betweenness"), loops = FALSE, normalized = FALSE, invert = TRUE)
  loops <- f(net, "Strength", loops = TRUE, normalized = FALSE, invert = TRUE)
  expect_false(isTRUE(all.equal(base$Strength, loops$Strength)))
  # with loops, strength gains twice the diagonal: an undirected self-loop
  # counts at both of its ends (the igraph convention cograph follows)
  w <- net$weights
  expect_equal(loops$Strength - base$Strength, 2 * unname(diag(w)[base$state]))
  nrm <- f(net, "Strength", loops = FALSE, normalized = TRUE, invert = TRUE)
  expect_true(all(nrm$Strength >= -1e-9 & nrm$Strength <= 1 + 1e-9))
  raw <- f(net, "Betweenness", loops = FALSE, normalized = FALSE, invert = FALSE)
  ref <- cograph::centrality(net, measures = "betweenness", loops = FALSE,
                             invert_weights = FALSE)
  expect_equal(raw$Betweenness[match(ref$node, raw$state)], ref$betweenness)
})

test_that("strength equals the off-diagonal row sums (independent check)", {
  net <- fixture_event_net()
  out <- cna("cna_centrality_df")(net, "Strength")
  w <- unclass(net$weights)
  diag(w) <- 0
  expect_equal(out$Strength, unname(rowSums(w)[out$state]))
})

test_that("group centralities equal per-network computations", {
  net <- fixture_group_net()
  out <- cna("cna_centrality_df")(net, "Strength")
  expect_setequal(unique(out$group), c("G1", "G2"))
  for (g in c("G1", "G2")) {
    ref <- cograph::centrality(net[[g]], measures = "strength", loops = FALSE)
    sub <- out[out$group == g, ]
    expect_equal(sub$Strength[match(ref$node, sub$state)], ref$strength_all)
  }
})

test_that("no measure selected is a classed error, not a silent default", {
  expect_error(cna("cna_centrality_df")(fixture_event_net(), character()),
               class = "cna_no_measures")
})

test_that("every centrality checkbox in every analysis maps to an engine measure", {
  map <- cna(".cna_measure_map")
  titles <- cna(".cna_measure_titles")
  expect_setequal(names(titles), names(map))
  for (a in ANALYSES) {
    opts <- grep("^centrality_[A-Z]", option_names(a), value = TRUE)
    expect_setequal(sub("^centrality_", "", opts), names(map))
  }
})

test_that("selected measures follow the checkboxes, in display order", {
  f <- cna("cna_selected_measures")
  expect_equal(f(list(centrality_Betweenness = TRUE, centrality_Strength = TRUE)),
               c("Strength", "Betweenness"))
  expect_length(f(list()), 0)
})

test_that("centrality table in the analysis == cograph (end to end)", {
  dat <- make_event_data()
  res <- run_event_single(dat, centrality_show_table = TRUE)
  tbl <- table_df(res, "centralityTable")
  net <- cna("cna_build_event_nest")(dat, opts_event())
  ref <- suppressMessages(cograph::centrality(net, measures = c("strength", "betweenness"),
                                              loops = FALSE, invert_weights = TRUE))
  ord <- match(ref$node, tbl$state)
  expect_equal(as.numeric(tbl$Strength[ord]), ref$strength_all)
  expect_equal(as.numeric(tbl$Betweenness[ord]), ref$betweenness)
})
