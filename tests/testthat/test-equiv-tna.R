skip_if_not_installed("tna")

# Cross-package equivalence with tna (Saqr, Tikka & Lopez-Pernas), whose
# definitions CNA says it follows: one-hot windows (import_onehot), the
# co-occurrence models (ctna / group_ctna), resampling on raw counts, and the
# case-dropping CS-coefficient (estimate_cs). Each block states whether the two
# packages agree exactly, up to a stated convention, or differ (and who is at
# fault). Written against tna 1.3.1; where tna >= 1.5.4 changed a definition the
# block checks the version (tna's NEWS for 1.5.4 describes each change).

tna_fixed_sliding <- utils::packageVersion("tna") >= "1.5.4"

srl_codes <- c("Plan", "Monitor", "Evaluate", "Reflect", "Read", "Summarize",
               "Practice", "Review")
fixture_codes <- paste0("code", seq_len(4L))
srl_onehot <- function() utils::read.csv(file.path(repo_root(), "data", "SRL_onehot.csv"),
                                         stringsAsFactors = FALSE)
srl_events <- function() utils::read.csv(file.path(repo_root(), "data", "SRL_events.csv"),
                                         stringsAsFactors = FALSE)

# A weight matrix as a plain double matrix, optionally reordered to `nodes`.
as_num <- function(w, nodes = rownames(w)) {
  w <- unclass(w)[nodes, nodes, drop = FALSE]
  matrix(as.numeric(w), nrow(w), dimnames = list(nodes, nodes))
}

# Nestimate warns that one-hot data with window 1 and no session give sparse
# networks, also when an Actor is set (docs/ONEHOT_AUDIT.md issue 16). Only that
# diagnosed warning is muffled; any other warning passes.
quiet_sparse <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (grepl("sparse networks", conditionMessage(w), fixed = TRUE))
    invokeRestart("muffleWarning")
})

# tna's one-hot import; import_onehot() tests missing(actor) / missing(session),
# so unused roles are left out of the call.
tna_onehot <- function(data, codes, window = 1L, type = "tumbling",
                       actor = NULL, session = NULL) {
  roles <- Filter(Negate(is.null), list(actor = actor, session = session))
  do.call(tna::import_onehot, c(list(data = data, cols = codes, window_size = window,
                                     window_type = type), roles))
}
tna_onehot_weights <- function(...) tna::ctna(tna_onehot(...))$weights

cna_onehot <- function(data, codes, window = 1L, type = "tumbling", actor = NULL,
                       session = NULL, group = NULL, scaling = "noScaling") {
  quiet_sparse(cna("cna_build_onehot_nest")(data, list(
    buildModel_variables_onehot = codes, buildModel_variables_actor = actor,
    buildModel_variables_session = session, buildModel_variables_group = group,
    buildModel_window = window, buildModel_window_type = type,
    buildModel_scaling = scaling)))
}

# tna's session split; prepare_data() reports progress with cli messages.
tna_sessions <- function(data, threshold, order = NULL) {
  args <- c(list(data = data, actor = "Actor", time = "Time", action = "Action",
                 time_threshold = threshold), if (!is.null(order)) list(order = order))
  suppressMessages(do.call(tna::prepare_data, args))
}
cna_event <- function(data, threshold, order = NULL, group = NULL, scaling = "noScaling")
  cna("cna_build_event_nest")(data, opts_event(
    buildModel_threshold = threshold, buildModel_variables_long_order = order,
    buildModel_variables_group = group, buildModel_scaling = scaling))

# One string per sequence (states joined in order, missing cells dropped).
seq_keys <- function(s) apply(as.matrix(s), 1L, \(r) paste(r[!is.na(r)], collapse = "-"))

# Independent one-hot window counts per unit: rows in blocks of `window`, a
# block adding tcrossprod(colSums(block)). Sliding first replaces the rows by one
# row per window marking the codes active in it: complete windows of `window`
# rows (tna >= 1.5.4, Nestimate a1754ab), or tna 1.3.1's widened windows, where
# row t >= 2 covers rows t - L, ..., t with L = max(1, window (window - 1) / 2)
# because each lag is applied to the already widened rows.
sliding_rows <- function(mat, window, widened = FALSE) {
  n <- nrow(mat)
  spans <- if (widened) {
    lag_span <- max(1L, (window * (window - 1L)) %/% 2L)
    lapply(seq_len(max(0L, n - 1L)) + 1L, \(t) max(1L, t - lag_span):t)
  } else if (n >= window) {
    lapply(seq_len(n - window + 1L), \(s) s:(s + window - 1L))
  } else list()
  active <- vapply(spans, \(r) as.numeric(colSums(mat[r, , drop = FALSE]) > 0),
                   numeric(ncol(mat)))
  matrix(active, ncol = ncol(mat), byrow = TRUE)
}
onehot_window_ref <- function(df, codes, unit, window, type = "tumbling", widened = FALSE) {
  per_unit <- function(mat) {
    if (identical(type, "sliding")) mat <- sliding_rows(mat, window, widened)
    n <- nrow(mat)
    blocks <- if (n) split(seq_len(n), (seq_len(n) - 1L) %/% window) else list()
    Reduce(`+`, lapply(blocks, \(b) tcrossprod(colSums(mat[b, , drop = FALSE]))),
           matrix(0, length(codes), length(codes)))
  }
  parts <- split(df, interaction(df[unit], drop = TRUE))
  out <- Reduce(`+`, lapply(parts, \(p) per_unit(as.matrix(p[codes]))))
  dimnames(out) <- list(codes, codes)
  out
}

# Group of each actor x session unit in first-seen order: the row order of
# import_onehot()'s output, which group_ctna() needs as a vector.
unit_groups <- function(df, units, group) {
  key <- do.call(paste, c(unname(as.list(df[units])), sep = "\r"))
  df[[group]][!duplicated(key)]
}
window_params <- function(oh) list(windowed = TRUE, window_size = attr(oh, "window_size"),
                                   window_span = attr(oh, "window_span"))

same_partition <- function(a, b) identical(outer(a, a, `==`), outer(b, b, `==`))

# tna's estimate_cs(), permutation_test() and communities() take no seed: run
# them under a seed and restore the caller's RNG state afterwards (CNA's own
# helper, also used for engine functions without a seed argument).
with_seed <- cna(".cna_with_seed")

# ---------------------------------------------------------------------------
# 1. One-hot networks
# ---------------------------------------------------------------------------

test_that("one-hot tumbling windows 1-4 and sliding window 2 equal tna::ctna(import_onehot())", {
  skip_on_cran()
  units <- list(actor = list(actor = "Actor"), session = list(session = "Session"),
                both = list(actor = "Actor", session = "Session"))
  grid <- rbind(
    expand.grid(window = seq_len(4L), type = "tumbling", unit = names(units), seed = c(7L, 11L, 23L),
                stringsAsFactors = FALSE),
    expand.grid(window = 2L, type = "sliding", unit = names(units), seed = c(7L, 11L, 23L),
                stringsAsFactors = FALSE))
  invisible(Map(function(window, type, unit, seed) {
    d <- make_onehot_data(n = 48, seed = seed)
    args <- c(list(data = d, codes = fixture_codes, window = window, type = type), units[[unit]])
    tw <- do.call(tna_onehot_weights, args)
    cw <- do.call(cna_onehot, args)$weights
    expect_equal(as_num(tw), as_num(cw, rownames(tw)),
                 info = sprintf("window %d %s, %s, seed %d", window, type, unit, seed))
  }, grid$window, grid$type, grid$unit, grid$seed))
})

test_that("SRL_onehot: tumbling 1-4 and sliding 2 equal tna, by actor, by session and by both", {
  skip_on_cran()
  srl <- srl_onehot()
  # Course is constant within actor: as Session alone it pools actors into two
  # units, as Session with Actor it splits nothing further.
  units <- list(actor = list(actor = "Actor"), session = list(session = "Course"),
                both = list(actor = "Actor", session = "Course"))
  grid <- rbind(expand.grid(window = seq_len(4L), type = "tumbling", unit = names(units),
                            stringsAsFactors = FALSE),
                expand.grid(window = 2L, type = "sliding", unit = names(units),
                            stringsAsFactors = FALSE))
  invisible(Map(function(window, type, unit) {
    args <- c(list(data = srl, codes = srl_codes, window = window, type = type), units[[unit]])
    tw <- do.call(tna_onehot_weights, args)
    cw <- do.call(cna_onehot, args)$weights
    expect_equal(as_num(tw), as_num(cw, rownames(tw)),
                 info = sprintf("window %d %s, %s", window, type, unit))
    # both follow the per-unit window definition
    expect_equal(as_num(cw, srl_codes),
                 onehot_window_ref(srl, srl_codes, unlist(units[[unit]]), window, type),
                 info = sprintf("reference, window %d %s, %s", window, type, unit))
  }, grid$window, grid$type, grid$unit))
})

test_that("sliding windows 1, 3, 4: CNA counts complete windows; tna 1.3.1 widens them (tna defect)", {
  # Hand-checked: actor a1 has rows A, AB, BC, C, AC; a2 has AB, BC.
  tiny <- data.frame(A = c(1, 1, 0, 0, 1, 1, 0), B = c(0, 1, 1, 0, 0, 1, 1),
                     C = c(0, 0, 1, 1, 1, 0, 1), Actor = c(rep("a1", 5), rep("a2", 2)))
  abc <- c("A", "B", "C")
  # CNA, window 3: a1's three complete windows are all ABC (one block of 3 rows,
  # 3 x 3 = 9 everywhere); a2 is shorter than the window and adds nothing.
  cw <- cna_onehot(tiny, abc, window = 3L, type = "sliding", actor = "Actor")$weights
  expect_equal(as_num(cw), matrix(9, 3, 3, dimnames = list(abc, abc)))
  tw <- tna_onehot_weights(tiny, abc, window = 3L, type = "sliding", actor = "Actor")
  if (tna_fixed_sliding) {
    expect_equal(as_num(tw), as_num(cw, rownames(tw)))
  } else {
    # tna 1.3.1: a1 rows 2-5 cover rows 1-2, 1-3, 1-4, 2-5 (AB, ABC, ABC, ABC),
    # blocks [AB, ABC, ABC] + [ABC]; a2 keeps its 2-row window ABC.
    expect_equal(as_num(tw), matrix(c(11, 11, 8, 11, 11, 8, 8, 8, 6), 3,
                                    dimnames = list(abc, abc)))
  }

  skip_on_cran()
  cases <- list(fixture = list(data = make_onehot_data(n = 48, seed = 7), codes = fixture_codes,
                               unit = c("Actor", "Session")),
                srl = list(data = srl_onehot(), codes = srl_codes, unit = "Actor"))
  grid <- expand.grid(window = c(1L, 3L, 4L), case = names(cases), stringsAsFactors = FALSE)
  invisible(Map(function(window, case) {
    cs <- cases[[case]]
    roles <- list(actor = "Actor", session = if (length(cs$unit) > 1L) "Session")
    cw <- do.call(cna_onehot, c(list(data = cs$data, codes = cs$codes, window = window,
                                     type = "sliding"), roles))$weights
    tw <- do.call(tna_onehot_weights, c(list(data = cs$data, codes = cs$codes, window = window,
                                             type = "sliding"), Filter(Negate(is.null), roles)))
    info <- sprintf("%s, sliding window %d", case, window)
    expect_equal(as_num(cw, cs$codes),
                 onehot_window_ref(cs$data, cs$codes, cs$unit, window, "sliding"), info = info)
    if (tna_fixed_sliding) {
      expect_equal(as_num(tw), as_num(cw, rownames(tw)), info = info)
    } else {
      widened <- onehot_window_ref(cs$data, cs$codes, cs$unit, window, "sliding", widened = TRUE)
      expect_equal(as_num(tw), widened[rownames(tw), rownames(tw)], info = info)
      expect_false(isTRUE(all.equal(as_num(tw), as_num(cw, rownames(tw)))), info = info)
    }
  }, grid$window, grid$case))
})

test_that("group_ctna() equals CNA's group networks only when given the window parameters (tna defect)", {
  d <- make_onehot_data(n = 96, seed = 3)
  units <- c("Actor", "Session")
  g <- unit_groups(d, units, "Group")
  configs <- data.frame(window = c(1L, 2L, 3L, 2L),
                        type = c("tumbling", "tumbling", "tumbling", "sliding"))
  invisible(Map(function(window, type) {
    oh <- tna_onehot(d, fixture_codes, window, type, actor = "Actor", session = "Session")
    gm <- tna::group_ctna(oh, group = g, params = window_params(oh))
    nets <- cna_onehot(d, fixture_codes, window, type, actor = "Actor", session = "Session",
                       group = "Group")
    expect_setequal(names(gm), names(nets))
    invisible(lapply(names(nets), \(k)
      expect_equal(as_num(gm[[k]]$weights), as_num(nets[[k]]$weights, rownames(gm[[k]]$weights)),
                   info = sprintf("%s, window %d %s", k, window, type))))
  }, configs$window, configs$type))

  # SRL: groups are fixed per actor
  srl <- srl_onehot()
  oh <- tna_onehot(srl, srl_codes, 2L, actor = "Actor")
  gm <- tna::group_ctna(oh, group = unit_groups(srl, "Actor", "Group"), params = window_params(oh))
  nets <- cna_onehot(srl, srl_codes, 2L, actor = "Actor", group = "Group")
  invisible(lapply(names(nets), \(k)
    expect_equal(as_num(gm[[k]]$weights), as_num(nets[[k]]$weights, rownames(gm[[k]]$weights)),
                 info = k)))

  # Without `params`, group_model() does not read import_onehot()'s window
  # attributes (build_model() does) and counts every pair of positions of the
  # whole actor x session row: off-diagonal n_c * n_d, diagonal choose(n_c, 2),
  # with n_c the rows in which code c is active.
  skip_if(utils::packageVersion("tna") > "1.5.5", "group_model() may read the window attributes")
  oh <- tna_onehot(d, fixture_codes, 2L, actor = "Actor", session = "Session")
  gm0 <- tna::group_ctna(oh, group = g)
  nets2 <- cna_onehot(d, fixture_codes, 2L, actor = "Actor", session = "Session", group = "Group")
  unwindowed <- function(rows) {
    parts <- split(rows, interaction(rows[units], drop = TRUE))
    Reduce(`+`, lapply(parts, function(p) {
      n_c <- colSums(as.matrix(p[fixture_codes]))
      m <- tcrossprod(n_c)
      diag(m) <- choose(n_c, 2)
      m
    }))
  }
  invisible(lapply(names(gm0), function(k) {
    ref <- unwindowed(d[d$Group == k, , drop = FALSE])
    dimnames(ref) <- list(fixture_codes, fixture_codes)
    nodes <- rownames(gm0[[k]]$weights)
    expect_equal(as_num(gm0[[k]]$weights), ref[nodes, nodes], info = k)
    expect_false(isTRUE(all.equal(as_num(gm0[[k]]$weights), as_num(nets2[[k]]$weights, nodes))),
                 info = k)
  }))
})

# ---------------------------------------------------------------------------
# 2. Sequence co-occurrence and sessions
# ---------------------------------------------------------------------------

test_that("sequence co-occurrence: Nestimate::build_network(co_occurrence) equals tna::ctna on wide sequences", {
  # Hand-checked: "A B A" has position pairs (1,2) A-B, (1,3) A-A, (2,3) B-A.
  one <- data.frame(T1 = "A", T2 = "B", T3 = "A", stringsAsFactors = FALSE)
  expect_equal(as_num(tna::ctna(one)$weights), matrix(c(1, 2, 2, 0), 2, dimnames = list(c("A", "B"), c("A", "B"))))
  expect_equal(as_num(Nestimate::build_network(one, method = "co_occurrence")$weights, c("A", "B")),
               matrix(c(1, 2, 2, 0), 2, dimnames = list(c("A", "B"), c("A", "B"))))
  invisible(lapply(c(1L, 2L, 3L), function(seed) {
    wide <- with_seed(seed, as.data.frame(t(replicate(15, {
      n <- sample(3:8, 1)
      c(sample(LETTERS[seq_len(5L)], n, TRUE), rep(NA, 8 - n))
    })), stringsAsFactors = FALSE))
    names(wide) <- paste0("T", seq_len(8L))
    tw <- tna::ctna(wide)$weights
    nw <- Nestimate::build_network(wide, method = "co_occurrence")$weights
    expect_equal(as_num(tw), as_num(nw, rownames(tw)), info = sprintf("seed %d", seed))
  }))
})

test_that("event sessions: CNA's sequences are tna::prepare_data()'s sessions, in the same order", {
  d <- make_event_data()
  invisible(lapply(c(60, 300, 900, 3000), function(thr) {
    cs <- cna("cna_event_sequences")(d, opts_event(buildModel_threshold = thr))
    pd <- tna_sessions(d, thr)
    expect_identical(seq_keys(cs), seq_keys(pd$sequence_data), info = sprintf("threshold %g", thr))
    tw <- tna::ctna(pd)$weights
    expect_equal(as_num(tw), as_num(cna_event(d, thr)$weights, rownames(tw)),
                 info = sprintf("network, threshold %g", thr))
  }))
  # SRL_events with an Order column for ties
  ev <- srl_events()
  invisible(lapply(c(300, 900), function(thr) {
    cs <- cna("cna_event_sequences")(ev, opts_event(buildModel_threshold = thr,
                                                    buildModel_variables_long_order = "Order"))
    pd <- tna_sessions(ev, thr, order = "Order")
    expect_identical(seq_keys(cs), seq_keys(pd$sequence_data), info = sprintf("SRL, threshold %g", thr))
    tw <- tna::ctna(pd)$weights
    expect_equal(as_num(tw), as_num(cna_event(ev, thr, order = "Order")$weights, rownames(tw)))
  }))
  # A gap equal to the threshold stays in the session in both (new session only
  # when the gap exceeds it); ISO 8601 strings parse to the same sessions.
  b <- data.frame(Actor = "x", Action = c("A", "B", "C", "A"),
                  Time = format(as.POSIXct("2021-01-01", tz = "UTC") + c(0, 100, 200, 301),
                                "%Y-%m-%dT%H:%M:%SZ"), stringsAsFactors = FALSE)
  cs <- cna("cna_event_sequences")(b, opts_event(buildModel_threshold = 100))
  expect_identical(seq_keys(cs), c("A-B-C", "A"))
  expect_identical(seq_keys(tna_sessions(b, 100)$sequence_data), c("A-B-C", "A"))
})

test_that("event group networks equal tna::group_ctna(<prepare_data>, group = 'Group')", {
  d <- make_event_data(n_actor = 30)
  pd <- tna_sessions(d, 3000)
  gm <- tna::group_ctna(pd, group = "Group")
  nets <- cna_event(d, 3000, group = "Group")
  expect_setequal(names(gm), names(nets))
  invisible(lapply(names(nets), \(k)
    expect_equal(as_num(gm[[k]]$weights), as_num(nets[[k]]$weights, rownames(gm[[k]]$weights)),
                 info = k)))
})

# ---------------------------------------------------------------------------
# 3. Centralities
# ---------------------------------------------------------------------------

test_that("centralities: strength, closeness, PageRank equal tna's; betweenness and RSP up to convention", {
  skip_on_cran()
  srl <- srl_onehot()
  ev <- srl_events()
  d <- make_event_data()
  cases <- list(
    onehot = list(tna = tna::ctna(tna_onehot(srl, srl_codes, 2L, actor = "Actor")),
                  cna = cna_onehot(srl, srl_codes, 2L, actor = "Actor")),
    srl_events = list(tna = tna::ctna(tna_sessions(ev, 900, order = "Order")),
                      cna = cna_event(ev, 900, order = "Order")),
    fixture = list(tna = tna::ctna(tna_sessions(d, 300)), cna = cna_event(d, 300)))
  invisible(Map(function(cs, nm) {
    ct <- as.data.frame(tna::centralities(cs$tna, measures = c(
      "InStrength", "OutStrength", "Closeness", "Betweenness", "PageRank", "BetweennessRSP")))
    states <- as.character(ct$state)
    cc <- cna("cna_centrality_df")(cs$cna, c("Strength", "Closeness", "Betweenness",
                                             "PageRank", "BetweennessRSP"))
    cc <- cc[match(states, cc$state), ]
    expect_equal(ct$OutStrength, ct$InStrength, info = nm)   # symmetric weights
    expect_equal(cc$Strength, ct$InStrength, info = nm)
    expect_equal(cc$Closeness, ct$Closeness, info = nm)
    expect_equal(cc$PageRank, ct$PageRank, info = nm)
    # tna builds a directed graph from the symmetric matrix, so every shortest
    # path is counted for (s, t) and (t, s): twice CNA's undirected betweenness.
    expect_equal(ct$Betweenness, 2 * cc$Betweenness, info = nm)
    # tna's RSP betweenness is rounded and shifted to a minimum of 1.
    rsp <- round(cc$BetweennessRSP)
    expect_equal(ct$BetweennessRSP, rsp - min(rsp) + 1, info = nm)

    # Loops: tna's directed graph counts w_ii once in the strength, CNA's
    # undirected graph twice (the igraph convention for a self-loop).
    w <- as_num(cs$cna$weights, states)
    ctl <- as.data.frame(tna::centralities(cs$tna, loops = TRUE, measures = "InStrength"))
    ccl <- cna("cna_centrality_df")(cs$cna, "Strength", loops = TRUE)
    expect_equal(ccl$Strength[match(states, ccl$state)], ctl$InStrength + unname(diag(w)), info = nm)

    # Normalisation: tna rescales to [0, 1] (min-max), cograph by a positive
    # factor (strength and PageRank / max, closeness x (n - 1)); min-max of
    # CNA's normalised values therefore equals tna's.
    minmax <- function(x) (x - min(x)) / (max(x) - min(x))
    ctn <- as.data.frame(tna::centralities(cs$tna, normalize = TRUE,
                                           measures = c("InStrength", "Closeness", "PageRank")))
    ccn <- cna("cna_centrality_df")(cs$cna, c("Strength", "Closeness", "PageRank"), normalized = TRUE)
    ccn <- ccn[match(states, ccn$state), ]
    expect_equal(minmax(ccn$Strength), ctn$InStrength, info = nm)
    expect_equal(minmax(ccn$Closeness), ctn$Closeness, info = nm)
    expect_equal(minmax(ccn$PageRank), ctn$PageRank, info = nm)
  }, cases, names(cases)))
})

test_that("Nestimate's built-in Diffusion, Clustering and RSP use tna's definitions", {
  # CNA's Diffusion (table and stability) is this one, range-normalised to 0-1
  # (next block); its Clustering and RSP are cograph's.
  net <- fixture_event_net()
  tm <- tna::ctna(tna_sessions(make_event_data(), 3000))
  ct <- as.data.frame(tna::centralities(tm, measures = c("Diffusion", "Clustering", "BetweennessRSP")))
  states <- as.character(ct$state)
  # what centrality_stability() correlates (its built-in measures)
  nc <- Nestimate::net_centrality(net, measures = c("Diffusion", "Clustering", "BetweennessRSP"),
                                  loops = FALSE, normalize_diffusion = FALSE)
  nc <- nc[match(states, as.character(nc$state)), ]
  expect_equal(nc$Diffusion, ct$Diffusion)
  expect_equal(nc$Clustering, ct$Clustering)
  expect_equal(nc$BetweennessRSP, ct$BetweennessRSP)
})

test_that("Diffusion and Clustering are tna's, on 0-1, and are what stability assesses", {
  # 2026-11-06: the table showed cograph's binary diffusion and unweighted
  # transitivity (both constant on a complete network) while stability
  # assessed tna's. Now both are tna's, on 0-1: Diffusion range-normalised
  # (tna's power series runs into the billions on counts), Clustering as
  # Zhang & Horvath on w / max(w) (tna's value / (2 max w); it grows with the
  # counts). 30 actors: enough sessions that every drop proportion drops a
  # case (tna skips a proportion that drops none). Two measures: tna 1.3.1
  # labels the states 1..n instead of by name when asked for one.
  d <- make_event_data(n_actor = 30)
  net <- cna_event(d, 3000)
  tm <- tna::ctna(tna_sessions(d, 3000))
  ct <- as.data.frame(tna::centralities(tm, measures = c("Diffusion", "Clustering")))
  range01 <- function(x) (x - min(x)) / (max(x) - min(x))
  w <- unclass(net$weights)
  diag(w) <- 0
  cc <- cna("cna_centrality_df")(net, c("Diffusion", "Clustering"))
  cc <- cc[match(as.character(ct$state), cc$state), ]
  expect_equal(cc$Diffusion, range01(ct$Diffusion))
  expect_equal(cc$Clustering, ct$Clustering / (2 * max(w)))
  expect_true(all(cc$Diffusion >= 0 & cc$Diffusion <= 1))
  expect_true(all(cc$Clustering > 0 & cc$Clustering <= 1))
  # the same draws as tna::estimate_cs() on tna's raw values: rescaling a
  # network's values by a constant does not change a correlation
  st <- cna("cna_stability_compute")(net, c("Diffusion", "Clustering"), iter = 50L,
                                     threshold = 0.7, certainty = 0.95, seed = 1L)
  ts <- with_seed(1L, tna::estimate_cs(tm, measures = c("Diffusion", "Clustering"), iter = 50L,
                                       threshold = 0.7, certainty = 0.95))
  invisible(lapply(c("Diffusion", "Clustering"), function(m) {
    expect_equal(unname(st$cs$correlations[[m]]), unname(ts[[m]]$correlations), info = m)
    expect_equal(unname(st$cs$cs[[m]]), ts[[m]]$cs_coefficient, info = m)
  }))
})

# tna's measure names -> the titles CNA's stability result carries.
tna_to_cna <- c(InStrength = "Strength", Closeness = "Closeness", Betweenness = "Betweenness")

# ---------------------------------------------------------------------------
# 4. Stability (case-dropping CS-coefficient)
# ---------------------------------------------------------------------------

test_that("stability: with the same seed CNA draws tna::estimate_cs()'s subsets and gets its CS-coefficients", {
  skip_on_cran()
  # Events: cases are sessions, in tna's order (checked above).
  d <- make_event_data(n_actor = 30)
  net <- cna_event(d, 3000)
  tm <- tna::ctna(tna_sessions(d, 3000))
  expect_equal(nrow(tm$data), nrow(net$data))
  # Strength and Closeness are defined alike in both; CNA's Diffusion and
  # Clustering are cograph's, not tna's (previous blocks), so they are not
  # compared. CNA titles the measures as its table does.
  st <- cna("cna_stability_compute")(net, c("Strength", "Closeness"),
                                     iter = 60L, threshold = 0.7, certainty = 0.95, seed = 3L)
  ts <- with_seed(3L, tna::estimate_cs(tm, measures = c("InStrength", "Closeness"),
                                       iter = 60L, threshold = 0.7, certainty = 0.95))
  expect_setequal(names(st$cs$cs), unname(tna_to_cna[names(ts)]))
  invisible(lapply(names(ts), function(m) {
    expect_equal(unname(st$cs$correlations[[tna_to_cna[[m]]]]), unname(ts[[m]]$correlations), info = m)
    expect_equal(unname(st$cs$cs[[tna_to_cna[[m]]]]), ts[[m]]$cs_coefficient, info = m)
  }))

  # One-hot: cases are actors (rows of import_onehot()'s output), windows 1-2.
  srl <- srl_onehot()
  invisible(lapply(c(1L, 2L), function(window) {
    neto <- cna_onehot(srl, srl_codes, window, actor = "Actor")
    tmo <- tna::ctna(tna_onehot(srl, srl_codes, window, actor = "Actor"))
    st <- cna("cna_stability_compute")(neto, c("Strength", "Betweenness", "Closeness"),
                                       iter = 40L, threshold = 0.7, certainty = 0.95, seed = 8L)
    ts <- with_seed(8L, tna::estimate_cs(tmo, measures = c("InStrength", "Betweenness", "Closeness"),
                                         iter = 40L, threshold = 0.7, certainty = 0.95))
    invisible(lapply(names(ts), function(m) {
      expect_equal(unname(st$cs$correlations[[tna_to_cna[[m]]]]), unname(ts[[m]]$correlations),
                   info = sprintf("window %d, %s", window, m))
      expect_equal(unname(st$cs$cs[[tna_to_cna[[m]]]]), ts[[m]]$cs_coefficient,
                   info = sprintf("window %d, %s", window, m))
    }))
  }))
})

test_that("grouped stability: CNA re-seeds every group, tna continues one stream", {
  skip_on_cran()
  srl <- srl_onehot()
  nets <- cna_onehot(srl, srl_codes, 2L, actor = "Actor", group = "Group")
  oh <- tna_onehot(srl, srl_codes, 2L, actor = "Actor")
  gm <- tna::group_ctna(oh, group = unit_groups(srl, "Actor", "Group"), params = window_params(oh))
  st <- cna("cna_stability_compute")(nets, c("Strength", "Closeness"), iter = 30L,
                                     threshold = 0.7, certainty = 0.95, seed = 9L)
  per_group <- lapply(names(gm), \(k)
    with_seed(9L, tna::estimate_cs(gm[[k]], measures = c("InStrength", "Closeness"), iter = 30L)))
  names(per_group) <- names(gm)
  invisible(lapply(names(gm), function(k) invisible(lapply(c("InStrength", "Closeness"), \(m)
    expect_equal(unname(st$cs[[k]]$correlations[[tna_to_cna[[m]]]]),
                 unname(per_group[[k]][[m]]$correlations), info = sprintf("%s, %s", k, m))))))
  # tna's grouped call: first group identical, the second continues the stream
  joint <- with_seed(9L, tna::estimate_cs(gm, measures = c("InStrength", "Closeness"), iter = 30L))
  first <- names(gm)[[1L]]
  second <- names(gm)[[2L]]
  expect_equal(unname(st$cs[[first]]$correlations$Strength), unname(joint[[first]]$InStrength$correlations))
  expect_false(isTRUE(all.equal(unname(st$cs[[second]]$correlations$Strength),
                                unname(joint[[second]]$InStrength$correlations))))
})

test_that("CNA stability runs where tna 1.3.1's estimate_cs() stops on an undefined closeness", {
  # Short sessions: some subsets leave a state unreachable, its closeness NaN.
  d <- make_event_data(n_actor = 12)
  net <- cna_event(d, 300)
  st <- cna("cna_stability_compute")(net, "Closeness", iter = 100L, threshold = 0.7,
                                     certainty = 0.95, seed = 3L)
  expect_length(cna("cna_stability_failures")(st), 0L)
  expect_true(any(is.na(st$cs$correlations$Closeness)))
  skip_if(utils::packageVersion("tna") >= "1.5.4", "tna 1.5.4 handles an undefined SD")
  tm <- tna::ctna(tna_sessions(d, 300))
  expect_error(with_seed(3L, tna::estimate_cs(tm, measures = "Closeness", iter = 100L)),
               class = "simpleError")
})

# ---------------------------------------------------------------------------
# 5. Bootstrap and permutation
# ---------------------------------------------------------------------------

test_that("bootstrap: same seed, identical p-values and intervals as tna::bootstrap()", {
  d <- make_event_data(n_actor = 30)
  net <- cna_event(d, 3000)
  tm <- tna::ctna(tna_sessions(d, 3000))
  w <- as_num(net$weights)
  bs <- cna("cna_bootstrap_compute")(net, iter = 200L, level = 0.05, seed = 4L)
  tb <- tna::bootstrap(tm, iter = 200L, level = 0.05, seed = 4L)
  s <- cna("cna_bootstrap_df")(bs, multi = FALSE, show_all = TRUE)
  t <- tb$summary
  # Convention: CNA lists each undirected pair once (diagonal included), tna
  # every non-zero cell, i.e. both directions.
  expect_equal(nrow(t), sum(w > 0))
  expect_equal(nrow(s), sum(w[upper.tri(w, diag = TRUE)] > 0))
  m <- match(paste(s$from, s$to), paste(t$from, t$to))
  expect_false(anyNA(m))
  expect_equal(s$weight, t$weight[m])
  expect_equal(s$p_value, t$p_value[m])
  expect_equal(s$ci_lower, t$ci_lower[m])
  expect_equal(s$ci_upper, t$ci_upper[m])
  expect_equal(s$cr_lower, t$cr_lower[m])
  expect_equal(s$cr_upper, t$cr_upper[m])
  expect_identical(s$significant == "Yes", t$sig[m])

  # Threshold inference: the default threshold is tna's 10% quantile of all cells.
  bt <- cna("cna_bootstrap_compute")(net, iter = 100L, level = 0.05, inference = "threshold", seed = 4L)
  tt <- tna::bootstrap(tm, iter = 100L, level = 0.05, method = "threshold", seed = 4L)
  expect_equal(unname(bt$edge_threshold), unname(stats::quantile(tm$weights, probs = 0.1)))
  mt <- match(paste(bt$summary$from, bt$summary$to), paste(tt$summary$from, tt$summary$to))
  expect_equal(bt$summary$p_value, tt$summary$p_value[mt])
  bt30 <- cna("cna_bootstrap_compute")(net, iter = 100L, level = 0.05, inference = "threshold",
                                       threshold = 30, seed = 4L)
  tt30 <- tna::bootstrap(tm, iter = 100L, level = 0.05, method = "threshold", threshold = 30, seed = 4L)
  m30 <- match(paste(bt30$summary$from, bt30$summary$to), paste(tt30$summary$from, tt30$summary$to))
  expect_equal(bt30$summary$p_value, tt30$summary$p_value[m30])
})

test_that("grouped bootstrap: each CNA group equals tna::bootstrap() of that group with the seed", {
  skip_on_cran()
  check_group <- function(bs, gm, info) invisible(lapply(names(gm), function(k) {
    tb <- tna::bootstrap(gm[[k]], iter = 100L, level = 0.05, seed = 2L)
    s <- bs[[k]]$summary
    m <- match(paste(s$from, s$to), paste(tb$summary$from, tb$summary$to))
    expect_equal(s$p_value, tb$summary$p_value[m], info = paste(info, k))
    expect_equal(s$ci_lower, tb$summary$ci_lower[m], info = paste(info, k))
    expect_equal(s$ci_upper, tb$summary$ci_upper[m], info = paste(info, k))
  }))
  # one-hot by actor, window 2
  srl <- srl_onehot()
  nets <- cna_onehot(srl, srl_codes, 2L, actor = "Actor", group = "Group")
  oh <- tna_onehot(srl, srl_codes, 2L, actor = "Actor")
  gm <- tna::group_ctna(oh, group = unit_groups(srl, "Actor", "Group"), params = window_params(oh))
  check_group(cna("cna_bootstrap_compute")(nets, iter = 100L, level = 0.05, seed = 2L), gm, "one-hot")
  # events; tna's grouped call seeds once and continues the stream, so only its
  # first group equals CNA's
  d <- make_event_data(n_actor = 30)
  ge <- tna::group_ctna(tna_sessions(d, 3000), group = "Group")
  be <- cna("cna_bootstrap_compute")(cna_event(d, 3000, group = "Group"), iter = 100L,
                                     level = 0.05, seed = 2L)
  check_group(be, ge, "events")
  joint <- tna::bootstrap(ge, iter = 100L, level = 0.05, seed = 2L)
  key <- \(s) paste(s$from, s$to)
  expect_equal(be$G1$summary$p_value, joint$G1$summary$p_value[match(key(be$G1$summary), key(joint$G1$summary))])
  expect_false(isTRUE(all.equal(be$G2$summary$p_value,
                                joint$G2$summary$p_value[match(key(be$G2$summary), key(joint$G2$summary))])))
})

test_that("permutation: identical differences; p-values and effect sizes agree within Monte Carlo error", {
  skip_on_cran()
  # Exact equality is not reachable: tna draws a full permutation per iteration
  # (sample(n)), Nestimate only the first group's indices (sample.int(n, n_x)),
  # so the streams part after the first draw. Bounds: with 2000 permutations the
  # SE of a difference of two independent p-value estimates is at most
  # sqrt(2 * 0.25 / 2000) = 0.016, so 0.06 is 3.8 SE; an effect size divides by
  # the SD of the permutation differences, whose relative SE is about
  # 1 / sqrt(2 * 2000) = 1.6% per package, so a 10% mean relative tolerance is
  # about 4.5 SE (Nestimate's population vs tna's sample SD differs by a factor
  # sqrt(2000 / 1999), negligible here).
  d <- make_event_data(n_actor = 30)
  srl <- srl_onehot()
  oh <- tna_onehot(srl, srl_codes, 2L, actor = "Actor")
  cases <- list(
    events = list(cna = cna_event(d, 3000, group = "Group"),
                  tna = tna::group_ctna(tna_sessions(d, 3000), group = "Group")),
    onehot = list(cna = cna_onehot(srl, srl_codes, 2L, actor = "Actor", group = "Group"),
                  tna = tna::group_ctna(oh, group = unit_groups(srl, "Actor", "Group"),
                                        params = window_params(oh))))
  invisible(Map(function(cs, nm) {
    pt <- cna("cna_permutation_compute")(cs$cna, iter = 2000L, seed = 5L)
    tp <- with_seed(5L, tna::permutation_test(cs$tna, iter = 2000L))
    s <- pt[[1]]$summary
    t <- tp[[1]]$edges$stats
    m <- match(paste(s$from, "->", s$to), t$edge_name)
    expect_false(anyNA(m), info = nm)
    expect_equal(s$diff, t$diff_true[m], info = nm)
    g1 <- names(cs$tna)[[1L]]
    expect_equal(s$weight_x, cs$tna[[g1]]$weights[cbind(s$from, s$to)], info = nm)
    expect_lte(max(abs(s$p_value - t$p_value[m])), 0.06)
    expect_equal(s$effect_size, t$effect_size[m], tolerance = 0.1, info = nm)
  }, cases, names(cases)))
})

test_that("permutation adjust: tna and Nestimate correct over all n x n cells, CNA over the unique edges", {
  d <- make_event_data(n_actor = 30)
  net <- cna_event(d, 3000, group = "Group")
  gm <- tna::group_ctna(tna_sessions(d, 3000), group = "Group")
  n <- nrow(gm[[1]]$weights)
  # tna: the n^2 cells, i.e. every undirected edge twice
  raw_t <- with_seed(5L, tna::permutation_test(gm, iter = 200L))[[1]]$edges$stats
  bh_t <- with_seed(5L, tna::permutation_test(gm, iter = 200L, adjust = "BH"))[[1]]$edges$stats
  bf_t <- with_seed(5L, tna::permutation_test(gm, iter = 200L, adjust = "bonferroni"))[[1]]$edges$stats
  expect_equal(nrow(raw_t), n^2)
  expect_equal(bh_t$p_value, stats::p.adjust(raw_t$p_value, "BH"))
  expect_equal(bf_t$p_value, pmin(1, raw_t$p_value * n^2))
  # Nestimate's own adjust follows tna (n^2 cells)
  raw_n <- Nestimate::permutation(net, iter = 200L, seed = 5L)[[1]]
  bh_n <- Nestimate::permutation(net, iter = 200L, seed = 5L, adjust = "BH")[[1]]
  cells <- raw_n$p_values
  cells[] <- stats::p.adjust(as.vector(raw_n$p_values), "BH")
  expect_equal(bh_n$summary$p_value, unname(cells[cbind(bh_n$summary$from, bh_n$summary$to)]))
  # CNA: the n (n + 1) / 2 unique edges (diagonal included)
  bh_c <- cna("cna_permutation_compute")(net, iter = 200L, seed = 5L, adjust = "BH")[[1]]$summary
  bf_c <- cna("cna_permutation_compute")(net, iter = 200L, seed = 5L, adjust = "bonferroni")[[1]]$summary
  expect_equal(nrow(bh_c), n * (n + 1) / 2)
  expect_equal(bh_c$p_value, stats::p.adjust(raw_n$summary$p_value, "BH"))
  expect_equal(bf_c$p_value, pmin(1, raw_n$summary$p_value * n * (n + 1) / 2))
})

# ---------------------------------------------------------------------------
# 6. Communities and cliques
# ---------------------------------------------------------------------------

test_that("communities: tna's partitions are CNA's on the A + t(A) graph tna builds (loops kept once)", {
  skip_on_cran()
  # tna::communities() runs igraph on graph_from_adjacency_matrix(mode = "plus")
  # (off-diagonal w_ij + w_ji = 2 w_ij, a loop keeps w_ii) or on the directed
  # graph (which these methods treat the same); cograph uses w_ij and w_ii. So
  # CNA given 2W with W's diagonal must reproduce tna. Louvain and Leiden are not
  # in tna 1.3.1.
  srl <- srl_onehot()
  ev <- srl_events()
  cases <- list(
    onehot1 = list(tna = tna::ctna(tna_onehot(srl, srl_codes, 1L, actor = "Actor")),
                   cna = cna_onehot(srl, srl_codes, 1L, actor = "Actor")),
    onehot2 = list(tna = tna::ctna(tna_onehot(srl, srl_codes, 2L, actor = "Actor")),
                   cna = cna_onehot(srl, srl_codes, 2L, actor = "Actor")),
    events = list(tna = tna::ctna(tna_sessions(ev, 900, order = "Order")),
                  cna = cna_event(ev, 900, order = "Order")))
  methods <- c(walktrap = "walktrap", fast_greedy = "fast_greedy",
               leading_eigen = "leading_eigenvector", spinglass = "spinglass",
               edge_betweenness = "edge_betweenness")
  # igraph's diagnosed note for edge betweenness (tna muffles it the same way)
  quiet_eb <- function(expr) withCallingHandlers(expr, warning = function(w) {
    if (grepl("Membership vector will be selected", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
  invisible(Map(function(cs, nm) {
    # in tna's node order: edge betweenness breaks ties by vertex order
    w <- as_num(cs$cna$weights, rownames(cs$tna$weights))
    plus <- w + t(w)
    diag(plus) <- diag(w)
    net_plus <- Nestimate::build_network(plus, method = "co_occurrence")
    invisible(Map(function(tna_m, cna_m) {
      ta <- with_seed(1L, tna::communities(cs$tna, methods = tna_m))$assignments
      cm <- quiet_eb(cna("cna_community_compute")(net_plus, cna_m, seed = 1L))[[1]]
      expect_true(same_partition(ta[[tna_m]], cm$community[match(ta$state, cm$node)]),
                  info = sprintf("%s, %s", nm, tna_m))
    }, names(methods), methods))
  }, cases, names(cases)))
  # On the network itself CNA's loops weigh twice as much relative to the
  # edges: with the large one-hot diagonal fast greedy splits differently.
  ta <- tna::communities(cases$onehot1$tna, methods = "fast_greedy")$assignments
  cm <- cna("cna_community_compute")(cases$onehot1$cna, "fast_greedy")[[1]]
  expect_false(same_partition(ta$fast_greedy, cm$community[match(ta$state, cm$node)]))
})

test_that("cliques: the same cliques as tna::cliques() for sizes 2-4, threshold inclusive", {
  ev <- srl_events()
  net <- cna_event(ev, 900, order = "Order")
  tm <- tna::ctna(tna_sessions(ev, 900, order = "Order"))
  w <- as_num(net$weights)
  boundary <- w["Evaluate", "Read"]   # an existing weight: kept by both (>=)
  grid <- expand.grid(size = 2:4, threshold = c(0, 400, boundary, 600))
  invisible(Map(function(size, threshold) {
    tc <- tna::cliques(tm, size = size, threshold = threshold)
    t_members <- sort(vapply(tc$weights, \(x) paste(sort(rownames(x)), collapse = ", "), ""))
    cd <- cna("cna_cliques_df")(cna("cna_cliques_compute")(net, size, threshold), net, size)
    c_members <- if (is.null(cd)) character() else
      sort(vapply(strsplit(cd$members, ", ", fixed = TRUE), \(s) paste(sort(s), collapse = ", "), ""))
    expect_identical(c_members, t_members, info = sprintf("size %d, threshold %g", size, threshold))
  }, grid$size, grid$threshold))
  # mean edge weight of each clique = mean of tna's clique matrix off the diagonal
  tc <- tna::cliques(tm, size = 3, threshold = boundary)
  cd <- cna("cna_cliques_df")(cna("cna_cliques_compute")(net, 3L, boundary), net, 3L)
  t_key <- vapply(tc$weights, \(x) paste(sort(rownames(x)), collapse = ", "), "")
  t_mean <- vapply(tc$weights, \(x) mean(x[upper.tri(x)]), numeric(1))
  c_key <- vapply(strsplit(cd$members, ", ", fixed = TRUE), \(s) paste(sort(s), collapse = ", "), "")
  expect_equal(cd$mean_weight, unname(t_mean[match(c_key, t_key)]))
})

# ---------------------------------------------------------------------------
# 7. Scaling and resampling of scaled models
# ---------------------------------------------------------------------------

test_that("scaling: one network equals tna; groups scale per group; rank keeps absent edges absent", {
  d <- make_event_data(n_actor = 12)
  pd <- tna_sessions(d, 300)
  invisible(lapply(c("minmax", "max", "rank"), function(sc) {
    tw <- tna::ctna(pd, scaling = sc)$weights
    expect_equal(as_num(tw), as_num(cna_event(d, 300, scaling = sc)$weights, rownames(tw)), info = sc)
  }))
  raw <- cna_event(d, 300, group = "Group")
  nodes <- rownames(raw$G1$weights)
  pooled <- unlist(lapply(raw, \(net) as.vector(as_num(net$weights, nodes))))
  global <- list(minmax = \(w) (w - min(pooled)) / (max(pooled) - min(pooled)),
                 max = \(w) w / max(pooled))
  invisible(lapply(c("minmax", "max"), function(sc) {
    nets <- cna_event(d, 300, group = "Group", scaling = sc)
    # CNA = tna's groupwise = TRUE; tna's default scales the pooled cells of all
    # groups at once
    gw <- tna::group_ctna(pd, group = "Group", scaling = sc, groupwise = TRUE)
    gg <- tna::group_ctna(pd, group = "Group", scaling = sc)
    named <- function(m) {
      w <- m$weights
      dimnames(w) <- list(m$labels, m$labels)
      as_num(w, nodes)
    }
    invisible(lapply(names(nets), function(k) {
      expect_equal(as_num(gw[[k]]$weights, nodes), as_num(nets[[k]]$weights, nodes),
                   info = paste(sc, k))
      # (tna's global rescaling returns the matrices without state names)
      expect_equal(named(gg[[k]]), global[[sc]](as_num(raw[[k]]$weights, nodes)),
                   info = paste("global", sc, k))
    }))
    expect_false(isTRUE(all.equal(unname(lapply(gg, named)),
                                  unname(lapply(nets, \(n) as_num(n$weights, nodes))))))
  }))
  # Rank with an absent edge (G2 has none of C-C): Nestimate ranks the non-zero
  # cells and keeps 0 (documented in its .apply_scaling); tna ranks every cell,
  # so the absent edge becomes the weakest edge.
  w <- as_num(raw$G2$weights)
  expect_true(any(w == 0))
  ref_nonzero <- w
  ref_nonzero[w != 0] <- rank(w[w != 0], ties.method = "average")
  ref_all <- w
  ref_all[] <- rank(w, ties.method = "average")
  expect_equal(as_num(cna_event(d, 300, group = "Group", scaling = "rank")$G2$weights, rownames(w)),
               ref_nonzero)
  gw <- tna::group_ctna(pd, group = "Group", scaling = "rank", groupwise = TRUE)
  expect_equal(as_num(gw$G2$weights), ref_all[rownames(gw$G2$weights), rownames(gw$G2$weights)])
})

test_that("tna resamples scaled single models and refuses only global group scaling; CNA refuses all", {
  # CNA's note said resampling runs "on raw counts, as in tna" (wording fixed
  # 2026-11-06): tna 1.3.1 in fact rescales every resample of a scaled model
  # and refuses only groups scaled globally (groupwise = FALSE).
  d <- make_event_data(n_actor = 30)
  pd <- tna_sessions(d, 3000)
  b <- tna::bootstrap(tna::ctna(pd, scaling = "minmax"), iter = 20L, seed = 1L)
  expect_s3_class(b, "tna_bootstrap")
  expect_equal(max(b$weights_orig), 1)
  expect_error(tna::bootstrap(tna::group_ctna(pd, group = "Group", scaling = "minmax"), iter = 20L),
               class = "rlang_error")
  expect_s3_class(with_seed(1L, tna::permutation_test(
    tna::group_ctna(pd, group = "Group", scaling = "minmax", groupwise = TRUE), iter = 20L)),
    "group_tna_permutation")
  ctx <- list(options = opts_event(buildModel_scaling = "minmax"),
              spec = cna(".cna_kinds")$event_single)
  expect_match(cna(".cna_resampling_problem")(ctx, "bootstrap"), "raw counts", fixed = TRUE)
})

# ---------------------------------------------------------------------------
# 8. Invariant and error path
# ---------------------------------------------------------------------------

test_that("invariant: relabelling and reordering actors changes neither tna's nor CNA's network", {
  d <- make_onehot_data(n = 48, seed = 11)
  actors <- unique(d$Actor)
  relabelled <- d
  relabelled$Actor <- sprintf("z%02d", rev(seq_along(actors)))[match(d$Actor, actors)]
  # actors in reverse order, each actor's rows kept in their order
  ord <- order(match(relabelled$Actor, rev(unique(relabelled$Actor))), seq_len(nrow(d)))
  moved <- relabelled[ord, , drop = FALSE]
  configs <- data.frame(window = c(1L, 2L, 3L, 2L),
                        type = c("tumbling", "tumbling", "tumbling", "sliding"))
  invisible(Map(function(window, type) {
    info <- sprintf("window %d %s", window, type)
    t0 <- tna_onehot_weights(d, fixture_codes, window, type, actor = "Actor", session = "Session")
    t1 <- tna_onehot_weights(moved, fixture_codes, window, type, actor = "Actor", session = "Session")
    c0 <- cna_onehot(d, fixture_codes, window, type, actor = "Actor", session = "Session")$weights
    c1 <- cna_onehot(moved, fixture_codes, window, type, actor = "Actor", session = "Session")$weights
    expect_equal(as_num(t1), as_num(t0, rownames(t1)), info = info)
    expect_equal(as_num(c1, rownames(t1)), as_num(c0, rownames(t1)), info = info)
    expect_equal(as_num(t1), as_num(c1, rownames(t1)), info = info)
  }, configs$window, configs$type))
})

test_that("a one-hot value other than 0/1 is a classed CNA error; tna counts it as inactive", {
  d <- make_onehot_data(n = 16)
  d$code1[3] <- 2
  expect_error(cna_onehot(d, fixture_codes, actor = "Actor"), class = "cna_not_binary")
  zero <- d
  zero$code1[3] <- 0
  expect_equal(tna_onehot_weights(d, fixture_codes, actor = "Actor"),
               tna_onehot_weights(zero, fixture_codes, actor = "Actor"))
})

test_that("a code that is never active: tna drops the node, CNA keeps it with no edges", {
  d <- make_onehot_data(n = 16)
  d$code4 <- 0
  tw <- tna_onehot_weights(d, fixture_codes, actor = "Actor")
  cw <- as_num(cna_onehot(d, fixture_codes, actor = "Actor")$weights)
  expect_setequal(rownames(tw), c("code1", "code2", "code3"))
  expect_setequal(rownames(cw), fixture_codes)
  expect_equal(unname(cw["code4", ]), rep(0, 4))
  expect_equal(as_num(tw), cw[rownames(tw), rownames(tw)])
})
