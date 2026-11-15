# Network-building FORMULAS: CNA's networks against independent base-R
# implementations written from the definitions (not against Nestimate calls;
# those live in test-build-equivalence.R). Engine: Nestimate 1.0.7 (a1754ab).
#
# Definitions established from Nestimate's source and checked by hand:
#
# Events (cna_build_event_nest -> prepare() -> .count_cooccurrence_wide):
#   * Sessions: per actor, events sorted by time (then Order, then row); a gap
#     STRICTLY greater than the threshold, in seconds, opens a new session
#     (a gap equal to it does not). Tied times never split. A missing time sorts
#     last and is a session of its own. Without a Time column (Order or not)
#     each actor is one session and the threshold is ignored.
#   * Missing actions (NA and the void markers "", "NA", "NaN", "%", "*") are
#     dropped from the counts but keep their place in time, so they can bridge
#     a gap (convention).
#   * Counts per session over all position pairs i < j: W[x, y] = n_x * n_y
#     (x != y) and W[x, x] = choose(n_x, 2), summed over sessions. Nodes are
#     every observed action, so an action seen only alone is an all-zero node.
#   * Groups: each group network is the network of that group's rows (sessions
#     re-derived within the subset) over the pooled action set.
#
# Scaling (.apply_scaling, include_zeros = TRUE for co_occurrence), over the
# whole matrix, diagonal included:
#   minmax (W - min W) / (max W - min W)  (unchanged when constant)
#   max    W / max |W|
#   rank   non-zero cells replaced by their average-tie rank among the n x n
#          non-zero cells (an undirected pair is two tied cells); zeros stay 0.
#
# Delimited (cooccurrence()): one transaction per row (items split on `sep`,
# trimmed, empty dropped, de-duplicated) or per Transaction ID; B = incidence,
# C = t(B) B, f = item frequency (transactions holding the item, = diag C).
# Similarities (van Eck & Waltman, 2009, with s_i = occurrences f_i):
#   jaccard C/(f_i+f_j-C); cosine C/sqrt(f_i f_j); inclusion C/min(f_i,f_j);
#   association C/(f_i f_j) (association strength WITHOUT the total-count
#   factor of the probabilistic affinity index); dice 2C/(f_i+f_j);
#   equivalence C^2/(f_i f_j); relative C / rowSums(C).
#   diagonal = FALSE zeroes diag(C) BEFORE normalising (f stays colSums(B));
#   min_occur drops items with f < min_occur first; threshold zeroes W < t
#   (strict, diagonal included); top_n zeroes off-diagonal cells below the
#   top_n-th largest upper-triangle value (ties at the cut-off are kept).
#
# One-hot (wtna co-occurrence): per actor x session unit, rows in data order;
#   window 1: crossprod(X); window w > 1: sum over tumbling blocks of w rows of
#   s s' (s = code totals of the block); sliding: first one presence row per
#   complete window of w rows, then the same blocks of w.

# --- comparison helpers --------------------------------------------------------

# Weight matrix as a plain numeric matrix (no class, names kept).
as_plain <- function(w) matrix(as.numeric(w), nrow(w), dimnames = dimnames(w))

# Same node set, same weights once both are aligned by node name (Nestimate
# sorts nodes alphabetically in places).
expect_weights <- function(w, ref, info = NULL) {
  w <- as_plain(w)
  testthat::expect_setequal(rownames(w), rownames(ref))
  # different node sets: reported above; indexing would abort the whole block
  if (!setequal(rownames(w), rownames(ref))) return(invisible(FALSE))
  nm <- rownames(ref)
  testthat::expect_equal(unname(w[nm, nm, drop = FALSE]), unname(ref[nm, nm, drop = FALSE]),
                         info = info)
}

# Symmetric 3 x 3 matrix over nodes a, b, c.
sym3 <- function(diag, ab, ac, bc, nodes = c("a", "b", "c")) {
  m <- diag(diag, 3)
  m[1, 2] <- m[2, 1] <- ab
  m[1, 3] <- m[3, 1] <- ac
  m[2, 3] <- m[3, 2] <- bc
  dimnames(m) <- list(nodes, nodes)
  m
}

# --- independent references ----------------------------------------------------

# Session number of each event of ONE actor: time order, a gap strictly above
# the threshold starts a new session, a missing time is a session of its own.
ref_session_index <- function(time, threshold) {
  o <- order(time)
  gap <- c(NA_real_, diff(as.numeric(time[o])))
  s <- cumsum(is.na(gap) | gap > threshold)
  s[order(o)]
}

# Event co-occurrence from the definition: per session, n_x * n_y between codes
# and choose(n_x, 2) on the diagonal, summed over sessions.
ref_event <- function(action, actor = NULL, time = NULL, threshold = Inf,
                      alphabet = NULL) {
  action <- as.character(action)
  n <- length(action)
  actor <- if (is.null(actor)) rep("all", n) else as.character(actor)
  session <- if (is.null(time)) rep(1L, n) else
    unsplit(lapply(split(time, actor), ref_session_index, threshold = threshold), actor)
  alphabet <- alphabet %||% sort(unique(action[!is.na(action)]))
  keep <- !is.na(action)
  key <- factor(paste(actor, session, sep = "\r")[keep])
  K <- matrix(as.numeric(table(factor(action[keep], levels = alphabet), key)),
              nrow = length(alphabet), dimnames = list(alphabet, NULL))
  W <- tcrossprod(K)
  diag(W) <- rowSums(K * (K - 1) / 2)
  W
}

ref_scale <- function(W, method) {
  switch(method,
    noScaling = W,
    minmax = {
      r <- range(W)
      if (r[2] > r[1]) (W - r[1]) / (r[2] - r[1]) else W
    },
    max = {
      m <- max(abs(W))
      if (m > 0) W / m else W
    },
    rank = {
      nz <- W != 0          # raw counts are whole numbers: exact comparison
      W[nz] <- rank(W[nz], ties.method = "average")
      W
    })
}

# Items of one record or transaction: trimmed, empty and missing dropped, unique.
ref_items <- function(x) {
  x <- trimws(as.character(x))
  unique(x[!is.na(x) & nzchar(x)])
}
ref_delim_transactions <- function(cells, sep)
  lapply(strsplit(as.character(cells), sep, fixed = TRUE), ref_items)
ref_long_transactions <- function(items, by) unname(lapply(split(items, by), ref_items))

ref_similarity <- function(C, f, similarity) {
  switch(similarity,
    none        = C,
    jaccard     = C / (outer(f, f, "+") - C),
    cosine      = C / sqrt(outer(f, f)),
    inclusion   = C / outer(f, f, pmin),
    association = C / outer(f, f),
    dice        = 2 * C / outer(f, f, "+"),
    equivalence = C^2 / outer(f, f),
    relative    = {
      rs <- rowSums(C)
      C / ifelse(rs > 0, rs, 1)
    })
}

ref_cooc <- function(transactions, similarity = "none", diagonal = TRUE,
                     min_occur = 1L, threshold = 0, top_n = 0L) {
  freq_all <- table(unlist(transactions))
  items <- sort(names(freq_all)[freq_all >= min_occur])
  B <- matrix(as.numeric(unlist(lapply(transactions, \(t) items %in% t))),
              ncol = length(items), byrow = TRUE, dimnames = list(NULL, items))
  C <- crossprod(B)
  f <- colSums(B)
  if (!diagonal) diag(C) <- 0
  W <- ref_similarity(C, f, similarity)
  W[W < threshold] <- 0
  upper <- W[upper.tri(W)]
  if (top_n > 0L && length(upper) > top_n) {
    cut <- sort(upper, decreasing = TRUE)[top_n]
    W[W < cut & row(W) != col(W)] <- 0
  }
  W
}

# One-hot co-occurrence from the window definition, per actor x session unit.
ref_onehot <- function(df, codes, units = NULL, window = 1L, sliding = FALSE) {
  X <- as.matrix(df[codes])
  X[is.na(X)] <- 0
  k <- length(codes)
  rows_by_unit <- if (length(units)) split(seq_len(nrow(df)), df[units], drop = TRUE)
                  else list(seq_len(nrow(df)))
  per_unit <- function(rows) {
    M <- X[rows, , drop = FALSE]
    if (window <= 1L) return(crossprod(M))
    if (sliding) {
      starts <- seq_len(max(nrow(M) - window + 1L, 0L))
      M <- matrix(as.numeric(unlist(lapply(starts, \(s)
        colSums(M[s:(s + window - 1L), , drop = FALSE]) > 0))),
        ncol = k, byrow = TRUE)
    }
    if (!nrow(M)) return(matrix(0, k, k))
    S <- rowsum(M, (seq_len(nrow(M)) - 1L) %/% window)
    crossprod(S)
  }
  W <- Reduce(`+`, lapply(rows_by_unit, per_unit), matrix(0, k, k))
  dimnames(W) <- list(codes, codes)
  W
}

# --- random data -----------------------------------------------------------------

# Events: 1-8 actors with 1-15 events each (single-event actors included), tied
# integer times, ~7% missing actions, shuffled rows; numeric (Unix seconds) or
# POSIXct times by seed. Group is per actor, RowGroup varies within actors.
rand_events <- function(seed) withr::with_seed(seed, {
  n_actor <- sample.int(8L, 1L)
  alphabet <- LETTERS[seq_len(sample(2:6, 1L))]
  actor <- rep(sprintf("u%02d", seq_len(n_actor)), sample.int(15L, n_actor, replace = TRUE))
  n <- length(actor)
  action <- sample(alphabet, n, replace = TRUE)
  action[stats::runif(n) < 0.07] <- NA
  action[1L] <- alphabet[1L]
  time <- sample(0:sample(c(10, 60, 600), 1L), n, replace = TRUE)
  if (seed %% 2L) time <- as.POSIXct(time, origin = "1970-01-01", tz = "UTC")
  d <- data.frame(Actor = actor, Action = action, Time = time,
                  Ord = stats::runif(n),
                  Group = ifelse(match(actor, unique(actor)) %% 2L, "G1", "G2"),
                  RowGroup = sample(c("R1", "R2"), n, replace = TRUE),
                  stringsAsFactors = FALSE)
  d[sample.int(n), , drop = FALSE]
})
rand_threshold <- function(seed) withr::with_seed(seed + 1000L, sample(c(1, 3, 10, 30, 1e7), 1L))

rand_onehot <- function(seed) withr::with_seed(seed, {
  n <- sample(10:60, 1L)
  codes <- paste0("c", seq_len(sample(3:5, 1L)))
  X <- matrix(stats::rbinom(n * length(codes), 1L, stats::runif(1L, 0.25, 0.6)), n)
  X[stats::runif(length(X)) < 0.03] <- NA
  d <- as.data.frame(X)
  names(d) <- codes
  d$Actor <- sample(sprintf("a%d", seq_len(sample(1:5, 1L))), n, replace = TRUE)
  d$Session <- sample(1:3, n, replace = TRUE)
  attr(d, "codes") <- codes
  d
})

rand_delimited <- function(seed) withr::with_seed(seed, {
  vocab <- c("network", "graph", "matrix", "model", "theory", "data", "time", "text")
  vocab <- vocab[seq_len(sample(3:8, 1L))]
  n <- sample(5:40, 1L)
  sep <- sample(c(";", ","), 1L)
  pad <- c("", " ", "  ")
  cells <- vapply(seq_len(n), \(i) {
    items <- sample(vocab, sample(0:5, 1L), replace = TRUE)          # repeats allowed
    items <- paste0(sample(pad, length(items), TRUE), items, sample(pad, length(items), TRUE))
    if (stats::runif(1L) < 0.2) items <- c(items, "", " ")            # empty items
    paste(items, collapse = sep)
  }, character(1))
  cells[stats::runif(n) < 0.05] <- NA
  data.frame(keywords = cells, stringsAsFactors = FALSE) |> structure(sep = sep)
})

rand_long <- function(seed) withr::with_seed(seed, {
  vocab <- c("network", "graph", "matrix", "model", "theory", "data")[seq_len(sample(3:6, 1L))]
  n <- sample(10:80, 1L)
  data.frame(paper = sample.int(sample(3:15, 1L), n, replace = TRUE),
             keyword = sample(vocab, n, replace = TRUE), stringsAsFactors = FALSE)
})

# A bijective relabelling of identifiers (new labels sort differently).
relabel <- function(x, seed) {
  ids <- unique(x)
  new <- withr::with_seed(seed, sample(sprintf("z%03d", seq_along(ids) * 7L)))
  new[match(x, ids)]
}

# CNA builders with the option lists of helper-cna.R.
build_event <- function(d, ...) cna("cna_build_event_nest")(d, opts_event(...))
build_delim <- function(d, ...) cna("cna_build_delimited_nest")(d, opts_delimited(...))
# Nestimate warns "One-hot data without windowing or sessions produces sparse
# networks" for window 1 without a Session (a documented advisory, diagnosed);
# only that warning is muffled, every other one propagates.
build_onehot <- function(d, ...) withCallingHandlers(
  cna("cna_build_onehot_nest")(d, opts_onehot(...)),
  warning = function(w) if (grepl("without windowing or sessions", conditionMessage(w)))
    invokeRestart("muffleWarning"))

# =============================================================================
# 1. Event co-occurrence
# =============================================================================

test_that("event pair counting: n_x * n_y between codes, choose(n, 2) on the diagonal", {
  t0 <- as.POSIXct("2021-01-01", tz = "UTC")
  # one session A A B C: AB = 2*1, AC = 2*1, BC = 1, AA = choose(2,2) = 1, BB = CC = 0
  w <- build_event(data.frame(Actor = "a", Action = c("A", "A", "B", "C"), Time = t0 + 0:3))$weights
  expect_weights(w, matrix(c(1, 2, 2, 2, 0, 1, 2, 1, 0), 3, dimnames = list(LETTERS[1:3], LETTERS[1:3])))
  # A B A A: AB = 3*1, AA = choose(3,2) = 3 (repeated codes count every position pair)
  w <- build_event(data.frame(Actor = "a", Action = c("A", "B", "A", "A"), Time = t0 + 0:3))$weights
  expect_weights(w, matrix(c(3, 3, 3, 0), 2, dimnames = list(c("A", "B"), c("A", "B"))))
  # two sessions add: (A B) and (A B B) -> AB = 1 + 2, BB = 0 + 1
  d <- data.frame(Actor = "a", Action = c("A", "B", "A", "B", "B"), Time = t0 + c(0, 1, 100, 101, 102))
  w <- build_event(d, buildModel_threshold = 10)$weights
  expect_weights(w, matrix(c(0, 3, 3, 1), 2, dimnames = list(c("A", "B"), c("A", "B"))))
})

test_that("session split: gap > threshold (seconds) splits, gap == threshold and ties do not", {
  t0 <- as.POSIXct("2021-01-01", tz = "UTC")
  ab <- function(time, thr) unname(build_event(
    data.frame(Actor = "a", Action = c("A", "B"), Time = time), buildModel_threshold = thr)$weights["A", "B"])
  expect_equal(ab(t0 + c(0, 10), 10), 1)     # equal: same session
  expect_equal(ab(t0 + c(0, 10), 9), 0)      # greater: new session
  expect_equal(ab(c(0, 10), 10), 1)          # numeric times are Unix seconds
  expect_equal(ab(c(0, 10), 9), 0)
  expect_equal(ab(t0 + c(5, 5), 1), 1)       # tied times never split
  # ISO strings with seconds: 5 s apart co-occur, 30 min apart do not
  d <- data.frame(Actor = "a", Action = c("A", "B", "C"),
                  Time = c("2025-01-01 08:00:00", "2025-01-01 08:00:05", "2025-01-01 08:30:00"))
  w <- build_event(d, buildModel_threshold = 60)$weights
  expect_weights(w, ref_event(d$Action, d$Actor, as.POSIXct(d$Time, tz = "UTC"), 60))
  expect_equal(unname(w["A", "C"]), 0)
})

test_that("Order without Time: one session per actor, threshold and order values irrelevant", {
  d <- data.frame(Actor = c("a", "a", "a", "b", "b"), Action = c("A", "B", "C", "A", "C"),
                  Ord = c(3, 1, 2, 1, 2))
  ref <- ref_event(d$Action, d$Actor)
  invisible(lapply(c(1, 900), \(thr) expect_weights(
    build_event(d, buildModel_variables_long_time = NULL, buildModel_variables_long_order = "Ord",
                buildModel_threshold = thr)$weights, ref, info = thr)))
  expect_weights(build_event(d, buildModel_variables_long_time = NULL)$weights, ref)
  # a threshold above every gap gives the same network as no Time at all
  d$Time <- c(0, 50, 100, 0, 7)
  expect_weights(build_event(d, buildModel_threshold = 1e6)$weights, ref)
})

test_that("edge cases: single events, singleton sessions, missing actions and times (conventions)", {
  t0 <- as.POSIXct("2021-01-01", tz = "UTC")
  # an actor with one event and an action seen only alone: a node with no edges
  d <- data.frame(Actor = c("a", "a", "b"), Action = c("A", "B", "Z"), Time = t0 + c(0, 1, 0))
  w <- build_event(d)$weights
  expect_weights(w, ref_event(d$Action, d$Actor, d$Time, 3000))
  expect_true("Z" %in% rownames(w))
  expect_equal(sum(w["Z", ]), 0)
  # every session of one event: all-zero matrix over the observed actions
  d <- data.frame(Actor = "a", Action = c("A", "B", "C"), Time = t0 + c(0, 100, 200))
  w <- build_event(d, buildModel_threshold = 10)$weights
  expect_setequal(rownames(w), c("A", "B", "C"))
  expect_equal(sum(w), 0)
  # a missing action is not counted but still bridges a time gap (convention)
  d <- data.frame(Actor = "a", Action = c("A", NA, "B"), Time = t0 + c(0, 8, 16))
  w <- build_event(d, buildModel_threshold = 10)$weights
  expect_setequal(rownames(w), c("A", "B"))
  expect_equal(unname(w["A", "B"]), 1)
  # void markers are missing actions (TraMineR/tna convention)
  d <- data.frame(Actor = "a", Action = c("A", "", "NA", "NaN", "%", "*", "B"), Time = t0 + 0:6)
  expect_weights(build_event(d)$weights, ref_event(c("A", NA, NA, NA, NA, NA, "B"), d$Actor, d$Time, 3000))
  # a missing time is a session of its own (sorted last, so it splits nothing)
  d <- data.frame(Actor = "a", Action = c("A", "B", "C", "D"), Time = t0 + c(0, 1, NA, 2))
  w <- build_event(d)$weights
  expect_weights(w, ref_event(d$Action, d$Actor, d$Time, 3000))
  expect_equal(sum(w["C", ]), 0)
})

test_that("Monte Carlo: event networks == the independent session/pair definition", {
  skip_on_cran()
  seeds <- seq_len(40L)
  invisible(lapply(seeds, \(s) {
    d <- rand_events(s)
    thr <- rand_threshold(s)
    w <- build_event(d, buildModel_threshold = thr)$weights
    expect_weights(w, ref_event(d$Action, d$Actor, d$Time, thr), info = sprintf("seed %d", s))
  }))
  # Order column, no Time: one session per actor whatever the threshold
  invisible(lapply(seeds[seeds <= 10L], \(s) {
    d <- rand_events(s)
    w <- build_event(d, buildModel_variables_long_time = NULL, buildModel_variables_long_order = "Ord",
                     buildModel_threshold = rand_threshold(s))$weights
    expect_weights(w, ref_event(d$Action, d$Actor), info = sprintf("order seed %d", s))
  }))
})

test_that("groups: each group network is its subset's network; actor-level groups sum to the pool", {
  # hand case: A(G1, 0) B(G2, 5) C(G1, 10), threshold 6. Pooled: one session.
  # G1's rows alone are 10 s apart, so its sessions split (subset convention).
  t0 <- as.POSIXct("2021-01-01", tz = "UTC")
  d <- data.frame(Actor = "a", Action = c("A", "B", "C"), Time = t0 + c(0, 5, 10),
                  Group = c("G1", "G2", "G1"))
  pooled <- build_event(d, buildModel_threshold = 6)$weights
  expect_equal(unname(pooled["A", "C"]), 1)
  grp <- build_event(d, buildModel_threshold = 6, buildModel_variables_group = "Group")
  expect_equal(unname(grp$G1$weights["A", "C"]), 0)
  expect_setequal(rownames(grp$G2$weights), c("A", "B", "C"))   # pooled node set

  skip_on_cran()
  invisible(lapply(seq_len(20L), \(s) {
    d <- rand_events(s)
    thr <- rand_threshold(s)
    alphabet <- sort(unique(d$Action[!is.na(d$Action)]))
    # row-level groups: network of the group's own rows, over the pooled actions
    rg <- build_event(d, buildModel_threshold = thr, buildModel_variables_group = "RowGroup")
    expect_setequal(names(rg), unique(d$RowGroup))
    invisible(lapply(names(rg), \(g) {
      sub <- subset(d, RowGroup == g)
      expect_weights(rg[[g]]$weights, ref_event(sub$Action, sub$Actor, sub$Time, thr, alphabet),
                     info = sprintf("seed %d group %s", s, g))
    }))
    # actor-level groups: sessions never cross groups, so the counts add up
    ag <- build_event(d, buildModel_threshold = thr, buildModel_variables_group = "Group")
    total <- Reduce(`+`, lapply(ag, \(n) as_plain(n$weights)[alphabet, alphabet]))
    expect_weights(total, as_plain(build_event(d, buildModel_threshold = thr)$weights),
                   info = sprintf("sum seed %d", s))
  }))
})

# =============================================================================
# 2. Scaling
# =============================================================================

test_that("scaling formulas on a hand-computed network", {
  t0 <- as.POSIXct("2021-01-01", tz = "UTC")
  # one session A A A B B C: AA = 3, AB = 6, AC = 3, BB = 1, BC = 2, CC = 0
  d <- data.frame(Actor = "a", Action = c("A", "A", "A", "B", "B", "C"), Time = t0 + 0:5)
  m3 <- function(aa, ab, ac, bb, bc, cc) matrix(c(aa, ab, ac, ab, bb, bc, ac, bc, cc), 3,
                                               dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  expect_weights(build_event(d)$weights, m3(3, 6, 3, 1, 2, 0))
  expect_weights(build_event(d, buildModel_scaling = "minmax")$weights, m3(.5, 1, .5, 1 / 6, 1 / 3, 0))
  expect_weights(build_event(d, buildModel_scaling = "max")$weights, m3(.5, 1, .5, 1 / 6, 1 / 3, 0))
  # ranks over the 7 non-zero cells: BB 1; BC,CB 2.5; AA,AC,CA 5; AB,BA 7.5
  expect_weights(build_event(d, buildModel_scaling = "rank")$weights, m3(5, 7.5, 5, 1, 2.5, 0))
  # no zero cell (A A B B): minmax maps the smallest observed cell to 0 (the
  # tna-compatible convention that includes zeros), max divides, rank ties
  d2 <- data.frame(Actor = "a", Action = c("A", "A", "B", "B"), Time = t0 + 0:3)
  m2 <- function(dg, ab) matrix(c(dg, ab, ab, dg), 2, dimnames = list(c("A", "B"), c("A", "B")))
  expect_weights(build_event(d2, buildModel_scaling = "minmax")$weights, m2(0, 1))
  expect_weights(build_event(d2, buildModel_scaling = "max")$weights, m2(.25, 1))
  expect_weights(build_event(d2, buildModel_scaling = "rank")$weights, m2(1.5, 3.5))
  # an all-zero network is left as it is
  d3 <- data.frame(Actor = c("a", "b"), Action = c("A", "B"), Time = t0 + 0:1)
  invisible(lapply(c("minmax", "max", "rank"), \(s)
    expect_equal(sum(abs(build_event(d3, buildModel_scaling = s)$weights)), 0, info = s)))
})

test_that("Monte Carlo: scaled event and one-hot networks == formula applied to the raw matrix", {
  skip_on_cran()
  methods <- c("minmax", "max", "rank")
  invisible(lapply(seq_len(15L), \(s) {
    d <- rand_events(s)
    thr <- rand_threshold(s)
    raw <- as_plain(build_event(d, buildModel_threshold = thr)$weights)
    invisible(lapply(methods, \(m) expect_weights(
      build_event(d, buildModel_threshold = thr, buildModel_scaling = m)$weights,
      ref_scale(raw, m), info = sprintf("event seed %d %s", s, m))))
    # groups: every group is scaled by its own range
    grp <- build_event(d, buildModel_threshold = thr, buildModel_variables_group = "RowGroup",
                       buildModel_scaling = "minmax")
    raw_grp <- build_event(d, buildModel_threshold = thr, buildModel_variables_group = "RowGroup")
    invisible(lapply(names(grp), \(g) expect_weights(
      grp[[g]]$weights, ref_scale(as_plain(raw_grp[[g]]$weights), "minmax"),
      info = sprintf("group seed %d %s", s, g))))

    oh <- rand_onehot(s)
    codes <- attr(oh, "codes")
    win <- 1L + s %% 3L
    raw <- ref_onehot(oh, codes, c("Actor", "Session"), win)
    invisible(lapply(methods, \(m) expect_weights(
      build_onehot(oh, buildModel_variables_onehot = codes, buildModel_window = win,
                   buildModel_scaling = m)$weights,
      ref_scale(raw, m), info = sprintf("onehot seed %d %s", s, m))))
  }))
})

test_that("delimited networks have no scaling option (convention)", {
  expect_false("buildModel_scaling" %in% option_names("CNADelimited"))
  dat <- make_delimited_data()
  expect_equal(build_delim(dat, buildModel_scaling = "max")$weights, build_delim(dat)$weights)
})

# =============================================================================
# 3. Delimited keywords
# =============================================================================

# records {a,b} {a,c} {a,b,c} {a}: f = (4, 2, 2); C_ab = 2, C_ac = 2, C_bc = 1
hand_records <- data.frame(keywords = c("a;b", "a;c", "a;b;c", "a"), stringsAsFactors = FALSE)

test_that("similarity formulas on a hand-computed example", {
  expected <- list(
    none        = sym3(c(4, 2, 2), 2, 2, 1),
    jaccard     = sym3(1, 2 / 4, 2 / 4, 1 / 3),
    cosine      = sym3(1, 2 / sqrt(8), 2 / sqrt(8), 1 / 2),
    inclusion   = sym3(1, 1, 1, 1 / 2),
    association = sym3(c(1 / 4, 1 / 2, 1 / 2), 2 / 8, 2 / 8, 1 / 4),
    dice        = sym3(1, 4 / 6, 4 / 6, 2 / 4),
    equivalence = sym3(1, 4 / 8, 4 / 8, 1 / 4),
    relative    = matrix(c(4 / 8, 2 / 5, 2 / 5, 2 / 8, 2 / 5, 1 / 5, 2 / 8, 1 / 5, 2 / 5), 3,
                         dimnames = list(letters[1:3], letters[1:3])))
  invisible(Map(\(sim, ref) expect_weights(build_delim(hand_records, buildModel_similarity = sim)$weights,
                                           ref, info = sim), names(expected), expected))
  # diagonal off: diag(C) zeroed before normalising, frequencies unchanged
  expect_weights(build_delim(hand_records, buildModel_similarity = "jaccard",
                             buildModel_diagonal = FALSE)$weights, sym3(0, 2 / 4, 2 / 4, 1 / 3))
  expect_weights(build_delim(hand_records, buildModel_similarity = "relative",
                             buildModel_diagonal = FALSE)$weights,
                 matrix(c(0, 2 / 3, 2 / 3, 2 / 4, 0, 1 / 3, 2 / 4, 1 / 3, 0), 3,
                        dimnames = list(letters[1:3], letters[1:3])))
})

test_that("threshold, min_occur and top_n follow their definitions", {
  # threshold: strictly-below zeroed, the diagonal included
  expect_weights(build_delim(hand_records, buildModel_edge_threshold = 2)$weights,
                 sym3(c(4, 2, 2), 2, 2, 0))
  expect_weights(build_delim(hand_records, buildModel_edge_threshold = 3)$weights,
                 sym3(c(4, 0, 0), 0, 0, 0))
  # min_occur = 2 keeps a, b, c (f >= 2); min_occur = 3 keeps a only
  expect_weights(build_delim(hand_records, buildModel_min_occur = 2L)$weights, sym3(c(4, 2, 2), 2, 2, 1))
  w <- build_delim(hand_records, buildModel_min_occur = 3L)$weights
  expect_equal(rownames(w), "a")
  expect_equal(unname(as_plain(w)), matrix(4))
  # top_n = 1 keeps the strongest edges; ties at the cut-off are all kept
  expect_weights(build_delim(hand_records, buildModel_top_n = 1L)$weights, sym3(c(4, 2, 2), 2, 2, 0))
  ties <- data.frame(keywords = c("a;b", "c;d", "a;c", "b;d"))
  w <- as_plain(build_delim(ties, buildModel_top_n = 2L)$weights)
  expect_equal(sum(w[upper.tri(w)] > 0), 4)
})

test_that("delimited parsing: whitespace and empty items dropped, repeats within a record once", {
  d <- data.frame(keywords = c("  network ;graph;; graph ;  ", "graph", "", NA, " ; ",
                               "network;network;matrix"), stringsAsFactors = FALSE)
  w <- build_delim(d)$weights
  # transactions: {network, graph}, {graph}, {network, matrix}
  expect_weights(w, matrix(c(2, 1, 1, 1, 2, 0, 1, 0, 1), 3,
                           dimnames = rep(list(c("network", "graph", "matrix")), 2)))
  expect_equal(build_delim(d)$params$n_transactions, 3L)
})

test_that("Monte Carlo: delimited networks == the transaction/similarity definition", {
  skip_on_cran()
  sims <- c("none", "jaccard", "cosine", "inclusion", "association", "dice", "equivalence", "relative")
  invisible(lapply(seq_len(30L), \(s) {
    d <- rand_delimited(s)
    sep <- attr(d, "sep")
    tr <- Filter(length, ref_delim_transactions(d$keywords, sep))
    max_f <- max(table(unlist(tr)))
    o <- withr::with_seed(s + 500L, list(
      similarity = sample(sims, 1L), diagonal = stats::runif(1L) < 0.7,
      min_occur = min(sample(1:3, 1L), max_f), top_n = sample(c(0L, 0L, 1L, 3L, 5L), 1L),
      thresholded = stats::runif(1L) < 0.5))
    o$threshold <- if (!o$thresholded) 0 else if (o$similarity == "none") 2 else 0.3
    o$thresholded <- NULL
    w <- build_delim(d, buildModel_sep = sep, buildModel_similarity = o$similarity,
                     buildModel_diagonal = o$diagonal, buildModel_min_occur = o$min_occur,
                     buildModel_edge_threshold = o$threshold, buildModel_top_n = o$top_n)$weights
    ref <- ref_cooc(tr, o$similarity, o$diagonal, o$min_occur, o$threshold, o$top_n)
    expect_weights(w, ref, info = sprintf("seed %d %s", s, paste(names(o), o, collapse = " ")))
  }))
})

test_that("Monte Carlo: Transaction-ID (one item per row) == the transaction definition", {
  skip_on_cran()
  invisible(lapply(seq_len(20L), \(s) {
    d <- rand_long(s)
    sim <- c("none", "jaccard", "association", "dice")[1L + s %% 4L]
    w <- build_delim(d, buildModel_variables_field = "keyword", buildModel_variables_by = "paper",
                     buildModel_similarity = sim)$weights
    expect_weights(w, ref_cooc(ref_long_transactions(d$keyword, d$paper), sim),
                   info = sprintf("seed %d %s", s, sim))
  }))
  # an item repeated within one transaction counts once
  d <- data.frame(paper = c(1, 1, 1, 2, 2), keyword = c("a", "a", "b", "b", "c"))
  w <- build_delim(d, buildModel_variables_field = "keyword", buildModel_variables_by = "paper")$weights
  expect_weights(w, sym3(c(1, 2, 1), 1, 0, 1))
})

test_that("Transaction-ID items are trimmed and blank items dropped, as in delimited mode", {
  # Fixed in CNA 2026-11-06 (.cna_transaction_items), pinned here. Was
  # (Nestimate .co_parse_long): Transaction-ID values are used verbatim,
  # so " a" and "a" are two nodes and a blank "  " becomes a node, while the
  # same keywords in delimited mode are trimmed and blanks dropped.
  d <- data.frame(paper = c(1, 1, 1, 2, 2, 3, 3), keyword = c("a", " a", "  ", "b", "", "a ", "b"),
                  stringsAsFactors = FALSE)
  w <- build_delim(d, buildModel_variables_field = "keyword", buildModel_variables_by = "paper")$weights
  expect_setequal(rownames(w), c("a", "b"))
  expect_weights(w, ref_cooc(ref_long_transactions(d$keyword, d$paper)))
})

test_that("several Item Columns with a Transaction ID pool the items of each transaction", {
  # Fixed in CNA 2026-11-06 (.cna_transaction_items pools the columns), pinned
  # here. Was (Nestimate cooccurrence(field = <2+ columns>, by =) / CNA passed the
  # vector unguarded): fails with "subscript out of bounds". The UI allows one
  # or more Item Columns together with a Transaction ID.
  d <- data.frame(paper = c(1, 1, 2), k1 = c("a", "b", "a"), k2 = c("c", "d", "c"),
                  stringsAsFactors = FALSE)
  w <- tryCatch(build_delim(d, buildModel_variables_field = c("k1", "k2"),
                            buildModel_variables_by = "paper")$weights,
                error = function(e) e)
  expect_false(inherits(w, "error"), info = if (inherits(w, "error")) conditionMessage(w))
  if (!inherits(w, "error"))
    expect_weights(w, ref_cooc(list(c("a", "b", "c", "d"), c("a", "c"))))
})

test_that("several Item Columns with a Separator pool each row's items", {
  d <- data.frame(k1 = c("a;b", "c", "a"), k2 = c("c", "c;d", "a;b"), stringsAsFactors = FALSE)
  w <- build_delim(d, buildModel_variables_field = c("k1", "k2"))$weights
  tr <- Map(\(x, y) ref_items(c(strsplit(x, ";")[[1]], strsplit(y, ";")[[1]])), d$k1, d$k2)
  expect_weights(w, ref_cooc(unname(tr)))
})

test_that("relative similarity is not offered: it is directed, Nestimate returns it undirected", {
  # Row normalisation is asymmetric by definition (b -> a = 2/5, a -> b = 2/8),
  # but Nestimate 1.0.7's cooccurrence() hard-codes directed = FALSE and its
  # edge list keeps only the upper triangle, silently dropping b -> a, c -> a,
  # c -> b. CNA's networks are undirected throughout, so the option is
  # withdrawn until Nestimate returns a directed network
  # (docs/EQUIVALENCE_AUDIT.md).
  sim <- option_def("CNADelimited", "buildModel_similarity")
  offered <- vapply(sim$options, function(o) if (is.list(o)) o$name else o, character(1))
  expect_false("relative" %in% offered)
  expect_setequal(offered, c("none", "jaccard", "cosine", "inclusion", "association",
                             "dice", "equivalence"))
})

test_that("a whitespace-only Separator is refused at input, not read as 'no separator' (convention)", {
  # The builder treats " " as no separator (each cell one item); the analysis
  # refuses it first, so space-separated keywords cannot be analysed.
  problem <- cna(".cna_input_problem")
  expect_match(problem(list(buildModel_variables_field = "k", buildModel_sep = " "),
                       list(data = "delimited")), "Separator")
  expect_null(problem(list(buildModel_variables_field = "k", buildModel_sep = ";"),
                      list(data = "delimited")))
})

# =============================================================================
# 4. One-hot sweep
# =============================================================================

test_that("Monte Carlo: one-hot == window definition (seeds x windows 1-4 x mode x units)", {
  skip_on_cran()
  grid <- expand.grid(seed = seq_len(8L), window = 1:4, type = c("tumbling", "sliding"),
                      units = c("actor", "actor_session"), stringsAsFactors = FALSE)
  invisible(Map(\(seed, window, type, units) {
    oh <- rand_onehot(seed)
    codes <- attr(oh, "codes")
    session <- if (units == "actor_session") "Session"
    w <- build_onehot(oh, buildModel_variables_onehot = codes, buildModel_variables_session = session,
                      buildModel_window = window, buildModel_window_type = type)$weights
    ref <- ref_onehot(oh, codes, c("Actor", session), window, sliding = type == "sliding")
    expect_weights(w, ref, info = sprintf("seed %d window %d %s %s", seed, window, type, units))
  }, grid$seed, grid$window, grid$type, grid$units))
})

# =============================================================================
# 5. Invariants
# =============================================================================

test_that("row-order permutation and actor relabelling leave event networks unchanged", {
  invisible(lapply(seq_len(10L), \(s) {
    d <- rand_events(s)
    thr <- rand_threshold(s)
    base <- build_event(d, buildModel_threshold = thr)$weights
    shuffled <- d[withr::with_seed(s, sample.int(nrow(d))), , drop = FALSE]
    expect_weights(build_event(shuffled, buildModel_threshold = thr)$weights, as_plain(base),
                   info = sprintf("permute seed %d", s))
    renamed <- within(d, Actor <- relabel(Actor, s))
    expect_weights(build_event(renamed, buildModel_threshold = thr)$weights, as_plain(base),
                   info = sprintf("relabel seed %d", s))
  }))
})

test_that("one-hot: unit order and actor labels do not matter; full row order only with window 1", {
  invisible(lapply(seq_len(6L), \(s) {
    oh <- rand_onehot(s)
    codes <- attr(oh, "codes")
    build <- function(d, window) build_onehot(d, buildModel_variables_onehot = codes,
                                              buildModel_window = window,
                                              buildModel_window_type = "sliding")$weights
    # reorder whole units (rows within a unit keep their order) and relabel actors
    units <- interaction(oh[c("Actor", "Session")], drop = TRUE)
    lv <- withr::with_seed(s, sample(levels(units)))
    moved <- oh[order(match(as.character(units), lv), seq_len(nrow(oh))), , drop = FALSE]
    moved$Actor <- relabel(moved$Actor, s)
    expect_weights(build(moved, 3L), as_plain(build(oh, 3L)), info = sprintf("units seed %d", s))
    shuffled <- oh[withr::with_seed(s, sample.int(nrow(oh))), , drop = FALSE]
    expect_weights(build(shuffled, 1L), as_plain(build(oh, 1L)), info = sprintf("rows seed %d", s))
  }))
})

test_that("delimited networks do not depend on record order", {
  invisible(lapply(seq_len(10L), \(s) {
    d <- rand_delimited(s)
    sim <- c("none", "jaccard", "cosine", "association", "relative")[1L + s %% 5L]
    base <- build_delim(d, buildModel_sep = attr(d, "sep"), buildModel_similarity = sim)$weights
    shuffled <- d[withr::with_seed(s, sample.int(nrow(d))), , drop = FALSE]
    expect_weights(build_delim(shuffled, buildModel_sep = attr(d, "sep"),
                               buildModel_similarity = sim)$weights,
                   as_plain(base), info = sprintf("seed %d %s", s, sim))
  }))
})

test_that("symmetric, non-negative, and whole counts when unscaled", {
  nets <- c(
    lapply(seq_len(10L), \(s) list(w = build_event(rand_events(s), buildModel_threshold = rand_threshold(s))$weights,
                                   counts = TRUE, label = sprintf("event %d", s))),
    lapply(seq_len(6L), \(s) {
      oh <- rand_onehot(s)
      list(w = build_onehot(oh, buildModel_variables_onehot = attr(oh, "codes"),
                            buildModel_window = 1L + s %% 4L)$weights,
           counts = TRUE, label = sprintf("onehot %d", s))
    }),
    lapply(c("none", "jaccard", "cosine", "inclusion", "association", "dice", "equivalence"), \(sim) {
      d <- rand_delimited(7L)
      list(w = build_delim(d, buildModel_sep = attr(d, "sep"), buildModel_similarity = sim)$weights,
           counts = sim == "none", label = sim)
    }),
    list(list(w = build_event(rand_events(3L), buildModel_scaling = "rank")$weights,
              counts = FALSE, label = "event rank")))
  invisible(lapply(nets, \(n) {
    w <- as_plain(n$w)
    expect_true(isSymmetric(w), info = n$label)
    expect_true(all(w >= 0), info = n$label)
    if (n$counts) expect_equal(w, round(w), info = n$label)
  }))
})

test_that("classed error paths: unreadable Time and a column given two roles", {
  bad_time <- data.frame(Actor = "a", Action = c("A", "B"),
                         Time = c("2021-01-01 00:00:00", "not a time"))
  expect_error(build_event(bad_time), class = "nestimate_bad_input")
  oh <- data.frame(c1 = c(1, 0), c2 = c(0, 1), Actor = c("a", "a"))
  expect_error(build_onehot(oh, buildModel_variables_onehot = c("c1", "c2", "Actor"),
                            buildModel_variables_session = NULL),
               class = "cna_bad_input")
})
