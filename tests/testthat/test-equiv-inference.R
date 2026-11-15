# Equivalence of CNA's inference layer with independent implementations:
# association rules (arules + a base-R enumeration), bootstrap, permutation,
# case-dropping stability and split-half reliability (base-R reimplementations
# that replay Nestimate's draws from the same seed), group comparison
# (closed-form metrics, igraph paths) and clustering (stringdist, cluster::pam,
# stats::hclust, base-R silhouettes).
#
# Conventions the references follow (read from Nestimate 1.0.7, a1754ab):
# * Event data: one resampling unit and one rule transaction per session
#   (actor x time gap), i.e. per row of net$data. A session's co-occurrence is
#   n n' with choose(n, 2) on the diagonal (n = state counts of the session).
# * One-hot data: one unit per actor x session, in the order of the levels of
#   interaction(Actor, Session) (what split() returns); one transaction per
#   window. Clustering takes the units in first-seen order.
# * Draws: bootstrap sample.int(n, n, TRUE) per iteration; permutation
#   sample.int(n_x + n_y, n_x) (paired: sample(c(TRUE, FALSE), n, TRUE));
#   stability / case-drop sample.int(n, n - floor(n * p)) per drop proportion;
#   split-half sample.int(n, round(n / 2)). Each procedure and each group
#   calls set.seed(seed) itself.

skip_if_not_installed("arules")
skip_if_not_installed("stringdist")
skip_if_not_installed("cluster")
skip_if_not_installed("igraph")

# ---------------------------------------------------------------------------
# Independent references
# ---------------------------------------------------------------------------

# Event sessions from the raw long data: a new session starts at a new actor
# or after a gap longer than `threshold` seconds. One transaction per session.
eqi_event_transactions <- function(d, threshold) {
  d <- d[order(d$Actor, d$Time), ]
  gap <- c(Inf, diff(as.numeric(d$Time)))
  new_actor <- c(TRUE, d$Actor[-1L] != d$Actor[-nrow(d)])
  session <- cumsum(new_actor | gap > threshold)
  unname(lapply(split(d$Action, session), \(a) sort(unique(a))))
}

# One-hot transactions: the codes active in each window of each unit.
eqi_onehot_transactions <- function(df, codes, units, window, mode) {
  rows_by_unit <- split(seq_len(nrow(df)), interaction(df[units], drop = TRUE))
  windows <- unlist(lapply(rows_by_unit, \(r) {
    n <- length(r)
    starts <- if (identical(mode, "overlapping")) {
      if (n >= window) seq_len(n - window + 1L) else integer()
    } else seq(1L, n, by = window)
    lapply(starts, \(s) r[s:min(n, s + window - 1L)])
  }), recursive = FALSE)
  X <- as.matrix(df[codes])
  unname(lapply(windows, \(w) codes[colSums(X[w, , drop = FALSE]) > 0]))
}

# Every rule A => B (B may hold several items) by exhaustive enumeration of
# the itemsets: support = n(AB)/N, confidence = n(AB)/n(A),
# lift = N n(AB) / (n(A) n(B)), conviction = (1 - P(B)) / (1 - confidence).
# Thresholds: support and confidence >= their minimum; lift in exact integer
# arithmetic (Agrawal & Srikant, 1994; Brin et al., 1997).
eqi_rules_ref <- function(tx, min_support, min_confidence, min_lift, max_length) {
  items <- sort(unique(unlist(tx)))
  n <- length(tx)
  X <- matrix(vapply(tx, \(t) items %in% t, logical(length(items))),
              ncol = length(items), byrow = TRUE)
  sets <- unlist(lapply(seq_len(min(max_length, length(items))),
                        \(k) utils::combn(items, k, simplify = FALSE)), recursive = FALSE)
  count <- vapply(sets, \(s) sum(rowSums(X[, match(s, items), drop = FALSE]) == length(s)),
                  numeric(1))
  key <- vapply(sets, paste, character(1), collapse = ", ")
  count_of <- stats::setNames(count, key)
  frequent <- which(count / n >= min_support & lengths(sets) >= 2L)
  rules <- do.call(rbind, lapply(frequent, \(i) {
    s <- sets[[i]]
    antecedents <- unlist(lapply(seq_len(length(s) - 1L),
                                 \(k) utils::combn(s, k, simplify = FALSE)), recursive = FALSE)
    data.frame(
      antecedent = vapply(antecedents, paste, character(1), collapse = ", "),
      consequent = vapply(antecedents, \(a) paste(setdiff(s, a), collapse = ", "), character(1)),
      count = count[[i]], stringsAsFactors = FALSE)
  }))
  if (is.null(rules)) return(data.frame())
  n_a <- unname(count_of[rules$antecedent])
  n_b <- unname(count_of[rules$consequent])
  rules$support <- rules$count / n
  rules$confidence <- rules$count / n_a
  rules$lift <- n * rules$count / (n_a * n_b)
  rules$conviction <- ifelse(rules$count == n_a, Inf, (1 - n_b / n) / (1 - rules$confidence))
  subset(rules, confidence >= min_confidence & n * count >= min_lift * n_a * n_b)
}

# arules' Apriori (single-item consequents, no lift filter). arules warns when
# maxlen stops the mining early; that one warning is expected here.
eqi_arules <- function(tx, min_support, min_confidence, max_length) {
  trans <- methods::as(tx, "transactions")
  r <- withCallingHandlers(
    arules::apriori(trans, parameter = list(supp = min_support, conf = min_confidence,
                                            minlen = 2L, maxlen = max_length,
                                            target = "rules"),
                    control = list(verbose = FALSE)),
    warning = function(w) {
      if (grepl("maxlen reached", conditionMessage(w), fixed = TRUE))
        invokeRestart("muffleWarning")
    })
  q <- arules::quality(r)
  side <- function(x) vapply(arules::LIST(x), \(v) paste(sort(v), collapse = ", "), character(1))
  data.frame(antecedent = side(arules::lhs(r)), consequent = side(arules::rhs(r)),
             support = q$support, confidence = q$confidence, lift = q$lift,
             count = q$count, stringsAsFactors = FALSE)
}

eqi_rule_key <- function(df) paste(df$antecedent, "=>", df$consequent)

# Compare a rules table with a reference on the same rule set.
eqi_expect_rules <- function(out, ref, cols, info) {
  expect_setequal(eqi_rule_key(out), eqi_rule_key(ref))
  at <- match(eqi_rule_key(ref), eqi_rule_key(out))
  invisible(lapply(cols, \(col)
    expect_equal(as.numeric(out[[col]][at]), as.numeric(ref[[col]]), info = paste(info, col))))
}

# Per-session co-occurrence counts: one row per session (row of the wide
# sequences), flattened column-major over `states`.
eqi_session_counts <- function(seqs, states) {
  matrix(vapply(seq_len(nrow(seqs)), \(i) {
    v <- unlist(seqs[i, ], use.names = FALSE)
    v <- v[!is.na(v)]
    n <- as.numeric(table(factor(v, levels = states)))
    C <- outer(n, n)
    diag(C) <- choose(n, 2)
    as.vector(C)
  }, numeric(length(states)^2)), ncol = length(states)^2, byrow = TRUE)
}

# The sessions an event network rests on.
eqi_net_counts <- function(net, states = rownames(net$weights))
  eqi_session_counts(net$data, states)

# Per-unit one-hot co-occurrence (tumbling windows): sum over windows of s s',
# s = code totals of the window. Units in interaction(Actor, Session) order.
eqi_onehot_unit_counts <- function(df, codes, units, window) {
  rows_by_unit <- split(seq_len(nrow(df)), interaction(df[units], drop = TRUE))
  X <- as.matrix(df[codes])
  matrix(vapply(rows_by_unit, \(r) {
    starts <- seq(1L, length(r), by = window)
    C <- Reduce(`+`, lapply(starts, \(s)
      tcrossprod(colSums(X[r[s:min(length(r), s + window - 1L)], , drop = FALSE]))))
    as.vector(C)
  }, numeric(length(codes)^2)), ncol = length(codes)^2, byrow = TRUE)
}

eqi_as_matrix <- function(v, states) matrix(v, length(states), dimnames = list(states, states))

# Evaluate `expr` from `seed`, restoring the caller's RNG state afterwards.
eqi_with_seed <- function(seed, expr) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv())
          else rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed)
  expr
}

# Bootstrap of unit counts U (rows = units): iter x cells.
eqi_bootstrap_draws <- function(U, iter, seed) eqi_with_seed(seed,
  matrix(vapply(seq_len(iter), \(i)
    colSums(U[sample.int(nrow(U), nrow(U), replace = TRUE), , drop = FALSE]),
    numeric(ncol(U))), nrow = iter, byrow = TRUE))

# Unpaired / paired permutation of pooled unit counts: iter x cells of x - y.
eqi_permutation_draws <- function(UX, UY, iter, seed, paired = FALSE) {
  pooled <- rbind(UX, UY)
  nx <- nrow(UX)
  nt <- nrow(pooled)
  eqi_with_seed(seed, matrix(vapply(seq_len(iter), \(i) {
    if (paired) {
      swap <- sample(c(TRUE, FALSE), nx, replace = TRUE)
      ix <- ifelse(swap, seq(nx + 1L, nt), seq_len(nx))
      iy <- ifelse(swap, seq_len(nx), seq(nx + 1L, nt))
    } else {
      ix <- sample.int(nt, nx)
      iy <- setdiff(seq_len(nt), ix)
    }
    colSums(pooled[ix, , drop = FALSE]) - colSums(pooled[iy, , drop = FALSE])
  }, numeric(ncol(pooled))), nrow = iter, byrow = TRUE))
}

# CS-coefficient (Epskamp, Borsboom & Fried, 2018): the largest drop
# proportion at which the correlation is >= threshold in at least `certainty`
# of the subsamples; 0 when none qualifies.
eqi_cs <- function(R, props, threshold = 0.7, certainty = 0.95) {
  above <- colMeans(R >= threshold, na.rm = TRUE)
  if (any(above >= certainty, na.rm = TRUE)) max(props[which(above >= certainty)]) else 0
}

# Case-dropping: for each drop proportion, `iter` subsamples of n - floor(n p)
# units; `stat(kept_counts)` gives the per-subsample statistic.
eqi_casedrop <- function(U, iter, seed, props, stat) {
  n <- nrow(U)
  eqi_with_seed(seed, vapply(props, \(p) vapply(seq_len(iter), \(i)
    stat(colSums(U[sample.int(n, n - floor(n * p)), , drop = FALSE])), numeric(1)),
    numeric(iter)))
}

# Silhouette widths (Rousseeuw, 1987); 0 for a singleton cluster.
eqi_silhouette <- function(dm, cl) {
  vapply(seq_along(cl), \(i) {
    own <- cl == cl[i]
    if (sum(own) == 1L) return(0)
    a <- sum(dm[i, own]) / (sum(own) - 1)
    b <- min(vapply(setdiff(unique(cl), cl[i]), \(k) mean(dm[i, cl == k]), numeric(1)))
    (b - a) / max(a, b)
  }, numeric(1))
}

# Mean within-cluster distance (0 for a singleton), per cluster 1..k.
eqi_within <- function(dm, cl, k = max(cl)) {
  vapply(seq_len(k), \(c) {
    m <- dm[cl == c, cl == c, drop = FALSE]
    if (nrow(m) > 1L) mean(m[lower.tri(m)]) else 0
  }, numeric(1))
}

eqi_same_partition <- function(a, b) {
  length(a) == length(b) &&
    nrow(unique(cbind(a, b))) == length(unique(a)) &&
    length(unique(a)) == length(unique(b))
}

# Hamming over the full width: a missing position is a symbol of its own, so
# a shorter sequence differs from a longer one at every extra position.
eqi_hamming <- function(M) {
  n <- nrow(M)
  matrix(vapply(seq_len(n * n), \(ij) {
    a <- M[(ij - 1L) %% n + 1L, ]
    b <- M[(ij - 1L) %/% n + 1L, ]
    sum(xor(is.na(a), is.na(b)) | (!is.na(a) & !is.na(b) & a != b))
  }, numeric(1)), n)
}

# Sequences as strings in an encoding of our own (first-seen order, letters
# from "Q"), truncated after the last observed state.
eqi_strings <- function(M) {
  states <- unique(M[!is.na(M)])
  chars <- intToUtf8(utf8ToInt("Q") + seq_along(states) - 1L, multiple = TRUE)
  vapply(seq_len(nrow(M)), \(i) {
    v <- M[i, ]
    last <- max(c(0L, which(!is.na(v))))
    enc <- ifelse(is.na(v), "#", chars[match(v, states)])
    paste0(enc[seq_len(last)], collapse = "")
  }, character(1))
}

eqi_cluster_opts <- function(dissimilarity = "hamming", method = "pam", k = 3L)
  list(clustering_k = k, clustering_dissimilarity = dissimilarity,
       clustering_method = method, buildModel_scaling = "noScaling",
       seed_set = TRUE, seed_value = 1L)

# One-hot data whose actors each hold two sessions of four rows; with
# `short = TRUE` actor a01's first session keeps a single row (shorter than
# any window above 1).
eqi_onehot_sessions <- function(short = TRUE) {
  d <- make_onehot_data(n = 96)
  d$Actor <- rep(sprintf("a%02d", seq_len(12)), each = 8)
  d$Session <- rep(rep(c(1L, 2L), each = 4), times = 12)
  if (short) d[-c(2L, 3L, 4L), ] else d
}

eqi_codes <- paste0("code", 1:4)

# Event fixture with several sessions per actor (400 s gaps split them).
eqi_event_data <- function() make_event_data(n_actor = 12, per = 10)
eqi_event_net <- function(...) cna("cna_build_event_nest")(
  eqi_event_data(), opts_event(buildModel_threshold = 400, ...))

# ---------------------------------------------------------------------------
# 1. Association rules
# ---------------------------------------------------------------------------

test_that("event rules: one transaction per session; same rule set as arules and a full enumeration", {
  net <- eqi_event_net()
  tx <- eqi_event_transactions(eqi_event_data(), threshold = 400)
  expect_gt(length(tx), length(unique(eqi_event_data()$Actor)))   # sessions, not actors
  cases <- list(c(supp = 0.1, conf = 0.5, len = 5), c(supp = 0.2, conf = 0.3, len = 2),
                c(supp = 0.15, conf = 0.6, len = 3), c(supp = 0.05, conf = 0.2, len = 4))
  invisible(lapply(cases, \(cs) {
    info <- paste(names(cs), cs, collapse = " ")
    # all rules, multi-item consequents included, against the enumeration
    rules_list <- cna("cna_rules_compute")(net, min_support = cs[["supp"]],
                                           min_confidence = cs[["conf"]],
                                           max_length = cs[["len"]])
    rules <- rules_list[[1]]
    expect_equal(rules$n_transactions, length(tx), info = info)
    ref <- eqi_rules_ref(tx, cs[["supp"]], cs[["conf"]], min_lift = 1, cs[["len"]])
    eqi_expect_rules(rules$rules, ref,
                     c("support", "confidence", "lift", "conviction", "count"), info)
    # the table keeps every rule when max_rows allows it
    out <- cna("cna_rules_df")(rules_list, max_rows = 10000L)
    expect_equal(nrow(out), nrow(ref), info = info)
    # single-item consequents, no lift filter, against arules' Apriori
    all_rules <- cna("cna_rules_compute")(net, min_support = cs[["supp"]],
                                          min_confidence = cs[["conf"]], min_lift = 0,
                                          max_length = cs[["len"]])[[1]]$rules
    single <- subset(all_rules, !grepl(",", consequent, fixed = TRUE))
    eqi_expect_rules(single, eqi_arules(tx, cs[["supp"]], cs[["conf"]], cs[["len"]]),
                     c("support", "confidence", "lift", "count"), info)
  }))
})

test_that("event rules per group come from each group's own sessions", {
  d <- eqi_event_data()
  net <- eqi_event_net(buildModel_variables_group = "Group")
  rules <- cna("cna_rules_compute")(net, min_support = 0.1, min_confidence = 0.3, min_lift = 0)
  expect_equal(names(rules), names(net))
  invisible(Map(\(r, g) {
    tx <- eqi_event_transactions(subset(d, Group == g), threshold = 400)
    expect_equal(r$n_transactions, length(tx), info = g)
    single <- subset(r$rules, !grepl(",", consequent, fixed = TRUE))
    eqi_expect_rules(single, eqi_arules(tx, 0.1, 0.3, 5L),
                     c("support", "confidence", "lift", "count"), g)
  }, rules, names(rules)))
})

test_that("one-hot rules: one transaction per window of each actor x session (arules + enumeration)", {
  d <- eqi_onehot_sessions()
  cases <- list(list(1L, "tumbling"), list(2L, "tumbling"), list(3L, "tumbling"),
                list(2L, "sliding"), list(3L, "sliding"))
  invisible(lapply(cases, \(case) {
    info <- paste(case, collapse = " ")
    mode <- if (case[[2]] == "sliding") "overlapping" else "non-overlapping"
    net <- cna("cna_build_onehot_nest")(d, opts_onehot(buildModel_window = case[[1]],
                                                      buildModel_window_type = case[[2]]))
    tx <- eqi_onehot_transactions(d, eqi_codes, c("Actor", "Session"), case[[1]], mode)
    rules <- cna("cna_rules_compute")(net, min_support = 0.2, min_confidence = 0.5,
                                      min_lift = 0, max_length = 4L)[[1]]
    expect_equal(rules$n_transactions, length(tx), info = info)
    expect_setequal(rules$items, eqi_codes)
    eqi_expect_rules(rules$rules, eqi_rules_ref(tx, 0.2, 0.5, min_lift = 0, 4L),
                     c("support", "confidence", "lift", "conviction", "count"), info)
    single <- subset(rules$rules, !grepl(",", consequent, fixed = TRUE))
    eqi_expect_rules(single, eqi_arules(tx, 0.2, 0.5, 4L),
                     c("support", "confidence", "lift", "count"), info)
  }))
})

test_that("invariant: shuffling the input rows leaves the rules unchanged", {
  d <- eqi_event_data()
  shuffled <- d[eqi_with_seed(31L, sample.int(nrow(d))), ]
  f <- function(x) cna("cna_rules_df")(cna("cna_rules_compute")(
    cna("cna_build_event_nest")(x, opts_event(buildModel_threshold = 400)),
    min_support = 0.1, min_confidence = 0.3, max_length = 4L), max_rows = 10000L)
  a <- f(d)
  b <- f(shuffled)
  expect_gt(nrow(a), 0)
  expect_identical(eqi_rule_key(a), eqi_rule_key(b))
  expect_equal(b, a)
})

# ---------------------------------------------------------------------------
# 2. Bootstrap
# ---------------------------------------------------------------------------

test_that("event bootstrap == independent resampling of sessions (same seed, same draws)", {
  net <- eqi_event_net()
  states <- rownames(net$weights)
  U <- eqi_net_counts(net)
  w <- unclass(net$weights)
  expect_equal(eqi_as_matrix(colSums(U), states), w)        # the units sum to the network
  iter <- 40L
  invisible(lapply(list(list(level = 0.05, seed = 7L), list(level = 0.1, seed = 8L)), \(cs) {
    B <- eqi_bootstrap_draws(U, iter, cs$seed)
    bs <- cna("cna_bootstrap_compute")(net, iter = iter, level = cs$level, seed = cs$seed)
    lo <- eqi_as_matrix(apply(B, 2, stats::quantile, probs = cs$level / 2, names = FALSE), states)
    hi <- eqi_as_matrix(apply(B, 2, stats::quantile, probs = 1 - cs$level / 2, names = FALSE), states)
    # stability inference: share of resamples outside [0.75 w, 1.25 w] (edges inclusive)
    outside <- B <= rep(0.75 * as.vector(w), each = iter) | B >= rep(1.25 * as.vector(w), each = iter)
    p <- eqi_as_matrix((colSums(outside) + 1) / (iter + 1), states)
    expect_equal(unclass(bs$mean), eqi_as_matrix(colMeans(B), states))
    expect_equal(unclass(bs$sd), eqi_as_matrix(apply(B, 2, stats::sd), states))
    expect_equal(unclass(bs$p_values), p)
    df <- cna("cna_bootstrap_df")(bs, multi = FALSE, show_all = TRUE)
    at <- cbind(df$from, df$to)
    # one row per unique edge (upper triangle, diagonal included, non-zero)
    expect_equal(nrow(df), sum(w[upper.tri(w, diag = TRUE)] != 0))
    expect_equal(df$p_value, p[at])
    expect_equal(df$ci_lower, lo[at])
    expect_equal(df$ci_upper, hi[at])
    expect_equal(df$cr_lower, 0.75 * w[at])
    expect_equal(df$cr_upper, 1.25 * w[at])
    expect_identical(df$significant, ifelse(p[at] < cs$level, "Yes", "No"))
    # basic interval: the percentile bounds reflected around the estimate
    basic <- cna("cna_bootstrap_compute")(net, iter = iter, level = cs$level,
                                          ci_method = "basic", seed = cs$seed)
    expect_equal(unclass(basic$ci_lower), 2 * w - hi)
    expect_equal(unclass(basic$ci_upper), 2 * w - lo)
    expect_equal(unclass(basic$p_values), p)                 # inference unchanged
  }))
})

test_that("one-hot bootstrap resamples actor x session units (independent reimplementation)", {
  d <- eqi_onehot_sessions()
  net <- cna("cna_build_onehot_nest")(d, opts_onehot(buildModel_window = 2L))
  states <- rownames(net$weights)
  U <- eqi_onehot_unit_counts(d, eqi_codes, c("Actor", "Session"), 2L)
  expect_equal(nrow(U), 24L)
  expect_equal(eqi_as_matrix(colSums(U), states), unclass(net$weights))
  B <- eqi_bootstrap_draws(U, 30L, 3L)
  bs <- cna("cna_bootstrap_compute")(net, iter = 30L, level = 0.05, seed = 3L)
  expect_equal(unclass(bs$mean), eqi_as_matrix(colMeans(B), states))
  expect_equal(unclass(bs$sd), eqi_as_matrix(apply(B, 2, stats::sd), states))
  expect_equal(unclass(bs$ci_upper),
               eqi_as_matrix(apply(B, 2, stats::quantile, probs = 0.975, names = FALSE), states))
})

test_that("grouped bootstrap: each group resampled from its own sessions with the seed", {
  net <- eqi_event_net(buildModel_variables_group = "Group")
  bs <- cna("cna_bootstrap_compute")(net, iter = 30L, level = 0.05, seed = 12L)
  df <- cna("cna_bootstrap_df")(bs, multi = TRUE, show_all = TRUE)
  invisible(Map(\(g_net, g) {
    states <- rownames(g_net$weights)
    w <- unclass(g_net$weights)
    B <- eqi_bootstrap_draws(eqi_net_counts(g_net), 30L, 12L)
    outside <- B <= rep(0.75 * as.vector(w), each = 30L) | B >= rep(1.25 * as.vector(w), each = 30L)
    p <- eqi_as_matrix((colSums(outside) + 1) / 31, states)
    rows <- subset(df, group == g)
    expect_equal(rows$p_value, p[cbind(rows$from, rows$to)], info = g)
    expect_equal(unclass(bs[[g]]$mean), eqi_as_matrix(colMeans(B), states), info = g)
  }, cna(".cna_nets")(net), names(net)))
})

test_that("threshold inference: p = share of resamples below the edge threshold", {
  net <- eqi_event_net()
  states <- rownames(net$weights)
  w <- unclass(net$weights)
  B <- eqi_bootstrap_draws(eqi_net_counts(net), 40L, 5L)
  # automatic threshold: the 10th percentile of the non-zero |weights|
  auto <- cna("cna_bootstrap_compute")(net, iter = 40L, level = 0.05, inference = "threshold",
                                       seed = 5L)
  thr <- stats::quantile(abs(w[w != 0]), 0.1, names = FALSE)
  expect_equal(unname(auto$edge_threshold), thr)
  expect_equal(unclass(auto$p_values),
               eqi_as_matrix((colSums(abs(B) < thr) + 1) / 41, states))
  # a user threshold
  fixed <- cna("cna_bootstrap_compute")(net, iter = 40L, level = 0.05, inference = "threshold",
                                        threshold = 15, seed = 5L)
  p15 <- eqi_as_matrix((colSums(abs(B) < 15) + 1) / 41, states)
  expect_equal(unclass(fixed$p_values), p15)
  expect_equal(unclass(fixed$significant), w * (p15 < 0.05))
})

test_that("bootstrap table works with threshold inference (one row per edge, its p-value)", {
  # Fixed 2026-11-06, pinned here. Was (CNA 1.1.0): cna_bootstrap_df() reads s$cr_lower / s$cr_upper,
  # which Nestimate omits for threshold inference, so data.frame() fails and
  # the jamovi section reports "This analysis could not be completed".
  res <- run_event_single(eqi_event_data(), buildModel_threshold = 400,
                          bootstrap_show_table = TRUE, bootstrap_iteration = 20L,
                          bootstrap_method = "threshold", seed_set = TRUE)
  expect_length(section_failures(res), 0)
  expect_gt(res$bootstrapTable$rowCount, 0)
  net <- eqi_event_net()
  bt <- cna("cna_bootstrap_compute")(net, iter = 20L, level = 0.05, inference = "threshold",
                                     seed = 2L)
  df <- cna("cna_bootstrap_df")(bt, multi = FALSE, show_all = TRUE)
  expect_equal(nrow(df), nrow(bt$summary))
  expect_equal(df$p_value, unclass(bt$p_values)[cbind(df$from, df$to)])
})

# ---------------------------------------------------------------------------
# 3. Permutation
# ---------------------------------------------------------------------------

test_that("unpaired permutation == independent shuffling of the sessions' group labels", {
  net <- eqi_event_net(buildModel_variables_group = "Group")
  states <- rownames(net$G1$weights)
  UX <- eqi_net_counts(net$G1, states)
  UY <- eqi_net_counts(net$G2, states)
  expect_false(nrow(UX) == nrow(UY))                       # unequal groups
  obs <- colSums(UX) - colSums(UY)
  D <- eqi_permutation_draws(UX, UY, 50L, 9L)
  p <- eqi_as_matrix((colSums(abs(D) >= rep(abs(obs), each = 50L)) + 1) / 51, states)
  null_sd <- sqrt(pmax(colMeans(D^2) - colMeans(D)^2, 0))     # population SD of the null
  es <- eqi_as_matrix(ifelse(null_sd > 0, obs / null_sd, 0), states)
  pt <- cna("cna_permutation_compute")(net, iter = 50L, seed = 9L)
  expect_identical(names(pt), "G1 vs G2")
  cmp <- pt[["G1 vs G2"]]
  expect_equal(unclass(cmp$p_values), p)
  expect_equal(unclass(cmp$effect_size), es)
  s <- cmp$summary
  at <- cbind(s$from, s$to)
  wx <- unclass(net$G1$weights)
  wy <- unclass(net$G2$weights)
  expect_equal(s$diff, (wx - wy)[at])
  expect_equal(s$weight_x, wx[at])
  expect_equal(s$weight_y, wy[at])
  # unique pairs present in either network, diagonal included
  present <- (wx != 0 | wy != 0) & upper.tri(wx, diag = TRUE)
  expect_equal(nrow(s), sum(present))
  df <- cna("cna_permutation_df")(pt, show_all = TRUE)
  expect_equal(sort(df$p_value), sort(p[at]))
})

test_that("paired permutation == independent within-pair swaps; unequal pairs are refused", {
  net <- fixture_group_net()
  states <- rownames(net$G1$weights)
  UX <- eqi_net_counts(net$G1, states)
  UY <- eqi_net_counts(net$G2, states)
  expect_equal(nrow(UX), nrow(UY))
  obs <- colSums(UX) - colSums(UY)
  D <- eqi_permutation_draws(UX, UY, 50L, 4L, paired = TRUE)
  pt <- cna("cna_permutation_compute")(net, iter = 50L, paired = TRUE, seed = 4L)
  expect_equal(unclass(pt[[1]]$p_values),
               eqi_as_matrix((colSums(abs(D) >= rep(abs(obs), each = 50L)) + 1) / 51, states))
  expect_error(cna("cna_permutation_compute")(eqi_event_net(buildModel_variables_group = "Group"),
                                              iter = 10L, paired = TRUE, seed = 1L),
               class = "nestimate_bad_input")
})

test_that("p-value correction == p.adjust over the unique edges of each comparison, every method", {
  # three groups; G1 is planted to favour A, so some p-values are small
  d <- make_event_data(n_actor = 15, per = 10)
  d$Group <- rep(c("G1", "G2", "G3"), each = 50)
  d$Action[d$Group == "G1"] <- eqi_with_seed(17L, sample(
    LETTERS[1:5], 50, replace = TRUE, prob = c(0.8, 0.05, 0.05, 0.05, 0.05)))
  net <- cna("cna_build_event_nest")(d, opts_event(buildModel_threshold = 400,
                                                   buildModel_variables_group = "Group"))
  expect_true(all(vapply(net, \(n) identical(rownames(n$weights), LETTERS[1:5]), logical(1))))
  raw <- cna("cna_permutation_compute")(net, iter = 60L, seed = 6L)
  expect_identical(names(raw), c("G1 vs G2", "G1 vs G3", "G2 vs G3"))
  expect_lt(min(raw[["G1 vs G2"]]$summary$p_value), 0.05)
  # every pair is its own two-group test drawn from the same seed
  states <- rownames(net$G1$weights)
  UX <- eqi_net_counts(net$G1, states)
  UZ <- eqi_net_counts(net$G3, states)
  D <- eqi_permutation_draws(UX, UZ, 60L, 6L)
  obs <- colSums(UX) - colSums(UZ)
  expect_equal(unclass(raw[["G1 vs G3"]]$p_values),
               eqi_as_matrix((colSums(abs(D) >= rep(abs(obs), each = 60L)) + 1) / 61, states))
  methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr")
  invisible(lapply(methods, \(m) {
    adj <- cna("cna_adjust_permutation")(raw, m, 0.05)
    invisible(Map(\(a, r, nm) {
      info <- paste(m, nm)
      ref <- stats::p.adjust(r$summary$p_value, method = m)
      expect_equal(a$summary$p_value, ref, info = info)
      expect_identical(a$summary$sig, ref < 0.05, info = info)
      # the matrix the plot reads carries the same values, symmetric
      at <- cbind(a$summary$from, a$summary$to)
      expect_equal(unclass(a$p_values)[at], ref, info = info)
      expect_equal(unclass(a$p_values)[at[, c(2L, 1L)]], ref, info = info)
      expect_equal(a$diff_sig, a$diff * (a$p_values < 0.05), info = info)
    }, adj, raw, names(raw)))
  }))
  expect_identical(cna("cna_adjust_permutation")(raw, "none", 0.05), raw)
  # Nestimate's own `adjust` corrects over all n x n cells (known, see LEARNINGS),
  # counting each undirected edge twice: with Holm the families differ.
  nest <- Nestimate::permutation(net, iter = 60L, adjust = "holm", seed = 6L)[["G1 vs G2"]]
  cells <- stats::p.adjust(as.vector(unclass(raw[["G1 vs G2"]]$p_values)), method = "holm")
  s <- nest$summary
  expect_equal(s$p_value, eqi_as_matrix(cells, states)[cbind(s$from, s$to)])
  expect_false(isTRUE(all.equal(
    s$p_value, cna("cna_adjust_permutation")(raw, "holm", 0.05)[["G1 vs G2"]]$summary$p_value)))
})

test_that("invariant: reversing the group order swaps x and y and negates the differences", {
  # groups follow the Group's factor levels (otherwise sorted), not the row
  # order (2026-11-06): reversing the levels reverses the comparison
  d <- eqi_event_data()
  reversed <- transform(d, Group = factor(Group, levels = c("G2", "G1")))
  f <- function(x, pair) cna("cna_permutation_compute")(cna("cna_build_event_nest")(
    x, opts_event(buildModel_threshold = 400, buildModel_variables_group = "Group")),
    iter = 20L, seed = 1L)[[pair]]$summary
  a <- f(d, "G1 vs G2")
  b <- f(reversed, "G2 vs G1")
  expect_identical(names(cna("cna_build_event_nest")(d[order(d$Group != "G2"), ],
    opts_event(buildModel_threshold = 400, buildModel_variables_group = "Group"))),
    c("G1", "G2"))
  expect_identical(paste(a$from, a$to), paste(b$from, b$to))
  expect_equal(b$diff, -a$diff)
  expect_equal(b$weight_x, a$weight_y)
  expect_equal(b$weight_y, a$weight_x)
})

# ---------------------------------------------------------------------------
# 4. Stability and reliability
# ---------------------------------------------------------------------------

test_that("CS-coefficients == independent case-dropping (pearson, spearman, kendall)", {
  net <- eqi_event_net()
  states <- rownames(net$weights)
  k <- length(states)
  U <- eqi_net_counts(net)
  props <- seq(0.1, 0.9, by = 0.1)
  strength <- function(v) { m <- matrix(v, k); diag(m) <- 0; colSums(m) }  # loops = FALSE
  w <- unclass(net$weights)
  edge_cells <- w != 0 & row(w) != col(w)          # existing off-diagonal edges
  orig_strength <- strength(colSums(U))
  invisible(lapply(c("pearson", "spearman", "kendall"), \(m) {
    st <- cna("cna_stability_compute")(net, "Strength", iter = 50L, threshold = 0.7,
                                       certainty = 0.95, method = m, seed = 3L)
    R <- eqi_casedrop(U, 50L, 3L, props, \(v) {
      s <- strength(v)
      if (stats::sd(s) == 0) NA_real_ else stats::cor(s, orig_strength, method = m)
    })
    expect_equal(unname(st$cs$correlations$Strength), R, info = m)
    expect_equal(unname(st$cs$cs[["Strength"]]), eqi_cs(R, props), info = m)
    # edge-weight case-drop: correlation of the existing edges with the original
    Rcd <- eqi_casedrop(U, 50L, 3L, props, \(v)
      stats::cor(v[edge_cells], w[edge_cells], method = m))
    expect_equal(unname(st$cd$correlations), Rcd, info = m)
    expect_equal(st$cd$cs, eqi_cs(Rcd, props), info = m)
    df <- cna("cna_stability_df")(st, multi = FALSE)
    expect_equal(df$measure, c("Strength", "Edge weights"), info = m)
    expect_equal(df$cs, c(eqi_cs(R, props), eqi_cs(Rcd, props)), info = m)
  }))
})

test_that("CS of a cograph-only measure follows the definition from its own correlations", {
  net <- eqi_event_net()
  st <- cna("cna_stability_compute")(net, c("Strength", "Eigenvector"), iter = 50L,
                                     threshold = 0.5, certainty = 0.9, seed = 11L)
  props <- st$cs$drop_prop
  cs <- vapply(st$cs$correlations, \(R) eqi_cs(R, props, 0.5, 0.9), numeric(1))
  expect_equal(unname(st$cs$cs), unname(cs))
  expect_equal(utils::head(cna("cna_stability_df")(st, multi = FALSE)$cs, 2L), unname(cs))
})

test_that("split-half reliability == independent random halves of the sessions", {
  net <- eqi_event_net()
  states <- rownames(net$weights)
  U <- eqi_net_counts(net)
  n <- nrow(U)
  metrics <- eqi_with_seed(3L, matrix(vapply(seq_len(50L), \(i) {
    a <- sample.int(n, round(n / 2))
    d <- colSums(U[a, , drop = FALSE]) - colSums(U[-a, , drop = FALSE])
    c(mean(abs(d)), stats::median(abs(d)),
      stats::cor(colSums(U[a, , drop = FALSE]), colSums(U[-a, , drop = FALSE])), max(abs(d)))
  }, numeric(4)), nrow = 4))
  st <- cna("cna_stability_compute")(net, "Strength", iter = 50L, threshold = 0.7,
                                     certainty = 0.95, seed = 3L)
  rel <- cna("cna_reliability_df")(st)
  expect_equal(rel$metric, c("mean_dev", "median_dev", "cor", "max_dev"))
  expect_equal(rel$mean, rowMeans(metrics))
  expect_equal(rel$sd, apply(metrics, 1, stats::sd))
  # grouped: one block of rows per group, labelled by the group
  gnet <- eqi_event_net(buildModel_variables_group = "Group")
  gst <- cna("cna_stability_compute")(gnet, "Strength", iter = 50L, threshold = 0.7,
                                      certainty = 0.95, seed = 3L)
  expect_equal(unique(cna("cna_reliability_df")(gst)$group), c("G1", "G2"))
})

test_that("reliability table of a single network uses the empty group label", {
  # Fixed 2026-11-06, pinned here. Was (cosmetic, CNA 1.1.0): the *_df() convention is group = "" for one
  # network; Nestimate labels an unnamed netobject by its method, so
  # cna_reliability_df() gives "co_occurrence" (the column is hidden in jamovi).
  st <- cna("cna_stability_compute")(eqi_event_net(), "Strength", iter = 50L,
                                     threshold = 0.7, certainty = 0.95, seed = 3L)
  expect_equal(unique(cna("cna_reliability_df")(st)$group),
               unique(cna("cna_stability_df")(st, multi = FALSE)$group))
})

# ---------------------------------------------------------------------------
# 5. Comparing two networks
# ---------------------------------------------------------------------------

eqi_scale <- function(W, scaling) {
  v <- as.vector(W)
  out <- switch(scaling,
    none = v,
    minmax = if (max(v) > min(v)) (v - min(v)) / (max(v) - min(v)) else v,
    max = if (max(v) != 0) v / max(v) else v,
    rank = { r <- rank(v); if (max(r) > min(r)) (r - min(r)) / (max(r) - min(r)) else r })
  matrix(out, nrow(W), dimnames = dimnames(W))
}

# Distance correlation as the V-statistic S1 + S2 - 2 S3 (Szekely, Rizzo &
# Bakirov, 2007, Theorem 1), not the double-centred form.
eqi_dcor <- function(x, y) {
  ax <- abs(outer(x, x, "-"))
  ay <- abs(outer(y, y, "-"))
  dcov2 <- function(a, b) mean(a * b) + mean(a) * mean(b) - 2 * mean(rowMeans(a) * rowMeans(b))
  sqrt(max(dcov2(ax, ay), 0) / sqrt(dcov2(ax, ax) * dcov2(ay, ay)))
}

# RV coefficient (Robert & Escoufier, 1976) as ||X'Y||^2 / (||X'X|| ||Y'Y||).
eqi_rv <- function(X, Y) {
  X <- sweep(X, 2, colMeans(X))
  Y <- sweep(Y, 2, colMeans(Y))
  sum(crossprod(X, Y)^2) / sqrt(sum(crossprod(X)^2) * sum(crossprod(Y)^2))
}

eqi_compare_metrics <- function(X, Y) {
  x <- as.vector(X); y <- as.vector(Y); d <- x - y; n <- nrow(X)
  either <- abs(x) + abs(y) > 0
  c("Mean Abs. Diff." = mean(abs(d)), "Median Abs. Diff." = stats::median(abs(d)),
    "RMS Diff." = sqrt(mean(d^2)), "Max Abs. Diff." = max(abs(d)),
    "Rel. Mean Abs. Diff." = mean(abs(d)) / mean(abs(y)),
    "CV Ratio" = (stats::sd(x) / mean(x)) / (stats::sd(y) / mean(y)),
    "Pearson" = stats::cor(x, y), "Spearman" = stats::cor(rank(x), rank(y)),
    "Kendall" = stats::cor(x, y, method = "kendall"), "Distance" = eqi_dcor(x, y),
    "Euclidean" = sqrt(sum(d^2)), "Manhattan" = sum(abs(d)),
    "Canberra" = sum(abs(d)[either] / (abs(x) + abs(y))[either]),
    "Bray-Curtis" = sum(abs(d)) / sum(abs(x) + abs(y)),
    "Frobenius" = sqrt(sum(d^2)) / sqrt(n / 2),   # tna's normalisation
    "Cosine" = sum(x * y) / sqrt(sum(x^2) * sum(y^2)),
    "Jaccard" = sum(pmin(abs(x), abs(y))) / sum(pmax(abs(x), abs(y))),
    "Dice" = 2 * sum(pmin(abs(x), abs(y))) / (sum(abs(x)) + sum(abs(y))),
    "Overlap" = sum(pmin(abs(x), abs(y))) / min(sum(abs(x)), sum(abs(y))),
    "RV" = eqi_rv(X, Y),
    "Rank Agreement" = mean(sign(X[-1L, ] - X[-n, ]) == sign(Y[-1L, ] - Y[-n, ])),
    "Sign Agreement" = mean(sign(x) == sign(y)))
}

# Side-by-side network metrics from first principles (igraph for paths and
# degree centralization; weights are path lengths, as in tna's summary).
eqi_network_metrics <- function(W) {
  g <- igraph::graph_from_adjacency_matrix(abs(W), mode = "undirected", weighted = TRUE,
                                           diag = FALSE)
  D <- igraph::distances(g)
  off <- D[row(D) != col(D) & is.finite(D)]
  strength <- rowSums(abs(W))
  degree <- rowSums(W != 0)
  c("Node Count" = nrow(W), "Edge Count" = sum(W[upper.tri(W, diag = TRUE)] != 0),
    "Mean Distance" = mean(off), "Mean Out-Strength" = mean(strength),
    "SD Out-Strength" = stats::sd(strength), "Mean In-Strength" = mean(strength),
    "SD In-Strength" = stats::sd(strength), "Mean Out-Degree" = mean(degree),
    "SD Out-Degree" = stats::sd(degree),
    "Centralization (Out-Degree)" = igraph::centr_degree(g, loops = FALSE)$centralization,
    "Reciprocity" = 1)
}

test_that("every comparison metric == an independent formula, under each compare scaling", {
  net <- fixture_group_net()
  invisible(lapply(c("none", "minmax", "max", "rank"), \(sc) {
    cmp <- cna("cna_compare_compute")(net, 1, 2, scaling = sc)
    X <- eqi_scale(unclass(net$G1$weights), sc)
    Y <- eqi_scale(unclass(net$G2$weights), sc)
    ref <- eqi_compare_metrics(X, Y)
    s <- cna("cna_compare_summary_df")(cmp)
    expect_identical(s$metric, names(ref), info = sc)
    expect_equal(s$value, unname(ref), info = sc)
    e <- cna("cna_compare_edges_df")(cmp, max_rows = 1000L)
    at <- cbind(e$from, e$to)
    expect_equal(nrow(e), nrow(X) * (nrow(X) + 1) / 2, info = sc)
    expect_equal(e$weight_i, X[at], info = sc)
    expect_equal(e$weight_j, Y[at], info = sc)
    expect_equal(e$difference, X[at] - Y[at], info = sc)
    expect_equal(e$relative_difference, abs(X[at] - Y[at]) / (X[at] + Y[at]), info = sc)
    nm <- cna("cna_compare_network_df")(cmp)
    rx <- eqi_network_metrics(X)
    ry <- eqi_network_metrics(Y)
    at_m <- match(names(rx), nm$metric)
    expect_false(anyNA(at_m))
    expect_equal(nm$group_i[at_m], unname(rx), info = sc)
    expect_equal(nm$group_j[at_m], unname(ry), info = sc)
  }))
})

test_that("invariant: swapping the two networks negates differences, keeps symmetric metrics", {
  net <- fixture_group_net()
  ab <- cna("cna_compare_compute")(net, 1, 2)
  ba <- cna("cna_compare_compute")(net, 2, 1)
  s_ab <- cna("cna_compare_summary_df")(ab)
  s_ba <- cna("cna_compare_summary_df")(ba)
  symmetric <- setdiff(s_ab$metric, c("Rel. Mean Abs. Diff.", "CV Ratio"))
  at <- match(symmetric, s_ab$metric)
  expect_equal(s_ba$value[at], s_ab$value[at])
  cv <- match("CV Ratio", s_ab$metric)
  expect_equal(s_ba$value[cv], 1 / s_ab$value[cv])
  e_ab <- cna("cna_compare_edges_df")(ab, max_rows = 1000L)
  e_ba <- cna("cna_compare_edges_df")(ba, max_rows = 1000L)
  key <- match(paste(e_ab$from, e_ab$to), paste(e_ba$from, e_ba$to))
  expect_equal(e_ba$difference[key], -e_ab$difference)
  expect_equal(e_ba$weight_i[key], e_ab$weight_j)
})

# ---------------------------------------------------------------------------
# 6. Clustering
# ---------------------------------------------------------------------------

eqi_event_seqs <- function() cna("cna_event_sequences")(
  make_event_data(n_actor = 14, per = 12), opts_event(buildModel_threshold = 600))

test_that("event dissimilarities == stringdist / base R on independently encoded sequences", {
  seqs <- eqi_event_seqs()
  M <- as.matrix(seqs)
  lengths_obs <- rowSums(!is.na(M))
  expect_gt(length(unique(lengths_obs)), 1L)           # unequal lengths are exercised
  strs <- eqi_strings(M)
  # sequences shorter than q = 2 have no q-grams: stringdist gives NaN against
  # a sequence that has some; Nestimate reads that as distance 1
  qgram_ref <- function(method) {
    m <- as.matrix(stringdist::stringdistmatrix(strs, method = method, q = 2))
    m[is.na(m)] <- 1
    m
  }
  refs <- list(
    hamming = eqi_hamming(M),
    osa = as.matrix(stringdist::stringdistmatrix(strs, method = "osa")),
    lv = as.matrix(stringdist::stringdistmatrix(strs, method = "lv")),
    dl = as.matrix(stringdist::stringdistmatrix(strs, method = "dl")),
    lcs = as.matrix(stringdist::stringdistmatrix(strs, method = "lcs")),
    qgram = qgram_ref("qgram"), cosine = qgram_ref("cosine"), jaccard = qgram_ref("jaccard"),
    jw = as.matrix(stringdist::stringdistmatrix(strs, method = "jw", p = 0.1)))
  invisible(Map(\(ref, diss) {
    cl <- cna(".cna_clusters")(seqs, eqi_cluster_opts(diss))
    expect_equal(unname(as.matrix(cl$distance)), unname(ref), info = diss)
    # PAM on the reference distances gives the same clusters
    expect_true(eqi_same_partition(
      cl$assignments, cluster::pam(stats::as.dist(ref), k = 3, diss = TRUE)$clustering),
      info = diss)
  }, refs, names(refs)))
  # Hamming of a shorter sequence counts every extra position of the longer one
  short <- which.min(lengths_obs)
  long <- which.max(lengths_obs)
  same_prefix <- sum(M[short, seq_len(lengths_obs[short])] != M[long, seq_len(lengths_obs[short])])
  expect_equal(refs$hamming[short, long], same_prefix + lengths_obs[long] - lengths_obs[short])
})

test_that("assignments == cluster::pam / hclust + cutree for every method; networks follow them", {
  seqs <- eqi_event_seqs()
  D <- stats::as.dist(eqi_hamming(as.matrix(seqs)))
  states <- sort(unique(stats::na.omit(unlist(seqs, use.names = FALSE))))
  U <- eqi_session_counts(seqs, states)
  methods <- c("pam", "ward.D", "ward.D2", "complete", "average", "single", "mcquitty",
               "median", "centroid")
  invisible(lapply(methods, \(m) {
    ref <- if (m == "pam") cluster::pam(D, k = 3, diss = TRUE)$clustering
           else stats::cutree(stats::hclust(D, method = m), k = 3)
    net <- cna("cna_cluster_networks")(seqs, eqi_cluster_opts(method = m))
    cl <- attr(net, "clustering")$assignments
    expect_true(eqi_same_partition(cl, ref), info = m)
    expect_equal(names(net), paste("Cluster", 1:3), info = m)
    # network "Cluster c" = co-occurrence of the sessions in that cluster
    invisible(Map(\(cnet, c) {
      W <- eqi_as_matrix(colSums(U[cl == c, , drop = FALSE]), states)
      rn <- rownames(cnet$weights)
      expect_equal(unclass(cnet$weights), W[rn, rn], info = paste(m, c))
    }, cna(".cna_nets")(net), seq_len(3L)))
  }))
})

test_that("cluster diagnostics: sizes, within-cluster distance, silhouettes == independent", {
  seqs <- eqi_event_seqs()
  dm <- eqi_hamming(as.matrix(seqs))
  invisible(lapply(c("pam", "average"), \(m) {
    net <- cna("cna_cluster_networks")(seqs, eqi_cluster_opts(method = m))
    cl <- attr(net, "clustering")$assignments
    dg <- cna("cna_cluster_diagnostics_df")(net)
    sil <- eqi_silhouette(dm, cl)
    expect_equal(sil, unname(cluster::silhouette(cl, stats::as.dist(dm))[, "sil_width"]), info = m)
    expect_equal(dg$cluster, paste("Cluster", 1:3), info = m)
    expect_equal(dg$size, tabulate(cl, 3L), info = m)
    expect_equal(dg$pct, 100 * tabulate(cl, 3L) / length(cl), info = m)
    expect_equal(dg$mean_within_dist, eqi_within(dm, cl), info = m)
    expect_equal(dg$silhouette, vapply(seq_len(3L), \(c) mean(sil[cl == c]), numeric(1)), info = m)
    expect_equal(attr(net, "clustering")$silhouette, mean(sil), info = m)
  }))
})

test_that("choosing k == independent fits for each k (event and one-hot units)", {
  choice_ref <- function(dm, ks, method) {
    D <- stats::as.dist(dm)
    do.call(rbind, lapply(ks, \(k) {
      cl <- if (method == "pam") cluster::pam(D, k = k, diss = TRUE)$clustering
            else stats::cutree(stats::hclust(D, method = method), k = k)
      sizes <- tabulate(cl, k)
      data.frame(k = k, silhouette = mean(eqi_silhouette(dm, cl)),
                 mean_within_dist = sum(sizes * eqi_within(dm, cl, k)) / sum(sizes),
                 min_size = min(sizes), max_size = max(sizes),
                 size_ratio = max(sizes) / min(sizes))
    }))
  }
  seqs <- eqi_event_seqs()
  dm <- eqi_hamming(as.matrix(seqs))
  invisible(lapply(c("pam", "ward.D2"), \(m) {
    cc <- cna("cna_cluster_choice_compute")(seqs, 2L, 5L, method = m, seed = 1L)
    expect_equal(cna("cna_cluster_choice_df")(cc), choice_ref(dm, 2:5, m), info = m)
  }))
  # one-hot: the units clustered are the actor x session windowed sequences
  d <- eqi_onehot_sessions()
  units <- cna("cna_onehot_units")(d, opts_onehot(buildModel_window = 2L))
  cc <- cna("cna_cluster_choice_compute")(units, 2L, 4L, seed = 1L)
  key <- paste(d$Actor, d$Session)
  X <- as.matrix(d[eqi_codes])
  flat <- lapply(unique(key), \(u) as.vector(t(X[key == u, , drop = FALSE])))
  L <- max(lengths(flat))
  F <- t(vapply(flat, \(f) c(f, numeric(L - length(f))), numeric(L)))
  expect_equal(cna("cna_cluster_choice_df")(cc),
               choice_ref(as.matrix(stats::dist(F, method = "manhattan")), 2:4, "pam"))
})

# The presence rows a unit is clustered on: its rows, or with sliding windows
# one row per complete window marking the codes active anywhere in it.
eqi_presence <- function(X, rows, window, sliding) {
  M <- X[rows, , drop = FALSE]
  if (!sliding || window <= 1L) return(M)
  if (nrow(M) < window) return(M[0, , drop = FALSE])
  t(vapply(seq_len(nrow(M) - window + 1L), \(s)
    as.numeric(colSums(M[s:(s + window - 1L), , drop = FALSE]) > 0), numeric(ncol(M))))
}

# L1 distances between units' flattened presence rows, zero-padded.
eqi_unit_l1 <- function(d, units, key, window, sliding) {
  X <- as.matrix(d[eqi_codes])
  flat <- lapply(units, \(u) as.vector(t(eqi_presence(X, which(key == u), window, sliding))))
  L <- max(lengths(flat))
  F <- t(vapply(flat, \(f) c(f, numeric(L - length(f))), numeric(L)))
  as.matrix(stats::dist(F, method = "manhattan"))
}

test_that("one-hot clustering: Hamming between units = L1 of their presence rows (tumbling, sliding)", {
  # every unit has a complete sliding window here (see the next test)
  cases <- list(list(1L, "tumbling", TRUE), list(3L, "tumbling", TRUE),
                list(2L, "sliding", FALSE), list(3L, "sliding", FALSE))
  invisible(lapply(cases, \(case) {
    info <- paste(case[[1]], case[[2]])
    d <- eqi_onehot_sessions(short = case[[3]])
    key <- paste(d$Actor, d$Session)
    H <- eqi_unit_l1(d, unique(key), key, case[[1]], case[[2]] == "sliding")
    o <- modifyList(opts_onehot(buildModel_window = case[[1]],
                                buildModel_window_type = case[[2]]),
                    c(eqi_cluster_opts(k = 2L), list(clustering_run = TRUE)))
    nets <- cna("cna_build_onehot_cluster_nest")(d, o)
    cl <- attr(nets, "clustering")
    expect_equal(unname(as.matrix(cl$distance)), unname(H), info = info)
    expect_true(eqi_same_partition(
      cl$assignments, cluster::pam(stats::as.dist(H), k = 2, diss = TRUE)$clustering), info = info)
    dg <- cna("cna_cluster_diagnostics_df")(nets)
    sil <- eqi_silhouette(H, cl$assignments)
    expect_equal(dg$size, tabulate(cl$assignments, 2L), info = info)
    expect_equal(dg$silhouette, vapply(seq_len(2L), \(c) mean(sil[cl$assignments == c]), numeric(1)),
                 info = info)
  }))
})

test_that("one-hot clustering runs when a unit is shorter than the sliding window", {
  # Fixed 2026-11-06, pinned here. Was (CNA 1.1.0): build_clusters() clusters only the units with a
  # complete sliding window (a shorter unit contributes nothing, as CNA's own
  # note says), but cna_build_onehot_cluster_nest() expects one assignment per
  # actor x session and stops: "one cluster assignment per actor/session".
  d <- eqi_onehot_sessions(short = TRUE)
  res <- run_cna("CNABinaryCluster", d, buildModel_window = 2L,
                 buildModel_window_type = "sliding", clustering_run = TRUE, clustering_k = 2L,
                 cluster_show_diagnostics = TRUE, seed_set = TRUE)
  expect_false(grepl("Could not build", res$errorText$content, fixed = TRUE))
  expect_equal(res$clusterDiagTable$rowCount, 2)
  key <- paste(d$Actor, d$Session)
  windowed <- Filter(\(u) sum(key == u) >= 2L, unique(key))
  o <- modifyList(opts_onehot(buildModel_window = 2L, buildModel_window_type = "sliding"),
                  c(eqi_cluster_opts(k = 2L), list(clustering_run = TRUE)))
  nets <- cna("cna_build_onehot_cluster_nest")(d, o)
  cl <- attr(nets, "clustering")
  expect_length(cl$assignments, length(windowed))
  expect_equal(unname(as.matrix(cl$distance)),
               unname(eqi_unit_l1(d, windowed, key, 2L, sliding = TRUE)))
})

# ---------------------------------------------------------------------------
# 7. Error paths
# ---------------------------------------------------------------------------

test_that("resampling a network that rests on one unit raises a classed error", {
  d <- eqi_event_data()
  no_actor <- cna("cna_build_event_nest")(d, opts_event(buildModel_variables_actor = NULL))
  expect_error(cna("cna_bootstrap_compute")(no_actor, iter = 20L, level = 0.05, seed = 1L),
               class = "nestimate_single_sequence")
  grouped <- cna("cna_build_event_nest")(d, opts_event(buildModel_variables_actor = NULL,
                                                       buildModel_variables_group = "Group"))
  expect_error(cna("cna_permutation_compute")(grouped, iter = 20L, seed = 1L),
               class = "nestimate_single_sequence")
  st <- cna("cna_stability_compute")(no_actor, "Strength", iter = 50L, threshold = 0.7,
                                     certainty = 0.95, seed = 1L)
  invisible(lapply(c("cs", "cd", "nr"), \(k)
    expect_s3_class(st[[k]], "nestimate_single_sequence")))
  expect_setequal(names(cna("cna_stability_failures")(st)), c("cs", "cd", "nr"))
  # a singleton cluster is one unit too, named in the message
  clusters <- cna("cna_cluster_networks")(eqi_event_seqs(), eqi_cluster_opts(method = "single"))
  expect_true(any(vapply(cna(".cna_nets")(clusters), \(n) nrow(n$data) == 1L, logical(1))))
  expect_error(cna("cna_bootstrap_compute")(clusters, iter = 20L, level = 0.05, seed = 1L),
               class = "nestimate_single_sequence")
})
