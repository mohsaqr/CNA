# One-hot (binary) input: what CNA hands Nestimate, and the results it gets
# back (Nestimate 1.0.7, a1754ab). Reference values are computed independently
# with crossprod() / tcrossprod(colSums(block)) per actor x session.

# The window definition of Nestimate >= 1.0.7 and tna::import_onehot(): with a
# window above 1, rows are counted in blocks of `window`; sliding mode first
# replaces the rows by one row per complete sliding window, marking the codes
# active anywhere in it (a unit shorter than the window contributes nothing).
onehot_ref <- function(df, codes, unit = NULL, window = 1L, mode = "non-overlapping") {
  sliding_rows <- function(M) {
    n <- nrow(M)
    if (n < window) return(M[0, , drop = FALSE])
    rows <- lapply(seq_len(n - window + 1L), \(i)
      as.numeric(colSums(M[i:(i + window - 1L), , drop = FALSE]) > 0))
    matrix(unlist(rows), ncol = ncol(M), byrow = TRUE)
  }
  windowed <- function(M) {
    if (window <= 1L) return(crossprod(M))
    if (mode == "overlapping") M <- sliding_rows(M)
    n <- nrow(M)
    blocks <- if (n) lapply(seq(1L, n, by = window), \(s) s:min(n, s + window - 1L)) else list()
    Reduce(`+`, lapply(blocks, \(b) tcrossprod(colSums(M[b, , drop = FALSE]))),
           matrix(0, ncol(M), ncol(M)))
  }
  parts <- if (is.null(unit)) list(df) else split(df, interaction(df[unit], drop = TRUE))
  Reduce(`+`, lapply(parts, \(p) windowed(as.matrix(p[codes]))))
}

# Independent one-hot transactions: the codes active in each window of each
# actor x session (a window of 1 is a row; sliding: every complete window).
onehot_transactions <- function(df, codes, unit, window = 1L, mode = "non-overlapping") {
  X <- as.matrix(df[codes])
  windows_of <- function(rows) {
    n <- length(rows)
    starts <- if (window <= 1L) seq_len(n)
              else if (mode == "non-overlapping") seq(1L, n, by = window)
              else if (n >= window) seq_len(n - window + 1L) else integer()
    lapply(starts, \(s) rows[s:min(n, s + window - 1L)])
  }
  by_unit <- split(seq_len(nrow(df)), interaction(df[unit], drop = TRUE))
  windows <- unlist(lapply(by_unit, windows_of), recursive = FALSE)
  unname(lapply(windows, \(r) codes[colSums(X[r, , drop = FALSE]) > 0]))
}

codes3 <- c("A", "B", "C")
tiny <- data.frame(A = c(1, 1, 0, 0, 1, 1, 0), B = c(0, 1, 1, 0, 0, 1, 1),
                   C = c(0, 0, 1, 1, 1, 0, 1), Actor = c(rep("a1", 5), rep("a2", 2)),
                   Session = c(1, 1, 1, 2, 2, 1, 1))
opts3 <- function(...) modifyList(list(
  buildModel_variables_onehot = codes3, buildModel_variables_actor = "Actor",
  buildModel_variables_session = NULL, buildModel_variables_group = NULL,
  buildModel_window = 1L, buildModel_window_type = "tumbling",
  buildModel_scaling = "noScaling"), list(...))

test_that("co-occurrence equals the per-actor / per-session window definition", {
  for (w in 1:3) for (type in c("tumbling", "sliding")) {
    mode <- if (type == "sliding") "overlapping" else "non-overlapping"
    net <- suppressWarnings(cna("cna_build_onehot_nest")(tiny, opts3(buildModel_window = w,
                                                                     buildModel_window_type = type)))
    expect_equal(unname(unclass(net$weights)),
                 unname(onehot_ref(tiny, codes3, "Actor", w, mode)), info = paste(w, type))
  }
  net <- suppressWarnings(cna("cna_build_onehot_nest")(tiny, opts3(buildModel_variables_session = "Session",
                                                                   buildModel_window = 2L)))
  expect_equal(unname(unclass(net$weights)),
               unname(onehot_ref(tiny, codes3, c("Actor", "Session"), 2L)))
})

test_that("a 'Session' column the user did not assign is never used", {
  # CNA passes only the assigned columns, so the network is the actor-only one
  net <- suppressWarnings(cna("cna_build_onehot_nest")(tiny, opts3(buildModel_window = 2L)))
  expect_equal(unname(unclass(net$weights)), unname(onehot_ref(tiny, codes3, "Actor", 2L)))
})

test_that("0/1 stored as numbers, logicals or labels all give the same network", {
  base <- suppressWarnings(cna("cna_build_onehot_nest")(tiny, opts3()))$weights
  as_factor <- tiny
  as_factor[codes3] <- lapply(tiny[codes3], factor)
  as_logical <- tiny
  as_logical[codes3] <- lapply(tiny[codes3], as.logical)
  as_text <- tiny
  as_text[codes3] <- lapply(tiny[codes3], function(x) ifelse(x == 1, "TRUE", "FALSE"))
  for (d in list(as_factor, as_logical, as_text))
    expect_equal(suppressWarnings(cna("cna_build_onehot_nest")(d, opts3()))$weights, base)
})

test_that("non-0/1 values are an error naming the column, not silent weights", {
  bad <- tiny
  bad$B[2] <- 2
  expect_error(cna("cna_build_onehot_nest")(bad, opts3()), class = "cna_not_binary")
  bad <- tiny
  bad$C <- factor(ifelse(tiny$C == 1, "Yes", "No"))
  err <- tryCatch(cna("cna_build_onehot_nest")(bad, opts3()), error = function(e) e)
  expect_s3_class(err, "cna_not_binary")
  expect_match(conditionMessage(err), "'C'.*Yes|No")
})

test_that("columns named Action, Session or session_id are read by their assigned role", {
  # Nestimate >= 1.0.7 gives a column a role only when it is declared
  renamed <- tiny
  names(renamed)[1:2] <- c("Action", "session_id")
  codes <- c("Action", "session_id", "C")
  net <- suppressWarnings(cna("cna_build_onehot_nest")(renamed, opts3(buildModel_variables_onehot = codes,
                                                                     buildModel_window = 2L)))
  expect_equal(rownames(net$weights), codes)
  expect_equal(unname(unclass(net$weights)), unname(onehot_ref(tiny, codes3, "Actor", 2L)))
  as_actor <- tiny
  names(as_actor)[names(as_actor) == "Actor"] <- "Session"
  net <- suppressWarnings(cna("cna_build_onehot_nest")(as_actor, opts3(buildModel_variables_actor = "Session")))
  expect_equal(unname(unclass(net$weights)), unname(onehot_ref(tiny, codes3, "Actor")))
})

test_that("rows with a missing Actor are dropped and counted; missing cells are noted", {
  d <- tiny
  d$Actor[2] <- NA
  d$C[7] <- NA
  net <- suppressWarnings(cna("cna_build_onehot_nest")(d, opts3()))
  kept <- d[!is.na(d$Actor), ]
  kept$C[is.na(kept$C)] <- 0
  expect_equal(unname(unclass(net$weights)), unname(onehot_ref(kept, codes3, "Actor")))
  notes <- attr(net, "cna_notes")
  expect_true(any(grepl("1 row with a missing Actor or Session was left out", notes)))
  expect_true(any(grepl("1 missing one-hot value was counted as 0", notes)))
})

test_that("sliding windows: actors shorter than the window are reported", {
  net <- suppressWarnings(cna("cna_build_onehot_nest")(tiny, opts3(buildModel_window = 3L,
                                                                   buildModel_window_type = "sliding")))
  expect_true(any(grepl("1 of 2 actors/sessions have fewer rows than the sliding window \\(3\\)",
                        attr(net, "cna_notes"))))
  tumbling <- suppressWarnings(cna("cna_build_onehot_nest")(tiny, opts3(buildModel_window = 3L)))
  expect_null(attr(tumbling, "cna_notes"))
})

test_that("the notes reach the analysis output", {
  d <- tiny
  d$Actor[2] <- NA
  res <- run_cna("CNABinarySingle", d, buildModel_variables_onehot = codes3,
                 buildModel_variables_actor = "Actor", buildModel_variables_session = NULL,
                 buildModel_show_summary = TRUE)
  expect_match(res$buildModelTitle$content, "1 row with a missing Actor or Session was left out")
})

# --- one-hot clustering ---------------------------------------------------------

srl_codes <- c("Plan", "Monitor", "Evaluate", "Reflect", "Read", "Summarize", "Practice", "Review")
srl_onehot <- function() utils::read.csv(file.path(repo_root(), "data", "SRL_onehot.csv"),
                                         stringsAsFactors = FALSE)
cluster_opts <- function(...) modifyList(list(
  buildModel_variables_onehot = srl_codes, buildModel_variables_actor = "Actor",
  buildModel_variables_session = NULL, buildModel_window = 1L,
  buildModel_window_type = "tumbling", buildModel_scaling = "noScaling",
  clustering_run = TRUE, clustering_k = 3L, clustering_dissimilarity = "hamming",
  clustering_method = "pam", seed_set = TRUE, seed_value = 1L), list(...))

unit_rows <- function(d, nets, k) {
  # the assignments follow the first-seen order of the actors
  unit_cluster <- stats::setNames(attr(nets, "clustering")$assignments, unique(d$Actor))
  d[unit_cluster[as.character(d$Actor)] == k, ]
}

test_that("each cluster's network is the one-hot network of its own rows (any window)", {
  d <- srl_onehot()
  cases <- list(list(1L, "tumbling"), list(6L, "tumbling"), list(3L, "sliding"))
  invisible(lapply(cases, \(case) {
    mode <- if (case[[2]] == "sliding") "overlapping" else "non-overlapping"
    nets <- cna("cna_build_onehot_cluster_nest")(d, cluster_opts(buildModel_window = case[[1]],
                                                                 buildModel_window_type = case[[2]]))
    expect_equal(names(nets), paste("Cluster", 1:3))
    expect_null(attr(nets, "cna_notes"))
    invisible(Map(function(net, k) expect_equal(
      unname(unclass(net$weights)),
      unname(onehot_ref(unit_rows(d, nets, k), srl_codes, "Actor", case[[1]], mode)),
      info = paste(case, collapse = " ")), nets, 1:3))
  }))
})

test_that("clusters come from Nestimate's windowed units; unscaled == build_network(<net_clustering>)", {
  d <- srl_onehot()
  o <- cluster_opts(buildModel_window = 3L, buildModel_window_type = "sliding")
  nets <- cna("cna_build_onehot_cluster_nest")(d, o)
  cl <- attr(nets, "clustering")
  expect_s3_class(cl, "net_clustering")
  expect_length(cl$assignments, length(unique(d$Actor)))
  ref <- Nestimate::build_network(cl)
  invisible(Map(function(net, r) {
    m <- unclass(r$weights)
    expect_equal(unclass(net$weights)[rownames(m), colnames(m)], m)
  }, nets, ref))
})

test_that("scaled cluster networks are scaled (build_network(<net_clustering>) ignores it)", {
  d <- srl_onehot()
  nets <- cna("cna_build_onehot_cluster_nest")(d, cluster_opts(buildModel_window = 2L,
                                                               buildModel_scaling = "minmax"))
  invisible(Map(function(net, k) {
    rows <- unit_rows(d, nets, k)
    ref <- Nestimate::build_network(rows[c(srl_codes, "Actor")], method = "co_occurrence",
                                    codes = srl_codes, format = "onehot", actor = "Actor",
                                    window_size = 2L, scaling = "minmax")
    expect_equal(unclass(net$weights), unclass(ref$weights))
    expect_lte(max(net$weights), 1)
  }, nets, seq_along(nets)))
})

test_that("Binary Clustering runs with diagnostics, choosing k, mixed view and resampling", {
  d <- srl_onehot()
  res <- run_cna("CNABinaryCluster", d, buildModel_variables_onehot = srl_codes,
                 buildModel_variables_actor = "Actor", buildModel_variables_session = NULL,
                 clustering_run = TRUE, clustering_k = 2L, seed_set = TRUE,
                 cluster_show_diagnostics = TRUE, cluster_show_choice = TRUE,
                 cluster_choice_kmax = 3L, buildModel_show_mixed = TRUE,
                 bootstrap_show_table = TRUE, bootstrap_iteration = 20L,
                 permutation_show_table = TRUE, permutation_iter = 20L)
  expect_length(section_failures(res), 0)
  expect_equal(sum(table_df(res, "clusterDiagTable")$size), length(unique(d$Actor)))
  expect_equal(res$clusterChoiceTable$rowCount, 2)
  expect_gt(res$bootstrapTable$rowCount, 0)
  expect_gt(res$permutationTable$rowCount, 0)
  expect_equal(render_image(res$buildModel_mixed), "ok")
})

# --- association rules on one-hot data -------------------------------------------

test_that("one-hot rules use the codes as items, one transaction per window", {
  d <- srl_onehot()
  key <- function(r) sort(paste(r$antecedent, r$consequent, sprintf("%.12f", r$support),
                                sprintf("%.12f", r$confidence)))
  cases <- list(list(1L, "tumbling"), list(2L, "tumbling"), list(3L, "sliding"))
  invisible(lapply(cases, \(case) {
    mode <- if (case[[2]] == "sliding") "overlapping" else "non-overlapping"
    net <- suppressWarnings(cna("cna_build_onehot_nest")(d, cluster_opts(
      buildModel_window = case[[1]], buildModel_window_type = case[[2]])))
    rules <- cna("cna_rules_compute")(net, min_support = 0.05, min_confidence = 0.3)[[1]]
    expect_setequal(rules$items, srl_codes)
    tx <- onehot_transactions(d, srl_codes, "Actor", case[[1]], mode)
    expect_equal(rules$n_transactions, length(tx))
    expect_equal(rules$frequent$support[rules$frequent$itemset == "Plan"],
                 mean(vapply(tx, \(t) "Plan" %in% t, logical(1))))
    ref <- Nestimate::association_rules(tx, min_support = 0.05, min_confidence = 0.3)
    expect_identical(key(rules$rules), key(ref$rules), info = paste(case, collapse = " "))
  }))
})

test_that("one-hot rules per group follow the networks' groups", {
  d <- srl_onehot()
  net <- suppressWarnings(cna("cna_build_onehot_nest")(d, cluster_opts(buildModel_variables_group = "Group")))
  rules <- cna("cna_rules_compute")(net, min_support = 0.05, min_confidence = 0.3)
  expect_equal(names(rules), names(net))
  expect_equal(rules$Surface$n_transactions, sum(d$Group == "Surface"))
})
