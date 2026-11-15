# Full backends through jmvcore: every analysis, every section switched on.
# No section may report a failure, every visible table must fill, and every
# visible image must render.

integration_case <- function(analysis) {
  ev <- regulation_data(n_actors = 250)
  # one-hot rows = sessions of Regulation_long, a 1 for every action in the
  # session (several codes per row, as real one-hot coding has)
  sess <- utils::read.csv(file.path(repo_root(), "data", "Regulation_sessions.csv"),
                          stringsAsFactors = FALSE)
  oh <- local({
    sess <- sess[sess$Actor %in% unique(ev$Actor), ]
    items <- strsplit(sess$Actions, "; ", fixed = TRUE)
    acts <- sort(unique(unlist(items)))
    m <- as.data.frame(t(vapply(items, function(x) as.integer(acts %in% x), integer(length(acts)))))
    names(m) <- acts
    m$Actor <- sess$Actor
    m$Group <- sess$Group
    m
  })
  switch(analysis,
    CNAEventSingle = list(ev, buildModel_variables_long_time = "Time_Posix",
                          buildModel_variables_group = "Group", compare_group_i = "A",
                          compare_group_j = "B", buildModel_threshold = 900),
    CNAEventCluster = list(ev, buildModel_variables_long_time = "Time_Posix",
                           buildModel_threshold = 900),
    # one-hot by Actor with windows of 2: resampling and stability take whole
    # actors and keep their windows
    CNABinarySingle = list(oh, buildModel_variables_onehot = setdiff(names(oh), c("Group", "Actor")),
                           buildModel_variables_actor = "Actor", buildModel_variables_session = NULL,
                           buildModel_variables_group = "Group", compare_group_i = "A",
                           compare_group_j = "B", buildModel_window = 2L),
    # clustering needs an Actor or a Session: they are the units clustered
    CNABinaryCluster = list(oh, buildModel_variables_onehot = setdiff(names(oh), c("Group", "Actor")),
                            buildModel_variables_actor = "Actor",
                            buildModel_variables_session = NULL, buildModel_window = 2L),
    CNADelimited = list(sess, buildModel_variables_field = "Actions"))
}


for (analysis in ANALYSES) {
  test_that(paste(analysis, "runs every section without a failure"), {
    case <- integration_case(analysis)
    extra <- list(rules_min_support = 0.05, cliques_threshold = 5)
    extra <- extra[intersect(names(extra), names(formals(cna(analysis))))]
    args <- c(list(analysis, case[[1]]), case[-1], all_sections_on(analysis), extra)
    res <- do.call(run_cna, args)
    expect_length(section_failures(res), 0)
    expect_false(isTRUE(res$errorText$visible))

    items <- res$itemNames
    tables <- Filter(function(nm) inherits(res[[nm]], "Table") && res[[nm]]$visible, items)
    expect_gt(length(tables), 5)
    # every visible table either has rows or a note saying why it has none
    for (nm in tables)
      expect_true(res[[nm]]$rowCount > 0 || length(res[[nm]]$notes) > 0,
                  label = paste(analysis, nm, "has rows or an explanatory note"))

    images <- Filter(function(nm) inherits(res[[nm]], "Image") && res[[nm]]$visible, items)
    expect_gt(length(images), 5)
    for (nm in images) expect_equal(render_image(res[[nm]]), "ok", label = paste(analysis, nm))
  })
}

test_that("nothing is shown until the required variables are assigned", {
  res <- run_cna("CNAEventSingle", make_event_data(), buildModel_variables_long_action = NULL,
                 centrality_show_table = TRUE)
  expect_false(res$centralityTable$visible)
  expect_false(res$buildModel_plot$visible)
  expect_false(isTRUE(res$errorText$visible))
})

test_that("an invalid input is explained, not raised", {
  res <- run_cna("CNABinarySingle", make_onehot_data(), buildModel_variables_onehot = "code1")
  expect_true(res$errorText$visible)
  expect_match(res$errorText$content, "at least 2 one-hot columns")
  res <- run_cna("CNABinaryCluster", make_onehot_data(), buildModel_variables_actor = NULL,
                 buildModel_variables_session = NULL, clustering_run = TRUE)
  expect_match(res$errorText$content, "needs an Actor or a Session")
})

test_that("clustering waits for 'Run clustering', but Choosing k works before it", {
  res <- run_cna("CNAEventCluster", make_event_data(n_actor = 12, per = 12),
                 clustering_run = FALSE, centrality_show_table = TRUE,
                 cluster_show_choice = TRUE, cluster_choice_kmax = 3L)
  expect_equal(res$centralityTable$rowCount, 0)
  expect_true(res$clusteringNote$visible)
  expect_equal(res$clusterChoiceTable$rowCount, 2)
})

test_that("group-only sections stay hidden without a group variable", {
  res <- run_event_single(make_event_data(), permutation_show_table = TRUE,
                          compare_show_summary = TRUE)
  expect_false(res$permutationTable$visible)
  expect_false(res$compareSummaryTable$visible)
})

test_that("a section clears its old message before running again (review fix)", {
  heading <- local({
    content <- "This analysis could not be completed: old error"
    list(setContent = function(x) content <<- x, get = function() content)
  })
  ctx <- list(results = list(h = heading), multi = FALSE)
  cna(".cna_run_section")(ctx, list(title = "h", run = function(ctx) NULL))
  expect_identical(heading$get(), "")
  cna(".cna_run_section")(ctx, list(title = "h", run = function(ctx) stop("boom")))
  expect_match(heading$get(), "could not be completed: boom")
})

test_that("one failing section does not stop the others", {
  # optimal community detection refuses > 50 nodes; centrality must still fill
  n <- 55
  dat <- data.frame(Actor = rep(seq_len(40), each = 6),
                    Action = sprintf("x%02d", (seq_len(240) * 7) %% n + 1),
                    stringsAsFactors = FALSE)
  res <- run_event_single(dat, buildModel_variables_long_time = NULL,
                          community_show_table = TRUE, community_methods = "optimal",
                          centrality_show_table = TRUE)
  expect_match(res$communityTitle$content, "could not be completed")
  expect_gt(res$centralityTable$rowCount, 0)
})
