# Bundled example datasets: each is listed in the manifest, reproducible from
# tools/make_datasets.R, runs in the analyses it is meant for, and the
# simulated ones show the structure that was planted in them.

read_dataset <- function(name)
  utils::read.csv(file.path(repo_root(), "data", name), stringsAsFactors = FALSE)

test_that("every bundled CSV is in the manifest and every manifest entry exists", {
  manifest <- yaml::read_yaml(file.path(repo_root(), "jamovi", "0000.yaml"))
  listed <- vapply(manifest$datasets, `[[`, "", "path")
  on_disk <- list.files(file.path(repo_root(), "data"), pattern = "\\.csv$")
  expect_setequal(listed, on_disk)
  descriptions <- vapply(manifest$datasets, `[[`, "", "description")
  expect_true(all(grepl("^(Real|Simulated)", descriptions)))   # provenance stated
})

test_that("0/1 datasets hold only 0 and 1 in their indicator columns", {
  srl <- read_dataset("SRL_onehot.csv")
  la <- read_dataset("LearningActivities.csv")
  hai <- read_dataset("HumanAI_turns.csv")
  ind <- list(srl[c("Plan", "Monitor", "Evaluate", "Reflect", "Read", "Summarize",
                    "Practice", "Review")],
              la[c("Reading", "Video", "Forum", "Quiz", "Coding", "Review")],
              hai[setdiff(names(hai), c("Session", "Project", "Date", "Speaker", "Turn", "Codes"))])
  invisible(lapply(ind, function(x) expect_true(all(unlist(x) %in% c(0, 1)))))
  # the 0/1 columns and the code list of a turn agree
  codes <- setdiff(names(hai), c("Session", "Project", "Date", "Speaker", "Turn", "Codes"))
  listed <- strsplit(hai$Codes, "; ", fixed = TRUE)
  from_columns <- apply(hai[codes] == 1, 1, function(r) sort(codes[r]))
  expect_identical(lapply(listed, sort), unname(from_columns))
})

test_that("SRL_events: the planted learning cycles are recovered", {
  srl <- read_dataset("SRL_events.csv")
  o <- list(buildModel_variables_long_action = "Action", buildModel_variables_actor = "Actor",
            buildModel_variables_long_time = "Time", buildModel_variables_long_order = "Order",
            buildModel_variables_group = "Group", buildModel_scaling = "noScaling",
            buildModel_threshold = 900)
  model <- cna("cna_build_event_nest")(srl, o)
  expect_setequal(names(model), c("Strategic", "Surface"))
  pairs <- cna("cna_build_mixed")(srl, o, "event_single", model)
  strongest <- function(tr, from) {
    w <- tr[from, ]
    names(w)[which.max(w)]
  }
  # Strategic: reflect -> plan, summarise -> monitor; Surface: read -> review
  expect_equal(strongest(pairs$Strategic$transition$weights, "Reflect"), "Plan")
  expect_equal(strongest(pairs$Strategic$transition$weights, "Summarize"), "Monitor")
  expect_equal(strongest(pairs$Surface$transition$weights, "Read"), "Review")
  # and the groups' co-occurrence networks differ significantly somewhere
  pt <- cna("cna_permutation_compute")(model, iter = 200L, seed = 1L)
  expect_gt(sum(pt[[1]]$summary$p_value < 0.05), 0)
})

test_that("each dataset runs in its intended analyses without a failing section", {
  srl <- read_dataset("SRL_events.csv")
  srl_oh <- read_dataset("SRL_onehot.csv")
  la <- read_dataset("LearningActivities.csv")
  hai <- read_dataset("HumanAI_turns.csv")
  reg_s <- read_dataset("Regulation_sessions.csv")
  srl_codes <- c("Plan", "Monitor", "Evaluate", "Reflect", "Read", "Summarize", "Practice", "Review")
  hai_codes <- setdiff(names(hai), c("Session", "Project", "Date", "Speaker", "Turn", "Codes"))
  shows <- list(centrality_show_table = TRUE, community_show_table = TRUE,
                buildModel_show_mixed = TRUE)
  runs <- list(
    run_cna("CNAEventSingle", srl, buildModel_variables_long_time = "Time",
            buildModel_variables_long_order = "Order", buildModel_variables_group = "Group",
            compare_group_i = "Strategic", compare_group_j = "Surface", buildModel_threshold = 900,
            centrality_show_table = TRUE, buildModel_show_mixed = TRUE,
            compare_show_edges = TRUE),
    run_cna("CNABinarySingle", srl_oh, buildModel_variables_onehot = srl_codes,
            buildModel_variables_actor = "Actor", buildModel_variables_session = NULL,
            buildModel_variables_group = "Group", compare_group_i = "Strategic",
            compare_group_j = "Surface", centrality_show_table = TRUE,
            buildModel_show_mixed = TRUE),
    run_cna("CNABinaryCluster", la, buildModel_variables_onehot = c("Reading", "Video", "Forum",
                                                                    "Quiz", "Coding", "Review"),
            buildModel_variables_actor = "Student", buildModel_variables_session = NULL,
            clustering_run = TRUE, centrality_show_table = TRUE, seed_set = TRUE),
    run_cna("CNABinarySingle", hai, buildModel_variables_onehot = hai_codes,
            buildModel_variables_actor = "Session", buildModel_variables_session = NULL,
            buildModel_variables_group = "Speaker", compare_group_i = "AI",
            compare_group_j = "Human", centrality_show_table = TRUE),
    run_cna("CNADelimited", hai, buildModel_variables_field = "Codes",
            centrality_show_table = TRUE, community_show_table = TRUE),
    run_cna("CNADelimited", reg_s, buildModel_variables_field = "Actions",
            centrality_show_table = TRUE))
  invisible(lapply(runs, function(res) {
    expect_length(section_failures(res), 0)
    expect_false(isTRUE(res$errorText$visible))
    expect_gt(res$centralityTable$rowCount, 0)
  }))
  # the mixed networks render for the SRL datasets
  invisible(lapply(runs[1:2], function(res)
    expect_equal(render_image(res$buildModel_mixed), "ok")))
})
