#!/usr/bin/env Rscript
# tools/make_datasets.R
#
# Builds the example datasets bundled with the module (data/*.csv). Every file
# is reproducible from this script; the manifest descriptions live in
# tools/generate_jamovi.R (write_manifest).
#
#   Regulation_sessions.csv  real    delimited: the Regulation_long sessions as
#                                    "a; b; c" action lists (one row per session)
#   SRL_events.csv           simulated  event (long) log of 8 self-regulated
#                                    learning (SRL) actions, two groups with
#                                    different learning cycles (simulab)
#   SRL_onehot.csv           simulated  the same log coded 0/1 per episode of
#                                    four consecutive events
#   LearningActivities.csv   simulated  0/1 activity indicators, 200 students x
#                                    30 time points (Nestimate::learning_activities)
#   HumanAI_turns.csv        real    human-AI pair-programming turns: codes as a
#                                    "; " list and as 0/1 columns
#                                    (Nestimate::human_long + ai_long)
#
# Usage, from the repository root:  Rscript tools/make_datasets.R

stopifnot("run from the repository root" = file.exists("DESCRIPTION"))
suppressPackageStartupMessages({
  library(Nestimate)
  library(simulab)
})

write_data <- function(df, name) {
  utils::write.csv(df, file.path("data", name), row.names = FALSE)
  cat(sprintf("Wrote data/%-24s %6d rows x %2d columns\n", name, nrow(df), ncol(df)))
}

# One 0/1 column per state: 1 when the state occurs among `values` of the row's
# key. Rows follow the order of first appearance of `key`.
indicator_columns <- function(key, values, states) {
  key <- factor(key, levels = unique(key))
  present <- lapply(states, function(s) as.integer(tapply(values == s, key, any)))
  stats::setNames(as.data.frame(present), states)
}

# ---------------------------------------------------------------------------
# Regulation_sessions.csv (real): sessions of Regulation_long as item lists
# ---------------------------------------------------------------------------

long <- utils::read.csv(file.path("data", "Regulation_long.csv"), stringsAsFactors = FALSE)
events <- Nestimate::prepare(long, actor = "Actor", action = "Action",
                             time = "Time_Posix", time_threshold = 900)$long_data
items <- tapply(events$Action, events$.session_id, function(a) paste(unique(a), collapse = "; "))
meta <- events[!duplicated(events$.session_id), c(".session_id", "Actor", "Achievers", "Group")]
sessions <- data.frame(Session = meta$.session_id, Actor = meta$Actor,
                       Achievers = meta$Achievers, Group = meta$Group,
                       Actions = as.character(items[meta$.session_id]),
                       stringsAsFactors = FALSE)
sessions <- sessions[order(sessions$Actor, sessions$Session), ]
stopifnot("one row per session" = !anyDuplicated(sessions$Session),
          "all sessions kept" = nrow(sessions) == length(unique(events$.session_id)))
write_data(sessions, "Regulation_sessions.csv")

# ---------------------------------------------------------------------------
# SRL_events.csv (simulated): two groups with different learning cycles
# ---------------------------------------------------------------------------
# Each state moves to one preferred next state with probability 0.5, stays with
# 0.15, and spreads the rest evenly. "Strategic" learners cycle through
# planning, reading, summarising, monitoring, evaluating and reflecting;
# "Surface" learners loop between reading, reviewing and practising.

srl_states <- c("Plan", "Monitor", "Evaluate", "Reflect", "Read", "Summarize",
                "Practice", "Review")

learning_cycle <- function(preferred, stay = 0.15, strength = 0.5) {
  n <- length(srl_states)
  m <- matrix((1 - stay - strength) / (n - 2), n, n, dimnames = list(srl_states, srl_states))
  diag(m) <- stay
  m[cbind(names(preferred), preferred)] <- strength
  stopifnot("each row is a probability distribution" = all(abs(rowSums(m) - 1) < 1e-12))
  m
}
strategic <- learning_cycle(c(Plan = "Read", Read = "Summarize", Summarize = "Monitor",
                              Monitor = "Evaluate", Evaluate = "Reflect", Reflect = "Plan",
                              Practice = "Monitor", Review = "Evaluate"))
surface <- learning_cycle(c(Plan = "Read", Read = "Review", Review = "Practice",
                            Practice = "Read", Monitor = "Read", Evaluate = "Review",
                            Summarize = "Read", Reflect = "Read"))

sim <- as.data.frame(simulab::simulate_event_log(
  groups = 2, actors = 100, courses = c("Course A", "Course B"), states = srl_states,
  sequence_length = c(15L, 40L), achievement_levels = c("Low", "High"),
  interval_range = c(20, 1200), transitions = list(strategic, surface), seed = 2026))
group_label <- c(`Group 1` = "Strategic", `Group 2` = "Surface")
srl <- data.frame(
  Actor = sim$id, Group = unname(group_label[sim$group]), Course = sim$course,
  Achievement = sim$achievement,
  Time = format(sim$timestamp, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
  Order = sim$period, Action = sim$state, stringsAsFactors = FALSE)
srl <- srl[order(srl$Group, srl$Actor, srl$Order), ]
stopifnot("no missing values" = !anyNA(srl),
          "both groups present" = setequal(unique(srl$Group), group_label))
write_data(srl, "SRL_events.csv")

# ---------------------------------------------------------------------------
# SRL_onehot.csv (simulated): the same log, 0/1 per episode of 4 events
# ---------------------------------------------------------------------------

episode_length <- 4L
srl$Episode <- stats::ave(srl$Order, srl$Actor,
                          FUN = function(o) (seq_along(o) - 1L) %/% episode_length + 1L)
key <- paste(srl$Actor, srl$Episode)
first <- !duplicated(key)
srl_onehot <- data.frame(srl[first, c("Actor", "Group", "Course", "Achievement", "Episode")],
                         indicator_columns(key, srl$Action, srl_states),
                         row.names = NULL, check.names = FALSE)
stopifnot("one row per actor episode" = nrow(srl_onehot) == length(unique(key)),
          "every episode has an action" = all(rowSums(srl_onehot[srl_states]) >= 1))
write_data(srl_onehot, "SRL_onehot.csv")

# ---------------------------------------------------------------------------
# LearningActivities.csv (simulated, Nestimate): 0/1 activities over time
# ---------------------------------------------------------------------------

la <- Nestimate::learning_activities
activities <- setdiff(names(la), "student")
learning <- data.frame(Student = la$student,
                       Time = stats::ave(la$student, la$student, FUN = seq_along),
                       la[activities], check.names = FALSE)
stopifnot("binary indicators" = all(unlist(learning[activities]) %in% c(0L, 1L)),
          "30 time points per student" = all(table(learning$Student) == 30L))
write_data(learning, "LearningActivities.csv")

# ---------------------------------------------------------------------------
# HumanAI_turns.csv (real, Nestimate): one row per conversation turn
# ---------------------------------------------------------------------------

turn_rows <- rbind(
  data.frame(Nestimate::human_long, Speaker = "Human", stringsAsFactors = FALSE),
  data.frame(Nestimate::ai_long, Speaker = "AI", stringsAsFactors = FALSE))
turn_rows <- turn_rows[order(turn_rows$session_id, turn_rows$timestamp,
                             turn_rows$message_id, turn_rows$code_order), ]
turn_key <- paste(turn_rows$session_id, turn_rows$message_id)
codes <- sort(unique(turn_rows$code))
first <- !duplicated(turn_key)
code_lists <- tapply(turn_rows$code, factor(turn_key, levels = unique(turn_key)),
                     function(x) paste(unique(x), collapse = "; "))
turns <- data.frame(
  Session = turn_rows$session_id[first], Project = turn_rows$project[first],
  Date = turn_rows$session_date[first], Speaker = turn_rows$Speaker[first],
  Turn = stats::ave(seq_len(sum(first)), turn_rows$session_id[first], FUN = seq_along),
  Codes = as.character(code_lists),
  indicator_columns(turn_key, turn_rows$code, codes),
  check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("one speaker per turn" = !anyDuplicated(turn_key[first]),
          "every turn coded" = all(rowSums(turns[codes]) >= 1),
          "17 codes" = length(codes) == 17L)
write_data(turns, "HumanAI_turns.csv")
