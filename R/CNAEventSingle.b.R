# Event Data (Single/Group) — Nestimate engine (CNA v2)

CNAEventSingleClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
  "CNAEventSingleClass",
  inherit = CNAEventSingleBase,
  private = list(
    .run = function() {

      # Shared helpers (functions live in the package namespace; source() is a
      # fallback for environments where they were not attached).
      if (!exists("cna_build_event_nest", mode = "function")) {
        source_file <- system.file("R", "cna_utils.R", package = "CNA")
        if (source_file == "") {
          pkg_dir <- find.package("CNA", quiet = TRUE)
          if (length(pkg_dir) > 0) {
            utils_path <- file.path(pkg_dir, "R", "cna_utils.R")
            if (file.exists(utils_path)) source(utils_path, local = TRUE)
          }
        } else {
          source(source_file, local = TRUE)
        }
      }

      is_multi <- !is.null(self$options$buildModel_variables_group)

      self$results$instructions$setContent(
        cna_instructions_html("event_single")
      )

      # Action is required
      if (is.null(self$options$buildModel_variables_long_action)) {
        self$results$errorText$setVisible(FALSE)
        return()
      }

      if (!self$options$buildModel_show_matrix &&
          !self$options$buildModel_show_plot &&
          !self$options$buildModel_show_histo) {
        self$results$buildModelTitle$setVisible(FALSE)
      }

      # ── Build the co-occurrence network (Nestimate) ──
      model <- NULL
      if (self$results$buildModelContent$isFilled()) {
        model <- self$results$buildModelContent$state
      } else if (!is.null(self$data) && ncol(self$data) >= 1) {
        tryCatch({
          model <- cna_build_event_nest(self$data, self$options)
        }, error = function(e) {
          error_msg <- tolower(as.character(e$message))
          if (grepl("time|date|posix|format", error_msg)) {
            self$results$errorText$setContent("Please enter an appropriate time format")
          } else {
            self$results$errorText$setContent(paste("Data preparation error:", e$message))
          }
          self$results$errorText$setVisible(TRUE)
        })
      }

      if (!is.null(model)) {
        if (!self$results$buildModelContent$isFilled()) {
          self$results$buildModelContent$setContent(model)
          self$results$buildModelContent$setState(model)
        }
        self$results$buildModelContent$setVisible(self$options$buildModel_show_matrix)
        self$results$buildModel_plot$setVisible(self$options$buildModel_show_plot)
        self$results$buildModel_histo$setVisible(self$options$buildModel_show_histo)
        self$results$buildModel_frequencies$setVisible(self$options$buildModel_show_frequencies)
        self$results$buildModel_mosaic$setVisible(self$options$buildModel_show_mosaic)
        self$results$buildModel_distribution$setVisible(self$options$buildModel_show_distribution)
      }

      # ── Analyses (all Nestimate / cograph) ──
      cna_run_centrality_nest(model, self$options, self$results, is_multi)
      cna_run_community_nest(model, self$options, self$results, is_multi)
      cna_run_bootstrap_nest(model, self$options, self$results, is_multi)
      cna_run_stability_nest(model, self$options, self$results, is_multi)
      if (is_multi) {
        cna_run_permutation_nest(model, self$options, self$results)
        cna_run_compare_nest(model, self$options, self$results, mode = "group")
      }
    },

    # ── Plot dispatchers (cograph::splot renders single + group grids) ──

    .showBuildModelPlot = function(image, ...) {
      cna_splot_model_nest(self$results$buildModelContent$state, self$options)
    },

    .showBuildModelHisto = function(image, ...) {
      if (!self$options$buildModel_show_histo) return(FALSE)
      is_multi <- !is.null(self$options$buildModel_variables_group)
      cna_splot_histo_nest(self$results$buildModelContent$state, is_multi)
    },

    .showBuildModelFrequencies = function(image, ...) {
      if (!self$options$buildModel_show_frequencies) return(FALSE)
      is_multi <- !is.null(self$options$buildModel_variables_group)
      cna_splot_freq_nest(self$results$buildModelContent$state, is_multi)
    },

    .showBuildModelMosaic = function(image, ...) {
      if (!self$options$buildModel_show_mosaic) return(FALSE)
      cna_splot_mosaic_nest(self$results$buildModelContent$state)
    },

    .showBuildModelDistribution = function(image, ...) {
      if (!self$options$buildModel_show_distribution) return(FALSE)
      cna_splot_distribution_nest(self$results$buildModelContent$state)
    },

    .showCentralityPlot = function(image, ...) {
      if (!self$options$centrality_show_plot) return(FALSE)
      is_multi <- !is.null(self$options$buildModel_variables_group)
      cna_splot_centrality_nest(self$results$buildModelContent$state,
                                self$options, is_multi)
    },

    .showCommunityPlot = function(image, ...) {
      if (!self$options$community_show_plot) return(FALSE)
      is_multi <- !is.null(self$options$buildModel_variables_group)
      cna_splot_community_nest(self$results$community_plot$state, self$options, is_multi)
    },

    .showBootstrapPlot = function(image, ...) {
      if (!self$options$bootstrap_show_plot) return(FALSE)
      cna_splot_bootstrap_nest(self$results$bootstrap_plot$state)
    },

    .showStabilityPlot = function(image, ...) {
      if (!self$options$stability_show_plot) return(FALSE)
      cna_splot_stability_nest(self$results$stability_plot$state)
    },

    .showPermutationPlot = function(image, ...) {
      if (!self$options$permutation_show_plot) return(FALSE)
      cna_splot_permutation_nest(self$results$permutation_plot$state)
    },

    .showComparePlot = function(image, ...) {
      if (!isTRUE(self$options$compare_show_plot)) return(FALSE)
      cna_splot_compare_nest(self$results$compare_plot$state, self$options)
    }
  )
)
