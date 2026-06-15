# Event Data (Clustering) — Nestimate engine (CNA v2)

CNAEventClusterClass <- if (requireNamespace('jmvcore', quietly=TRUE)) R6::R6Class(
  "CNAEventClusterClass",
  inherit = CNAEventClusterBase,
  private = list(
    .run = function() {

      if (!exists("cna_build_event_cluster_nest", mode = "function")) {
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

      is_multi <- TRUE

      self$results$instructions$setContent(
        cna_instructions_html("event_cluster")
      )

      if (is.null(self$options$buildModel_variables_long_action)) {
        self$results$errorText$setVisible(FALSE)
        return()
      }
      if (is.null(self$options$buildModel_variables_actor)) {
        self$results$errorText$setContent("Clustering requires an Actor variable.")
        self$results$errorText$setVisible(TRUE)
        return()
      }

      # Clustering gate
      if (!isTRUE(self$options$clustering_run)) {
        self$results$clusteringNote$setContent(
          '<div style="border: 2px solid #e6f4fe; border-radius: 15px; padding: 15px; background-color: #e6f4fe; margin-top: 10px;">
          <p style="text-align:center; font-weight: bold; color: #1a5276; font-size: 14px;">
          Check "Run Clustering Analysis" to start clustering.
          </p></div>')
        self$results$clusteringNote$setVisible(TRUE)
        return()
      }

      if (!self$options$buildModel_show_matrix &&
          !self$options$buildModel_show_plot &&
          !self$options$buildModel_show_histo) {
        self$results$buildModelTitle$setVisible(FALSE)
      }

      # ── Cluster sequences + build per-cluster co-occurrence networks ──
      model <- NULL
      if (self$results$buildModelContent$isFilled()) {
        model <- self$results$buildModelContent$state
      } else if (!is.null(self$data) && ncol(self$data) >= 1) {
        tryCatch({
          model <- cna_build_event_cluster_nest(self$data, self$options)
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
      }

      # ── Analyses (clusters behave exactly like groups) ──
      cna_run_centrality_nest(model, self$options, self$results, is_multi)
      cna_run_community_nest(model, self$options, self$results, is_multi)
      cna_run_bootstrap_nest(model, self$options, self$results, is_multi)
      cna_run_permutation_nest(model, self$options, self$results)
      cna_run_compare_nest(model, self$options, self$results, mode = "cluster")
    },

    # ── Plot dispatchers ──

    .showBuildModelPlot = function(image, ...) {
      cna_splot_model_nest(self$results$buildModelContent$state, self$options)
    },

    .showBuildModelHisto = function(image, ...) {
      if (!self$options$buildModel_show_histo) return(FALSE)
      cna_splot_histo_nest(self$results$buildModelContent$state, TRUE)
    },

    .showBuildModelFrequencies = function(image, ...) {
      if (!self$options$buildModel_show_frequencies) return(FALSE)
      cna_splot_freq_nest(self$results$buildModelContent$state, TRUE)
    },

    .showCentralityPlot = function(image, ...) {
      if (!self$options$centrality_show_plot) return(FALSE)
      cna_splot_centrality_nest(self$results$buildModelContent$state, self$options, TRUE)
    },

    .showCommunityPlot = function(image, ...) {
      if (!self$options$community_show_plot) return(FALSE)
      cna_splot_community_nest(self$results$community_plot$state, self$options, TRUE)
    },

    .showBootstrapPlot = function(image, ...) {
      if (!self$options$bootstrap_show_plot) return(FALSE)
      cna_splot_bootstrap_nest(self$results$bootstrap_plot$state)
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
