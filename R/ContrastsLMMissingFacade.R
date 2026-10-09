#' LM + group-mean imputation contrast analysis facade
#'
#' Encapsulates the pipeline: \code{\link[prolfqua]{strategy_lm}} ->
#' \code{\link[prolfqua]{build_model}} -> \code{\link[prolfqua]{Contrasts}} ->
#' \code{\link[prolfqua]{ContrastsModerated}} -> merge with
#' \code{\link{ContrastsMissing}}.
#'
#' Proteins without a fitted model get their contrasts filled in from the
#' group-mean imputation method (\code{\link{ContrastsMissing}}).
#' Registered as \code{method = "lm_missing"} for
#' \code{\link[prolfqua]{build_contrast_analysis}} when the package is loaded.
#'
#' @return An R6 class generator.
#' @export
#' @examples
#' istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 30)
#' lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
#' lfqdata$rename_response("transformedIntensity")
#' contrasts <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
#' fa <- ContrastsLMMissingFacade$new(lfqdata, "~ group_", contrasts)
#' head(fa$get_contrasts())
#' fa$to_wide()
#' fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm_missing")
ContrastsLMMissingFacade <- R6::R6Class(
  "ContrastsLMMissingFacade",
  inherit = prolfqua::ContrastsFacadeBase,
  public = list(
    #' @description
    #' initialize
    #' @param lfqdata LFQData object
    #' @param modelstr model formula string (e.g. "~ group_")
    #' @param contrasts named character vector of contrasts
    #' @param weights column name for per-observation weights (default:
    #'   \code{lfqdata$nr_children_col()}). Pass \code{NULL} for unweighted.
    #' @param ... passed to \code{\link[prolfqua]{strategy_lm}}
    initialize = function(lfqdata, modelstr, contrasts, weights = lfqdata$nr_children_col(), ...) {
      full_formula <- private$.setup(lfqdata, contrasts, "lm_missing", modelstr)
      strat <- prolfqua::strategy_lm(full_formula, weights = weights, ...)
      self$model <- prolfqua::build_model(lfqdata, strat)
      base_contrast <- prolfqua::ContrastsModerated$new(prolfqua::Contrasts$new(self$model, contrasts))
      missing_contrast <- ContrastsMissing$new(lfqdata, contrasts = contrasts)
      self$contrast <- prolfqua::merge_contrasts_results(base_contrast, missing_contrast)$merged
    },
    #' @description get \code{\link[prolfqua]{ContrastsPlotter}} built from the stamped
    #'   facade output (so modelName is the facade key, not the merged leg names)
    #' @param fc_threshold fold change threshold
    #' @param fdr_threshold FDR threshold
    get_Plotter = function(fc_threshold = 1, fdr_threshold = 0.1) {
      prolfqua::ContrastsPlotter$new(
        self$get_contrasts(),
        subject_id = self$contrast$subject_id,
        fcthresh = fc_threshold,
        volcano = list(
          list(score = "p.value", thresh = fdr_threshold),
          list(score = "FDR", thresh = fdr_threshold)
        ),
        score = list(list(score = "statistic", thresh = 5))
      )
    }
  )
)
