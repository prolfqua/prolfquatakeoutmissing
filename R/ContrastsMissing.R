#' Compute contrasts with group mean imputation
#'
#' If there are no observations in one of the groups for some of the
#' proteins, the group mean cannot be estimated. Therefore, assuming
#' that the observation is missing because the protein abundance is
#' below the detection limit, we substitute the unobserved group with
#' the median of protein abundances observed only in one sample of the
#' group. The variance of a protein is estimated using the pooled
#' variance of all observations of all groups.
#'
#' This is a group-mean substitution, not a per-protein model fit. prolfqua
#' offers \code{\link[prolfqua]{build_model_impute}} and the \code{lm_impute}
#' facade, which refit failed proteins with LOD-imputed values.
#'
#' @seealso \code{\link{MissingHelpersContrasts}}, \code{\link{ContrastsLMMissingFacade}}
#' @return An R6 class generator.
#' @export
#' @examples
#' istar <- prolfqua::sim_lfq_data_protein_config(Nprot = 120, weight_missing = .4)
#' lProt <- prolfqua::LFQData$new(istar$data, istar$config)
#' lProt$rename_response("transformedIntensity")
#'
#' Contr <- c("dil.b_vs_a" = "group_A - group_Ctrl")
#' csi <- ContrastsMissing$new(lProt, contrasts = Contr)
#' csi$get_contrast_sides()
#'
#' res <- csi$get_contrasts()
#' stopifnot(nrow(res) == length(unique(istar$data$protein_Id)))
#' stopifnot(sum(is.na(res$p.value)) == 0)
#' csi$column_description()
#' p <- csi$get_Plotter()$volcano()
#'
#' dd <- prolfqua::sim_lfq_data_2factor_config(Nprot = 100, weight_missing = 0.1)
#' Contrasts <- c("c1" = "TreatmentA - TreatmentB",
#'                "C2" = "BackgroundX- BackgroundZ",
#'                "c3" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
#'                "c4" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`"
#'                )
#' lProt <- prolfqua::LFQData$new(dd$data, dd$config)
#' lProt$rename_response("transformedIntensity")
#' csi <- ContrastsMissing$new(lProt, contrasts = Contrasts)
#' res <- csi$get_contrasts()
ContrastsMissing <- R6::R6Class(
  "ContrastsMissing",
  inherit = prolfqua::ContrastsInterface,
  public = list(
    #' @field subject_id subject_id e.g. protein_ID column
    subject_id = character(),
    #' @field contrasts array with contrasts (see example)
    contrasts = character(),
    #' @field model_name model name
    model_name = character(),
    #' @field contrast_result data frame with results of contrast computation
    contrast_result = NULL,
    #' @field lfqdata data frame
    lfqdata = NULL,
    #' @field confint confidence interval
    confint = 0.95,
    #' @field p.adjust function to adjust p-values
    p.adjust = NULL,
    #' @description
    #' initialize
    #' @param lfqdata LFQData
    #' @param contrasts array of contrasts (see example)
    #' @param confint confidence interval
    #' @param p.adjust method for p-value adjustment - default Benjamini Hochberg
    #' @param model_name default "groupAverage"
    initialize = function(
      lfqdata,
      contrasts,
      confint = 0.95,
      p.adjust = prolfqua::adjust_p_values,
      model_name = "groupAverage"
    ) {
      self$subject_id <- lfqdata$relevant_hierarchy_keys()
      self$contrasts <- contrasts
      self$model_name <- model_name
      self$lfqdata <- lfqdata
      self$confint <- confint
      self$p.adjust <- p.adjust
    },
    #' @description
    #' get contrasts sides
    #'
    get_contrast_sides = function() {
      tt <- self$contrasts[grep("-", self$contrasts)]
      tt <- tibble::tibble(contrast = names(tt), rhs = tt)
      tt |>
        dplyr::mutate(rhs = gsub("[` ]", "", .data$rhs)) |>
        tidyr::separate("rhs", c("group_1", "group_2"), sep = "-")
    },
    #' @description
    #' table with results of contrast computation
    #' @param all FALSE, do not show all columns (default)
    get_contrasts = function(all = FALSE) {
      if (is.null(self$contrast_result)) {
        config <- self$lfqdata$get_config()
        if (config$hierarchy_depth < length(config$hierarchy_keys())) {
          rlang::abort(c(
            "Invalid argument.",
            x = "`lfqdata` must be aggregated first (hierarchy_depth must equal the number of hierarchy keys).",
            i = sprintf(
              "Actual value: hierarchy_depth = %s; hierarchy_keys = %s.",
              config$hierarchy_depth,
              length(config$hierarchy_keys())
            )
          ))
        }
        mh1 <- MissingHelpersContrasts$new(self$lfqdata$data_long(), config, prob = 0.5, weighted = TRUE)
        result <- mh1$get_contrasts(Contrasts = self$contrasts, confint = self$confint, all = all)
        result <- self$p.adjust(result, column = "p.value", group_by_col = "contrast", newname = "FDR")
        result <- result |> dplyr::rename(diff = "estimate", sigma = "sd", std.error = "sdT")
        result <- result |>
          dplyr::mutate(
            std.error.unmoderated = .data$std.error,
            df.unmoderated = .data$df
          )
        # group-mean substitution, not a model fit: every row is a fallback
        result$estimate_type <- "missing_fallback"
        result <- dplyr::mutate(result, modelName = self$model_name, .before = 1)
        result <- dplyr::relocate(result, "estimate_type", .after = "modelName")
        self$contrast_result <- dplyr::ungroup(result)
      }
      res <- self$contrast_result
      stopifnot(all(super$column_description()$column_name %in% colnames(res)))
      invisible(res)
    },
    #' @description
    #' get ContrastsPlotter
    #' @return Contrast_Plotter
    get_Plotter = function() {
      prolfqua::ContrastsPlotter$new(
        self$get_contrasts(),
        subject_id = self$subject_id,
        volcano = list(list(score = "p.value", thresh = 0.1), list(score = "FDR", thresh = 0.1))
      )
    }
  )
)
