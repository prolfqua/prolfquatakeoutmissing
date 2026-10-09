#' Group-mean imputation and pooled-variance contrasts
#'
#' Extends \code{\link[prolfqua]{MissingHelpers}} (statistics per group and the
#' limit of detection) with the group-mean imputation and pooled-variance
#' t-test used by \code{\link{ContrastsMissing}}. Group means are weighted by
#' the LOD: \eqn{(LOD \cdot nrNAs + meanAbundance \cdot nrMeasured) / nrReplicates}.
#'
#' @return An R6 class generator.
#' @export
#' @examples
#' Contrasts <- c("group.b-a" = "group_A - group_B", "group.a-ctrl" = "group_A - group_Ctrl")
#' dd <- prolfqua::sim_lfq_data_protein_config(Nprot = 100, weight_missing = 2)
#' mh <- MissingHelpersContrasts$new(dd$data, dd$config, prob = 0.8, weighted = TRUE)
#' xx <- mh$get_stats()
#' xx <- mh$get_lod()
#' xx <- mh$impute_weighted_lod()
#' xx <- mh$impute_lod()
#' xx <- mh$get_poolvar()
#' bb <- mh$get_contrast_estimates(Contrasts)
#' mh$get_contrasts(Contrasts)
#'
#' dd <- prolfqua::sim_lfq_data_2factor_config(Nprot = 100, weight_missing = 0.1)
#'
#' Contrasts <- c("c1" = "TreatmentA - TreatmentB",
#'                "C2" = "BackgroundX- BackgroundZ",
#'                "c3" = "`TreatmentA:BackgroundX` - `TreatmentA:BackgroundZ`",
#'                "c4" = "`TreatmentB:BackgroundX` - `TreatmentB:BackgroundZ`"
#'                )
#' mh <- MissingHelpersContrasts$new(dd$data, dd$config, prob = 0.8, weighted = TRUE)
#' mh$get_stats()$interaction |> table()
#' mh$get_contrast_estimates(Contrasts)
#'
MissingHelpersContrasts <- R6::R6Class(
  "MissingHelpersContrasts",
  inherit = prolfqua::MissingHelpers,
  public = list(
    #' @field weighted should we weight the LOD
    weighted = TRUE,
    #' @description
    #' initialize
    #' @param data data
    #' @param config config
    #' @param prob default 0.5, median of groups with one observed value
    #' @param weighted should group average be computed used weighting, default TRUE.
    initialize = function(data, config, prob = 0.5, weighted = TRUE) {
      super$initialize(data, config, prob = prob)
      self$weighted <- weighted
    },
    #' @description
    #' compute group averages using weighted lod
    #'
    impute_weighted_lod = function() {
      toimp <- self$get_stats()
      toimp$meanAbundanceZero <- ifelse(is.na(toimp$meanAbundance), 0, toimp$meanAbundance)
      imputed_data <- toimp |>
        dplyr::mutate(
          meanAbundanceImp = (.data$nrMeasured * .data$meanAbundanceZero + .data$nrNAs * self$get_lod()) /
            .data$nrReplicates
        )
      return(imputed_data)
    },
    #' @description
    #' if group average absent substitute with LOD
    #'
    impute_lod = function() {
      toimp <- self$get_stats()
      toimp$meanAbundanceImp <- ifelse(is.na(toimp$meanAbundance), self$get_lod(), toimp$meanAbundance)
      return(toimp)
    },
    #' @description
    #' compute pooled var per protein
    #' @param prob prob of sd from proteins where it can be computed
    get_poolvar = function(prob = 0.75) {
      imputed_data <- private$.imputed()
      pooled <- prolfqua::poolvar(imputed_data, self$config)
      pooled <- dplyr::select(pooled, -dplyr::all_of(c(self$config$factor_keys_depth()[1], "var")))

      pooled_zero <- pooled[pooled$df > 0 & pooled$sd > 0, ]
      meandf <- pooled_zero |>
        dplyr::summarize(
          n = 1,
          df = 1,
          sd = stats::quantile(.data$sd, prob = prob, na.rm = TRUE),
          sdT = stats::quantile(.data$sdT, prob = prob, na.rm = TRUE)
        )

      minsd <- 1
      meandf$sd <- ifelse(meandf$sd > 0, meandf$sd, minsd)
      meandf$sdT <- ifelse(meandf$sdT > 0, meandf$sdT, minsd)

      pooled <- pooled |> dplyr::mutate(sd = ifelse(is.na(.data$sd) | .data$sd == 0, meandf$sd, .data$sd))
      pooled <- pooled |> dplyr::mutate(sdT = ifelse(is.na(.data$sdT) | .data$sdT == 0, meandf$sdT, .data$sdT))
      pooled <- pooled |> dplyr::mutate(df = ifelse(.data$df == 0, 1, .data$df))
      return(pooled)
    },
    #' @description
    #' get contrast estimates
    #' @param Contrasts named array with contrasts
    get_contrast_estimates = function(Contrasts) {
      lt <- private$.imputed()
      lt$is_missing <- ifelse(lt$nrNAs == lt$nrReplicates, 1, 0)
      hierarchy_keys <- self$config$hierarchy_keys()
      # contrasts of one per-group statistic, computed on its protein x group wide table
      contrast_of <- function(values_from) {
        wide <- tidyr::pivot_wider(
          lt,
          id_cols = dplyr::all_of(hierarchy_keys),
          names_from = "interaction",
          values_from = dplyr::all_of(values_from)
        )
        prolfqua::get_contrast(dplyr::ungroup(wide), hierarchy_keys, Contrasts)
      }

      imputed <- contrast_of("meanAbundanceImp")
      imputed$avgAbd <- (imputed$group_1 + imputed$group_2) / 2
      imputed <- imputed |>
        dplyr::rename(meanAbundanceImp_group_1 = "group_1", meanAbundanceImp_group_2 = "group_2")
      nrs <- contrast_of("is_missing") |>
        dplyr::select(dplyr::all_of(c(hierarchy_keys, "contrast", indic = "estimate")))
      nr_measured <- contrast_of("nrMeasured") |>
        dplyr::select(dplyr::all_of(c(
          hierarchy_keys,
          "contrast",
          nrMeasured_group_1 = "group_1",
          nrMeasured_group_2 = "group_2"
        )))
      imputed <- imputed |>
        dplyr::inner_join(nrs, by = c(hierarchy_keys, "contrast")) |>
        dplyr::inner_join(nr_measured, by = c(hierarchy_keys, "contrast"))

      # an estimate saying that a completely missing group is the more abundant one is set to 0
      imputed |>
        dplyr::mutate(
          estimate = ifelse(
            (.data$indic < 0 & .data$estimate < 0) | (.data$indic > 0 & .data$estimate > 0),
            0,
            .data$estimate
          )
        )
    },
    #' @description
    #' compute contrasts
    #' @param Contrasts vector with contrasts
    #' @param confint compute confint
    #' @param all return all columns, default FALSE
    get_contrasts = function(Contrasts, confint = 0.95, all = FALSE) {
      imputed <- self$get_contrast_estimates(Contrasts = Contrasts)
      pooled <- self$get_poolvar()
      result <- dplyr::inner_join(imputed, pooled, by = self$config$hierarchy_keys())
      result <- private$add_p_values(result, confint = confint)
      if (!all) {
        result <- dplyr::select(result, -dplyr::all_of(c("nrMeasured", "mean", "n.groups", "n", "meanAll")))
      }
      return(result)
    }
  ),
  private = list(
    .imputed = function() {
      if (self$weighted) self$impute_weighted_lod() else self$impute_lod()
    },
    add_p_values = function(result, confint = 0.95) {
      result <- dplyr::mutate(
        result,
        statistic = .data$estimate / .data$sdT,
        p.value = 2 * stats::pt(abs(.data$statistic), df = .data$df, lower.tail = FALSE)
      )
      prqt <- -stats::qt((1 - confint) / 2, df = result$df)
      result$conf.low <- result$estimate - prqt * (result$sdT)
      result$conf.high <- result$estimate + prqt * (result$sdT)
      return(result)
    }
  )
)
