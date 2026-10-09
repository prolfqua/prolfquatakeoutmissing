make_protein_lfqdata <- function(Nprot = 30, weight_missing = 0.2, seed = 42) {
  istar <- prolfqua::sim_lfq_data_protein_config(Nprot = Nprot, weight_missing = weight_missing, seed = seed)
  lfqdata <- prolfqua::LFQData$new(istar$data, istar$config)
  lfqdata$rename_response("transformedIntensity")
  lfqdata
}

CONTRASTS <- c("A_vs_Ctrl" = "group_A - group_Ctrl")
MODELSTR <- "~ group_"

test_that("MissingHelpersContrasts computes imputed contrasts", {
  dd <- prolfqua::sim_lfq_data_protein_config(Nprot = 100, weight_missing = 2)
  mh <- MissingHelpersContrasts$new(dd$data, dd$config, prob = 0.8, weighted = TRUE)
  expect_true(is.numeric(mh$get_lod()))
  expect_true("meanAbundanceImp" %in% colnames(mh$impute_weighted_lod()))
  expect_true("meanAbundanceImp" %in% colnames(mh$impute_lod()))
  res <- mh$get_contrasts(c("group.a-ctrl" = "group_A - group_Ctrl"))
  expect_true(all(c("estimate", "sdT", "p.value", "conf.low", "conf.high") %in% colnames(res)))
})

test_that("ContrastsMissing", {
  lProt <- make_protein_lfqdata(Nprot = 20, weight_missing = 0.4)
  csi <- ContrastsMissing$new(lProt, contrasts = c("AvsCtrl" = "group_A - group_Ctrl"))

  res <- csi$get_contrasts()
  expect_s3_class(res, "data.frame")
  expect_true(all(c("diff", "p.value", "FDR") %in% colnames(res)))
  expect_equal(res$std.error.unmoderated, res$std.error)
  expect_equal(res$df.unmoderated, res$df)
  expect_equal(sum(is.na(res$p.value)), 0)
  expect_true(all(res$estimate_type == "missing_fallback"))

  expect_s3_class(csi$get_contrast_sides(), "data.frame")
  expect_s3_class(csi$get_Plotter(), "ContrastsPlotter")
})

test_that("ContrastsMissing refuses non-aggregated data", {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, seed = 42)
  lpep <- prolfqua::LFQData$new(istar$data, istar$config)
  expect_error(ContrastsMissing$new(lpep, contrasts = CONTRASTS)$get_contrasts(), "aggregated")
})

test_that("ContrastsMissing merges with ContrastsLimma", {
  lProt <- make_protein_lfqdata(Nprot = 50, weight_missing = 0.3, seed = 1234)
  mod <- prolfqua::build_model_limma(lProt, prolfqua::strategy_limma("transformedIntensity ~ group_"))
  contr_limma <- prolfqua::ContrastsLimma$new(mod, CONTRASTS)
  csi <- ContrastsMissing$new(lProt, contrasts = CONTRASTS)

  merged <- prolfqua::merge_contrasts_results(contr_limma, csi)
  expect_true(inherits(merged$merged, "ContrastsTable"))
  expect_true(nrow(merged$merged$get_contrasts()) > 0)
})

test_that("ContrastsLMMissingFacade initialises and returns correct structure", {
  fa <- ContrastsLMMissingFacade$new(make_protein_lfqdata(), MODELSTR, CONTRASTS)
  expect_true(inherits(fa, "ContrastsLMMissingFacade"))
  expect_false(is.null(fa$model))
  expect_false(is.null(fa$contrast))
  expect_s3_class(fa$get_contrasts(), "data.frame")
  expect_s3_class(fa$to_wide(), "data.frame")
  expect_s3_class(fa$get_Plotter(), "ContrastsPlotter")
  expect_true(all(fa$get_contrasts()$modelName == "lm_missing"))
})

test_that("loading the package registers lm_missing for build_contrast_analysis", {
  entry <- prolfqua::lookup_facade("lm_missing")
  expect_equal(entry$package, "prolfquatakeoutmissing")
  fa <- prolfqua::build_contrast_analysis(make_protein_lfqdata(), MODELSTR, CONTRASTS, method = "lm_missing")
  expect_true(inherits(fa, "ContrastsLMMissingFacade"))
})

test_that("ContrastsLMMissingFacade requires aggregated data", {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 10, seed = 42)
  lpep <- prolfqua::LFQData$new(istar$data, istar$config)
  expect_error(ContrastsLMMissingFacade$new(lpep, MODELSTR, CONTRASTS), "requires aggregated LFQData")
})

test_that("lm_missing leaves fewer proteins missing than lm", {
  istar <- prolfqua::sim_lfq_data_peptide_config(Nprot = 80, seed = 42)
  lfq_peptide <- prolfqua::LFQData$new(istar$data, istar$config)
  lfq_peptide <- lfq_peptide$get_Transformer()$log2()$lfq
  lfq_protein <- prolfqua::AggregateMedpolish$new(lfq_peptide, "protein")$aggregate()
  contrasts2 <- c("A_vs_Ctrl" = "group_A - group_Ctrl", "B_vs_Ctrl" = "group_B - group_Ctrl")

  fa_lm <- prolfqua::ContrastsLMFacade$new(lfq_protein, MODELSTR, contrasts2)
  fa_lm_missing <- ContrastsLMMissingFacade$new(lfq_protein, MODELSTR, contrasts2)
  expect_true(nrow(fa_lm_missing$get_missing()) <= nrow(fa_lm$get_missing()))

  n_expected <- dplyr::n_distinct(lfq_protein$data_long()$protein_Id) * length(contrasts2)
  expect_equal(nrow(fa_lm_missing$get_contrasts()) + nrow(fa_lm_missing$get_missing()), n_expected)
})
