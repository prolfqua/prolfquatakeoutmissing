# prolfquatakeoutmissing 0.1.0

* First version. Takes over `ContrastsMissing`, the `lm_missing` facade (`ContrastsLMMissingFacade`) and the group-mean imputation methods of `MissingHelpers` (now `MissingHelpersContrasts`), which were removed from prolfqua 1.9.0. Loading the package registers `method = "lm_missing"` for `prolfqua::build_contrast_analysis()`.
