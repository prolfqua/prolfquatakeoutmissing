# prolfquatakeoutmissing 0.1.1

* `lm_missing` weights by the child count times the identification confidence, like prolfqua 1.9.11's `lm` facade; with confidence 1 the results are unchanged.
* `ContrastsLMMissingFacade$new()` (`method = "lm_missing"`) takes `weights` and no further arguments, like prolfqua's `lm` facade; an argument it does not take is an error instead of being passed to `prolfqua::strategy_lm()`.
* Computes the contrasts of group means and the pooled variances with its own copies of `get_contrast()` and `poolvar()`, which prolfqua 1.9.5 deprecates.
* Requires prolfqua >= 1.9.5; examples and tests use its renamed simulator argument `n_proteins`.
* `lm_missing` registers with `needs = "aggregated"`, the name prolfqua 1.9.1 uses for facades on aggregated input; prolfqua 1.9.1 no longer accepts `"same"`. Earlier prolfqua versions treat both values alike.

# prolfquatakeoutmissing 0.1.0

* First version. Takes over `ContrastsMissing`, the `lm_missing` facade (`ContrastsLMMissingFacade`) and the group-mean imputation methods of `MissingHelpers` (now `MissingHelpersContrasts`), which were removed from prolfqua 1.9.0. Loading the package registers `method = "lm_missing"` for `prolfqua::build_contrast_analysis()`.
