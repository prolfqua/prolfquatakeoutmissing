# prolfquatakeoutmissing

Group-mean imputation contrasts taken out of [prolfqua](https://github.com/fgcz/prolfqua) 1.9.0.

- `ContrastsMissing`: group-mean substitution with a pooled-variance t-test
- `ContrastsLMMissingFacade`: `lm` contrasts, with proteins lacking a fit filled in from `ContrastsMissing`
- `MissingHelpersContrasts`: the imputation and pooled-variance helpers behind them

Loading the package registers `method = "lm_missing"`:

```r
library(prolfquatakeoutmissing)
fa <- prolfqua::build_contrast_analysis(lfqdata, "~ group_", contrasts, method = "lm_missing")
```

For new analyses prefer prolfqua's `lm_impute` facade, which refits failed proteins with LOD-imputed values.

## Install

```bash
make install
```
