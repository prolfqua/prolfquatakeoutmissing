# Helpers copied from prolfqua 1.9.5, which deprecates them; only this package used them.

#' Split a contrast "LHS - RHS" into its left and right side expressions
#' @param contrast one contrast expression string
#' @return list with `lhs`, `rhs` and `full` expressions, or NULL if the contrast is not a top-level difference
#' @noRd
.contrast_sides_expr <- function(contrast) {
  expr <- rlang::parse_expr(contrast)
  if (rlang::is_call(expr, "-") && length(expr) == 3L) {
    list(lhs = expr[[2]], rhs = expr[[3]], full = expr)
  } else {
    NULL
  }
}

#' Compute fold changes of named contrasts on a wide table with one column per group
#'
#' A contrast must be a difference "LHS - RHS"; `group_1` and `group_2` are the evaluated sides and `estimate` is
#' their difference. A later contrast may refer to an earlier one by name.
#' @param data wide data.frame, one column per group
#' @param hierarchy_keys hierarchy key columns
#' @param contrasts named character vector of contrasts
#' @return tibble with the `hierarchy_keys`, `group_1`, `group_2`, `estimate`, `group_1_name`, `group_2_name` and
#'   `contrast`
#' @noRd
.get_contrast <- function(data, hierarchy_keys, contrasts) {
  res <- vector(mode = "list", length(contrasts))
  names(res) <- names(contrasts)
  for (i in seq_along(contrasts)) {
    cname <- names(contrasts)[i]
    sx <- .contrast_sides_expr(contrasts[i])
    if (is.null(sx)) {
      rlang::abort(paste0(
        "get_contrast: contrast '",
        cname,
        "' (",
        contrasts[i],
        ") is not of the required form 'LHS - RHS'. ",
        "A contrast must be a difference of two group expressions, e.g. 'group_A - group_B' or ",
        "'(group_A + group_B)/2 - group_Ctrl'."
      ))
    }
    # a contrast naming a level absent from `data` errors here, once per contrast
    dd <- dplyr::mutate(data, group_1 = !!sx$lhs, group_2 = !!sx$rhs, estimate = !!sx$full)
    # keep this contrast as a column so a later contrast can refer to it
    data[[cname]] <- dd[["estimate"]]
    df <- dplyr::select(dd, dplyr::all_of(hierarchy_keys), "group_1", "group_2", "estimate")
    df$group_1_name <- rlang::as_label(sx$lhs)
    df$group_2_name <- rlang::as_label(sx$rhs)
    df$contrast <- cname
    res[[cname]] <- df
  }
  dplyr::ungroup(dplyr::bind_rows(res))
}

#' Pool the per-group variances of one protein
#' @param x per-group statistics of one protein (`nrMeasured`, `meanAbundance`, `var`)
#' @return one-row data.frame with the pooled variance, its degrees of freedom and means
#' @noRd
.compute_pooled <- function(x) {
  xm <- x |> dplyr::filter(.data$nrMeasured > 0)
  mean_all <- sum(xm$meanAbundance * xm$nrMeasured) / sum(xm$nrMeasured)
  x <- x |> dplyr::filter(.data$nrMeasured > 1)
  n <- x$nrMeasured
  pool.n <- sum(n)
  n.groups <- length(x$var)
  pool.var <- sum((n - 1) * x$var) / (pool.n - n.groups)
  pool.mean <- sum(x$meanAbundance * n) / pool.n
  data.frame(
    n.groups = n.groups,
    n = pool.n,
    df = pool.n - n.groups,
    sd = sqrt(pool.var),
    sdT = sqrt(pool.var * 2 / (pool.n / n.groups)),
    var = pool.var,
    mean = if (is.na(pool.mean)) mean_all else pool.mean,
    meanAll = mean_all,
    nrMeasured = sum(xm$nrMeasured)
  )
}

#' Pooled variance per protein from per-group statistics
#' @param res1 per-group statistics, one row per protein and group
#' @param config AnalysisConfiguration
#' @return one row per protein with the pooled statistics; the first factor column is set to "pooled"
#' @noRd
.poolvar <- function(res1, config) {
  resp <- tidyr::nest(res1, data = -dplyr::all_of(config$hierarchy_keys()))
  pooled <- dplyr::bind_rows(lapply(resp$data, .compute_pooled))
  resp$data <- NULL
  resp <- dplyr::bind_cols(resp, pooled)
  resp[[config$factor_keys()[1]]] <- "pooled"
  resp
}
