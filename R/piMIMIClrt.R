#' @title PI-MIMIC with Effect Sizes (Delta R²) using Score Tests
#'
#' @description
#' Implements DIF detection using the PI-MIMIC framework. It provides significance tests
#' via Score tests (as in `piMIMIC`) and additionally calculates effect sizes (ΔR²)
#' for uniform and non-uniform DIF by comparing nested models.
#'
#' @inheritParams piMIMIC
#'
#' @return
#' A list with:
#' \itemize{
#'   \item \code{DIF.Global} - data.frame: global.chi2, df, p.value (same as `piMIMIC`)
#'   \item \code{DIF.Uniforme} - data.frame: item, uDIF.chi2, df, p.value
#'   \item \code{DIF.NoUniforme} - data.frame: item, nuDIF.chi2, df, p.value
#'   \item \code{DeltaR2.uDIF} - data.frame: item, delta_R2
#'   \item \code{DeltaR2.nuDIF} - data.frame: item, delta_R2
#'   \item \code{fit} - unrestricted model (M3) for further use
#' }
#'
#' @importFrom lavaan cfa lavTestScore parameterEstimates inspect
#' @importFrom scripty prods mimicparam
#' @export
piMIMIClrt <- function(data, items, cov, lvname = "LatFact",
                       est = "MLM", Oort.adj = FALSE, p.crit = 0.05) {

  # ---- Basic checks ----
  if (!is.character(est) || nchar(est) == 0) {
    stop("Error: Estimator must be a non-empty string (see lavaan documentation).")
  }
  if (!all(c(items, cov) %in% colnames(data))) {
    stop("Error: Some item names or the covariate are not present in the data frame.")
  }
  if (!is.character(lvname) || lvname == "") {
    stop("Error: The latent variable name ('lvname') must be a non-empty string.")
  }

  # Convert factor covariate to numeric if needed
  if (is.factor(data[[cov]])) {
    data[[cov]] <- as.numeric(data[[cov]])
  } else if (!is.numeric(data[[cov]])) {
    stop(paste("Error: The covariate", cov, "must be either a factor or numeric."))
  }

  # ---- Prepare data with product indicators ----
  df_pi <- scripty::prods(df = data[, c(items, cov)],
                          covname = cov,
                          match.op = FALSE,
                          doubleMC.op = TRUE)

  cov_lat <- paste0(cov, "lat")
  int_fac <- paste0("LFacX", cov)

  # ---- Base model syntax (without DIF effects for a specific item) ----
  base_syntax <- paste0(
    lvname, " =~ ", paste(items, collapse = " + "), "\n",
    cov_lat, " =~ ", cov, "\n",
    int_fac, " =~ ", paste0(items, ".", cov, collapse = " + "), "\n",
    paste(paste0(items, " ~~ ", items, ".", cov), collapse = "\n")
  )

  # ---- Fit unrestricted model (M3) ----
  fit_m3_full <- tryCatch(
    lavaan::cfa(base_syntax,
                data = df_pi,
                estimator = est,
                meanstructure = TRUE),
    error = function(e) NULL
  )

  if (is.null(fit_m3_full)) {
    stop("Unrestricted model (M3) did not converge.")
  }

  # ---- Extract R² from M3 ----
  r2_m3 <- lavaan::inspect(fit_m3_full, "rsquare")
  if (is.matrix(r2_m3)) {
    r2_m3 <- as.vector(r2_m3)
    names(r2_m3) <- items
  }

  # ---- Helper function to fit a model with constraints for a given item ----
  fit_item_model <- function(item, type = c("direct", "none")) {
    type <- match.arg(type)
    mod_syntax <- base_syntax

    if (type == "none") {
      constr <- paste0(item, " ~ 0*", cov_lat, " + 0*", int_fac)
    } else { # "direct"
      constr <- paste0(item, " ~ 0*", int_fac)
    }

    mod_syntax <- paste(mod_syntax, constr, sep = "\n")

    fit <- tryCatch(
      lavaan::cfa(mod_syntax,
                  data = df_pi,
                  estimator = est,
                  meanstructure = TRUE),
      error = function(e) NULL
    )
    return(fit)
  }

  # ---- Use piMIMIC logic for significance tests (Score tests from unrestricted model) ----
  # Get parameters to release (using scripty::mimicparam)
  mimic_param <- scripty::mimicparam(fit_m3_full)

  # Global test (both effects for each item)
  any.out <- do.call(rbind, lapply(mimic_param, function(x) {
    test <- tryCatch(lavaan::lavTestScore(fit_m3_full, add = x)$test,
                     error = function(e) NULL)
    if (!is.null(test) && nrow(test) > 0) {
      test[1, c("X2", "df", "p")]
    } else {
      data.frame(X2 = NA_real_, df = NA_integer_, p = NA_real_)
    }
  }))

  # Univariate tests (each parameter separately)
  sep.out <- tryCatch(
    lavaan::lavTestScore(fit_m3_full, add = as.character(mimic_param)),
    error = function(e) NULL
  )

  # Build DIF tables (same as piMIMIC)
  if (!is.null(sep.out) && !is.null(sep.out$uni)) {
    oddnum <- seq(1, length(sep.out$uni$lhs), 2)
    evennum <- seq(2, length(sep.out$uni$lhs), 2)

    df_dif_uniforme <- data.frame(
      Item = sep.out$uni$lhs[oddnum],
      uDIF.Chi2 = round(sep.out$uni$X2[oddnum], 2),
      df = sep.out$uni$df[oddnum],
      p.value = round(sep.out$uni$p.value[oddnum], 3)
    )

    df_dif_nouniforme <- data.frame(
      Item = sep.out$uni$lhs[evennum],
      nuDIF.Chi2 = round(sep.out$uni$X2[evennum], 2),
      df = sep.out$uni$df[evennum],
      p.value = round(sep.out$uni$p.value[evennum], 3)
    )
  } else {
    # Fallback if sep.out fails
    df_dif_uniforme <- data.frame(Item = items, uDIF.Chi2 = NA_real_, df = NA_integer_, p.value = NA_real_)
    df_dif_nouniforme <- data.frame(Item = items, nuDIF.Chi2 = NA_real_, df = NA_integer_, p.value = NA_real_)
  }

  # Global DIF table
  df_dif_global <- data.frame(
    DIF.Global.Chi2 = round(any.out$X2, 2),
    df = any.out$df,
    p.value = round(any.out$p, 3)
  )

  # ---- Calculate Delta R² by fitting restricted models ----
  delta_r2_u <- data.frame(Item = items, delta_R2 = NA_real_)
  delta_r2_nu <- data.frame(Item = items, delta_R2 = NA_real_)

  for (i in seq_along(items)) {
    item <- items[i]

    # Fit M1 (no effects) and M2 (only direct)
    fit_m1 <- fit_item_model(item, type = "none")
    fit_m2 <- fit_item_model(item, type = "direct")

    if (!is.null(fit_m1) && !is.null(fit_m2)) {
      # Extract R² for this item
      r2_m1_vec <- lavaan::inspect(fit_m1, "rsquare")
      r2_m2_vec <- lavaan::inspect(fit_m2, "rsquare")

      # Ensure scalar
      if (is.matrix(r2_m1_vec)) {
        r2_m1 <- r2_m1_vec[item, item]
      } else {
        r2_m1 <- r2_m1_vec[item]
      }
      if (is.matrix(r2_m2_vec)) {
        r2_m2 <- r2_m2_vec[item, item]
      } else {
        r2_m2 <- r2_m2_vec[item]
      }

      # Compute ΔR²
      delta_r2_u[i, "delta_R2"] <- as.numeric(r2_m2 - r2_m1)
      delta_r2_nu[i, "delta_R2"] <- as.numeric(r2_m3[item] - r2_m2)
    } else {
      warning(paste("Restricted models for item", item, "did not converge. Delta R² set to NA."))
    }
  }

  # ---- Apply Oort adjustment if requested ----
  if (Oort.adj) {
    baseline <- fit_m3_full@test$standard
    chi0 <- as.numeric(baseline[["stat"]])
    df0 <- as.numeric(baseline[["df"]])

    add_oort <- function(df_table, df_col = "df") {
      if (is.null(df_table) || nrow(df_table) == 0) return(df_table)
      df_vals <- df_table[[df_col]]
      K <- qchisq(1 - p.crit, df_vals)
      crit.Oort <- (chi0 / (K + df0 - 1)) * K
      df_table$crit.Oort <- round(crit.Oort, 3)
      df_table
    }

    df_dif_global <- add_oort(df_dif_global, "df")
    df_dif_uniforme <- add_oort(df_dif_uniforme, "df")
    df_dif_nouniforme <- add_oort(df_dif_nouniforme, "df")
  }

  # ---- Prepare output ----
  out <- list(
    DIF.Global = df_dif_global,
    DIF.Uniforme = df_dif_uniforme,
    DIF.NoUniforme = df_dif_nouniforme,
    DeltaR2.uDIF = delta_r2_u,
    DeltaR2.nuDIF = delta_r2_nu,
    fit = fit_m3_full
  )

  class(out) <- "piMIMIClrt"
  return(out)
}

#' @export
print.piMIMIClrt <- function(x, ...) {
  cat("PI-MIMIC LRT Results (Score tests + Delta R²)\n")
  cat("=============================================\n\n")
  cat("Global DIF (both effects):\n")
  print(x$DIF.Global)
  cat("\nUniform DIF (direct effect):\n")
  print(x$DIF.Uniforme)
  cat("\nNon-uniform DIF (interaction):\n")
  print(x$DIF.NoUniforme)
  cat("\nDelta R² for Uniform DIF:\n")
  print(x$DeltaR2.uDIF)
  cat("\nDelta R² for Non-uniform DIF:\n")
  print(x$DeltaR2.nuDIF)
  invisible(x)
}
