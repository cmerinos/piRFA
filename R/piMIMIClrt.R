#' @title PI-MIMIC with Likelihood Ratio Tests (LRT) for DIF Detection
#'
#' @description
#' Implements Differential Item Functioning (DIF) analysis using the Product of Indicators (PI)
#' approach within a Multiple-Indicators Multiple-Causes (MIMIC) framework, based on likelihood
#' ratio tests (LRT) comparing nested models. This function complements `piMIMICscore` (which
#' uses score tests) by providing LRT-based inference and additionally calculating effect sizes
#' (ΔR²) for uniform and non-uniform DIF.
#'
#' @param data DataFrame containing items and the covariate.
#' @param items Vector of item names within `data`.
#' @param cov Name of the covariate in `data` (can be categorical or numeric). If categorical, it must be a factor.
#' @param lvname Name for the latent variable in the model (default is `"LatFact"`).
#' @param est Abbreviation of the estimator to use (see lavaan documentation). Default is `"MLM"`.
#' @param Oort.adj Logical. If `TRUE`, applies Oort's critical value adjustment to the chi-square tests.
#'                 Default is `FALSE`.
#' @param p.crit Numeric. Significance level used to compute the chi-square critical value.
#'               Only used if `Oort.adj = TRUE`. Default is `0.05`.
#'
#' @return
#' A list with the following components:
#' \itemize{
#'   \item \code{LRT.Global} - data.frame with columns: Item, global.chi2, df, p.value.
#'   \item \code{LRT.Uniforme} - data.frame with columns: Item, uniforme.chi2, df, p.value.
#'   \item \code{LRT.NoUniforme} - data.frame with columns: Item, nouniforme.chi2, df, p.value.
#'   \item \code{DeltaR2.uDIF} - data.frame with columns: Item, delta_R2.
#'   \item \code{DeltaR2.nuDIF} - data.frame with columns: Item, delta_R2.
#'   \item \code{fit} - the unrestricted (Model 3) lavaan object, for further inspection or plotting.
#' }
#'
#' @details
#' The function fits three nested models for each item:
#' \itemize{
#'   \item \strong{Model 1 (M1)}: No effects (direct and interaction) are estimated for the studied item.
#'   \item \strong{Model 2 (M2)}: Only the direct effect (uniform DIF) is estimated; the interaction is fixed to zero.
#'   \item \strong{Model 3 (M3)}: Both direct and interaction effects are freely estimated (unrestricted model).
#' }
#' Then, likelihood ratio tests are performed:
#' \itemize{
#'   \item \strong{Global DIF}: M1 vs M3 (2 df)
#'   \item \strong{Uniform DIF}: M1 vs M2 (1 df)
#'   \item \strong{Non-uniform DIF}: M2 vs M3 (1 df)
#' }
#' Additionally, ΔR² is computed for uniform and non-uniform DIF:
#' \itemize{
#'   \item ΔR²\_uniforme = R²(M2) - R²(M1)   (variance explained by the direct effect)
#'   \item ΔR²\_nouniforme = R²(M3) - R²(M2) (variance explained by the interaction)
#' }
#'
#' If \code{Oort.adj = TRUE}, an adjusted critical value is reported alongside the p-values,
#' following the correction proposed by Oort (1992, 1998) to control Type I error inflation.
#'
#' @importFrom lavaan cfa lavTestLRT lavTestScore parameterEstimates inspect
#' @importFrom scripty prods
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

  # ---- Base model syntax (without DIF effects) ----
  base_syntax <- paste0(
    lvname, " =~ ", paste(items, collapse = " + "), "\n",
    cov_lat, " =~ ", cov, "\n",
    int_fac, " =~ ", paste0(items, ".", cov, collapse = " + "), "\n",
    paste(paste0(items, " ~~ ", items, ".", cov), collapse = "\n")
  )

  # ---- Function to fit a model with constraints for a given item ----
  fit_item_model <- function(item, type = c("full", "direct", "none")) {
    type <- match.arg(type)

    mod_syntax <- base_syntax

    if (type == "none") {
      constr <- paste0(item, " ~ 0*", cov_lat, " + 0*", int_fac)
    } else if (type == "direct") {
      constr <- paste0(item, " ~ 0*", int_fac)
    } else {
      constr <- NULL
    }

    if (!is.null(constr)) {
      mod_syntax <- paste(mod_syntax, constr, sep = "\n")
    }

    fit <- tryCatch(
      lavaan::cfa(mod_syntax,
                  data = df_pi,
                  estimator = est,
                  meanstructure = TRUE),
      error = function(e) NULL
    )
    return(fit)
  }

  # ---- Fit unrestricted model (M3) once ----
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

  # Extract R² from M3 (unrestricted) - ensure it's a named vector
  r2_m3 <- lavaan::inspect(fit_m3_full, "rsquare")
  # If it's a matrix or data.frame, convert to vector
  if (is.matrix(r2_m3) || is.data.frame(r2_m3)) {
    r2_m3 <- as.vector(r2_m3)
    names(r2_m3) <- items
  }

  # ---- Initialize output containers ----
  n_items <- length(items)

  # LRT tables (with correct column names)
  lrt_global <- data.frame(Item = items, global.chi2 = NA_real_, df = 2L, p.value = NA_real_)
  lrt_uniforme <- data.frame(Item = items, uniforme.chi2 = NA_real_, df = 1L, p.value = NA_real_)
  lrt_nouniforme <- data.frame(Item = items, nouniforme.chi2 = NA_real_, df = 1L, p.value = NA_real_)

  # Delta R² tables (single column)
  delta_r2_u <- data.frame(Item = items, delta_R2 = NA_real_)
  delta_r2_nu <- data.frame(Item = items, delta_R2 = NA_real_)

  # ---- Determine if we can use LRT or need Score tests ----
  # For non-ML estimators, lavTestLRT may not work; use lavTestScore instead
  use_score <- !grepl("^ML$", est, ignore.case = TRUE)

  # ---- Loop over items ----
  for (i in seq_along(items)) {
    item <- items[i]

    # Fit M1 (no effects) and M2 (only direct)
    fit_m1 <- fit_item_model(item, type = "none")
    fit_m2 <- fit_item_model(item, type = "direct")

    if (is.null(fit_m1) || is.null(fit_m2)) {
      warning(paste("Models for item", item, "did not converge. Skipping."))
      next
    }

    # Extract R² from M1 and M2 (ensure scalar)
    r2_m1_val <- lavaan::inspect(fit_m1, "rsquare")
    r2_m2_val <- lavaan::inspect(fit_m2, "rsquare")

    # Extract the value for the specific item
    if (is.matrix(r2_m1_val) || is.data.frame(r2_m1_val)) {
      r2_m1 <- r2_m1_val[item, item]  # Adjust if needed
    } else {
      r2_m1 <- r2_m1_val[item]
    }
    if (is.matrix(r2_m2_val) || is.data.frame(r2_m2_val)) {
      r2_m2 <- r2_m2_val[item, item]
    } else {
      r2_m2 <- r2_m2_val[item]
    }

    # Compute ΔR²
    delta_r2_u[i, "delta_R2"] <- r2_m2 - r2_m1
    delta_r2_nu[i, "delta_R2"] <- r2_m3[item] - r2_m2

    # ---- Perform tests ----
    if (use_score) {
      # Use Score tests (like piMIMICscore)
      # Global test: add both parameters
      params_global <- c(paste0(item, "~", cov_lat), paste0(item, "~", int_fac))
      score_global <- tryCatch(
        lavaan::lavTestScore(fit_m1, add = params_global)$test,
        error = function(e) NULL
      )
      if (!is.null(score_global) && nrow(score_global) > 0) {
        lrt_global[i, "global.chi2"] <- score_global[1, "X2"]
        lrt_global[i, "df"] <- score_global[1, "df"]
        lrt_global[i, "p.value"] <- score_global[1, "p"]
      }

      # Uniform test: add only direct effect (from M1)
      score_unif <- tryCatch(
        lavaan::lavTestScore(fit_m1, add = paste0(item, "~", cov_lat))$test,
        error = function(e) NULL
      )
      if (!is.null(score_unif) && nrow(score_unif) > 0) {
        lrt_uniforme[i, "uniforme.chi2"] <- score_unif[1, "X2"]
        lrt_uniforme[i, "df"] <- score_unif[1, "df"]
        lrt_uniforme[i, "p.value"] <- score_unif[1, "p"]
      }

      # Non-uniform test: add interaction (from M2)
      score_nu <- tryCatch(
        lavaan::lavTestScore(fit_m2, add = paste0(item, "~", int_fac))$test,
        error = function(e) NULL
      )
      if (!is.null(score_nu) && nrow(score_nu) > 0) {
        lrt_nouniforme[i, "nouniforme.chi2"] <- score_nu[1, "X2"]
        lrt_nouniforme[i, "df"] <- score_nu[1, "df"]
        lrt_nouniforme[i, "p.value"] <- score_nu[1, "p"]
      }

    } else {
      # Use LRT (only for ML estimator)
      # Global: M1 vs M3
      lrt_global_i <- tryCatch(
        lavaan::lavTestLRT(fit_m1, fit_m3_full, method = "satorra.bentler.2001"),
        error = function(e) NULL
      )
      if (!is.null(lrt_global_i) && nrow(lrt_global_i) >= 2) {
        lrt_global[i, "global.chi2"] <- lrt_global_i[2, "Chisq diff"]
        lrt_global[i, "df"] <- lrt_global_i[2, "Df diff"]
        lrt_global[i, "p.value"] <- lrt_global_i[2, "Pr(>Chisq)"]
      }

      # Uniform: M1 vs M2
      lrt_unif_i <- tryCatch(
        lavaan::lavTestLRT(fit_m1, fit_m2, method = "satorra.bentler.2001"),
        error = function(e) NULL
      )
      if (!is.null(lrt_unif_i) && nrow(lrt_unif_i) >= 2) {
        lrt_uniforme[i, "uniforme.chi2"] <- lrt_unif_i[2, "Chisq diff"]
        lrt_uniforme[i, "df"] <- lrt_unif_i[2, "Df diff"]
        lrt_uniforme[i, "p.value"] <- lrt_unif_i[2, "Pr(>Chisq)"]
      }

      # Non-uniform: M2 vs M3
      lrt_nounif_i <- tryCatch(
        lavaan::lavTestLRT(fit_m2, fit_m3_full, method = "satorra.bentler.2001"),
        error = function(e) NULL
      )
      if (!is.null(lrt_nounif_i) && nrow(lrt_nounif_i) >= 2) {
        lrt_nouniforme[i, "nouniforme.chi2"] <- lrt_nounif_i[2, "Chisq diff"]
        lrt_nouniforme[i, "df"] <- lrt_nounif_i[2, "Df diff"]
        lrt_nouniforme[i, "p.value"] <- lrt_nounif_i[2, "Pr(>Chisq)"]
      }
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

    lrt_global <- add_oort(lrt_global, "df")
    lrt_uniforme <- add_oort(lrt_uniforme, "df")
    lrt_nouniforme <- add_oort(lrt_nouniforme, "df")
  }

  # ---- Prepare output ----
  out <- list(
    LRT.Global = lrt_global,
    LRT.Uniforme = lrt_uniforme,
    LRT.NoUniforme = lrt_nouniforme,
    DeltaR2.uDIF = delta_r2_u,
    DeltaR2.nuDIF = delta_r2_nu,
    fit = fit_m3_full
  )

  class(out) <- "piMIMIClrt"
  return(out)
}

#' @export
print.piMIMIClrt <- function(x, ...) {
  cat("PI-MIMIC LRT Results\n")
  cat("=====================\n\n")
  cat("Global DIF (M1 vs M3):\n")
  print(x$LRT.Global)
  cat("\nUniform DIF (M1 vs M2):\n")
  print(x$LRT.Uniforme)
  cat("\nNon-uniform DIF (M2 vs M3):\n")
  print(x$LRT.NoUniforme)
  cat("\nDelta R² for Uniform DIF:\n")
  print(x$DeltaR2.uDIF)
  cat("\nDelta R² for Non-uniform DIF:\n")
  print(x$DeltaR2.nuDIF)
  invisible(x)
}
