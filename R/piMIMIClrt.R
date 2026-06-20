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
#'   \item \code{DeltaR2.uDIF} - data.frame with columns: Item, delta_R2_uniforme.
#'   \item \code{DeltaR2.nuDIF} - data.frame with columns: Item, delta_R2_nouniforme.
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
#' @references
#' Kolbe, L., & Jorgensen, T. D. (2018). Using product indicators in restricted factor analysis
#' models to detect nonuniform measurement bias. In *Quantitative Psychology* (pp. 235–245). Springer.
#'
#' Oort, F. J. (1992). Using restricted factor analysis to detect item bias. *Methodika*, 6, 150–160.
#'
#' Oort, F. J. (1998). Simulation study of item bias detection with restricted factor analysis.
#' *Structural Equation Modeling*, 5, 107–124.
#'
#' @examples
#' \dontrun{
#' library(psych)
#' data("bfi")
#' data.bfi <- bfi[, c("N1","N2","N3","N4","N5","gender")]
#' data.bfi <- data.bfi[complete.cases(data.bfi), ]
#' data.bfi$gender <- as.factor(data.bfi$gender)
#'
#' res.lrt <- piMIMIClrt(data = data.bfi,
#'                       items = c("N1","N2","N3","N4","N5"),
#'                       cov = "gender",
#'                       lvname = "Neuroticism",
#'                       est = "MLM")
#'
#' res.lrt$LRT.Global
#' res.lrt$DeltaR2.uDIF
#' }
#'
#' @importFrom lavaan cfa lavTestLRT parameterEstimates
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
  # This is the skeleton we will modify for each item.
  base_syntax <- paste0(
    lvname, " =~ ", paste(items, collapse = " + "), "\n",
    cov_lat, " =~ ", cov, "\n",
    int_fac, " =~ ", paste0(items, ".", cov, collapse = " + "), "\n",
    paste(paste0(items, " ~~ ", items, ".", cov), collapse = "\n")
  )

  # ---- Function to fit a model with constraints for a given item ----
  #   type: "none"    -> both direct and interaction fixed to 0 (M1)
  #         "direct"  -> interaction fixed to 0 (M2)
  #         "full"    -> no constraints (M3)
  fit_item_model <- function(item, type = c("full", "direct", "none")) {
    type <- match.arg(type)

    # Start with base syntax
    mod_syntax <- base_syntax

    # Add constraints for the studied item
    if (type == "none") {
      # Fix both direct and interaction to 0
      constr <- paste0(item, " ~ 0*", cov_lat, " + 0*", int_fac)
    } else if (type == "direct") {
      # Only interaction fixed to 0
      constr <- paste0(item, " ~ 0*", int_fac)
    } else { # "full"
      # No constraints (free both)
      constr <- NULL
    }

    if (!is.null(constr)) {
      mod_syntax <- paste(mod_syntax, constr, sep = "\n")
    }

    # Fit the model
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
  fit_m3 <- fit_item_model(items[1], type = "full")  # placeholder, will be refitted per item
  # But we need the unrestricted model without per-item constraints; actually it's the same for all items.
  # Fit it once with the base syntax (no constraints).
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

  # Extract R² from M3 (unrestricted)
  r2_m3 <- lavaan::inspect(fit_m3_full, "rsquare")[items]

  # ---- Initialize output containers ----
  n_items <- length(items)
  lrt_global <- data.frame(Item = items, global.chi2 = NA_real_, df = 2L, p.value = NA_real_)
  lrt_uniforme <- data.frame(Item = items, uniforme.chi2 = NA_real_, df = 1L, p.value = NA_real_)
  lrt_nouniforme <- data.frame(Item = items, nouniforme.chi2 = NA_real_, df = 1L, p.value = NA_real_)
  delta_r2_u <- data.frame(Item = items, delta_R2_uniforme = NA_real_)
  delta_r2_nu <- data.frame(Item = items, delta_R2_nouniforme = NA_real_)

  # ---- Loop over items ----
  for (i in seq_along(items)) {
    item <- items[i]

    # Fit M1 (no effects) and M2 (only direct)
    fit_m1 <- fit_item_model(item, type = "none")
    fit_m2 <- fit_item_model(item, type = "direct")

    # If either fails, skip to next item
    if (is.null(fit_m1) || is.null(fit_m2)) {
      warning(paste("Models for item", item, "did not converge. Skipping."))
      next
    }

    # Extract R² from M1 and M2
    r2_m1 <- lavaan::inspect(fit_m1, "rsquare")[item]
    r2_m2 <- lavaan::inspect(fit_m2, "rsquare")[item]

    # Compute ΔR²
    delta_r2_u[i] <- r2_m2 - r2_m1          # uniform: direct effect added
    delta_r2_nu[i] <- r2_m3[item] - r2_m2   # non-uniform: interaction added

    # Likelihood Ratio Tests
    # Global: M1 vs M3 (2 df)
    lrt_global_i <- tryCatch(
      lavaan::lavTestLRT(fit_m1, fit_m3_full, method = "default"),
      error = function(e) NULL
    )
    if (!is.null(lrt_global_i) && nrow(lrt_global_i) >= 2) {
      # The second row corresponds to the difference
      lrt_global$global.chi2[i] <- lrt_global_i[2, "Chisq diff"]
      lrt_global$df[i] <- lrt_global_i[2, "Df diff"]
      lrt_global$p.value[i] <- lrt_global_i[2, "Pr(>Chisq)"]
    }

    # Uniform: M1 vs M2 (1 df)
    lrt_unif_i <- tryCatch(
      lavaan::lavTestLRT(fit_m1, fit_m2, method = "default"),
      error = function(e) NULL
    )
    if (!is.null(lrt_unif_i) && nrow(lrt_unif_i) >= 2) {
      lrt_uniforme$uniforme.chi2[i] <- lrt_unif_i[2, "Chisq diff"]
      lrt_uniforme$df[i] <- lrt_unif_i[2, "Df diff"]
      lrt_uniforme$p.value[i] <- lrt_unif_i[2, "Pr(>Chisq)"]
    }

    # Non-uniform: M2 vs M3 (1 df)
    lrt_nounif_i <- tryCatch(
      lavaan::lavTestLRT(fit_m2, fit_m3_full, method = "default"),
      error = function(e) NULL
    )
    if (!is.null(lrt_nounif_i) && nrow(lrt_nounif_i) >= 2) {
      lrt_nouniforme$nouniforme.chi2[i] <- lrt_nounif_i[2, "Chisq diff"]
      lrt_nouniforme$df[i] <- lrt_nounif_i[2, "Df diff"]
      lrt_nouniforme$p.value[i] <- lrt_nounif_i[2, "Pr(>Chisq)"]
    }
  }

  # ---- Apply Oort adjustment if requested ----
  if (Oort.adj) {
    # Obtain baseline chi2 and df from M3
    baseline <- fit_m3_full@test$standard
    chi0 <- as.numeric(baseline[["stat"]])
    df0 <- as.numeric(baseline[["df"]])

    # Function to add adjusted critical value to a table
    add_oort <- function(df_table, df_col = "df") {
      if (is.null(df_table) || nrow(df_table) == 0) return(df_table)
      df_vals <- df_table[[df_col]]
      # K is the original critical value at p.crit
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

# S3 print method (optional, but nice)
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
