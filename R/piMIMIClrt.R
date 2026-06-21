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
#'   \item \code{DIF.Global} - data.frame: global.chi2, df, p.value
#'   \item \code{DIF.Uniforme} - data.frame: item, uDIF.chi2, df, p.value
#'   \item \code{DIF.NoUniforme} - data.frame: item, nuDIF.chi2, df, p.value
#'   \item \code{DeltaR2.uDIF} - data.frame: item, delta_R2
#'   \item \code{DeltaR2.nuDIF} - data.frame: item, delta_R2
#'   \item \code{fit} - unrestricted model (M3) for further use
#' }
#'
#' @examples
#' \dontrun{
#' library(psych)
#' data("bfi")
#' data.bfi <- bfi[, c("N1","N2","N3","N4","N5","gender")]
#' data.bfi <- data.bfi[complete.cases(data.bfi), ]
#' data.bfi$gender <- as.factor(data.bfi$gender)
#'
#' res <- piMIMIClrt(data = data.bfi,
#'                   items = c("N1","N2","N3","N4","N5"),
#'                   cov = "gender",
#'                   lvname = "Neuroticism",
#'                   est = "MLM")
#'
#' res$DIF.Global
#' res$DeltaR2.uDIF
#' }
#'
#' @importFrom lavaan cfa lavTestScore parameterEstimates inspect
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

  # ---- Base model syntax (without DIF effects for a specific item) ----
  base_syntax <- paste0(
    lvname, " =~ ", paste(items, collapse = " + "), "\n",
    cov_lat, " =~ ", cov, "\n",
    int_fac, " =~ ", paste0(items, ".", cov, collapse = " + "), "\n",
    paste(paste0(items, " ~~ ", items, ".", cov), collapse = "\n")
  )

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

  # Extract R² from M3 (unrestricted) - ensure named vector
  r2_m3 <- lavaan::inspect(fit_m3_full, "rsquare")
  if (is.matrix(r2_m3)) {
    r2_m3 <- as.vector(r2_m3)
    names(r2_m3) <- items
  }

  # ---- Helper function to fit a model with constraints for a given item ----
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

  # ---- Initialize output containers ----
  n_items <- length(items)

  # Tables for Score tests (same structure as piMIMIC)
  # We'll build them incrementally
  lrt_global <- data.frame(DIF.Global.Chi2 = NA_real_, df = NA_integer_, p.value = NA_real_)
  df_dif_uniforme <- data.frame(Item = items, uDIF.Chi2 = NA_real_, df = NA_integer_, p.value = NA_real_)
  df_dif_nouniforme <- data.frame(Item = items, nuDIF.Chi2 = NA_real_, df = NA_integer_, p.value = NA_real_)

  # Delta R² tables
  delta_r2_u <- data.frame(Item = items, delta_R2 = NA_real_)
  delta_r2_nu <- data.frame(Item = items, delta_R2 = NA_real_)

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

    # ---- Extract R² for this item ----
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

    # ---- Score tests (same logic as piMIMIC) ----
    # Global: add both parameters to M1
    params_global <- c(paste0(item, "~", cov_lat), paste0(item, "~", int_fac))
    score_global <- tryCatch(
      lavaan::lavTestScore(fit_m1, add = params_global)$test,
      error = function(e) NULL
    )
    if (!is.null(score_global) && nrow(score_global) > 0) {
      lrt_global[i, "DIF.Global.Chi2"] <- round(score_global[1, "X2"], 2)
      lrt_global[i, "df"] <- score_global[1, "df"]
      lrt_global[i, "p.value"] <- round(score_global[1, "p"], 3)
    }

    # Uniform: add only direct effect to M1
    score_unif <- tryCatch(
      lavaan::lavTestScore(fit_m1, add = paste0(item, "~", cov_lat))$test,
      error = function(e) NULL
    )
    if (!is.null(score_unif) && nrow(score_unif) > 0) {
      df_dif_uniforme[i, "uDIF.Chi2"] <- round(score_unif[1, "X2"], 2)
      df_dif_uniforme[i, "df"] <- score_unif[1, "df"]
      df_dif_uniforme[i, "p.value"] <- round(score_unif[1, "p"], 3)
    }

    # Non-uniform: add interaction to M2
    score_nu <- tryCatch(
      lavaan::lavTestScore(fit_m2, add = paste0(item, "~", int_fac))$test,
      error = function(e) NULL
    )
    if (!is.null(score_nu) && nrow(score_nu) > 0) {
      df_dif_nouniforme[i, "nuDIF.Chi2"] <- round(score_nu[1, "X2"], 2)
      df_dif_nouniforme[i, "df"] <- score_nu[1, "df"]
      df_dif_nouniforme[i, "p.value"] <- round(score_nu[1, "p"], 3)
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
    df_dif_uniforme <- add_oort(df_dif_uniforme, "df")
    df_dif_nouniforme <- add_oort(df_dif_nouniforme, "df")
  }

  # ---- Prepare output ----
  out <- list(
    DIF.Global = lrt_global,
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
  cat("PI-MIMIC LRT Results (Score tests with effect sizes)\n")
  cat("===================================================\n\n")
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
