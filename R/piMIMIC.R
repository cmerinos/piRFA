#' @title DIF analysis using PI-MIMIC with Score Test (Oort adjustment optional)
#'
#' @description
#' Implements the product indicator (PI) approach for MIMIC models to detect
#' uniform and non‑uniform DIF using the score test (modification indices).
#' Optionally applies Oort's critical value adjustment to control Type I error.
#'
#' @param data Data frame containing items and the covariate.
#' @param items Character vector of item names.
#' @param cov Name of the covariate (numeric or factor).
#' @param lvname Name of the latent variable (default "LatFact").
#' @param est Estimator for lavaan (default "MLM").
#' @param Oort.adj Logical; if `TRUE`, applies Oort's adjustment to the critical value.
#' @param p.crit Numeric; significance level for the Oort adjustment.
#'
#' @return A list with:
#' \item{DIF.Global}{Data frame with global DIF test (2 df): Chi², df, p-value, and Oort critical value if requested.}
#' \item{DIF.Uniforme}{Data frame with uniform DIF test (1 df).}
#' \item{DIF.NoUniforme}{Data frame with non‑uniform DIF test (1 df).}
#' \item{SEPC.uDIF}{Standardized expected parameter change for uniform DIF.}
#' \item{SEPC.nuDIF}{Standardized expected parameter change for non‑uniform DIF.}
#' \item{fit}{The fitted lavaan object.}
#'
#' @details
#' This function implements Differential Item Functioning (DIF) analysis using
#' the Product of Indicators approach (PI; Kolbe & Jorgensen, 2018) within the
#' Restricted Factor Analysis (RFA; Oort, 1998) framework. This method operates under a
#' MIMIC scheme (Finch, 2005), incorporating latent variable interactions using PI (Kolbe et al., 2018, 2019;
#' Kolbe, Jorgensen, & Molenaar, 2020). It allows for the evaluation of uniform (uDIF) and non-uniform DIF (nuDIF)
#' with covariates that can be categorical (e.g., sex) or continuous (e.g., self-esteem, conscientiousness).
#'
#' Estimation is performed via `lavaan::cfa`, and DIF statistical tests are based
#' on the Score test (LRT). By default, chi-square tests compare an unrestricted model
#' (with covariate effect and interaction parameters freely estimated) and a restricted model
#' (with these parameters fixed to zero), using the standard chi-square distribution.
#'
#' Simulation studies (Oort, 1992, 1998; Kim, Yoon & Lee, 2011) have shown that
#' the LR test under MIMIC may suffer from inflated Type I error rates. To address this, Oort proposed
#' a correction of the chi-square critical value:
#'
#' \deqn{K' = (χ²₀ / (K + df₀ - 1)) * K}
#' where χ²₀ and df₀ are from the baseline (full invariance) model, and K is
#' the original critical value. This adjustment is recommended when the baseline
#' model shows evidence of misfit (χ²₀/df₀ > 1), as it helps control Type I error.
#'
#' where:
#' \itemize{
#'   \item \eqn{K'} is the adjusted critical value,
#'   \item \eqn{K} is the original critical value from the chi-square distribution at significance level \eqn{p.crit},
#'   \item \eqn{\chi^2_0} is the chi-square statistic of the baseline model,
#'   \item \eqn{df_0} is the corresponding degrees of freedom of the baseline model.
#' }
#'
#' Interpretation:
#' \itemize{
#'   \item The reported \code{p.value} corresponds to the standard chi-square test.
#'   \item When \code{Oort.adj = TRUE}, the function also reports the adjusted critical value (\code{crit.Oort}).
#'   \item Users can compare the observed chi-square statistic against both thresholds:
#'         the conventional critical value (via \code{p.value}) and the Oort-adjusted critical value.
#'   \item This comparison allows assessment of how conclusions may differ when controlling for
#'         potential inflation of Type I error in MIMIC-PI models.
#'
#' @examples
#' ### Example 1: simulated data -------------
#' set.seed(123)
#' Exmp1.data <- data.frame(
#'   grp = sample(0:1, 100, replace = TRUE),  # Group variable
#'   item1 = sample(1:5, 100, replace = TRUE),
#'   item2 = sample(1:5, 100, replace = TRUE),
#'   item3 = sample(1:5, 100, replace = TRUE),
#'   item4 = sample(1:5, 100, replace = TRUE))
#'
#' res1 <- piMIMIC(data = Exmp1.data , items = c("item1","item2","item3"), cov = "grp")
#' res1$DIF.Global
#'
#' ### Example 2: Using the 'bfi' dataset from the 'psych' package -------------
#' library(psych)
#' data("bfi")
#'
#' data.bfi <- bfi[, c("N1","N2","N3","N4","N5","gender")]
#' data.bfi <- data.bfi[complete.cases(data.bfi), ]
#' data.bfi$gender <- as.factor(data.bfi$gender)
#'
#' neuro.items <- c("N1","N2","N3","N4","N5")
#'
#' # Run DIF analysis with Oort adjustment
#' res.bfi <- piMIMIC(data = data.bfi, items = neuro.items, cov = "gender",
#'                  lvname = "Neuroticism", est = "MLM",
#'                  Oort.adj = TRUE, p.crit = 0.05)
#' res.bfi$DIF.Global
#'
#' @references
#' Kim, E. S., Yoon, M., & Lee, T. (2011). Testing Measurement Invariance Using MIMIC:
#' Likelihood Ratio Test With a Critical Value Adjustment. *Educational and Psychological Measurement, 72*(3), 469–492.
#' https://doi.org/10.1177/0013164411427395
#'
#' Oort, F. J. (1992). Using restricted factor analysis to detect item bias. *Psychological Methods*, 37, 547–567.
#'
#' Oort, F. J. (1998). Simulation study of item bias detection with restricted factor analysis.
#' *Structural Equation Modeling*, 5, 107–124.
#'
#' French, B. F., & Finch, W. H. (2008). Multigroup confirmatory factor analysis: Locating the invariant referent variables.
#' *Structural Equation Modeling*, 15(1), 96–113.
#'
#' Stark, S., Chernyshenko, O. S., & Drasgow, F. (2006). Detecting differential item functioning with confirmatory factor analysis and item response theory: Toward a unified strategy.
#' *Journal of Applied Psychology*, 91(6), 1292–1306.
#'
#' Kolbe, L., & Jorgensen, T. D. (2018). Using product indicators in restricted factor analysis models to detect nonuniform measurement bias.
#' In *Quantitative Psychology* (pp. 235–245). Springer.
#'
#' Kolbe, L., & Jorgensen, T. D. (2019). Using restricted factor analysis to select anchor items and detect differential item functioning.
#' *Behavior Research Methods*, 51, 138–151.
#'
#' Kolbe, L., Jorgensen, T. D., & Molenaar, D. (2020). The Impact of Unmodeled Heteroskedasticity on Assessing Measurement Invariance in Single-group Models.
#' *Structural Equation Modeling*, 28(1), 82–98.
#'
#' Whittaker, T. A. (2012). Estimation of Standardized Expected Parameter Change for DIF Detection.
#' *Educational and Psychological Measurement*, 72(3), 342-357.
#'
#' Garnier-Villarreal, M., & Jorgensen, T. D. (2024). Evaluating Local Model Misspecification with Modification Indices in Bayesian Structural Equation Modeling.
#' *Structural Equation Modeling*, 1–15.
#'
#' @importFrom lavaan cfa lavTestScore parameterestimates
#' @importFrom scripty prods mimicparam
#'
#' @export
piMIMIC <- function(data, items, cov, lvname = "LatFact", est = "MLM",
                    Oort.adj = FALSE, p.crit = 0.05) {

  if (!is.character(est) || nchar(est) == 0) {
    stop("Error: Estimator must be a non-empty string (see lavaan documentation).")
  }

  cov_lat <- paste0(cov, "lat")

  if (!all(c(items, cov) %in% colnames(data))) {
    stop("Error: Some item names or the covariate are not present in the data frame.")
  }
  if (!is.character(lvname) || lvname == "") {
    stop("Error: The latent variable name ('lvname') must be a non-empty string.")
  }

  if (is.factor(data[[cov]])) {
    data[[cov]] <- as.numeric(data[[cov]])
  } else if (!is.numeric(data[[cov]])) {
    stop(paste("Error: The covariate", cov, "must be either a factor or numeric."))
  }

  df_pi <- scripty::prods(df = data[, c(items, cov)],
                          covname = cov,
                          match.op = FALSE,
                          doubleMC.op = TRUE)

  model_mimic <- paste0(
    lvname, " =~ ", paste(items, collapse = " + "), "\n",
    cov_lat, " =~ ", cov, "\n",
    paste0("LFacX", cov, " =~ ", paste0(items, ".", cov, collapse = " + "), "\n"),
    paste(paste0(items, " ~~ ", items, ".", cov), collapse = "\n")
  )

  fit <- lavaan::cfa(model_mimic,
                     data = df_pi,
                     estimator = est,
                     meanstructure = TRUE)

  mimic_param <- scripty::mimicparam(fit)

  # --- Función interna para extraer resultados ---
  mimicout_modificado <- function(fit.mimic, mimic.param, cov, Oort.adj, p.crit) {
    ests <- as.data.frame(lavaan::parameterestimates(fit.mimic))
    uniqnames <- unique(ests$lhs)
    lvname <- uniqnames[1]

    any.out <- do.call(rbind, lapply(mimic.param, function(x)
      lavaan::lavTestScore(fit.mimic, add = x)$test))

    sep.out <- lavaan::lavTestScore(fit.mimic, add = as.character(mimic.param))

    # Baseline chi2 and df
    baseline <- fit.mimic@test$standard
    chi0 <- as.numeric(baseline[["stat"]])
    df0  <- as.numeric(baseline[["df"]])

    # Oort adjustment (critical values)
    if (Oort.adj) {
      K_global <- qchisq(1 - p.crit, 2)
      K_uniform <- qchisq(1 - p.crit, 1)
      crit.global <- (chi0 / (K_global + df0 - 1)) * K_global
      crit.uniform <- (chi0 / (K_uniform + df0 - 1)) * K_uniform
    }

    # DIF.Global
    df_dif_global <- data.frame(
      Item = items,
      Chi2 = round(any.out$X2, 3),
      df = any.out$df,
      p.value = round(any.out$p, 4)
    )
    if (Oort.adj) df_dif_global$crit.Oort <- round(crit.global, 3)

    # DIF.Uniforme
    oddnum <- seq(1, length(sep.out$uni$lhs), 2)
    df_dif_uniforme <- data.frame(
      Item = gsub("~.*$", "", sep.out$uni$lhs[oddnum]),
      Chi2 = round(sep.out$uni$X2[oddnum], 3),
      df = sep.out$uni$df[oddnum],
      p.value = round(sep.out$uni$p.value[oddnum], 4)
    )
    if (Oort.adj) df_dif_uniforme$crit.Oort <- round(crit.uniform, 3)

    # DIF.NoUniforme
    evennum <- seq(2, length(sep.out$uni$lhs), 2)
    df_dif_nouniforme <- data.frame(
      Item = gsub("~.*$", "", sep.out$uni$lhs[evennum]),
      Chi2 = round(sep.out$uni$X2[evennum], 3),
      df = sep.out$uni$df[evennum],
      p.value = round(sep.out$uni$p.value[evennum], 4)
    )
    if (Oort.adj) df_dif_nouniforme$crit.Oort <- round(crit.uniform, 3)

    # SEPC
    sepc_values <- lavaan::lavTestScore(fit.mimic,
                                        add = as.character(mimic.param),
                                        univariate = TRUE,
                                        standardized = TRUE,
                                        cov.std = TRUE,
                                        epc = TRUE)$epc

    df_sepc <- sepc_values[sepc_values$lhs %in% items &
                             sepc_values$rhs %in% c(paste0(cov, "lat"), paste0("LFacX", cov)),
                           c("lhs", "op", "rhs", "epc", "sepc.all")]

    colnames(df_sepc) <- c("Item", "Operator", "Effect", "EPC", "SEPC.ALL")
    df_sepc_u <- df_sepc[grepl("lat$", df_sepc$Effect), ]
    df_sepc_nu <- df_sepc[grepl("^LFacX", df_sepc$Effect), ]

    return(list(
      DIF.Global = df_dif_global,
      DIF.Uniforme = df_dif_uniforme,
      DIF.NoUniforme = df_dif_nouniforme,
      SEPC.uDIF = df_sepc_u,
      SEPC.nuDIF = df_sepc_nu
    ))
  }

  resultados_DIF <- mimicout_modificado(fit, mimic_param, cov, Oort.adj, p.crit)
  resultados_DIF$fit <- fit
  return(resultados_DIF)
}
