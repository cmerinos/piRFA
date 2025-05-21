#' @title Function to Analyze DIF with PI-MIMIC (Product of Indicators)
#'
#' @description
#' The `piRFA` function analyzes Differential Item Functioning (DIF)
#' using the multiple-indicators multiple-causes (MIMIC) framework with product of indicators (PI).
#' It relies on `lavaan` and `scripty`. Uniform and non-uniform DIF can be evaluated in a measurement scale,
#' through statistical tests and an effect size approximation.
#'
#' @param data DataFrame containing items and the covariate.
#' @param items Vector of item names within `data`.
#' @param cov Name of the covariate in `data` (can be categorical or numeric). If categorical, it must be a factor.
#' @param lvname Name for the latent variable in the model (default is `"LatFact"`).
#' @param est Abbreviation of the estimator to use; options are: "MLM", "MLR", "ULSMV" (default is `"MLM"`).
#'
#' @return
#' The function returns a list with the following DataFrames:
#'  \itemize{
#'    \item \code{DIF_Global} - Global DIF results.
#'    \item \code{DIF_Uniforme} - Uniform DIF results.
#'    \item \code{DIF_NoUniforme} - Non-uniform DIF results.
#'    \item \code{SEPC} - SEPC coefficients for uniform and non-uniform DIF.
#'  }
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
#' on the likelihood ratio test (LRT), comparing an unrestricted model
#' (with covariate effect and interaction parameters freely estimated) and a restricted model
#' (with these parameters fixed to zero), using a chi-square test with 2 degrees of freedom (`df = 2`).
#'
#' The implementation relies on `lavaan::lavTestScore` and `lavaan::cfa`, integrating modified functionalities
#' from `scripty::mimic` and `scripty::prods`. Unlike `scripty`, this function
#' avoids creating external files, allowing all output to be accessed in the R session.
#'
#' As an effect size indicator, `piRFA()` uses the releasing fixed-to-zero
#' parameters approach (Whittaker, 2012; Garnier-Villarreal & Jorgensen, 2024). The expected parameter change
#' (`epc`) and the standardized expected parameter change (`sepc.all`) are estimated following Chou & Bentler (1993).
#' `sepc.all` represents the standardized unit change upon releasing the zero restriction
#' in the analyzed item parameters. Values `sepc.all ≥ .20` are considered indicative of
#' possible model misspecification (Whittaker, 2012), while values `≥ .10`
#' (Kaplan, 1989) are a threshold for a moderate effect.
#'
#' The user can select the appropriate estimator among `"MLM"`, `"MLR"` and `"ULSMV"`,
#' depending on normality and robustness assumptions.
#'
#' @examples
#' ### Example 1 -------------
#'
#' set.seed(123)
#' Exmp1.data <- data.frame(
#'   grp = sample(0:1, 100, replace = TRUE),  # Group variable
#'   item1 = sample(1:5, 100, replace = TRUE),
#'   item2 = sample(1:5, 100, replace = TRUE),
#'   item3 = sample(1:5, 100, replace = TRUE),
#'   item4 = sample(1:5, 100, replace = TRUE))
#'
#' # Run DIF analysis, and full output
#' piRFA(data = Exmp1.data , items = c("item1", "item2", "item3"), cov = "grp")
#'
#' # Specific output: Uniform DIF (uDIF)
#' Exmp1.output <- piRFA(data = Exmp1.data , items = c("item1", "item2", "item3"), cov = "grp")
#'
#' Exmp1.output$DIF_Uniforme
#'
#' # Specific output: full SEPC
#' piRFA(data = Exmp1.data , items = c("item1", "item2", "item3"), cov = "grp")$SEPC
#'
#' # Specific output: full SEPC for uniform DIF
#' piRFA(data = Exmp1.data , items = c("item1", "item2", "item3"), cov = "grp")$SEPC[c(1,3), ]
#'
#'
#'
#' @references
#' Chou, C.-P., & Bentler, P. (1993). Invariant standardized estimated parameter change for model modification in covariance structure analysis. *Multivariate Behavioral Research*, 28, 97–110. https://doi.org/10.1207/s15327906mbr2801_6
#'
#' Finch, H. (2005). The MIMIC model as a method for detecting DIF: Comparison with Mantel–Haenszel, SIBTEST, and the IRT likelihood ratio. *Applied Psychological Measurement*, 29, 278–295. https://doi.org/10.1177/0146621605275728
#'
#' Garnier-Villarreal, M., & Jorgensen, T. D. (2024). Evaluating Local Model Misspecification with Modification Indices in Bayesian Structural Equation Modeling. *Structural Equation Modeling: A Multidisciplinary Journal*, 1–15. https://doi.org/10.1080/10705511.2024.2413128
#'
#' Kaplan D. (1989). Model Modification in Covariance Structure Analysis: Application of the Expected Parameter Change Statistic. *Multivariate Behavioral Research*, 24(3), 285–305. https://doi.org/10.1207/s15327906mbr2403_2
#'
#' Kolbe, L., & Jorgensen, T. D. (2018). Using product indicators in restricted factor analysis models to detect nonuniform measurement bias. In M. Wiberg, S. A. Culpepper, R. Janssen, J. González, & D. Molenaar (Eds.), Quantitative psychology: The 82nd Annual Meeting of the Psychometric Society, Zurich, Switzerland, 2017 (pp. 235–245). New York, NY: Springer. https://doi.org/10.1007/978-3-319-77249-3_20
#'
#' Kolbe, L., & Jorgensen, T. D. (2019). Using restricted factor analysis to select anchor items and detect differential item functioning. *Behavior Research Methods*, 51, 138–151. https://doi.org/10.3758/s13428-018-1151-3
#'
#' Kolbe, L., Jorgensen, T. D., & Molenaar, D. (2020). The Impact of Unmodeled Heteroskedasticity on Assessing Measurement Invariance in Single-group Models. *Structural Equation Modeling: A Multidisciplinary Journal*, 28(1), 82–98. https://doi.org/10.1080/10705511.2020.1766357
#'
#' Whittaker, T. A. (2012). Estimation of Standardized Expected Parameter Change for DIF Detection. *Educational and Psychological Measurement*, 72(3), 342-357.
#'
#' @seealso
#' \code{\link[lavaan]{cfa}}, \code{\link[scripty]{mimic}}, \code{\link[scripty]{prods}}, \code{\link{piRFA.plot}}
#'
#' @importFrom lavaan cfa lavTestScore parameterestimates
#' @importFrom scripty prods mimicparam
#' @export
piRFA <- function(data, items, cov, lvname = "LatFact", est = "MLM") {

  # Check if estimator is valid
  if (!est %in% c("MLM", "MLR", "ULSMV")) {
    stop("Error: Estimator must be 'MLM', 'MLR' or 'ULSMV'.")
  }

  mimicout_modificado <- function(fit.mimic, mimic.param, cov) {
    ests <- as.data.frame(lavaan::parameterestimates(fit.mimic))
    uniqnames <- unique(ests$lhs)
    lvname <- uniqnames[1]
    covlvname <- uniqnames[2]

    any.out <- do.call(rbind, lapply(mimic.param, function(x)
      lavaan::lavTestScore(fit.mimic, add = x)$test))

    sep.out <- lavaan::lavTestScore(fit.mimic, add = as.character(mimic.param))

    # Create `$DIF_Global` with "any.chi2", "any.df", "any.p"
    df_dif_global <- data.frame(
      DIF_Global_Chi2 = round(any.out$X2, 2),
      df = any.out$df,
      p_value = round(any.out$p, 3)
    )

    # Keep `$DIF_Uniforme` and `$DIF_NoUniforme`
    oddnum <- seq(1, length(sep.out$uni$lhs), 2)
    evennum <- seq(2, length(sep.out$uni$lhs), 2)

    df_dif_uniforme <- data.frame(
      Item = sep.out$uni$lhs[oddnum],
      uDIF_Chi2 = round(sep.out$uni$X2[oddnum], 2),
      df = sep.out$uni$df[oddnum],
      p_value = sep.out$uni$p.value[oddnum]
    )

    df_dif_nouniforme <- data.frame(
      Item = sep.out$uni$lhs[evennum],
      nuDIF_Chi2 = round(sep.out$uni$X2[evennum], 2),
      df = sep.out$uni$df[evennum],
      p_value = sep.out$uni$p.value[evennum]
    )

    # Extract "EPC" and "SEPC.ALL" for `$SEPC`
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

    return(list(
      DIF_Global = df_dif_global,
      DIF_Uniforme = df_dif_uniforme,
      DIF_NoUniforme = df_dif_nouniforme,
      SEPC = df_sepc
    ))
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

  fit <- cfa(model_mimic,
             data = df_pi,
             estimator = est,
             meanstructure = TRUE)

  mimic_param <- scripty::mimicparam(fit)

  resultados_DIF <- mimicout_modificado(fit, mimic_param, cov)
  resultados_DIF$fit <- fit  # <-- Añade el objeto lavaan
  return(resultados_DIF)
}
