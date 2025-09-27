#' @title Function to Analyze DIF with PI-RFA/MIMIC (PI: Product of Indicators)
#'
#' @description
#' The `piRFA` function analyzes Differential Item Functioning (DIF)
#' using the Restricted Factor Analysis (RFA) framework with product of indicators (PI), and
#' modeled like multiple-indicators multiple-causes (MIMIC) .
#' It relies on `lavaan` and `scripty`. Uniform and non-uniform DIF can be evaluated in a measurement scale,
#' through statistical tests and an effect size approximation.
#'
#' @param data DataFrame containing items and the covariate.
#' @param items Vector of item names within `data`.
#' @param cov Name of the covariate in `data` (can be categorical or numeric). If categorical, it must be a factor.
#' @param lvname Name for the latent variable in the model (default is `"LatFact"`).
#' @param est Abbreviation of the estimator to use. Must be a non-empty string (see lavaan documentation).
#'
#' @return
#' The function returns a list with the following DataFrames:
#'  \itemize{
#'    \item \code{DIF.Global} - Global DIF results.
#'    \item \code{DIF.Uniforme} - Uniform DIF results.
#'    \item \code{DIF.NoUniforme} - Non-uniform DIF results.
#'    \item \code{SEPC.uDIF} - SEPC coefficients for uniform DIF.
#'    \item \code{SEPC.nuDIF} - SEPC coefficients for non-uniform DIF.
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
#' The user can select the appropriate estimator depending on normality and robustness assumptions.
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
#' Exmp1.output$DIF.Uniforme
#'
#' # Specific output: SEPC for uniform DIF
#' Exmp1.output$SEPC.uDIF
#'
#' # Specific output: SEPC for non-uniform DIF
#' Exmp1.output$SEPC.nuDIF
#'
#' ### Example 2: Using the 'bfi' dataset from the 'psych' package
#' library(psych)
#' data("bfi")
#'
#' # Select Neuroticism items + gender as covariate
#' data.bfi <- bfi[, c("N1","N2","N3","N4","N5","gender")]
#' data.bfi <- data.bfi[complete.cases(data.bfi), ]
#' data.bfi$gender <- as.factor(data.bfi$gender)
#'
#' neuro.items <- c("N1","N2","N3","N4","N5")
#'
#' # Run DIF analysis
#' res.bfi <- piRFA(data = data.bfi,
#'                  items = neuro.items,
#'                  cov = "gender",
#'                  lvname = "Neuroticism",
#'                  est = "MLM")
#'
#' # Global DIF results
#' res.bfi$DIF.Global
#'
#' # SEPC for uniform DIF
#' res.bfi$SEPC.uDIF
#'
#' @references
#' (mismos que ya tenías…)
#'
#' @seealso
#' \code{\link[lavaan]{cfa}}, \code{\link[scripty]{mimic}}, \code{\link[scripty]{prods}}, \code{\link{piRFA.plot}}
#'
#' @importFrom lavaan cfa lavTestScore parameterestimates
#' @importFrom scripty prods mimicparam
#' @export
piRFA <- function(data, items, cov, lvname = "LatFact", est = "MLM") {

  # Check estimator: allow any non-empty string
  if (!is.character(est) || nchar(est) == 0) {
    stop("Error: Estimator must be a non-empty string (see lavaan documentation).")
  }

  mimicout_modificado <- function(fit.mimic, mimic.param, cov) {
    ests <- as.data.frame(lavaan::parameterestimates(fit.mimic))
    uniqnames <- unique(ests$lhs)
    lvname <- uniqnames[1]

    any.out <- do.call(rbind, lapply(mimic.param, function(x)
      lavaan::lavTestScore(fit.mimic, add = x)$test))

    sep.out <- lavaan::lavTestScore(fit.mimic, add = as.character(mimic.param))

    df_dif_global <- data.frame(
      DIF.Global.Chi2 = round(any.out$X2, 2),
      df = any.out$df,
      p.value = round(any.out$p, 3)
    )

    oddnum <- seq(1, length(sep.out$uni$lhs), 2)
    df_dif_uniforme <- data.frame(
      Item = sep.out$uni$lhs[oddnum],
      uDIF.Chi2 = round(sep.out$uni$X2[oddnum], 2),
      df = sep.out$uni$df[oddnum],
      p.value = sep.out$uni$p.value[oddnum]
    )

    evennum <- seq(2, length(sep.out$uni$lhs), 2)
    df_dif_nouniforme <- data.frame(
      Item = sep.out$uni$lhs[evennum],
      nuDIF.Chi2 = round(sep.out$uni$X2[evennum], 2),
      df = sep.out$uni$df[evennum],
      p.value = sep.out$uni$p.value[evennum]
    )

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
  resultados_DIF$fit <- fit
  return(resultados_DIF)
}
