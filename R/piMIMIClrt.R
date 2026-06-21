#' @title Function to Analyze DIF with PI-MIMIC using Likelihood Ratio Tests
#'
#' @description
#' This function implements Differential Item Functioning (DIF) analysis using
#' the Product of Indicators (PI) approach within the MIMIC framework.
#' Unlike the original piMIMIC, this version uses Likelihood Ratio Tests (LRT)
#' between unrestricted and restricted models (fixing DIF parameters to zero)
#' and reports the change in R² (delta R²) as an effect size measure.
#'
#' @param data Data frame containing items and the covariate.
#' @param items Character vector of item names.
#' @param cov Name of the covariate (must be numeric or factor; if factor, it
#'   will be converted to numeric).
#' @param lvname Name for the latent variable (default "LatFact").
#' @param est Estimator to use in lavaan (e.g., "MLM", "ML", "ULS"). Must be a
#'   non‑empty character string.
#'
#' @return A list with the following components:
#' \item{DIF.Global}{Data frame with global DIF test (2 df) for each item:
#'   Item, Chi2, df, p.value.}
#' \item{DIF.Uniforme}{Data frame with uniform DIF test (1 df) for each item.}
#' \item{DIF.NoUniforme}{Data frame with non‑uniform DIF test (1 df).}
#' \item{DeltaR2.uDIF}{Data frame with delta R² for uniform DIF (unrestricted R²
#'   minus restricted R²).}
#' \item{DeltaR2.nuDIF}{Data frame with delta R² for non‑uniform DIF.}
#' \item{fit}{The unrestricted lavaan model object.}
#'
#' @details
#' The function constructs a MIMIC model with latent interaction using the
#' double‑mean‑centering product indicator method (scripty::prods).
#' The unrestricted model includes direct effects of the covariate latent
#' variable and the interaction latent variable on each item.
#' For each item, three restricted models are fitted:
#' 1. Global: both direct effects fixed to 0 (2 df test).
#' 2. Uniform: only the covariate effect fixed to 0 (1 df).
#' 3. Non‑uniform: only the interaction effect fixed to 0 (1 df).
#' LRT compares each restricted model against the unrestricted model.
#' Additionally, the R² for each item is extracted from each model, and the
#' difference (unrestricted - restricted) is reported as an effect size.
#'
#' @examples
#' \dontrun{
#' library(psych)
#' data(bfi)
#' data.bfi <- bfi[, c("N1","N2","N3","N4","N5","gender")]
#' data.bfi <- data.bfi[complete.cases(data.bfi), ]
#' data.bfi$gender <- as.factor(data.bfi$gender)
#' neuro.items <- c("N1","N2","N3","N4","N5")
#' res <- piMIMIClrt(data = data.bfi, items = neuro.items,
#'                    cov = "gender", lvname = "Neuroticism", est = "MLM")
#' res$DIF.Global
#' res$DeltaR2.uDIF
#' }
#'
#' @importFrom lavaan cfa lavTestLRT lavInspect
#' @importFrom scripty prods
#' @export
piMIMIClrt <- function(data, items, cov, lvname = "LatFact", est = "MLM") {

  # ---- Argument checks ----
  if (!is.character(est) || nchar(est) == 0) {
    stop("Error: Estimator must be a non-empty string (see lavaan documentation).")
  }
  if (!all(c(items, cov) %in% colnames(data))) {
    stop("Error: Some item names or the covariate are not present in the data frame.")
  }
  if (!is.character(lvname) || lvname == "") {
    stop("Error: The latent variable name ('lvname') must be a non-empty string.")
  }

  # Convert factor covariate to numeric
  if (is.factor(data[[cov]])) {
    data[[cov]] <- as.numeric(data[[cov]])
  } else if (!is.numeric(data[[cov]])) {
    stop(paste("Error: The covariate", cov, "must be either a factor or numeric."))
  }

  # ---- Prepare product indicators (double mean centering) ----
  df_pi <- scripty::prods(df = data[, c(items, cov)],
                          covname = cov,
                          match.op = FALSE,
                          doubleMC.op = TRUE)

  # Names of product indicators (assumed to be "item.cov")
  prod_names <- paste0(items, ".", cov)
  cov_lat <- paste0(cov, "lat")          # latent covariate factor name
  int_name <- paste0("LFacX", cov)       # latent interaction factor name

  # ---- Build unrestricted model syntax ----
  # Base part: latent factors, residual covariances, and fix residual variance of cov indicator
  base_syntax <- paste0(
    lvname, " =~ ", paste(items, collapse = " + "), "\n",
    cov_lat, " =~ ", cov, "\n",
    int_name, " =~ ", paste(prod_names, collapse = " + "), "\n",
    paste(paste0(items, " ~~ ", prod_names), collapse = "\n"), "\n",
    paste0(cov, " ~~ 0*", cov)   # fix residual variance of covariate indicator
  )

  # Regressions of each item on cov_lat and int_name (direct DIF effects)
  reg_lines <- character()
  for (it in items) {
    reg_lines <- c(reg_lines,
                   paste0(it, " ~ ", cov_lat),          # uniform DIF
                   paste0(it, " ~ ", int_name))        # non-uniform DIF
  }
  unrestricted_syntax <- paste0(base_syntax, "\n", paste(reg_lines, collapse = "\n"))

  # ---- Fit unrestricted model ----
  fit_un <- lavaan::cfa(unrestricted_syntax,
                        data = df_pi,
                        estimator = est,
                        meanstructure = TRUE)

  # R² for all items in unrestricted model
  r2_un <- lavaan::lavInspect(fit_un, "rsquare")[items]

  # ---- Initialize output data frames ----
  n_items <- length(items)
  out_global <- data.frame(Item = items, Chi2 = NA_real_, df = NA_integer_, p.value = NA_real_,
                           stringsAsFactors = FALSE)
  out_uniform <- data.frame(Item = items, Chi2 = NA_real_, df = NA_integer_, p.value = NA_real_,
                            stringsAsFactors = FALSE)
  out_nouniform <- data.frame(Item = items, Chi2 = NA_real_, df = NA_integer_, p.value = NA_real_,
                              stringsAsFactors = FALSE)
  out_delta_u <- data.frame(Item = items, delta.R2 = NA_real_, stringsAsFactors = FALSE)
  out_delta_nu <- data.frame(Item = items, delta.R2 = NA_real_, stringsAsFactors = FALSE)

  # For each item, fit restricted models and perform LRT
  for (i in seq_along(items)) {
    it <- items[i]
    # positions of the two regression lines for this item in reg_lines
    pos_uniform <- 2 * (i - 1) + 1
    pos_nouniform <- 2 * (i - 1) + 2

    # ---- 1. Global restriction (remove both lines) ----
    kept_global <- reg_lines[-c(pos_uniform, pos_nouniform)]
    syntax_global <- paste0(base_syntax, "\n", paste(kept_global, collapse = "\n"))
    fit_global <- lavaan::cfa(syntax_global,
                              data = df_pi,
                              estimator = est,
                              meanstructure = TRUE)
    # LRT global (2 df)
    lrt_global <- lavaan::lavTestLRT(fit_global, fit_un, method = "default")
    # Extract difference test statistics (second row)
    out_global$Chi2[i] <- lrt_global[2, "Chisq diff"]
    out_global$df[i]   <- lrt_global[2, "Df diff"]
    out_global$p.value[i] <- lrt_global[2, "Pr(>Chisq)"]

    # R² for this item in global restricted model
    r2_global_it <- lavaan::lavInspect(fit_global, "rsquare")[it]
    # delta R² for global? We don't need to report it separately; but we compute
    # for uniform and non-uniform below.

    # ---- 2. Uniform restriction (remove cov_lat line only) ----
    kept_uniform <- reg_lines[-pos_uniform]
    syntax_uniform <- paste0(base_syntax, "\n", paste(kept_uniform, collapse = "\n"))
    fit_uniform <- lavaan::cfa(syntax_uniform,
                               data = df_pi,
                               estimator = est,
                               meanstructure = TRUE)
    lrt_uniform <- lavaan::lavTestLRT(fit_uniform, fit_un, method = "default")
    out_uniform$Chi2[i] <- lrt_uniform[2, "Chisq diff"]
    out_uniform$df[i]   <- lrt_uniform[2, "Df diff"]
    out_uniform$p.value[i] <- lrt_uniform[2, "Pr(>Chisq)"]

    r2_uniform_it <- lavaan::lavInspect(fit_uniform, "rsquare")[it]
    out_delta_u$delta.R2[i] <- r2_un[it] - r2_uniform_it

    # ---- 3. Non-uniform restriction (remove int_name line only) ----
    kept_nouniform <- reg_lines[-pos_nouniform]
    syntax_nouniform <- paste0(base_syntax, "\n", paste(kept_nouniform, collapse = "\n"))
    fit_nouniform <- lavaan::cfa(syntax_nouniform,
                                 data = df_pi,
                                 estimator = est,
                                 meanstructure = TRUE)
    lrt_nouniform <- lavaan::lavTestLRT(fit_nouniform, fit_un, method = "default")
    out_nouniform$Chi2[i] <- lrt_nouniform[2, "Chisq diff"]
    out_nouniform$df[i]   <- lrt_nouniform[2, "Df diff"]
    out_nouniform$p.value[i] <- lrt_nouniform[2, "Pr(>Chisq)"]

    r2_nouniform_it <- lavaan::lavInspect(fit_nouniform, "rsquare")[it]
    out_delta_nu$delta.R2[i] <- r2_un[it] - r2_nouniform_it
  }

  # ---- Return results ----
  list(
    DIF.Global = out_global,
    DIF.Uniforme = out_uniform,
    DIF.NoUniforme = out_nouniform,
    DeltaR2.uDIF = out_delta_u,
    DeltaR2.nuDIF = out_delta_nu,
    fit = fit_un
  )
}
