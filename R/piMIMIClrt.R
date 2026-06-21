#' @title DIF analysis with PI-MIMIC using Likelihood Ratio Tests (LRT)
#'
#' @description
#' Implements the product indicator (PI) approach for MIMIC models to detect
#' uniform and non-uniform DIF. Uses LRT between unrestricted and restricted
#' models (fixing DIF parameters to zero) and reports change in R² as effect size.
#' Follows the syntax and logic of Kolbe & Jorgensen (2018), Table 2.
#'
#' @param data Data frame containing items and the covariate.
#' @param items Character vector of item names.
#' @param cov Name of the covariate (numeric or factor; factors are converted).
#' @param lvname Name of the latent variable (default "LatFact").
#' @param est Estimator for lavaan (default "MLM"; can be "ML", "ULS", etc.).
#'
#' @return A list with:
#' \item{DIF.Global}{Global test (2 df) for each item: Chi2, df, p.value.}
#' \item{DIF.Uniforme}{Uniform DIF test (1 df).}
#' \item{DIF.NoUniforme}{Non-uniform DIF test (1 df).}
#' \item{DeltaR2.uDIF}{Change in R² (unrestricted - restricted) for uniform DIF.}
#' \item{DeltaR2.nuDIF}{Change in R² for non-uniform DIF.}
#' \item{fit}{Unrestricted lavaan model object.}
#'
#' @examples
#' \dontrun{
#' library(psych)
#' data(bfi)
#' data.bfi <- bfi[, c("N1","N2","N3","N4","N5","gender")]
#' data.bfi <- data.bfi[complete.cases(data.bfi), ]
#' data.bfi$gender <- as.factor(data.bfi$gender)
#' neuro.items <- c("N1","N2","N3","N4","N5")
#' res <- piMIMIC_lrt(data = data.bfi, items = neuro.items,
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
    stop("Error: Estimator must be a non-empty string.")
  }
  if (!all(c(items, cov) %in% colnames(data))) {
    stop("Error: Some item names or the covariate are not present in the data frame.")
  }
  if (!is.character(lvname) || lvname == "") {
    stop("Error: The latent variable name ('lvname') must be a non-empty string.")
  }

  # Convert factor covariate to numeric (as in piMIMIC)
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

  # Names of product indicators (item.cov)
  prod_names <- paste0(items, ".", cov)
  cov_lat <- paste0(cov, "lat")          # latent covariate
  int_name <- paste0("LFacX", cov)       # latent interaction

  # ---- Build unrestricted model syntax ----
  # Base part: latent factors, residual covariances, fix covariate indicator variance
  base_syntax <- paste0(
    lvname, " =~ ", paste(items, collapse = " + "), "\n",
    cov_lat, " =~ ", cov, "\n",
    int_name, " =~ ", paste(prod_names, collapse = " + "), "\n",
    paste(paste0(items, " ~~ ", prod_names), collapse = "\n"), "\n",
    paste0(cov, " ~~ 0*", cov)
  )

  # Regressions: for each item, both cov_lat and int_name predict the item
  # Using the compact notation "cov_lat + int_name =~ item"
  reg_lines <- paste0(cov_lat, " + ", int_name, " =~ ", items)
  unrestricted_syntax <- paste0(base_syntax, "\n", paste(reg_lines, collapse = "\n"))

  # ---- Fit unrestricted model ----
  fit_un <- tryCatch(
    lavaan::cfa(unrestricted_syntax,
                data = df_pi,
                estimator = est,
                meanstructure = TRUE),
    error = function(e) {
      # If fails, try with se = "none" to avoid inversion issues
      lavaan::cfa(unrestricted_syntax,
                  data = df_pi,
                  estimator = est,
                  meanstructure = TRUE,
                  se = "none")
    }
  )

  # R² for items in unrestricted model
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

  # ---- For each item, build and compare restricted models ----
  for (i in seq_along(items)) {
    it <- items[i]

    # --- Global restriction: remove the entire regression line for this item ---
    kept_global <- reg_lines[-i]
    syntax_global <- paste0(base_syntax, "\n", paste(kept_global, collapse = "\n"))
    fit_global <- tryCatch(
      lavaan::cfa(syntax_global,
                  data = df_pi,
                  estimator = est,
                  meanstructure = TRUE),
      error = function(e) {
        lavaan::cfa(syntax_global,
                    data = df_pi,
                    estimator = est,
                    meanstructure = TRUE,
                    se = "none")
      }
    )
    # LRT global (2 df)
    lrt_global <- tryCatch(
      lavaan::lavTestLRT(fit_global, fit_un, method = "default"),
      error = function(e) {
        # fallback to anova if lavTestLRT fails
        stats::anova(fit_global, fit_un)
      }
    )
    out_global$Chi2[i] <- lrt_global[2, "Chisq diff"]
    out_global$df[i]   <- lrt_global[2, "Df diff"]
    out_global$p.value[i] <- lrt_global[2, "Pr(>Chisq)"]

    # R² for this item in global restricted model (only for this item)
    r2_global_it <- lavaan::lavInspect(fit_global, "rsquare")[it]

    # --- Uniform restriction: remove cov_lat effect, keep int_name ---
    # Replace the line with only "int_name =~ item" (i.e., item ~ int_name)
    reg_uniform <- reg_lines
    reg_uniform[i] <- paste0(int_name, " =~ ", it)
    syntax_uniform <- paste0(base_syntax, "\n", paste(reg_uniform, collapse = "\n"))
    fit_uniform <- tryCatch(
      lavaan::cfa(syntax_uniform,
                  data = df_pi,
                  estimator = est,
                  meanstructure = TRUE),
      error = function(e) {
        lavaan::cfa(syntax_uniform,
                    data = df_pi,
                    estimator = est,
                    meanstructure = TRUE,
                    se = "none")
      }
    )
    lrt_uniform <- tryCatch(
      lavaan::lavTestLRT(fit_uniform, fit_un, method = "default"),
      error = function(e) stats::anova(fit_uniform, fit_un)
    )
    out_uniform$Chi2[i] <- lrt_uniform[2, "Chisq diff"]
    out_uniform$df[i]   <- lrt_uniform[2, "Df diff"]
    out_uniform$p.value[i] <- lrt_uniform[2, "Pr(>Chisq)"]

    r2_uniform_it <- lavaan::lavInspect(fit_uniform, "rsquare")[it]
    out_delta_u$delta.R2[i] <- r2_un[it] - r2_uniform_it

    # --- Non-uniform restriction: remove int_name effect, keep cov_lat ---
    reg_nouniform <- reg_lines
    reg_nouniform[i] <- paste0(cov_lat, " =~ ", it)
    syntax_nouniform <- paste0(base_syntax, "\n", paste(reg_nouniform, collapse = "\n"))
    fit_nouniform <- tryCatch(
      lavaan::cfa(syntax_nouniform,
                  data = df_pi,
                  estimator = est,
                  meanstructure = TRUE),
      error = function(e) {
        lavaan::cfa(syntax_nouniform,
                    data = df_pi,
                    estimator = est,
                    meanstructure = TRUE,
                    se = "none")
      }
    )
    lrt_nouniform <- tryCatch(
      lavaan::lavTestLRT(fit_nouniform, fit_un, method = "default"),
      error = function(e) stats::anova(fit_nouniform, fit_un)
    )
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
