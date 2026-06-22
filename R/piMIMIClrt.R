#' @title DIF analysis with PI-MIMIC using Likelihood Ratio Tests (LRT)
#'
#' @description
#' Implements the product indicator (PI) approach for MIMIC models to detect
#' uniform and non‑uniform DIF. Uses LRT between unrestricted and restricted
#' models and reports change in R² as effect size. Optionally applies Oort's
#' critical value adjustment to control Type I error inflation.
#'
#' @param data Data frame containing items and the covariate.
#' @param items Character vector of item names.
#' @param cov Name of the covariate (numeric or factor; factors are converted).
#' @param lvname Name of the latent variable (default "LatFact").
#' @param est Estimator for lavaan (default "MLM"; can be "ML", "ULS", etc.).
#' @param anchor Must be one of:
#'        * `"rest"`: iterative "un‑contra‑todos" approach (each item tested against all others as anchors).
#'        * `"none"`: no fixed anchors (all items are tested against a model without effects).
#'        * A character vector with specific item names to use as fixed anchors (at least two recommended).
#' @param Oort.adj Logical; if `TRUE`, applies Oort's critical value adjustment.
#' @param p.crit Numeric; significance level for the Oort adjustment (default 0.05).
#' @param return_models Logical; if `TRUE`, returns the fitted model objects (only for the last iteration if `anchor = "rest"`).
#' @param adjust Character; p-value adjustment method passed to `p.adjust`
#'        (e.g., "bonferroni", "holm", "fdr"). Default "none".
#' @param ... Additional arguments passed to `lavaan::cfa`.
#'
#' @return A list with:
#' \item{DIF.Global}{Data frame: Item, Chi2 (2 df), p.value, and crit.Oort if adjusted.}
#' \item{DIF.Uniforme}{Data frame: Item, Chi2 (1 df), p.value, and crit.Oort if adjusted.}
#' \item{DIF.NoUniforme}{Data frame: Item, Chi2 (1 df), p.value, and crit.Oort if adjusted.}
#' \item{DeltaR2.Global}{Data frame: Item, DeltaR² (full vs no‑DIF model).}
#' \item{DeltaR2.uDIF}{Data frame: Item, DeltaR² (full vs b=0 model).}
#' \item{DeltaR2.nuDIF}{Data frame: Item, DeltaR² (full vs c=0 model).}
#' \item{fit}{The fitted unrestricted (full) lavaan object (only for fixed anchors or `"none"`; `NULL` for `"rest"`).}
#' \item{constrained_fits}{If `return_models=TRUE`, a nested list of fitted models (only for fixed anchors or `"none"`; `NULL` for `"rest"`).}
#'
#' @details
#' When `anchor = "rest"`, the function iterates over each item, treating all other
#' items as anchors. This is the approach described in Oort (1998) and commonly used
#' in MIMIC DIF detection. It does not require pre‑selected anchor items and is
#' suitable when no clear anchor set is available.
#'
#' When `anchor = "none"`, all items are evaluated against a model with no effects
#' (i.e., the baseline model is fully invariant). This is equivalent to a single‑step
#' score test, but using LRT.
#'
#' When a character vector of item names is provided, those items are used as fixed
#' anchors (they never have DIF effects). The remaining items are tested.
#'
#' The Oort adjustment modifies the critical chi-square value:
#' \deqn{K' = (χ²₀ / (K + df₀ - 1)) * K}
#' where χ²₀ and df₀ are from the baseline (full invariance) model, and K is
#' the original critical value. This adjustment is recommended when the baseline
#' model shows evidence of misfit (χ²₀/df₀ > 1), as it helps control Type I error.
#'
#' @references
#' Oort, F. J. (1998). Simulation study of item bias detection with restricted
#' factor analysis. *Structural Equation Modeling, 5*, 107–124.
#'
#' Kolbe, L., & Jorgensen, T. D. (2018). Using product indicators in restricted
#' factor analysis models to detect nonuniform measurement bias.
#' In *Quantitative Psychology* (pp. 235–245). Springer.
#'
#' @importFrom lavaan cfa lavTestLRT lavInspect
#' @importFrom scripty prods
#' @importFrom stats p.adjust qchisq
#' @export
piMIMIClrt <- function(data, items, cov, lvname = "LatFact", est = "MLM",
                       anchor, Oort.adj = FALSE, p.crit = 0.05,
                       return_models = FALSE, adjust = "none", ...) {

  # ---- Verificar que anchor esté especificado ----
  if (missing(anchor)) {
    stop("You must specify 'anchor'. Use 'rest' for iterative approach, 'none' for no anchors, or a vector of item names.")
  }

  # ---- Si anchor == "rest" (iterativo "un-contra-todos") ----
  if (length(anchor) == 1 && anchor == "rest") {
    message("Using 'rest' approach: each item is tested against all other items as anchors.")
    message("Consider using at least two fixed anchors for better identification, if available.")

    n_items <- length(items)
    # Inicializar tablas de resultados
    res_global <- data.frame(Item = items, Chi2 = NA, df = 2, p.value = NA)
    res_uniform <- data.frame(Item = items, Chi2 = NA, df = 1, p.value = NA)
    res_nonuniform <- data.frame(Item = items, Chi2 = NA, df = 1, p.value = NA)
    delta_global <- data.frame(Item = items, DeltaR2 = NA)
    delta_uniform <- data.frame(Item = items, DeltaR2 = NA)
    delta_nonuniform <- data.frame(Item = items, DeltaR2 = NA)

    if (Oort.adj) {
      res_global$crit.Oort <- NA
      res_uniform$crit.Oort <- NA
      res_nonuniform$crit.Oort <- NA
    }

    # Bucle sobre cada ítem
    for (i in seq_along(items)) {
      item_eval <- items[i]
      anchors_rest <- items[-i]  # todos los demás como anclas

      # Llamada recursiva con los anclas específicos
      sub_result <- piMIMIClrt(
        data = data,
        items = items,
        cov = cov,
        lvname = lvname,
        est = est,
        anchor = anchors_rest,
        Oort.adj = Oort.adj,
        p.crit = p.crit,
        return_models = FALSE,
        adjust = "none",
        ...
      )

      # Extraer resultados para el ítem evaluado (única fila en cada tabla)
      res_global[i, c("Chi2", "p.value")] <- sub_result$DIF.Global[1, c("Chi2", "p.value")]
      res_uniform[i, c("Chi2", "p.value")] <- sub_result$DIF.Uniforme[1, c("Chi2", "p.value")]
      res_nonuniform[i, c("Chi2", "p.value")] <- sub_result$DIF.NoUniforme[1, c("Chi2", "p.value")]

      if (Oort.adj) {
        res_global[i, "crit.Oort"] <- sub_result$DIF.Global[1, "crit.Oort"]
        res_uniform[i, "crit.Oort"] <- sub_result$DIF.Uniforme[1, "crit.Oort"]
        res_nonuniform[i, "crit.Oort"] <- sub_result$DIF.NoUniforme[1, "crit.Oort"]
      }

      delta_global[i, "DeltaR2"] <- sub_result$DeltaR2.Global[1, "DeltaR2"]
      delta_uniform[i, "DeltaR2"] <- sub_result$DeltaR2.uDIF[1, "DeltaR2"]
      delta_nonuniform[i, "DeltaR2"] <- sub_result$DeltaR2.nuDIF[1, "DeltaR2"]
    }

    # Ajuste de p-valores si se solicita
    if (adjust != "none") {
      res_global$p.value <- stats::p.adjust(res_global$p.value, method = adjust)
      res_uniform$p.value <- stats::p.adjust(res_uniform$p.value, method = adjust)
      res_nonuniform$p.value <- stats::p.adjust(res_nonuniform$p.value, method = adjust)
    }

    out <- list(
      DIF.Global = res_global,
      DIF.Uniforme = res_uniform,
      DIF.NoUniforme = res_nonuniform,
      DeltaR2.Global = delta_global,
      DeltaR2.uDIF = delta_uniform,
      DeltaR2.nuDIF = delta_nonuniform,
      fit = NULL,
      constrained_fits = NULL
    )
    class(out) <- "piMIMIC"
    return(out)
  }

  # ---- Si anchor == "none" (sin anclas fijas) ----
  if (length(anchor) == 1 && anchor == "none") {
    anchor_items <- character(0)
    message("No anchor items specified. All items will be tested for DIF.")
  } else {
    # ---- Vector de nombres de ítems como anclas fijas ----
    anchor_items <- anchor
    if (!all(anchor_items %in% items)) {
      stop("All anchor items must be in 'items'.")
    }
    if (length(anchor_items) < 2) {
      warning("Using only one anchor item may lead to identification issues.",
              " Consider using at least two anchors.")
    }
    message("Using fixed anchors: ", paste(anchor_items, collapse = ", "))
  }

  tested_items <- setdiff(items, anchor_items)
  if (length(tested_items) == 0) stop("No items left to test after removing anchors.")

  # ---- Preparar covariate centrada ----
  cov_orig <- data[[cov]]
  if (is.factor(cov_orig)) {
    cov_num <- as.numeric(cov_orig) - 1
  } else {
    cov_num <- as.numeric(cov_orig)
  }
  cov_centered <- cov_num - mean(cov_num, na.rm = TRUE)
  data[[paste0(cov, "_cent")]] <- cov_centered

  # ---- Crear productos indicadores (double‑mean‑centering) usando scripty ----
  if (!requireNamespace("scripty", quietly = TRUE)) {
    stop("Package 'scripty' is required for product indicator creation.")
  }
  prod_data <- scripty::prods(
    df = data[, c(items, paste0(cov, "_cent"))],
    covname = paste0(cov, "_cent"),
    match.op = FALSE,
    doubleMC.op = TRUE
  )
  prod_names <- paste0(items, ".", cov, "_cent")
  if (!all(prod_names %in% names(prod_data))) {
    stop("Product indicator columns were not created correctly.")
  }

  # ---- Construir modelo completo ----
  cov_fac <- paste0(cov, "_fac")
  int_fac <- paste0(lvname, "_x_", cov)

  syntax_lv <- paste0(lvname, " =~ ", paste(items, collapse = " + "))
  syntax_cov <- paste0(cov_fac, " =~ 1*", paste0(cov, "_cent"), "\n",
                       paste0(cov, "_cent"), " ~~ 0*", paste0(cov, "_cent"))
  syntax_int <- paste0(int_fac, " =~ ", paste(prod_names, collapse = " + "))

  # Regresiones solo para los ítems no ancla
  reg_lines <- character(length(tested_items))
  for (i in seq_along(tested_items)) {
    item <- tested_items[i]
    b_label <- paste0("b", i)
    c_label <- paste0("c", i)
    reg_lines[i] <- paste0(item, " ~ ", b_label, "*", cov_fac, " + ",
                           c_label, "*", int_fac)
  }
  reg_part <- paste(reg_lines, collapse = "\n")

  covariances <- paste0(lvname, " ~~ ", cov_fac, "\n",
                        lvname, " ~~ ", int_fac, "\n",
                        cov_fac, " ~~ ", int_fac)

  model_full <- paste(syntax_lv, syntax_cov, syntax_int,
                      reg_part, covariances, sep = "\n")

  # ---- Ajustar modelo completo ----
  fit_full <- lavaan::cfa(model = model_full,
                          data = prod_data,
                          estimator = est,
                          ...)
  if (!lavaan::lavInspect(fit_full, "converged")) {
    warning("Full model did not converge. Results may be unreliable.")
  }

  # ---- Extraer R² del modelo completo ----
  rsq_full <- lavaan::lavInspect(fit_full, "rsquare")
  rsq_full_tested <- rsq_full[tested_items]

  # ---- Modelo base para el ajuste Oort (sin efectos) ----
  if (Oort.adj) {
    base_model <- model_full
    for (i in seq_along(tested_items)) {
      base_model <- gsub(paste0("\\+ b", i, "\\*", cov_fac), "", base_model)
      base_model <- gsub(paste0("\\+ c", i, "\\*", int_fac), "", base_model)
    }
    base_model <- gsub("\n\n", "\n", base_model)
    fit_base <- lavaan::cfa(model = base_model, data = prod_data,
                            estimator = est, ...)
    chi0 <- lavaan::lavInspect(fit_base, "fit")["chisq"]
    df0  <- lavaan::lavInspect(fit_base, "fit")["df"]
    K_global <- qchisq(1 - p.crit, 2)
    K_uniform <- qchisq(1 - p.crit, 1)
    crit.global <- (chi0 / (K_global + df0 - 1)) * K_global
    crit.uniform <- (chi0 / (K_uniform + df0 - 1)) * K_uniform
  }

  # ---- Preparar contenedores de resultados ----
  n_items <- length(tested_items)
  results <- list(
    DIF.Global = data.frame(Item = tested_items, Chi2 = NA, df = 2, p.value = NA),
    DIF.Uniforme = data.frame(Item = tested_items, Chi2 = NA, df = 1, p.value = NA),
    DIF.NoUniforme = data.frame(Item = tested_items, Chi2 = NA, df = 1, p.value = NA),
    DeltaR2.Global = data.frame(Item = tested_items, DeltaR2 = NA),
    DeltaR2.uDIF = data.frame(Item = tested_items, DeltaR2 = NA),
    DeltaR2.nuDIF = data.frame(Item = tested_items, DeltaR2 = NA)
  )
  if (Oort.adj) {
    results$DIF.Global$crit.Oort <- NA
    results$DIF.Uniforme$crit.Oort <- NA
    results$DIF.NoUniforme$crit.Oort <- NA
  }
  if (return_models) {
    constrained_fits <- vector("list", n_items)
    names(constrained_fits) <- tested_items
  }

  # ---- Bucle sobre cada ítem no ancla ----
  for (i in seq_along(tested_items)) {
    item <- tested_items[i]
    b_label <- paste0("b", i)
    c_label <- paste0("c", i)

    # Modelo M0: b=0 y c=0
    syntax_m0 <- gsub(paste0(b_label, "\\*", cov_fac), paste0("0*", cov_fac), model_full)
    syntax_m0 <- gsub(paste0(c_label, "\\*", int_fac), paste0("0*", int_fac), syntax_m0)
    fit_m0 <- lavaan::cfa(model = syntax_m0, data = prod_data, estimator = est, ...)

    # Modelo Mb: b=0, c libre
    syntax_mb <- gsub(paste0(b_label, "\\*", cov_fac), paste0("0*", cov_fac), model_full)
    fit_mb <- lavaan::cfa(model = syntax_mb, data = prod_data, estimator = est, ...)

    # Modelo Mc: c=0, b libre
    syntax_mc <- gsub(paste0(c_label, "\\*", int_fac), paste0("0*", int_fac), model_full)
    fit_mc <- lavaan::cfa(model = syntax_mc, data = prod_data, estimator = est, ...)

    # R² de cada modelo
    rsq_m0 <- lavaan::lavInspect(fit_m0, "rsquare")[item]
    rsq_mb <- lavaan::lavInspect(fit_mb, "rsquare")[item]
    rsq_mc <- lavaan::lavInspect(fit_mc, "rsquare")[item]

    # Pruebas LRT
    lrt_method <- if (est == "MLM") "satorra.bentler.2001" else "default"

    # Global (M0 vs full)
    lrt_global <- lavaan::lavTestLRT(fit_m0, fit_full, method = lrt_method)
    chisq_global <- if (nrow(lrt_global) == 2) lrt_global[2, "Chisq diff"] else lrt_global[2, "Chisq"]
    p_global <- if (nrow(lrt_global) == 2) lrt_global[2, "Pr(>Chisq)"] else lrt_global[2, "P"]

    # Uniforme (Mb vs full)
    lrt_ub <- lavaan::lavTestLRT(fit_mb, fit_full, method = lrt_method)
    chisq_ub <- if (nrow(lrt_ub) == 2) lrt_ub[2, "Chisq diff"] else lrt_ub[2, "Chisq"]
    p_ub <- if (nrow(lrt_ub) == 2) lrt_ub[2, "Pr(>Chisq)"] else lrt_ub[2, "P"]

    # No uniforme (Mc vs full)
    lrt_nu <- lavaan::lavTestLRT(fit_mc, fit_full, method = lrt_method)
    chisq_nu <- if (nrow(lrt_nu) == 2) lrt_nu[2, "Chisq diff"] else lrt_nu[2, "Chisq"]
    p_nu <- if (nrow(lrt_nu) == 2) lrt_nu[2, "Pr(>Chisq)"] else lrt_nu[2, "P"]

    # Guardar resultados
    results$DIF.Global[i, c("Chi2", "p.value")] <- c(chisq_global, p_global)
    results$DIF.Uniforme[i, c("Chi2", "p.value")] <- c(chisq_ub, p_ub)
    results$DIF.NoUniforme[i, c("Chi2", "p.value")] <- c(chisq_nu, p_nu)

    if (Oort.adj) {
      results$DIF.Global[i, "crit.Oort"] <- crit.global
      results$DIF.Uniforme[i, "crit.Oort"] <- crit.uniform
      results$DIF.NoUniforme[i, "crit.Oort"] <- crit.uniform
    }

    results$DeltaR2.Global[i, "DeltaR2"] <- rsq_full_tested[i] - rsq_m0
    results$DeltaR2.uDIF[i, "DeltaR2"] <- rsq_full_tested[i] - rsq_mb
    results$DeltaR2.nuDIF[i, "DeltaR2"] <- rsq_full_tested[i] - rsq_mc

    if (return_models) {
      constrained_fits[[i]] <- list(M0 = fit_m0, Mb = fit_mb, Mc = fit_mc)
    }
  }

  # ---- Ajuste de p‑valores si se solicita ----
  if (adjust != "none") {
    results$DIF.Global$p.value <- stats::p.adjust(results$DIF.Global$p.value, method = adjust)
    results$DIF.Uniforme$p.value <- stats::p.adjust(results$DIF.Uniforme$p.value, method = adjust)
    results$DIF.NoUniforme$p.value <- stats::p.adjust(results$DIF.NoUniforme$p.value, method = adjust)
  }

  # ---- Salida ----
  out <- list(
    DIF.Global = results$DIF.Global,
    DIF.Uniforme = results$DIF.Uniforme,
    DIF.NoUniforme = results$DIF.NoUniforme,
    DeltaR2.Global = results$DeltaR2.Global,
    DeltaR2.uDIF = results$DeltaR2.uDIF,
    DeltaR2.nuDIF = results$DeltaR2.nuDIF,
    fit = fit_full
  )
  if (return_models) out$constrained_fits <- constrained_fits

  class(out) <- "piMIMIC"
  return(out)
}

#' @export
print.piMIMIC <- function(x, ...) {
  cat("PI-MIMIC DIF results\n")
  cat("====================\n\n")
  cat("Global DIF (2 df):\n")
  print(x$DIF.Global, row.names = FALSE)
  cat("\nUniform DIF (1 df):\n")
  print(x$DIF.Uniforme, row.names = FALSE)
  cat("\nNon-uniform DIF (1 df):\n")
  print(x$DIF.NoUniforme, row.names = FALSE)
  cat("\nDelta R² (effect sizes):\n")
  cat("Global:\n")
  print(x$DeltaR2.Global, row.names = FALSE)
  cat("Uniform:\n")
  print(x$DeltaR2.uDIF, row.names = FALSE)
  cat("Non-uniform:\n")
  print(x$DeltaR2.nuDIF, row.names = FALSE)
  invisible(x)
}
