#' Profile Plot of DIF for PI-MIMIC Output (Product of Indicators within a Multiple-Indicators Multiple-Causes Framework)
#'
#' This function plots the expected response curves (linear prediction) across levels of the latent trait (θ)
#' for each group defined by a covariate. It is useful to visually inspect DIF—especially non-uniform DIF.
#'
#' @param resultados Output object from the `piRFA()` function.
#' @param data The dataset used in the `piRFA()` analysis.
#' @param item The name of the item to be plotted (as a character string).
#' @param cov The name of the grouping covariate used in `piRFA()` (as a character string).
#' @param thetaRange A numeric vector of length 2 specifying the latent trait range (default is -3 to 3).
#' @param nPoints Number of evaluation points along the latent trait (default = 200).
#' @param themeOption  A `ggplot2` theme object to customize the plot appearance.
#' Common options include:
#' `theme_minimal()` (default), `theme_bw()`, `theme_classic()`,
#' `theme_light()`, `theme_dark()`, `theme_void()`, `theme_gray()`, `theme_linedraw()`,
#' `theme_test()` and others.
#' @param labels Optional vector of labels for the groups. If NULL, the original group names are used.
#' @param parType Type of parameters to be used when computing the expected curves:
#' \itemize{
#'   \item `"none"`: Ignores parameter estimates from `lavaan` and fixes `Intercept = 0`, `Loading = 1`.
#'   \item `"unstd"`: Uses unstandardized estimates (`est`) from the fitted `lavaan` model (original scale of the item).
#'   \item `"std"`: Uses standardized estimates (`std.all`) from the fitted `lavaan` model (standardized metric).
#' }
#'
#' @details
#' The function uses the parameter estimates from `piRFA()` to compute predicted responses
#' for each group along a latent variable continuum. It is especially helpful in visualizing
#' non-uniform DIF, where the slope or pattern of response varies by group.
#'
#' The formula used is a linear prediction:
#' \deqn{y = intercept + loading * theta + direct_effect * group + interaction * (theta × group)}
#'
#' The choice of `parType` determines how the parameters are applied:
#' \itemize{
#'   \item `"none"`: Parameters are fixed at `Intercept = 0`, `Loading = 1`, so the plot only
#'   reflects the *direct* and *interaction* effects. This isolates the DIF component.
#'   \item `"unstd"`: Uses unstandardized estimates (`est`) from the fitted `lavaan` model,
#'   corresponding to the original observed scale of the item.
#'   \item `"std"`: Uses standardized estimates (`std.all`) from the fitted `lavaan` model,
#'   making the prediction curves directly comparable across items in a standardized metric.
#' }
#'
#' This function assumes a linear model and is best interpreted when DIF is suspected or detected.
#'
#' @return A `ggplot2` object showing the DIF profile for PI-MIMIC output
#' (Multiple-Indicators Multiple-Causes with Product of Indicators).
#'
#' @examples
#' \dontrun{
#' # Run piRFA
#' results <- piRFA(data = mydata, items = c("item1", "item2"), cov = "group")
#'
#' # DIF profile for item1
#' piRFA.profile(resultados = results, data = mydata, item = "item1", cov = "group")
#'
#' # DIF profile with standardized parameters
#' piRFA.profile(resultados = results, data = mydata, item = "item1", cov = "group",
#'                parType = "std")
#' }
#'
#' @export
piRFA.profile <- function(resultados, data, item, cov,
                          thetaRange = c(-3, 3),
                          nPoints = 200,
                          themeOption = ggplot2::theme_minimal(),
                          labels = NULL,
                          parType = c("none", "unstd", "std")) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  parType <- match.arg(parType)

  groups <- unique(data[[cov]])
  if (length(groups) < 2) stop("`cov` debe tener al menos 2 niveles.")

  # -------- DIF effects (uDIF / nuDIF) --------
  if (is.null(resultados$SEPC.uDIF) && is.null(resultados$SEPC.nuDIF)) {
    stop("No se encontraron `SEPC.uDIF` ni `SEPC.nuDIF` en `resultados`.")
  }
  sepc_data <- rbind(resultados$SEPC.uDIF, resultados$SEPC.nuDIF)
  coefs_item <- sepc_data[sepc_data$Item == item, , drop = FALSE]

  # Elegir columna para los efectos DIF según parType
  sepc_col <- if (parType == "std" && "SEPC.ALL" %in% names(coefs_item)) "SEPC.ALL" else "EPC"
  if (parType == "std" && sepc_col != "SEPC.ALL") {
    warning("SEPC.ALL no disponible; usando EPC (no estandarizado) para efectos DIF.")
  }

  direct      <- coefs_item[[sepc_col]][coefs_item$Effect == paste0(cov, "lat")]
  interaction <- coefs_item[[sepc_col]][grepl(paste0("LFacX", cov), coefs_item$Effect)]
  if (length(direct) == 0)      direct <- 0
  if (length(interaction) == 0) interaction <- 0

  # -------- Intercept & Loading --------
  intercept <- 0
  loading   <- 1

  if (parType != "none" && !is.null(resultados$fit) && requireNamespace("lavaan", quietly = TRUE)) {
    pe <- try(lavaan::parameterEstimates(resultados$fit, standardized = TRUE), silent = TRUE)
    if (!inherits(pe, "try-error")) {

      get1 <- function(x) if (length(x) > 0) x[1] else NA_real_

      # loading
      lrow <- pe[pe$op == "=~" & pe$rhs == item, , drop = FALSE]
      if (nrow(lrow) > 0) {
        if (parType == "std") {
          val_std <- get1(lrow$std.all)
          if (!is.na(val_std)) {
            loading <- val_std
          } else {
            val_un  <- get1(lrow$est)
            if (!is.na(val_un)) {
              loading <- val_un
              warning("Loading std.all es NA; se usó 'est' (no estandarizado).")
            }
          }
        } else { # unstd
          val_un <- get1(lrow$est)
          if (!is.na(val_un)) loading <- val_un
        }
      }

      # intercept
      irow <- pe[pe$op == "~1" & pe$lhs == item, , drop = FALSE]
      if (nrow(irow) > 0) {
        if (parType == "std") {
          # Por convención, intercepto estandarizado = 0 si std.all es NA
          val_std <- get1(irow$std.all)
          if (!is.na(val_std)) {
            intercept <- val_std
          } else {
            intercept <- 0
            message("Intercepto std.all no definido; se fijó en 0 (métrica estandarizada).")
          }
        } else { # unstd
          val_un <- get1(irow$est)
          if (!is.na(val_un)) intercept <- val_un
        }
      }
    }
  }

  # -------- Debug --------
  message("Item: ", item,
          "\n  Intercept   = ", round(intercept, 4),
          "\n  Loading     = ", round(loading, 4),
          "\n  Direct      = ", round(direct, 4),
          "\n  Interaction = ", round(interaction, 4),
          "\n  parType     = ", parType)

  # -------- Predicciones --------
  theta_seq <- seq(thetaRange[1], thetaRange[2], length.out = nPoints)
  pred_df <- do.call(rbind, lapply(groups, function(g) {
    group_indicator <- ifelse(g == groups[1], 0, 1)
    y_pred <- intercept + loading * theta_seq +
      direct * group_indicator +
      interaction * theta_seq * group_indicator
    data.frame(theta = theta_seq, predicted = y_pred, group = as.character(g))
  }))

  # Etiquetas
  if (!is.null(labels)) {
    if (length(labels) != length(groups)) stop("Length de `labels` debe igualar # de grupos.")
    pred_df$group <- factor(pred_df$group, levels = groups, labels = labels)
  } else {
    pred_df$group <- factor(pred_df$group, levels = groups)
  }

  # -------- Plot --------
  p <- ggplot2::ggplot(pred_df, ggplot2::aes(x = theta, y = predicted, color = group)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = paste("DIF Profile -", item),
      x = "Latent Trait (θ)",
      y = "Expected Response",
      color = "Group"
    ) +
    themeOption

  return(p)
}
