#' Profile Plot of DIF for PI-MIMIC Output (Product of Indicators within a Multiple-Indicators Multiple-Causes Framework)
#'
#' This function plots the expected response curves (linear prediction) across levels of the latent trait (θ)
#' for each group defined by a covariate. It is useful to visually inspect DIF—especially non-uniform DIF.
#'
#' @param resultados Output object from the `piRFA()` function.
#' @param data The dataset used in the `piRFA()` analysis.
#' @param item The name of the item to be plotted (as a character string).
#' @param cov The name of the grouping covariate used in `piRFA()` (as a character string).
#' @param theta_range A numeric vector of length 2 specifying the latent trait range (default is -3 to 3).
#' @param n_points Number of evaluation points along the latent trait (default = 200).
#' @param theme_option A `ggplot2` theme object to customize the plot appearance.
#' Common options include: #' theme_minimal() (default), theme_bw(), theme_classic(),
#' theme_light(), theme_dark(), theme_void(), theme_gray(), theme_linedraw(), theme_test()
#' and others.
#' @param labels Optional character vector to rename groups in the plot legend.
#' If `NULL` (default), the original group names from `data[[cov]]` are used.
#'
#' @details
#' The function uses the parameter estimates from `piRFA()` to compute predicted responses
#' for each group along a latent variable continuum.
#'
#' The formula used is a linear prediction:
#' \deqn{y = intercept + loading * theta + direct_effect * group + interaction * (theta × group)}
#'
#' @return A `ggplot2` object showing the DIF profile for PI-MIMIC output.
#'
#' @export
piRFA.profile <- function(resultados, data, item, cov,
                          theta_range = c(-3, 3),
                          n_points = 200,
                          theme_option = ggplot2::theme_minimal(),
                          labels = NULL) {

  # Load required package
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required but not installed.")
  }

  # Build theta sequence
  theta_seq <- seq(from = theta_range[1], to = theta_range[2], length.out = n_points)

  # Get coefficients
  sepc_data <- resultados$SEPC
  coefs_item <- sepc_data[sepc_data$Item == item, ]

  # Extract parameters
  intercept   <- coefs_item$EPC[grepl("Intercept", coefs_item$Effect)]
  loading     <- coefs_item$EPC[grepl("LFac", coefs_item$Effect)]
  direct      <- coefs_item$EPC[grepl(paste0("~", cov), coefs_item$Effect)]
  interaction <- coefs_item$EPC[grepl(paste0("LFacX", cov), coefs_item$Effect)]

  # Defaults
  if (length(intercept) == 0)   intercept <- 0
  if (length(loading) == 0)     loading <- 1
  if (length(direct) == 0)      direct <- 0
  if (length(interaction) == 0) interaction <- 0

  # Predictions by group
  groups <- unique(data[[cov]])
  pred_df <- do.call(rbind, lapply(groups, function(g) {
    group_indicator <- ifelse(g == groups[1], 0, 1)
    theta <- theta_seq
    y_pred <- intercept + loading * theta +
      direct * group_indicator +
      interaction * theta * group_indicator
    data.frame(theta = theta, predicted = y_pred, group = as.character(g))
  }))

  # Rename groups if labels are provided
  if (!is.null(labels)) {
    if (length(labels) != length(unique(pred_df$group))) {
      stop("Length of 'labels' must match the number of groups.")
    }
    pred_df$group <- factor(pred_df$group,
                            levels = groups,
                            labels = labels)
  } else {
    pred_df$group <- factor(pred_df$group, levels = groups)
  }

  # Create the plot
  p <- ggplot2::ggplot(pred_df, ggplot2::aes(x = theta, y = predicted, color = group)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::labs(
      title = paste("DIF Profile -", item),
      x = "Latent Trait (θ)",
      y = "Expected Response",
      color = "Group"
    ) +
    theme_option

  return(p)
}
