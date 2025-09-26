#' DIF Profile Plot for PI-MIMIC output (Multiple-Indicatros Multiple-Causes with Product of Indicators)
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
#' @param theme_option  A `ggplot2` theme object to customize the plot appearance.
#' Common options include:
#' `theme_minimal()` (default), `theme_bw()`, `theme_classic()`,
#' `theme_light()`, `theme_dark()`, `theme_void()`, `theme_gray()`, `theme_linedraw()`,
#' `theme_test()` and others..
#'
#' @details
#' The function uses the parameter estimates from `piRFA()` to compute predicted responses
#' for each group along a latent variable continuum. It is especially helpful in visualizing
#' non-uniform DIF, where the slope or pattern of response varies by group.
#'
#' The formula used is a linear prediction:
#' \deqn{y = intercept + loading * theta + direct_effect * group + interaction * (theta × group)}
#'
#' This function assumes a linear model and is best interpreted when DIF is suspected or detected.
#' 
#' @return A `ggplot2` object showing the DIF profile for PI-MIMIC output 
#' (Multiple-Indicatros Multiple-Causes with Product of Indicators).
#' 
#'
#' @examples
#' \dontrun{
#' # Run piRFA
#' results <- piRFA(data = mydata, items = c("item1", "item2"), cov = "group")
#'
#' # DIF profile for item1
#' piRFA.profile(resultados = results, data = mydata, item = "item1", cov = "group")
#' }
#' 
#' @export
piRFA.profile <- function(resultados, data, item, cov,
                          theta_range = c(-3, 3),
                          n_points = 200,
                          theme_option = ggplot2::theme_minimal()) {
  
  # Load required package
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required but not installed.")
  }
  
  # Build theta sequence
  theta_seq <- seq(from = theta_range[1], to = theta_range[2], length.out = n_points)
  
  # Get coefficients
  sepc_data <- resultados$SEPC
  coefs_item <- sepc_data[sepc_data$Item == item, ]
  
  # Extract direct effect (uniform DIF)
  uniform_effect <- coefs_item$EPC[grepl(cov, coefs_item$Effect)]
  
  # Extract interaction effect (non-uniform DIF)
  interaction_effect <- coefs_item$EPC[grepl(paste0("LFacX", cov), coefs_item$Effect)]
  
  # Simulate predicted responses
  pred_df <- data.frame()
  
  for (g in unique(data[[cov]])) {
    theta <- theta_seq
    direct <- ifelse(length(uniform_effect) > 0, uniform_effect, 0)
    interaction <- ifelse(length(interaction_effect) > 0, interaction_effect, 0)
    
    y_pred <- theta + direct * (data[[cov]][1] == g) + interaction * theta * (data[[cov]][1] == g)
    
    temp_df <- data.frame(
      theta = theta,
      predicted = y_pred,
      group = as.character(g)
    )
    
    pred_df <- rbind(pred_df, temp_df)
  }
  
  # Rename factor levels if binary
  if (length(unique(data[[cov]])) == 2) {
    levels <- unique(data[[cov]])
    pred_df$group <- factor(pred_df$group,
                            levels = levels,
                            labels = c("Group 1", "Group 2"))
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
  
  print(p)
}
