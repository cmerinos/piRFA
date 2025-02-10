#' @title PI-MIMIC graphics: EPC and SEPC
#' 
#' @description
#' La función `piRFA.plot()` genera gráficos de los coeficientes EPC (EPC y SEPC.ALL) para DIF uniforme y no uniforme, 
#' facilitando la interpretación de los efectos DIF en los ítems evaluados por la función \code{\link{piRFA}}.
#' 
#' @param resultados Lista de resultados generados por \code{\link{piRFA}}.
#' @param cov Nombre de la covariable usada en \code{\link{piRFA}}.
#' @param theme_option Tema para `ggplot2` (por defecto `theme_minimal()`).
#' @param bar_colors Colores para los tipos de DIF (por defecto azul y rojo).
#' @param punto_corte Línea vertical que indica el umbral de DIF (por defecto `0.10`). Puede incluirse un vector, e.g., `c(-0.10, 0.10)`.
#' 
#' @return 
#' Un gráfico de barras con los coeficientes SEPC para DIF uniforme y no uniforme.
#' 
#' @details
#' \code{\link{piRFA.plot}} toma la salida de \code{\link{piRFA}}, centrada en el unstandardized (epc) o standardized
#' metric of the expected parameter change (spec.all). El usuario debe elegir cual quiere graficar.
#' 
#' 
#' @examples
#' ### Example 1 -------------
#' 
#' # Example 1 data
#' set.seed(123)
#' Exmp1.data <- data.frame(
#'  grp = sample(0:1, 100, replace = TRUE),  # Variable de grupo
#'  item1 = sample(1:5, 100, replace = TRUE),
#'  item2 = sample(1:5, 100, replace = TRUE),
#'  item3 = sample(1:5, 100, replace = TRUE))
#' 
#' # Output and save
#' Exmp1.output <- piRFA(data = Exmp1.data , items = c("Item1", "Item2", "Item3"), cov = "grp")
#' 
#' # Create graph
#' piRFA.plot(Exmp1.output, cov = "grp")
#' 
#' # Add cutoff range
#' piRFA.plot(Exmp1.output, cov = "grp", punto_corte = c(-.10, .10))
#' 
#' 
#' @seealso
#' \code{\link{ggplot2}}, \code{\link{piRFA}}
#'
#' @export

piRFA.plot <- function(resultados, cov, 
                       metric = "SEPC.ALL",  # ¡Aquí cambiamos el nombre de la columna!
                       theme_option = theme_minimal(),  # Tema de ggplot
                       bar_colors = c("DIF Uniforme" = "blue", "DIF No Uniforme" = "red"),
                       punto_corte = 0.10) {  # Línea vertical en SEPC para DIF alto
  
  library(ggplot2)
  
  # **Validación del argumento `metric`**
  if (!metric %in% c("EPC", "SEPC.ALL")) {  # ¡También corregimos aquí!
    stop("Error: El argumento 'metric' debe ser 'EPC' o 'SEPC.ALL'.")
  }
  
  # **Extraer datos de SEPC según la métrica seleccionada**
  df_sepc <- resultados$SEPC
  
  # Validar si SEPC tiene datos
  if (nrow(df_sepc) > 0) {
    
    # Crear nombres dinámicos para la covariable en SEPC
    cov_lat <- paste0(cov, "lat")
    fac_interaccion <- paste0("LFacX", cov)
    
    # Filtrar filas de SEPC dinámicamente según la covariable utilizada
    df_sepc_filtrado <- df_sepc[grepl(paste0(cov_lat, "|", fac_interaccion), df_sepc$Effect), ]
    
    # Corregir nombres de DIF
    df_sepc_filtrado$Tipo <- ifelse(grepl(cov_lat, df_sepc_filtrado$Effect), 
                                    "DIF Uniforme", "DIF No Uniforme")
    
    # Asegurar que la variable `Tipo` tiene los niveles correctos
    df_sepc_filtrado$Tipo <- factor(df_sepc_filtrado$Tipo, levels = c("DIF Uniforme", "DIF No Uniforme"))
    
    # Generar gráfico con la métrica seleccionada
    p1 <- ggplot(df_sepc_filtrado, aes(x = Item, y = !!sym(metric), fill = Tipo)) +
      geom_bar(stat = "identity", position = "dodge") +
      geom_hline(yintercept = punto_corte, linetype = "dashed", color = "black") +  # Línea de corte
      coord_flip() +
      labs(title = paste(metric, "por Ítem - Covariable:", cov), 
           x = "Ítems",
           y = paste(metric, "Estandarizado")) +
      scale_fill_manual(values = bar_colors) +  # Colores personalizables
      theme_option  # Tema de ggplot
    
    print(p1)
    
  } else {
    cat("No hay datos suficientes para generar el gráfico SEPC.\n")
  }
}
