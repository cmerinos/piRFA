#' @title Función para Analizar DIF con PI-MIMIC (Product of Indicators)
#'
#' @description
#' La función `piRFA` analiza el funcionamiento diferencial de ítems (DIF)
#' utilizando el multiple-indicators multiple-causes (MIMIC) framework with product of indicators (PI). Se apoya en el `lavaan` y `scripty`.
#' Permite evaluar DIF uniforme y no uniforme en una escala de medición, mediante pruebas estadísticas y una aproximación del tamaño del efecto.
#'
#' @param data DataFrame que contiene los ítems y la covariable.
#' @param items Vector de nombres de los ítems dentro de `data`.
#' @param cov Nombre de la covariable en `data` (puede ser categórica o numérica). Si la covariable es categórica, entonces debe estar en tipo factor.
#' @param lvname Nombre creado para la variable latente en el modelo (por defecto `"LatFact"`).
#' @param est Abreviatura del estimador elegido; opciones: "MLM", "MLR", "ULSMV" (por defecto `"MLM"`).
#'
#' @return
#' La función devuelve una lista con los siguientes DataFrames:
#'  \itemize{
#'    \item \code{DIF_Global} - Resultados global de DIF.
#'    \item \code{DIF_Uniforme} - Resultados DIF Uniforme.
#'    \item \code{DIF_NoUniforme} - Resultados DIF No Uniforme.
#'    \item \code{SEPC} - Coeficientes SEPC para DIF uniforme y no uniforme.
#'  }
#'
#' @details
#' Esta función implementa el análisis de Funcionamiento Diferencial de Ítems (DIF) utilizando
#' el enfoque de Producto de Indicadores (PI; Kolbe & Jorgensen, 2018) dentro del marco
#' de Restricted Factor Analysis (RFA; Oort, 1998). Este método opera bajo un esquema
#' MIMIC (Finch, 2005), incorporando interacciones entre variables latentes mediante PI (Kolbe et al., 2018, 2019;
#' Kolbe, Jorgensen, & Molenaar, 2020). Permite evaluar DIF uniforme (uDIF) y DIF no uniforme (nuDIF)
#' con covariables tanto categóricas (ej., sexo) como continuas (ej., autoestima, conciencia).
#'
#' La estimación se realiza mediante `lavaan::cfa`, y las pruebas estadísticas de DIF se basan
#' en el likelihood ratio test (LRT), utilizando una comparación entre un modelo no restringido
#' (con parámetros de efecto de la covariable y su interacción estimados libremente) y un modelo
#' restringido (con estos parámetros fijados en cero), evaluado mediante una prueba chi-cuadrado
#' con 2 grados de libertad (`df = 2`).
#'
#' La implementación se basa en `lavaan::lavTestScore` y `lavaan::cfa`, integrando funcionalidades
#' modificadas de `scripty::mimic` y `scripty::prods`. A diferencia de `scripty`, esta función
#' evita la creación de archivos externos, permitiendo que toda la salida sea accesible en la sesión de R.
#'
#' Como indicador de tamaño del efecto, `piRFA()` adopta el enfoque de releasing fixed-to-zero
#' parameters (Whittaker, 2012; Garnier-Villarreal & Jorgensen, 2024). Se estiman el expected parameter change
#' (`epc`) y el standardized expected parameter change (`sepc.all`), siguiendo la metodología de Chou & Bentler (1993).
#' `sepc.all` representa el cambio en unidades estandarizadas al liberar la restricción del valor cero
#' en los parámetros analizados de los ítems. Valores `sepc.all ≥ .20` se consideran indicativos de
#' una posible mala especificación del modelo (Whittaker, 2012), mientras que valores `≥ .10`
#' (Kaplan, 1989) representan un límite para un efecto moderado.
#'
#' El usuario puede seleccionar el estimador apropiado entre `"MLM"`, `"MLR"` y `"ULSMV"`,
#' dependiendo del supuesto de normalidad y robustez deseado.
#'
#' @examples
#' ### Example 1 -------------
#'
#' set.seed(123)
#' Exmp1.data <- data.frame(
#'   grp = sample(0:1, 100, replace = TRUE),  # Variable de grupo
#'   item1 = sample(1:5, 100, replace = TRUE),
#'   item2 = sample(1:5, 100, replace = TRUE),
#'   item3 = sample(1:5, 100, replace = TRUE),
#'   item4 = sample(1:5, 100, replace = TRUE))
#'
#'
#' Run DIF analysis, and full output
#' piRFA(data = Exmp1.data , items = c("Item1", "Item2", "Item3"), cov = "grp")
#'
#' Specific output: Uniform DIF (uDIF)
#' Exmp1.output <- piRFA(data = Exmp1.data , items = c("Item1", "Item2", "Item3"), cov = "grp")
#'
#' Exmp1.output$DIF_Uniforme
#'
#' Specific output: full SPEC
#' piRFA(data = Exmp1.data , items = c("Item1", "Item2", "Item3"), cov = "grp")$SPEC
#'
#' Specific output: full SPEC for uniform DIF
#' piRFA(data = Exmp1.data , items = c("Item1", "Item2", "Item3"), cov = "grp")$SPEC[c(1,3), ]
#'
#'
#'
#' @references
#' Chou, C.-P., & Bentler, P. (1993). Invariant standardized estimated parameter change for model modification in covariance structure analysis. \emph{Multivariate Behavioral #' Research}, 28, 97–110. https://doi.org/10.1207/s15327906mbr2801_6
#'
#' Finch, H. (2005). The MIMIC model as a method for detecting DIF: Comparison with Mantel–Haenszel, SIBTEST, and the IRT likelihood ratio. \emph{Applied Psychological Measurement}, #' 29, 278–295. https://doi.org/10.1177/0146621605275728
#'
#' Garnier-Villarreal, M., & Jorgensen, T. D. (2024). Evaluating Local Model Misspecification with Modification Indices in Bayesian Structural Equation Modeling. \emph{Structural #' Equation Modeling: A Multidisciplinary Journal}, 1–15. https://doi.org/10.1080/10705511.2024.2413128
#'
#' Kaplan D. (1989). Model Modification in Covariance Structure Analysis: Application of the Expected Parameter Change Statistic. \emph{Multivariate Behavioral Research}, 24(3), 285–#' 305. https://doi.org/10.1207/s15327906mbr2403_2
#'
#' Kolbe, L., & Jorgensen, T. D. (2018). Using product indicators in restricted factor analysis models to detect nonuniform measurement bias. In M. Wiberg, S. A. Culpepper, R. Janssen, #' J. González, & D. Molenaar (Eds.), Quantitative psychology: The 82nd Annual Meeting of the Psychometric Society, Zurich, Switzerland, 2017 (pp. 235–245). New York, NY: Springer. https://doi.org/10.1007/978-3-319-77249-3_20
#'
#' Kolbe, L., & Jorgensen, T. D. (2019). Using restricted factor analysis to select anchor items and detect differential item functioning. \emph{Behavior Research Methods}, 51, 138–#' 151. https://doi.org/10.3758/s13428-018-1151-3
#'
#' Kolbe, L., Jorgensen, T. D., & Molenaar, D. (2020). The Impact of Unmodeled Heteroskedasticity on Assessing Measurement Invariance in Single-group Models. Structural Equation Modeling: A Multidisciplinary Journal, 28(1), 82–98. https://doi.org/10.1080/10705511.2020.1766357
#'
#' Whittaker, T. A. (2012). Estimation of Standardized Expected Parameter Change for DIF Detection. \emph{Educational and Psychological Measurement}, 72(3), 342-357.
#'
#'
#'
#' @seealso
#' \code{\link[lavaan]{cfa}}, \code{\link[scripty]{mimic}}, \code{\link[scripty]{prods}}
#'
#'@export

piRFA <- function(data, items, cov, lvname = "LatFact", est = "MLM") {

  # Verificar si el estimador ingresado es válido
  if (!est %in% c("MLM", "MLR", "ULSMV")) {
    stop("Error: El estimador debe ser 'MLM', 'MLR' o 'ULSMV'.")
  }

  library(lavaan)

  mimicout_modificado <- function(fit.mimic, mimic.param, cov) {
    ests <- as.data.frame(lavaan::parameterestimates(fit.mimic))
    uniqnames <- unique(ests$lhs)
    lvname <- uniqnames[1]
    covlvname <- uniqnames[2]

    any.out <- do.call(rbind, lapply(mimic.param, function(x)
      lavaan::lavTestScore(fit.mimic, add = x)$test))

    sep.out <- lavaan::lavTestScore(fit.mimic, add = as.character(mimic.param))

    # **Nuevo: Crear `$DIF_Global` con "any.chi2", "any.df", "any.p"**
    df_dif_global <- data.frame(
      DIF_Global_Chi2 = round(any.out$X2, 2),
      df = any.out$df,
      p_value = round(any.out$p, 3)
    )

    # **Mantener `$DIF_Uniforme` y `$DIF_NoUniforme`**
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

    # **Extraer "EPC" y "SEPC_ALL" para `$SEPC`**
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
    stop("Error: Algunos nombres de ítems o la covariable no están en la base de datos.")
  }

  if (!is.character(lvname) || lvname == "") {
    stop("Error: El nombre de la variable latente ('lvname') debe ser un string no vacío.")
  }

  if (is.factor(data[[cov]])) {
    data[[cov]] <- as.numeric(data[[cov]])
  } else if (!is.numeric(data[[cov]])) {
    stop(paste("Error: La covariable", cov, "debe ser un factor o numérica."))
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

  return(resultados_DIF)
}
