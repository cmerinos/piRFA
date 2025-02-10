# PIRFA: Product of Indicators (PI) for MIMIC/RFA Models in DIF detection

<!-- badges: start -->

<!-- badges: end -->

**piRFA** is an R package for detecting **Differential Item Functioning (DIF)** using the **Product of Indicators (PI)** approach within a **MIMIC**/**RFA** framework (**Multiple Indicators Multiple Causes**/**Restricted Factor Analysis**).

## Installation

You can install the development version of PIRFA like so:

``` r
# Install from GitHub
devtools::install_github("cmerinos/piRFA")
library(piRFA)
```

## Example

Basic example:

``` r
library(piRFA)

# Load data
set.seed(123)
example_data <- data.frame(
  group = sample(0:1, 100, replace = TRUE), 
  item1 = sample(1:5, 100, replace = TRUE),
  item2 = sample(1:5, 100, replace = TRUE),
  item3 = sample(1:5, 100, replace = TRUE))

# Run DIF analysis
results <- piRFA(data = example_data, items = c("item1", "item2", "item3"), cov = "group")

results

# View specific results
print(results$DIF_Global)
print(results$SEPC)


# Plot results
piRFA.plot(results, cov = "group")
```

## References

Kolbe, L., & Jorgensen, T. D. (2018). Using product indicators in restricted factor analysis models to detect nonuniform measurement bias. In M. Wiberg, S. A. Culpepper, R. Janssen, \#' J. González, & D. Molenaar (Eds.), Quantitative psychology: The 82nd Annual Meeting of the Psychometric Society, Zurich, Switzerland, 2017 (pp. 235–245). New York, NY: Springer. [https://doi.org/10.1007/978-3-319-77249-3_20](doi:https://doi.org/10.1007/978-3-319-77249-3_20){.uri}

Kolbe, L., & Jorgensen, T. D. (2019). Using restricted factor analysis to select anchor items and detect differential item functioning. Behavior Research Methods, 51, 138–151. <https://doi.org/10.3758/s13428-018-1151-3>

Kolbe, L., Jorgensen, T. D., & Molenaar, D. (2020). The Impact of Unmodeled Heteroskedasticity on Assessing Measurement Invariance in Single-group Models. Structural Equation Modeling: A Multidisciplinary Journal, 28(1), 82–98. <https://doi.org/10.1080/10705511.2020.1766357>
