#' Factors Zoo Dataset
#'
#' Monthly values of portfolio returns and factors from July 1976 to December
#' 2017, originally studied in Wan et al. (2024). The dataset includes 18
#' classical, publicly available factors such as Fama-French five factors,
#' q-factors, and intermediary asset pricing factors. The dataset also contains 
#' 202 excess portfolio returns, constructed as the difference between raw 
#' portfolio returns and the risk-free rate.
#'
#' @format A data frame with 498 rows and 221 variables:
#' \describe{
#'   \item{date}{Monthly date index, from July 1976 to December 2017.}
#'   \item{MktRF, SMB, HML, RMW, CMA}{Fama-French five factors.}
#'   \item{HML_Devil, gma, orgcap, BAB, QMJ, HXZ_IA, HXZ_ROE, convind, 
#'   cash, hire, gad, ala}{13 additional classical asset pricing factors.}
#'   \item{portfolio1, ..., portfolio202}{Excess returns of 202 portfolios,
#'     computed as raw portfolio return minus the risk-free rate.}
#' }
#'
#' @references
#' Wan, R., Li, Y., Lu, W. & Song, R. (2024). Mining the factor zoo: Estimation 
#' of  latent factor models with sufficient proxies. \emph{Journal of 
#' Econometrics}, 239(2), 105386.
#'
#' Fama, E. F. and French, K. R. (2015). A five-factor asset pricing model.
#' \emph{Journal of Financial Economics}, 116(1), 1-22.
#'
#' Hou, K., Xue, C. and Zhang, L. (2015). Digesting anomalies: An investment 
#' approach. \emph{The Review of Financial Studies}, 28(3), 650-705.
#'
#' He, Z., Kelly, B. and Manela, A. (2017). Intermediary asset pricing: New 
#' evidence from many asset classes.
#' \emph{Journal of Financial Economics}, 126(1), 1-35.
"asset"


#' Boston Housing Dataset
#'
#' A classical benchmark dataset for regression analysis, containing housing
#' values and socioeconomic variables across census tracts in the Boston
#' Standard Metropolitan Statistical Areas.
#'
#' @format A data frame with 506 rows and 14 variables:
#' \describe{
#'   \item{CRIM}{Per capita crime rate by town.}
#'   \item{ZN}{Proportion of residential land zoned for lots over 25,000
#'     square feet.}
#'   \item{INDUS}{Proportion of non-retail business acres per town.}
#'   \item{CHAS}{Charles River dummy variable.}
#'   \item{NOX}{Nitric oxides concentration.}
#'   \item{RM}{Average number of rooms per dwelling.}
#'   \item{AGE}{Proportion of owner-occupied units built prior to 1940.}
#'   \item{DIS}{Weighted distances to five Boston employment centres.}
#'   \item{RAD}{Index of accessibility to radial highways.}
#'   \item{TAX}{Full-value property-tax rate per USD 10,000.}
#'   \item{PTRATIO}{Pupil-teacher ratio by town.}
#'   \item{B}{Proportion of African Americans by town}
#'   \item{LSTAT}{Percentage of lower status of the population.}
#'   \item{MEDV}{Median value of owner-occupied homes in USD 1000s.}
#' }
#'
#' @references
#' Harrison, D. and Rubinfeld, D. L. (1978). Hedonic housing prices and the
#' demand for clean air. \emph{Journal of Environmental Economics and
#' Management}, 5(1), 81-102.
"house"