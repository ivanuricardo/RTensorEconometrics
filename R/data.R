#' GVAR Traditional Data
#' 
#' A flattened version of the Global VAR data. Captures GDP, inflation, and 
#' interest rate of 32 countries. For more information, go to
#' https://www.mohaddes.org/gvar. 
#' 
#' @format: A matrix with 161 observations and 96 variables. 
#'  
#' \describe{
#'  \item{y.ARGENTINA}{Argentina GDP}
#'  \item{Dp.ARGENTINA}{Argentina Inflation Rate}
#'  \item{r.ARGENTINA}{Argentina Interest Rate}
#' }
#' 
#' @source: <https://www.mohaddes.org/gvar>
"traditional_data"

#' GVAR Tensor Data
#' 
#' Global VAR data in tensor form. Data format is a 161x32x3 array with 161
#' observations, 32 countries, and 3 economic indicators per country.
#' 
#' @format An array with 161 observations, 32 countries, and 3 economic variables
#' 
#' @source https://www.mohaddes.org/gvar.
"tensor_data"