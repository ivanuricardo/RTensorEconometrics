#' Global Vector Autoregression Data Set
#'
#' This data set contains a global vector autoregression (VAR) model with 5
#' countries and 3 economic variables. Each column corresponds to a specific
#' country and economic variable combination.
#'
#' @format A time series data frame with 120 observations of 15 variables:
#'   \describe{
#'     \item{USA_GDP}{Gross Domestic Product (GDP) for the United States}
#'     \item{USA_CPI}{Consumer Price Index (CPI) for the United States}
#'     \item{USA_IR}{Interest rate for the United States}
#'     \item{CAN_GDP}{Gross Domestic Product (GDP) for Canada}
#'     \item{CAN_CPI}{Consumer Price Index (CPI) for Canada}
#'     \item{CAN_IR}{Interest rate for Canada}
#'     \item{UK_GDP}{Gross Domestic Product (GDP) for the United Kingdom}
#'     \item{UK_CPI}{Consumer Price Index (CPI) for the United Kingdom}
#'     \item{UK_IR}{Interest rate for the United Kingdom}
#'     \item{GER_GDP}{Gross Domestic Product (GDP) for Germany}
#'     \item{GER_CPI}{Consumer Price Index (CPI) for Germany}
#'     \item{GER_IR}{Interest rate for Germany}
#'     \item{JPN_GDP}{Gross Domestic Product (GDP) for Japan}
#'     \item{JPN_CPI}{Consumer Price Index (CPI) for Japan}
#'     \item{JPN_IR}{Interest rate for Japan}
#'   }
#'
#' @source This data set is simulated for demonstration purposes.
#'
#' @examples
#' data(global_var)
#' head(global_var)
#'
#' @keywords data
#' @docType data
#' @name traditional_data
"traditional_data"

#' Global Vector Autoregression Data Set (Array Format)
#'
#' This data set contains a global vector autoregression (VAR) model with 5
#' countries and 3 economic variables. The data is rearranged into an array
#' with dimensions 161 x 32 x 5, where the first dimension represents the time
#' index, the second dimension represents the country index, and the third
#' dimension represents the economic variable index.
#'
#' @format An array with dimensions 161 x 32 x 5.
#'
#' @dimnames The array has the following dimension names:
#'   \describe{
#'     \item{Time_Index}{The time index (1 to 161)}
#'     \item{Country_Index}{The country index (1 to 32)}
#'     \item{Economic_Variable_Index}{The economic variable index (1 to 5)}
#'   }
#'
#' @dim The array has 161 rows, 32 columns, and 5 economic variables.
#'
#' @source This data set is simulated for demonstration purposes.
#'
#' @examples
#' data(global_var_array)
#' dim(global_var_array)
#' dimnames(global_var_array)
#'
#' @keywords data
#' @docType data
#' @name tensor_data
"tensor_data"

#' Global Vector Autoregression Data Set LEVELS
#'
#' This data set contains a global vector autoregression (VAR) model with 5
#' countries and 3 economic variables. Each column corresponds to a specific
#' country and economic variable combination.
#'
#' @format A time series data frame with 120 observations of 15 variables:
#'   \describe{
#'     \item{USA_GDP}{Gross Domestic Product (GDP) for the United States}
#'     \item{USA_CPI}{Consumer Price Index (CPI) for the United States}
#'     \item{USA_IR}{Interest rate for the United States}
#'     \item{CAN_GDP}{Gross Domestic Product (GDP) for Canada}
#'     \item{CAN_CPI}{Consumer Price Index (CPI) for Canada}
#'     \item{CAN_IR}{Interest rate for Canada}
#'     \item{UK_GDP}{Gross Domestic Product (GDP) for the United Kingdom}
#'     \item{UK_CPI}{Consumer Price Index (CPI) for the United Kingdom}
#'     \item{UK_IR}{Interest rate for the United Kingdom}
#'     \item{GER_GDP}{Gross Domestic Product (GDP) for Germany}
#'     \item{GER_CPI}{Consumer Price Index (CPI) for Germany}
#'     \item{GER_IR}{Interest rate for Germany}
#'     \item{JPN_GDP}{Gross Domestic Product (GDP) for Japan}
#'     \item{JPN_CPI}{Consumer Price Index (CPI) for Japan}
#'     \item{JPN_IR}{Interest rate for Japan}
#'   }
#'
#' @source This data set is simulated for demonstration purposes.
#'
#' @examples
#' data(global_var)
#' head(global_var)
#'
#' @keywords data
#' @docType data
#' @name traditional_data
"traditional_data_levels"

#' Global Vector Autoregression Data Set (Array Format, LEVELS)
#'
#' This data set contains a global vector autoregression (VAR) model with 5
#' countries and 3 economic variables. The data is rearranged into an array
#' with dimensions 161 x 32 x 5, where the first dimension represents the time
#' index, the second dimension represents the country index, and the third
#' dimension represents the economic variable index.
#'
#' @format An array with dimensions 161 x 32 x 5.
#'
#' @dimnames The array has the following dimension names:
#'   \describe{
#'     \item{Time_Index}{The time index (1 to 161)}
#'     \item{Country_Index}{The country index (1 to 32)}
#'     \item{Economic_Variable_Index}{The economic variable index (1 to 5)}
#'   }
#'
#' @dim The array has 161 rows, 32 columns, and 5 economic variables.
#'
#' @source This data set is simulated for demonstration purposes.
#'
#' @examples
#' data(global_var_array)
#' dim(global_var_array)
#' dimnames(global_var_array)
#'
#' @keywords data
#' @docType data
#' @name tensor_data
"tensor_data_levels"