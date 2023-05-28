test_that("data_aligns", {
  # Compare the levels time series and the differenced time series, make sure they 
  # correspond to the correct dimension
  library(readxl)
  library(dplyr)
  
  ###### Constructing Data
  pathname <- "/home/ivan/Desktop/Projects/Research/TensorEconometrics/data-raw/GVARDatabase(1979Q2-2019Q4)/CountryData.xls"
  sheet_names <- excel_sheets(pathname)
  data <- read_excel(pathname, sheet = 1)
  gvar_data <- data %>% 
    mutate(country = sheet_names[1])
  
  for (sheet in (2:(length(sheet_names)-1))) {
    appending_data <- read_excel(pathname, sheet = sheet)
    appending_data <- appending_data%>%mutate(country = sheet_names[sheet])
    new_data <- plyr::rbind.fill(gvar_data, appending_data)
    gvar_data <- new_data
  }
  
  data("tensor_data")
  
  # Take differences of gdp for Argentina and remove first observation
  diff_Arg_y <- diff(gvar_data[gvar_data$country == "ARGENTINA", "y"])[2:162]
  
  # Verify this matches with first fiber of the tensor
  expect_equal(tensor_data[,1,1], diff_Arg_y)
  
  # Verify diff of gdp for Belgium also aligns
  diff_Belgium_y <- diff(gvar_data[gvar_data$country == "BELGIUM", "y"])[2:162]
  expect_equal(tensor_data[,4,1], diff_Belgium_y)
  
  # Verify levels of Dp for Spain matches. Note we remove first two obs for differences
  Spain_Dp <- gvar_data[gvar_data$country == "SPAIN", "Dp"][3:163]
  expect_equal(tensor_data[,26,2], Spain_Dp)
  
  # Verify levels of r for Turkey
  Turkey_r <- gvar_data[gvar_data$country == "TURKEY", "r"][3:163]
  expect_equal(tensor_data[,30, 3], Turkey_r)
})
