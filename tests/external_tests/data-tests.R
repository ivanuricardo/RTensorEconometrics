# Test that tensor_data gives the correct output
# Compare the levels time series and the differenced time series, make sure they 
# correspond to the correct dimension
library(readxl)
library(dplyr)

###### Constructing Data
pathname <- "data/GVARDatabase(1979Q2-2019Q4)/CountryData(1979Q2-2019Q4).xls"
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
tensor_data[,1,1] == diff_Arg_y

# Verify diff of gdp for Belgium also aligns
diff_Belgium_y <- diff(gvar_data[gvar_data$country == "BELGIUM", "y"])[2:162]
tensor_data[,4,1] == diff_Belgium_y

