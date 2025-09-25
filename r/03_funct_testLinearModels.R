# Function: testLinearModels
# Fits two linear regression models to assess the relationship between lagged AQI and COPD case counts.
# Tests for a specified county name, the dataframe containinf AQI and COPD data, and a gapPeriod value to determine how much to shift AQI values.
# Outputs a list containing linear model summaries for OverallCount and OverallCount/OverallPopulation according to gapPeriod

testLinearModels <- function(countyName, data, gapPeriod) {
  
  countyData <- subset(data, County == countyName)
  
  countyData$LaggedAQI <- c(rep(NA, gapPeriod), countyData$AQI[1:(nrow(countyData) - gapPeriod)])
  
  countyData <- na.omit(countyData)
  
  model1 <- lm(OverallCount ~ LaggedAQI, data = countyData)
  model2 <- lm((OverallCount / OverallPopulation) ~ LaggedAQI, data = countyData)
  
  return(list(CountModel = summary(model1), CountPopModel = summary(model2)))
}