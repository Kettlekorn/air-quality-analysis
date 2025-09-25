# Function: extractPValuesForCounties
# Uses testLinearModels for multiple counties across different lag periods of 1-9 years. 
# Takes a vector containing the names of counties to perform the function on and the merged dataframe containing AQI and COPD data
# Creates an empty dataframe for storing results and extracts p-values for hypothesis testing.
# It then outputs the resulting dataframe with added entries to County, GapPeriod, Model, and PValue

extractPValuesForCounties <- function(counties, data) {
  
  results <- data.frame(County = character(), GapPeriod = integer(), Model = character(), PValue = numeric())
  
  for (county in counties) {
    
    for (gap in 1:9) {
      
      modelResults <- testLinearModels(county, data, gap)
      
      pValueCount <- getPValue(modelResults$CountModel)
      pValueCountPop <- getPValue(modelResults$CountPopModel)
      
      results <- rbind(results, data.frame(County = county, GapPeriod = gap, Model = "Count", PValue = pValueCount))
      
      results <- rbind(results, data.frame(County = county, GapPeriod = gap, Model = "Count/Population", PValue = pValueCountPop))
    }
  }
  
  return(results)
}