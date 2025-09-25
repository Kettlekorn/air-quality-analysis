# Function: processAqiYearly
# Takes a CSV file as an input and filters to only include rows with State.Name containing California
# Extracts the year from the Date.Local column and stores it in a new column, Year
# It then calculates the yearly average AQI per county before returning it
# It also prints the output for a call of countUniqueCounties to each processed file.
processAqiYearly <- function(fileName) {
  
  aqiData <- read.csv(fileName, stringsAsFactors = FALSE)
  
  colnames(aqiData) <- gsub(" ", ".", colnames(aqiData))
  
  caData <- aqiData[aqiData$State.Name == "California", ]
  
  caData$Year <- substr(caData$Date.Local, 1, 4)
  
  yearlyAqi <- aggregate(AQI ~ Year + County.Name, data = caData, FUN = mean, na.rm = TRUE)
  
  return(yearlyAqi)
}