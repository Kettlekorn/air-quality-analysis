# Function: bhProcedure
# To perform the Benjamini-Hochberg Procedure and identify statistically significant p-values while controlling the false discovery rate
# Takes a dataframe containing the p-values extracted from multiple tests for Count and Count/Population with AQI and the defaulted fdrLevel of 0.05 for statistical significance.
# Then outputs the dataframe with added columns for ranking, BH threshold, and significance determined by a boolean statement.
# Follows the formula Threshold = ( Rank / Total P-Values ) × fdrLevel, sorting the dataframe by ascending p-values and assigning ranks beforehand.

bhProcedure <- function(dataFrame, fdrLevel = 0.05) {
  
  dataFrame <- dataFrame[order(dataFrame$PValue), ]
  
  totalPValues <- nrow(dataFrame)
  
  dataFrame$Rank <- seq_len(totalPValues)

  dataFrame$BHThreshold <- (dataFrame$Rank / totalPValues) * fdrLevel
  
  significantIndex <- max(which(dataFrame$PValue <= dataFrame$BHThreshold), na.rm = TRUE)
  
  dataFrame$Significant <- dataFrame$Rank <= significantIndex
  
  return(dataFrame)
}