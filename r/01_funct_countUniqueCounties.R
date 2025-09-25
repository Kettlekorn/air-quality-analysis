# Function: countUniqueCounties
# Takes a CSV file assumed to only contain entries with State.Name California as an input and returns the integer count of unique California counties.

countUniqueCounties <- function(caData) {
  uniqueCountyCount <- length(unique(caData$County.Name))
  return(uniqueCountyCount)
}