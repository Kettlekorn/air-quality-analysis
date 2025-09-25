# Function: getPValue
# Takes a linear model object (from testLinearModels output),  and returns the p-value for the independent variable (Lagged AQI)

getPValue <- function(model) {
  coefs <- model$coefficients
  return(coefs[2, 4])
}
