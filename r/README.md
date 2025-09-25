# README FOR P-VALUE ANALYSIS IN COPD AND AIR QUALITY DATA
# By Brandon Wu

This project includes data files, R scripts, R Markdown files, and the
final report. Below is a breakdown of each file’s role:


## Packages:

"00_requirements.R" - Installs usmap, ggplot2, and maps for converting 
county names to FIPS codes before displaying frequency of p-values per county 
on a county separated map of California.

## Functions:

"01_funct_countUniqueCounties.R" - Counts the number of unique counties
in AQI data sets for validation during cleaning. Takes the AQI data
frame as an input and outputs an integer value.

"02_funct_processAqiYearly.R" - Takes the file name of an AQI data set,
then filters, cleans, and aggregate AQI data for each year, returning
annual PM10 concentrations per county. Outputs a processed data set that
is compiled into a singular comprehensive table.

"03_funct_testLinearModels.R" - Performs linear regression to model the
relationship between PM10 exposure and COPD hospitalizations,
incorporating lag periods. Takes county name (string), merged data set
(data frame), and gap period (integer) to output two linear regression
model objects.

"04_funct_getPValue.R" - Extracts p-value from regression model objects.

"05_funct_extractPValuesForCounties.R" - Given a vector list of county
names, the function iterates through counties and uses getPValue across
various lag periods (1-9 years), outputting the data frame with
extracted p-values for different gap periods.

"06_funct_bhProcedure.R" - Given a data frame with p-values in a column,
the Benjamini-Hochberg correction is implemented, adjusting for false
positives in multiple hypothesis testing before outputting a data frame
with adjusted significance threshold and marked p-values for
significance.

## R Markdown Files:

"AQICOPDProject.Rmd" - The main analysis script, containing all data
cleaning, statistical modeling, and visualization steps, along with
generating the merged dataset and p-values.

"FinalReport_WU_BRANDON.Rmd" - The final report, summarizing the
findings along with explanations, regression results, and
visualizations.

## Data Files:

"AllYearsAqiCA.csv" - Aggregated AQI data set containing PM10
concentration per county per year.

"cleanedAqiData.csv" - Processed AQI data set with California-only data,
standardized county names, and yearly averages.

"cleanedCOPDdata.csv" - Filtered COPD hospitalization data set,
including only COPD cases for individuals aged 40+ and removing missing
data.

"mergedData.csv" - Final data set combining AQI and COPD data by county
and year, used for regression modeling.

"countyEntryCounts.csv" - Summary file counting data availability per
county, showing how many years of COPD and AQI data overlap.

"pValuesDF.csv" - Extracted p-values from regression models, before
applying multiple testing correction.

"pValuesCount.csv" - Processed p-values for the model analyzing COPD
hospitalization counts.

"pValuesCountPop.csv" - Processed p-values for the model analyzing COPD
hospitalization rates per population.

## AQI Datasets:

"daily_81102_1990.csv" "daily_81102_1991.csv" "daily_81102_1992.csv"
"daily_81102_1993.csv" "daily_81102_1994.csv" "daily_81102_1995.csv"
"daily_81102_1996.csv" "daily_81102_1997.csv" "daily_81102_1998.csv"
"daily_81102_1999.csv" "daily_81102_2000.csv" "daily_81102_2001.csv"
"daily_81102_2002.csv" "daily_81102_2003.csv" "daily_81102_2004.csv"
"daily_81102_2005.csv" "daily_81102_2006.csv" "daily_81102_2007.csv"
"daily_81102_2008.csv" "daily_81102_2009.csv" "daily_81102_2010.csv"
"daily_81102_2011.csv" "daily_81102_2012.csv" "daily_81102_2013.csv"
"daily_81102_2014.csv"

-   These files contain daily PM10 concentration levels and are used to
    generate yearly average AQI per county before merging with COPD
    hospitalization data.

## Hospitalization Dataset:

"CAHospitalizations.csv" - Original dataset downloaded and renamed from the
California Health and Human Services Agency archive

## Plot Images:
(For visualization in the final report)

"Count_Barplot" - Bar plot of significant P-values for count

"Count_Population_Barplot" - Bar plot of significant P-values for count/population

"CAPPlot" - County-separated map of California, following gradient from white 
to red, depending on frequency of significant P-values.

To reproduce the analysis, only the hospitilization data (CAHospitalizations.csv) and  daily AQI datasets
(daily_81102_1990.csv to daily_81102_2014.csv) are required when the
code is run in AQICOPDProject.Rmd. These daily AQI files serve as input,
and the project processes them to generate cleanedAqiData.csv, cleanedCOPDdata.csv, 
and mergedData.csv, which are then saved for direct loading in the final report. If not already 
installed, run the script in 00_requirements.R to download packages usmap, ggplot2, and maps.
