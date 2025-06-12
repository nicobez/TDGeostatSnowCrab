################################
### TD Geostatistique 
### nicolas.bez@ird.fr
### octobre 2024
################################

rm(list=ls())

source("Script/loadingLibrariesAndAddsOnFunctions.R")

source("Script/readingAndPreparingDataFrame.R")

### Setting the y-axis to be labelled horizontally
### and saving the default graphical environment
par(las=1)
defPar <- par(no.readonly = TRUE)

source("Script/geostatisticAnalyses.R")

