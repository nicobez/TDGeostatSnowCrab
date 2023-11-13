################################
### Reading and preparing data 
################################

################################
### Snow crab 
################################
myDf <- read.table(file="Data/sc_2008_2010.txt",header=T) 

################################
# Transforming coordinates : centred + correction by the cosine of the latitude
# The new coordinates are in KILOMETERS
################################
myDf$lonKM <- (myDf$Lon-mean(myDf$Lon))*60*1852*cos(myDf$Lat/180*pi)/1000
myDf$latKM <- (myDf$Lat-mean(myDf$Lat))*60*1852/1000

################################
# Adding a density of sc in g/m2
# Weight are in kg 
# Ares_swept are in m2
################################
myDf$density <- myDf$Commercial.weight*1000/myDf$Area_swept


################################
### Polygon
################################
temp <- read.table(file="Data/file03_Polygon_gulf.txt", header=F,sep="\t") 
myPolyLine <- list(x=(-temp[,2]-mean(myDf$Lon))*60*1852*cos(temp[,1]/180*pi)/1000,
                   y=(temp[,1]-mean(myDf$Lat))*60*1852/1000)
