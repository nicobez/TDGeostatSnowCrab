
###############################################################
plot(myDf$Lon,myDf$Lat)
map("worldHires",add=T)
names(myDf)

# Defining the working DB (gstlearn)
myDb <- Db_create()
myDb$addColumns(myDf$lonKM,'longitude',ELoc_X(),0)
myDb$addColumns(myDf$latKM,'latitude',ELoc_X(),1)
myDb$addColumns(myDf$Year,'year')
myDb$addColumns(myDf$density,'density')
myDb$addColumns(myDf$Depth,'depth')
myDb$addColumns(myDf$Temperature,'temperature')
myDb

### Mapping the data
#ratio <- 1/cos(47/180*pi)
symbols(myDb['longitude'],myDb['latitude'],circles=sqrt(myDb['density']),
        inches=0.2,bg=rgb(0,0,1,0.5),xlab="Easting",ylab="Norting")
points(myDb['longitude'][myDb['density']==0],
       myDb['latitude'][myDb['density']==0],
       pch=3,col=2)

temp <- map("worldHires",regions="Canada",plot=F)
myCoastLine <- list(x=(temp$x-mean(myDf$Lon))*60*1852*cos(temp$y/180*pi)/1000,
                    y=(temp$y-mean(myDf$Lat))*60*1852/1000)
rm(temp)
lines(myCoastLine)

### Creating a sub-db for year 2010 cleaned from the un_useful data
sel2010 <- myDb['year']==2010
myDb$addSelection(sel2010)
myDb2010 <- Db_clone(myDb)
myDb2010 = Db_createReduce(myDb2010)
myDb2010$deleteColumns('NewSel')

symbols(myDb2010['longitude'],myDb2010['latitude'],circles=sqrt(myDb2010['density']),
        inches=0.15,bg=rgb(0,0,1,0.5),asp=1,xlab="Easting",ylab="Norting")
points(myDb2010['longitude'][myDb2010['density']==0],
       myDb2010['latitude'][myDb2010['density']==0],
       pch=3,col=2)
lines(myCoastLine)

### Polygon : create a gstlearn polygon 
myPoly = Polygons()
poly1 = PolyElem(x=myPolyLine$x, y = myPolyLine$y)
myPoly$addPolyElem(poly1)
rm(poly1)

### Defining the working variable  
myDb2010$setLocator('density',ELoc_Z(),0)

# ### and selecting the inner data points
# db_polygon(myDb2010, myPoly)

##########################
### Variography
##########################

# Create an empty object containing computation's characteristics
myVarioParamOmni = VarioParam()
# Define computation's characteristics
temp = DirParam_create(npas=55,dpas=10)
# Add them to the dedicated object
myVarioParamOmni$addDir(temp)
# Clean
rm(temp)
# Create an empty object containing the output 
myVarioOmni = Vario(myVarioParamOmni)
# Compute the variogram
myVarioOmni$compute(myDb2010,ECalcVario_VARIOGRAM())
# Plot the result
plotVario(myVarioOmni)

### variogramme directionnel
myVarioParamDirectional = VarioParam()
temp = DirParam_createMultiple(ndir=4, npas=40, dpas=10)
myVarioParamDirectional$addMultiDirs(temp)
myVarioDirectional = Vario(myVarioParamDirectional)
myVarioDirectional$compute(myDb2010,ECalcVario_VARIOGRAM())
plotVario(myVarioDirectional)

### variogramme directionnel avec pas et npas variables par direction
myVarioParamDirectional = VarioParam()
temp = DirParam_create(npas=35,dpas=10,angle2D = 0,tolang=25)
myVarioParamDirectional$addDir(temp)
temp = DirParam_create(npas=15,dpas=10,angle2D = 45,tolang=25)
myVarioParamDirectional$addDir(temp)
temp = DirParam_create(npas=15,20,dpas=10,angle2D = 90,tolang=25)
myVarioParamDirectional$addDir(temp)
temp = DirParam_create(npas=35,dpas=10,angle2D = 135,tolang=25)
myVarioParamDirectional$addDir(temp)

myVarioDirectional2 = Vario(myVarioParamDirectional)
myVarioDirectional2$compute(myDb2010,ECalcVario_VARIOGRAM())

plotVario(myVarioDirectional2,inches=0.1)

##########################
### model du vario 2010
##########################

# Printing all possible models
ECov_printAll()

myModel = Model_createFromDb(myDb2010)
myModel$fit(vario=myVarioOmni, types=ECov_fromKeys(c("NUGGET","SPHERICAL","LINEAR")))
#optvar = Option_VarioFit(flag_noreduce = TRUE)) #mauto = Option_AutoFit$setWmode(wmode=2)
plotModel(myVarioOmni,myModel)

##########################
### Krigeage local = cartographie 
##########################

### grille de krigeage
myGrid = DbGrid_createCoveringDb(dbin=myDb2010,dx=c(1,1))
db_polygon(myGrid, myPoly)

# choix d'un voisinage
myNeigh = NeighMoving_create(flag_xvalid=FALSE,nmaxi=50,radius=100)

kriging(myDb2010,myGrid,myModel,myNeigh)

plotGrid(myGrid,varName="Kriging.density.estim",polyName = 'Polygon',
         coastLine = myCoastLine,title="Kriging density (g/m²)")

plotGrid(myGrid,varName="Kriging.density.stdev",polyName = 'Polygon',
         coastLine = myCoastLine,title="Standard deviation of the kriging density (g/m²)")

plotGrid(myGrid,varName="Kriging.density.stdev",polyName = 'Polygon',
         coastLine = myCoastLine,title="Standard deviation of the kriging density (g/m²)",legend=F)
points(myDb2010['longitude'],myDb2010['latitude'])

##########################
### Global estimation by kriging
##########################

# Cvv is discretised in a consistent manner wrt to the kriging grid 
# A high def grid is thus not recommended
# Users have to make a compromise between the numerical precision of the computation of Cvv and 
# the computing time.
# To this end, start with a very crude grid and refine it up to the convergence of the estimation std

# mode = CovCalcMode()
# mode$setAsVario(T)
# myModel$evalCvv(ext=myPoly,ndisc=c(10,10),angles=c(0,0),mode=mode)

myGridCrude = DbGrid_createCoveringDb(dbin=myDb2010,dx=c(20,20))
db_polygon(myGridCrude, myPoly)
res <- global_kriging(dbin=myDb2010,dbout=myGridCrude,model=myModel,ivar0=0,verbose=T)

geoEstimStdTest <- res$sse

myGridCrude = DbGrid_createCoveringDb(dbin=myDb2010,dx=c(10,10))
db_polygon(myGridCrude, myPoly)
res <- global_kriging(dbin=myDb2010,dbout=myGridCrude,model=myModel,ivar0=0,verbose=T)

geoEstimStdTest[2] <- res$sse

myGridCrude = DbGrid_createCoveringDb(dbin=myDb2010,dx=c(5,5))
db_polygon(myGridCrude, myPoly)
res <- global_kriging(dbin=myDb2010,dbout=myGridCrude,model=myModel,ivar0=0,verbose=T)

geoEstimStdTest[3] <- res$sse

myGridCrude = DbGrid_createCoveringDb(dbin=myDb2010,dx=c(2.5,2.5))
db_polygon(myGridCrude, myPoly)
res <- global_kriging(dbin=myDb2010,dbout=myGridCrude,model=myModel,ivar0=0,verbose=T)

geoEstimStdTest[4] <- res$sse

plot(c(20,10,5,2.5),geoEstimStdTest,type="b")
#  ==> convergence OK

geoEstimStd <- geoEstimStdTest[4]
geoEstim <- res$zest
geoCV <- res$cvgeo

### Computing CViid

db_polygon(myDb2010, myPoly)
myDb2010$getActiveSampleNumber()
# ==> N = 279

iidEstimVariance <- myDb2010$getVariance('density',useSel = TRUE)/myDb2010$getActiveSampleNumber()
iidEstimStd <- sqrt(iidEstimVariance)
iidEstim <- myDb2010$getMean('density',useSel = TRUE)
iidCV <- iidEstimStd/iidEstim

V <- myPoly$getSurface() # in square kilometers
iidQ <- iidEstim * V * 10^6 / 10^6 # x 10^6 to get the surface in square meters
# / 10^6 to get tonnes instead of grams

geoQ <- geoEstim * V * 10^6 / 10^6

# the biomass in the polygon is estimated 
#  34 757 tonnes by statistical method and
#  32 164 tonnes using geostatistics, i.e. 7.5% less.

# The CV is estimated
#  6.75% by statistical method and
#  4.63% by geostatistical approach

# The confidence intervals, under the assumption that the errors are Gaussian, 
# associated with a 5% risk, are 
#  [30 065 ; 39 448] tonnes by statistical method and
#  [29 187 ; 35 141] tonnes by geostatistical approach


plot(seq(25000, 50000,100),dnorm(x=seq(25000, 50000,100),mean = geoQ,sd = geoEstimStd*V),
     type='l',col=2,xlab = 'Biomass (tonnes)',ylab="",yaxs='i',ylim=c(0,0.0003))
lines(seq(25000, 50000,100),dnorm(x=seq(25000, 50000,100),mean= iidQ,sd=iidEstimStd*V),col=1)
legend("topright",legend=c("Statistical approach (i.i.d.)","Geostatistical approach"),text.col = c(1,2))

############################## FIN


