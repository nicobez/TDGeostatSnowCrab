plot(myDf$Lon,myDf$Lat,xlab="Easting",ylab="Norting")
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


# # Defining a relevant distance lag for computing the variogram
# # CAUTIOUS : the digisDist function must be applied when the zoom option of RStudio is 100%
# # See : Tools/Global Options/Appearance
# digitDist()

# Create an empty object containing computation's characteristics
myVarioParamOmni = VarioParam()
# Define computation's characteristics
temp = DirParam_create(nlag=55,dlag=10)
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
temp = DirParam_createMultiple(ndir=4, nlag=40, dlag=10)
myVarioParamDirectional$addMultiDirs(temp)
myVarioDirectional = Vario(myVarioParamDirectional)
myVarioDirectional$compute(myDb2010,ECalcVario_VARIOGRAM())
plotVario(myVarioDirectional)

### variogramme directionnel avec lag et nlag variables par direction
myVarioParamDirectional = VarioParam()
temp = DirParam_create(nlag=35,dlag=10,angle2D = 0,tolang=25)
myVarioParamDirectional$addDir(temp)
temp = DirParam_create(nlag=15,dlag=10,angle2D = 45,tolang=25)
myVarioParamDirectional$addDir(temp)
temp = DirParam_create(nlag=15,20,dlag=10,angle2D = 90,tolang=25)
myVarioParamDirectional$addDir(temp)
temp = DirParam_create(nlag=35,dlag=10,angle2D = 135,tolang=25)
myVarioParamDirectional$addDir(temp)

myVarioDirectional2 = Vario(myVarioParamDirectional)
myVarioDirectional2$compute(myDb2010,ECalcVario_VARIOGRAM())

plotVario(myVarioDirectional2,inches=0.1)

##########################
### model du vario 2010
##########################

# Printing all possible models
ECov_printAll()

# Automatic fit of a variogram model
myModel0 = Model_create() #Model_createFromDb(myDb2010)
myModel0$fit(vario=myVarioOmni, types=ECov_fromKeys(c("NUGGET","SPHERICAL","SPHERICAL")))
#optvar = Option_VarioFit(flag_noreduce = TRUE)) #mauto = Option_AutoFit$setWmode(wmode=2)
plotVarioModel(myVarioOmni,myModel0)

# Manual copy of the variogram model to overcome a small bug in global_kriging() 
# that will be fixed in the next GSTLEARN version
myModel <- Model_createFromParam(type = ECov_fromKey("NUGGET"), sill=myModel0$getSill(0,0,0))
myModel$addCovFromParam(type=ECov_fromKey("SPHERICAL"),
                         range=myModel0$getRange(1),
                         sill=myModel0$getSill(1,0,0))
myModel$addCovFromParam(type=ECov_fromKey("SPHERICAL"),
                         range=myModel0$getRange(2),
                         sill=myModel0$getSill(2,0,0))

# IMPORTANT : 
# Add a mean drift to the model
# This will filter out the mean when kriging (i.e. to get kriging weights that sum to one) 
myModel$addDrift(DriftM())

##########################
### Krigeage local = cartographie 
##########################

### kriging grid
myGrid = DbGrid_createCoveringDb(dbin=myDb2010,dx=c(1,1))
db_polygon(myGrid, myPoly)

# choosing the neighbourhood
myNeigh = NeighMoving_create(flag_xvalid=FALSE,nmaxi=50,radius=100)

kriging(myDb2010,myGrid,myModel,myNeigh)

plotGrid(myGrid,varName="Kriging.density.estim",polyName = 'Polygon',
         coastLine = myCoastLine,title="Kriging density (g/m²)")

plotGrid(myGrid,varName="Kriging.density.stdev",polyName = 'Polygon',
         coastLine = myCoastLine,title="Standard deviation of the kriging density (g/m²)")

plotGrid(myGrid,varName="Kriging.density.stdev",polyName = 'Polygon',
         coastLine = myCoastLine,title="Standard deviation of the kriging density (g/m²)",legend=F)
points(myDb2010['longitude'],myDb2010['latitude'])
lines(myPoly$getX(0),myPoly$getY(0))

##########################
### Global estimation by kriging
##########################

# Civ is discretised in a consistent manner wrt to the kriging grid 
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

kriEstimStdTest <- res$sse

myGridCrude = DbGrid_createCoveringDb(dbin=myDb2010,dx=c(10,10))
db_polygon(myGridCrude, myPoly)
res <- global_kriging(dbin=myDb2010,dbout=myGridCrude,model=myModel,ivar0=0,verbose=T)

kriEstimStdTest[2] <- res$sse

myGridCrude = DbGrid_createCoveringDb(dbin=myDb2010,dx=c(5,5))
db_polygon(myGridCrude, myPoly)
res <- global_kriging(dbin=myDb2010,dbout=myGridCrude,model=myModel,ivar0=0,verbose=T)

kriEstimStdTest[3] <- res$sse

myGridCrude = DbGrid_createCoveringDb(dbin=myDb2010,dx=c(2.5,2.5))
db_polygon(myGridCrude, myPoly)
res <- global_kriging(dbin=myDb2010,dbout=myGridCrude,model=myModel,ivar0=0,verbose=T)

kriEstimStdTest[4] <- res$sse

par(mai=c(1,1,0.2,0.2))
plot(c(20,10,5,2.5),kriEstimStdTest,xlab="Pixels' size of the discretization grid (geographical units)",ylab="")
mtext("Global kriging standard deviation",2,line=4,las=0)
arrows(x0=c(20,10,5),y0=kriEstimStdTest[-4],x1=c(10,5,2.5),y1=kriEstimStdTest[-1])
par(defPar)
#  ==> convergence OK for a 2.5 x 2.5 grid

kriEstimStd <- res$sse
kriEstim <- res$zest
kriCV <- res$cvgeo

V <- myPoly$getSurface() # in square kilometers i.e. the units of the coordinates


# Computing the biomass over the polygon
kriQ <- kriEstim * V * 10^6 / 10^6  # x 10^6 to get the surface in square meters
                                    # / 10^6 to get tonnes instead of grams
kriQ

# Confidence interval 
# Under the hypothesis of a Gaussian distribution

qMin <- kriQ*0.8 ; qMax <- kriQ*1.2
tempX <- seq(qMin, qMax, length=100)
tempY <- dnorm(x=tempX,mean = kriQ,sd = kriEstimStd*V)
plot(tempX,tempY,
     type='l',col=2,xlab = 'Biomass (tonnes)',ylab="",yaxs='i',ylim=c(0,1.2*max(tempY)),
     las=1,main="Probability distribution of the biomass estimator")
q5 <- qnorm(0.05,mean = kriQ,sd = kriEstimStd*V)
q95 <- qnorm(0.95,mean = kriQ,sd = kriEstimStd*V)

segments(q5,0,q5,dnorm(q5,mean = kriQ,sd = kriEstimStd*V),lty=2)
segments(q95,0,q95,dnorm(q95,mean = kriQ,sd = kriEstimStd*V),lty=2)
segments(kriQ,0,kriQ,dnorm(kriQ,mean = kriQ,sd = kriEstimStd*V))
text(kriQ,1.05*max(tempY),paste0((round(kriQ,0))," tonnes"))
text(q5,1.05*dnorm(x=q5,mean = kriQ,sd = kriEstimStd*V),paste0((round(q5,0))," tonnes"),pos=2)
text(q95,1.05*dnorm(x=q95,mean = kriQ,sd = kriEstimStd*V),paste0((round(q95,0))," tonnes"),pos=4)


# the estimated biomass in the polygon is kriQ tonnes.
# The confidence interval, under the assumption that the errors are Gaussian, 
# associated with a 10% risk, is [q5 tonnes ; q95 tonnes]  

############################## FIN


