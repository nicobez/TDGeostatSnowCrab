install.packages("gstlearn",repos="https://soft.mines-paristech.fr/cran")
library(gstlearn)

myDb <- Db_create()
myDb$addColumns(c(0,1),'longitude',ELoc_X(),0)
myDb$addColumns(c(0,1),'latitude',ELoc_X(),1)
myDb$addColumns(c(3,7),'density',ELoc_Z(),0)

myGrid = DbGrid_create(nx=c(10,10))

### Pour fitter un vario il me faut un vario expérimental
### Pour cela je passe par une simu non cond d'un sphérique

mySimu = DbGrid_create(nx=c(50,50))
modelSpherical = Model_createFromParam(type = ECov_SPHERICAL(), range = 10)
simtub(NULL, mySimu, modelSpherical, nbtuba=1000, seed = 43431)

### Calcul du varioexp
varioParamOmni = VarioParam_createOmniDirection(nlag=30, dlag=1) 
varioexp = Vario(varioParamOmni)
varioexp$compute(mySimu)
plotVario(varioexp)

#####################################################
### cas#1 : UN SPHERIQUE
#####################################################
# Fit du varioexp 
myModelFitted = Model_create() #Model_createFromDb(myDb2010)
myModelFitted$fit(vario=varioexp, types=ECov_fromKeys(c("SPHERICAL")))
plotVarioModel(varioexp,myModelFitted)

# Définition à la main du même modèle
myModelDefined <- Model_createFromParam(type = ECov_fromKey("SPHERICAL"), 
                                        range = myModelFitted$getRange(0),
                                        sill=myModelFitted$getSill(0,0,0))

resFitted <- global_kriging(dbin=myDb,dbout=myGrid,model=myModelFitted,ivar0=0,verbose=T)
resDefined <- global_kriging(dbin=myDb,dbout=myGrid,model=myModelDefined,ivar0=0,verbose=T)
resFitted$weights
resDefined$weights

# ==> les deux situations sont OK mais les weights ne sont pas fournis

#####################################################
### cas#2 : UN LINEAR
#####################################################
# Fit du varioexp 
myModelFitted = Model_create() #Model_createFromDb(myDb2010)
myModelFitted$fit(vario=varioexp, types=ECov_fromKeys(c("LINEAR")))
plotVarioModel(varioexp,myModelFitted)

# Définition à la main du même modèle
myModelDefined <- Model_createFromParam(type = ECov_fromKey("LINEAR"), 
                                        range = myModelFitted$getRange(0),
                                        sill=myModelFitted$getSill(0,0,0))

# global kriging
resFitted <- global_kriging(dbin=myDb,dbout=myGrid,model=myModelFitted,ivar0=0,verbose=T)
resDefined <- global_kriging(dbin=myDb,dbout=myGrid,model=myModelDefined,ivar0=0,verbose=T)

# ==> dans les deux cas, on récupère zest mais sse=0

#####################################################
### cas#3 : ajout d'une structure (e.g. nugget)
#####################################################
mySimu['Simu'] <- mySimu['Simu'] + rnorm(mySimu$getNSample(),0,0.6)

### Calcul du varioexp
varioParamOmni = VarioParam_createOmniDirection(nlag=30, dlag=1) 
varioexp = Vario(varioParamOmni)
varioexp$compute(mySimu)
plotVario(varioexp)

# Fit du varioexp 
myModelFitted = Model_create() #Model_createFromDb(myDb2010)
myModelFitted$fit(vario=varioexp, types=ECov_fromKeys(c("NUGGET","SPHERICAL")))
plotVarioModel(varioexp,myModelFitted)

# Définition à la main du même modèle
myModelDefined <- Model_createFromParam(type = ECov_fromKey("NUGGET"), 
                                        range = myModelFitted$getRange(0),
                                        sill=myModelFitted$getSill(0,0,0))
myModelDefined$addCovFromParam(type=ECov_fromKey("SPHERICAL"),
                               range=myModelFitted$getRange(1),
                               sill=myModelFitted$getSill(1,0,0))

# global kriging
resFitted <- global_kriging(dbin=myDb,dbout=myGrid,model=myModelFitted,ivar0=0,verbose=T)
resDefined <- global_kriging(dbin=myDb,dbout=myGrid,model=myModelDefined,ivar0=0,verbose=T)

# ==> ok avec le modèle défini A LA MAIN MAIS PAS AVEC LE MODELE FITTED
resFitted$sse
resDefined$sse

# et pourtant estimation est ok dans les deux cas
resFitted$zest
resDefined$zest
