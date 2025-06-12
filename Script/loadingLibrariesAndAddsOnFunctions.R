######################
# Loading libraries
######################
# Installing (if needed) and opening required CRAN libraries
CRAN_package_list <- c('maps','mapdata','mapproj','wesanderson')
not_in_list <- CRAN_package_list[!(CRAN_package_list %in% installed.packages())]
lapply(not_in_list,install.packages,dependencies=TRUE)
#
lapply(CRAN_package_list, require, character.only = TRUE)
# clean
rm(CRAN_package_list,not_in_list)

# Installing the up-to-date gstlearn library
install.packages("gstlearn",repos="https://soft.mines-paristech.fr/cran")
library(gstlearn)

# Choosing a color palette
myPalette = wes_palette("Zissou1", 16, type = "continuous")

# default par
defPar <- par(no.readonly = TRUE)

######################
### Personnal graphical functions
# plotVario
# plotModel
# plotGrid
######################
plotVario <- function(vario,col=myPalette[c(1,6,13,16)],opacity="80",inches=0.15,legend=T){
  ### Plot univariate directional or omnidirectional variograms  
  # to get NON transparent colour set opacity="" 
  # Turn the palette color into rgb
  colRgb <- rgb(t(col2rgb(col)),maxColorValue = 255)
  ndir <- vario$getNDir() # 2025/05/26 vario$getDirectionNumber()
  dir <- NULL 
  for(i in 0:(ndir-1)){
    h <- vario$getAllHh(idir=i)
    g <- vario$getAllGg(idir=i)
    w <- vario$getAllSw(idir=i)
    symbols(h,g,circles=sqrt(w),inches=inches,bg=paste0(colRgb[i+1],opacity),add=i>0,
            xlim=c(0,1.1*vario$getHmax()),ylim=c(0,1.1*vario$getGmax()),xaxs="i",yaxs="i",las=1,
            xlab='Distances',ylab='Variogram',cex.lab=1.5,cex.axis=1.5)
    lines(h,g,col=colRgb[i+1])
    codir <- vario$getDirParam(i)$getCodirs()
    dir[i+1] <- atan(codir[2]/codir[1])*180/pi
  }
  abline(h=vario$getVar(ivar=0,jvar=0),lty=2,col=1,lwd=1.5)
  if(legend) {
    legend("bottomright",
           paste0(as.character(dir),"°"),text.col=col[1:ndir],
           lty=1,lwd=2,col=col[1:ndir],
           title='Direction',title.col=1)
  }
}

####################
plotVarioModel <- function(vario,model,col=myPalette[c(1,6,13,16)],asCov=F){
  ### Draw model on an existing empirical variogram
  mode = CovCalcMode()
  mode$setAsVario(!asCov)
  nDir <- vario$getNDir()
  plotVario(vario)
  hh <- seq(0,par("usr")[2],length=200)
  for(i in 0:(nDir-1)){
    codir <- vario$getDirParam(i)$getCodirs()
    gg <- model$sample(hh, codir=codir, mode=mode)
    if(model$getCovName(0) == "Nugget Effect") gg[1] <- model$getSill(0,0,0)
    lines(hh,gg,col=col[i+1],lwd=2)
  }
}

####################
plotModel <- function(model,hMax=1,dir=t(c(1,0)),col=NULL,asCov=F){
  ### Draw model 
  if(is.null(col)) col <- 1:NROW(dir)
  mode = CovCalcMode()
  mode$setAsVario(!asCov)
  nDir <- NROW(dir)
  plot(c(0,hMax),c(0,model$getTotalSill()),type="n",xaxs="i",yaxs="i",las=1,
            xlab='Distances',ylab='Variogram',cex.lab=1.5,cex.axis=1.5)
  hh <- seq(0,hMax,length=200)
  for(i in 1:nDir){
    codir <- dir[i,]
    gg <- model$sample(hh, codir=codir, mode=mode)
    if(model$getCovName(0) == "Nugget Effect") gg[1] <- model$getSill(0,0,0)
    lines(hh,gg,col=col[i],lwd=2)
  }
}


######################
plotGrid <- function(grid,varName,col=myPalette,polyName=NULL,title=NULL,coastLine=NULL,legend = T,legendTitle=NULL,probs=NULL,...){
  # Assumes that par(las=1)
  nx <- grid$getNX(0)
  z <- grid[varName]
  if(!is.null(probs)){
    z[z > quantile(z,probs=probs,na.rm=T)] <- quantile(z,probs=probs,na.rm=T)
  }
  if(legend){
    zLegend <- seq(min(z,na.rm=T),max(z,na.rm=T),length=length(col))
    layout(matrix(c(1,1,1,0,2,0),ncol=2,byrow=F), widths=c(10,1),heights=rep(5,3))
    par(mar=c(5,5,4,1)+0.1,cex=1)#.axis=1.5,cex.lab=1.5,cex.main=1.5)
    image(unique(grid['x1']), unique(grid['x2']),
          if(is.null(polyName)) matrix(z,nx) 
          else matrix(ifelse(grid['Polygon']==1,z,NA),nx),col=col,
          xlab="Longitude (km)",ylab="Latitude (km)",
          main=title,asp=1,...)
    lines(coastLine)
    # Draw the color legend
    par(mai=c(0.1,0.5,0.1,0.1))
    image(1, zLegend, t(zLegend), col= col, axes=FALSE, xlab="",ylab="",main=legendTitle)
    axis(2)
    par(defPar)
  } else {
    image(unique(grid['x1']), unique(grid['x2']),
          if(is.null(polyName)) matrix(z,nx) 
          else matrix(ifelse(grid['Polygon']==1,z,NA),nx),col=col,
          xlab="Longitude (km)",ylab="Latitude (km)",
          las=1,main=title,asp=1,...)
    lines(coastLine)
  }
}

######################
digitDist <- function(){
  #       Estim the distance between two points "                       
  #       pointed on the graph with the mouse."                          
  #       The unit of the distance is those of the graph."                                  
  cat("point the two points \n")
  toto <- locator(2)
  d <- sqrt(diff(toto$y)^2 + diff(toto$x)^2)
  cat(round(d, 2), "\n")
  d
}

