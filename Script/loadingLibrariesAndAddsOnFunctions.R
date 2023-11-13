######################
# Loading libraries
######################
install.packages("gstlearn",repos="https://soft.mines-paristech.fr/cran")
library(gstlearn)
library(maps)
library(mapdata)
library(mapproj)
library(wesanderson) ; myPalette = wes_palette("Zissou1", 16, type = "continuous")
library(ggplot2)

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
  ndir <- vario$getDirectionNumber()
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
plotModel <- function(vario,model,col=myPalette[c(1,6,13,16)],opacity="80",asCov=F){
  ### Draw model on an existing empirical variogram
  mode = CovCalcMode()
  mode$setAsVario(!asCov)
  ndir <- vario$getDirectionNumber()
  plotVario(vario)
  hh <- seq(0,par("usr")[2],length=200)
  for(i in 0:(ndir-1)){
    codir <- vario$getDirParam(i)$getCodirs()
    gg <- model$sample(hh, codir=codir, mode=mode)
    if(model$getCovName(0) == "Nugget Effect") gg[1] <- model$getSill(0,0,0)
    lines(hh,gg,col=col[i+1],lwd=2)
  }
}



######################
plotGrid <- function(grid,varName,col=myPalette,polyName=NULL,title=NULL,coastLine=NULL,legend = T,legendTitle=NULL,...){
  # Assumes that par(las=1)
  nx <- grid$getNX(0)
  z <- grid[varName]
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

