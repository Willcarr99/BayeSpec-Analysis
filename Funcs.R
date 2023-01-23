## Defined Functions
##---------------------------------------------------------------

## Help calls another file
h <- function()source("Help.R")

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Let the log scale
logy <- function(){
  
  if(logscale=="")
    logscale<<-"y"
  else
    logscale<<-""

  replot()
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Load the dat file
load.dat <- function(name="spec.dat"){
  filename <<- name
  spec <<- read.table(name,header=FALSE,skip=1)
  Energy <<- spec[,1] ##1:length(spec[,1]) 
  dEnergy <<- rep(0,length.out=length(Energy))
  spec[,1] <<- 1:length(spec[,1])
  Energybkup <<- Energy
  spec[spec[,2]<0,2] <<- 0
  draw.spec()
##  if(exists("RefLines"))rm(RefLines)

  ## Load the Field
  fieldfile <<- sub(".dat",".field",name,fixed=TRUE)
  if(file.exists(fieldfile))
    Field <<- scan(fieldfile,nmax=1,quiet=TRUE)

  
  ## Counts in the zero channel
  ##cat("There are ",spec[1,2]," counts in 0 out of ",sum(spec[-1,2]),"\n",
  ##    "So correct peaks by factor: ",1+spec[1,2]/sum(spec[-1,2]),
  ##    "\n\n",sep="")
}
q
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Quit
quit <- function(){
  tkdestroy(tt)
  par(oldpar)
  dev.off()
  rm(list=ls(pos=1),inherits=TRUE)
  source("~/.Rprofile")
  source("Spec.R")
  cat("Done cleaning up!\n")
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Open X11 window and draw the spectrum
draw.spec <- function(){

  ## make real histrogram step output
  ##  make.steps()
  
  ## Define the maxes and mins
  orig.spec <<- spec
  
  xrange <<- c(min(spec[,1]),max(spec[,1]))
  yrange <<- c(0,max(spec[,2])*1.05)

  X11(width=11,height=5,xpos=0,ypos=0)
  oldpar <<- par(newpar)
  doplot(xrange,yrange)
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Print the spectrum to a file
print.spec <- function(){

  ## Enter the filename
  filename <- scan(nlines=1,quiet=TRUE,what=character())

  outputpdf(filename,width=7,height=4)

}
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot
doplot <- function(xrange,yrange){

  if(logscale=="y" && yrange[1]==0)yrange[1]=0.1

  tspec <- spec
  ##tspec[tspec[,2]<1,2] <- 0.01
  
  plot(tspec,type="s",
       xlab="BRho (T.m)",ylab="Counts",xlim=xrange,
       ylim=yrange,xaxt="n",yaxs="i",log=logscale,
       main=filename)

## If the comparison spectrum exists, plot it
  if(exists("c.spec")){
    lines(c.spec,type="s",lty=1,col="darkred")
  }

  
  xat <- axTicks(3)+0.5
  ind <- sapply(axTicks(3),function(x)which.min(abs(x-spec[,1])))
  ylabels <- spec[ind,1]
  if(axTicks(3)[1]==0){
    xat <- xat[2:length(xat)]
    ylabels <- ylabels[2:length(ylabels)]
  }
  axis(3,at=xat,labels=ylabels)

  ERange <- Energy[xrange]
  xlabels <- pretty(ERange)
  xlabels <- xlabels[2:(length(xlabels)-1)]
  ind <- sapply(xlabels,function(x)which.min(abs(x-Energy)))
  axis(1,at=ind,labels=xlabels)
  
  peakfind()

  if(exists("RefLines")){
    abline(v=RefLines,lty=2,col="red")
  }

  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Get a point
getpoint <- function(){

  while(!is.null(test <- locator(1))){
    cat(Energy[floor(test$x)]," +- ",dEnergy[floor(test$x)],"\n")
  }
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## zoom to whole spectrum
za <- function(){

  ## Define the maxes and mins
  xrange <<- c(min(spec[,1]),max(spec[,1]))
  yrange <<- c(0,max(spec[,2])*1.05)

  doplot(xrange,yrange)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## zoom in x direction (click to zoom)
zx <- function(){

  ## Get user to point on plot
  cat("Click on region to zoom\n")
  positions <- locator(2)
  positions$x[1] <- max(positions$x[1],1)
  
  ## Need to floor the clicks to make integers
  ## for now, work with uncalibrated
  xrange <<- c(floor(positions$x[1]),floor(positions$x[2]))
  yrange <<- c(0,max(spec[xrange[1]:xrange[2],2]*1.05))

  doplot(xrange,yrange)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Zoom out in the x-range
zo <- function(){

  xmin <- floor(max(spec[1,1],xrange[1]-diff(xrange)/4))
  xmax <- floor(min(spec[length(spec[,1]),1],xrange[2]+diff(xrange)/4))
  xrange <<- c(xmin,xmax)
  yrange <<- c(0,max(spec[xrange[1]:xrange[2],2]*1.05))
  
  doplot(xrange,yrange)
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## zoom in y direction (click to zoom)
zy <- function(){

  ## Get user to point on plot
  cat("Click on region to zoom\n")
  positions <- locator(2)

  ## Need to floor the clicks to make integers
  yrange <<- c(max(0, floor(min(positions$y))),floor(max(positions$y)))

  doplot(xrange,yrange)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# auto scale y-direction
ay <- function(){

  yrange <<- c(0,max(spec[xrange[1]:xrange[2],2]*1.05))
  doplot(xrange,yrange)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate gross area of a peak
gross <- function(){

  ## Call functin to select a single area
  region <- select.single()
  E <- region[,3]
  chn <- region[,1]
  counts <- region[,2]

  ## Find the area
  area <- sum(counts)
  area.unc <- sqrt(area)

  cat("\n--------------------------\n")
  cat("From E=",E[1],"keV to ",E[length(E)],"keV\n",sep="")
  pos <- peakpos(E,counts)
  cat("Peak at E = ",pos[1]," +- ",pos[2]," keV\n",sep="")
  cat("Chan  =",chn[1]," to ",chn[length(chn)],"\n",sep="")
  cat("Gross Area = ",area,"+-",area.unc,"\n")

}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## replot current view (remove other decorations)
replot <- function(){
  doplot(xrange,yrange)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## replot current view in inverted colors
invert <- function(){
  par(list(bg="white",fg="black",col.axis="black",
               col.lab="black",col.main="black",
               col.sub="black"))
  replot()
  par(newpar)
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Rebin the data (save the old data)
rebin <- function(){

  cat("Enter rebin factor")
  rfactor <<- scan(nlines=1,quiet=TRUE)

  if(rfactor < 0){
    spec <<- orig.spec
    if(exists("c.spec")) c.spec <<- orig.c.spec
  } else {
    ##spec[,2] <<- filter(spec[,2],rep(1/rfactor,rfactor),sides=2)
    for(i in seq(from=1,to=(length(spec[,1])-1),by=rfactor)){
      j <- seq(from=i,to=i+rfactor-1)
      spec[j,2] <- sum(spec[j,2])/rfactor
      if(exists("c.spec"))c.spec[j,2] <- sum(c.spec[j,2])/rfactor
##      spec[i+1,2] <- spec[i,2]
    }
    spec <<- spec
    if(exists("c.spec"))c.spec <<- c.spec
  }
  replot()

}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Overlay another spectrum
overlay <- function(name){

  if(name != ""){
    c.filename <<- name
    c.spec <<- read.table(c.filename,header=FALSE)
    c.Energy <<- spec[,1]
    c.spec[,1] <- 1:length(spec[,1])
    c.spec[c.spec[,2]<0,2] <<- 0
    orig.c.spec <<- c.spec
    replot()
  } else {
    if(exists("c.spec"))rm(c.spec,orig.c.spec)
    replot()
  }

}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calibrate the spectrum using a linear calibration
cal <- function(){

  ## Get the calibration coefficients
  cat("Enter Peak 1 (Energy, Channel)")
  peak1 <- scan(nlines=1,quiet=TRUE)
  cat("Enter Peak 2 (Energy, Channel)")
  peak2 <- scan(nlines=1,quiet=TRUE)

  E <- c(peak1[1],peak2[1])
  chn <- c(peak1[2],peak2[2])

  ## Do the fit with a simple linear model
  fit <- lm(E~chn)

  ## The calibration coefficients
  Cal <<- as.double(fit$coefficients)

  ## Change the Energy
  Energy <<- spec[,1]*Cal[2] + Cal[1]

  replot()
  
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate Net Area
net <- function(){

  ## Select the side bands and peak
  region <- select.triple()
  chn <- region[,1]
  peak.counts <- region[,2]
  bg.counts <- region[,3]
  bg.counts.unc <- region[,4]
  E <- region[,5]

  ## Calculate where the peak is
  pos <- peakpos(E,peak.counts-bg.counts)
  pos[2] <- pos[2]*diff(Energy)[2]
  chn.pos <- which.min(abs(pos[1]-Energy))
  
  ## Calculate the posterior and get important points
  peak.sum <- sum(peak.counts)
  bg.sum <- sum(bg.counts)
  bg.sum.unc <- sqrt(bg.sum+sum(bg.counts.unc^2))

  points <- upperlim(peak.sum,bg.sum,bg.sum.unc)

  cat("\n--------------------------\n",
      "Peak at E = ",pos[1]," +- ",sqrt(pos[2]^2+dEnergy[chn.pos]^2),
      " keV\n",sep="")
  
  cat(points,"\n")
  cat("Net Counts = ",points[3]," +- ",(points[4]-points[2])/2,"\n")
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Calculate upper limit
ul <- function(){
  region <- select.triple()
  chn <- region[,1]
  peak.counts <- region[,2]
  bg.counts <- region[,3]
  bg.counts.unc <- region[,4]
  E <- region[,5]

  ## Calculate the posterior and get important points
  peak.sum <- sum(peak.counts)
  bg.sum <- sum(bg.counts)
  bg.sum.unc <- sqrt(bg.sum+sum(bg.counts.unc^2))

  points <- upperlim(peak.sum,bg.sum,bg.sum.unc)
  
  cat(100*ConfUL,"% Upper limit = ",points[5],"\n\n",sep="")
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
upperlim <- function(peak.sum,bg.sum,bg.sum.unc){
  
  ## The prior to use
  m <- BayesNorm
  Prior <- function(s,b)1/(s+b)^m

  ## Background probability density function
  bg.prob <- function(b)exp(-(bg.sum-b)^2/(2*bg.sum.unc^2))
  
  ## liklihood function use Poisson for low counts
  if(peak.sum<100){
    cat("Using Poisson Statistics\n")
    liklihood <- function(s,b){
      (exp(-(s+b))*(s+b)^peak.sum/factorial(peak.sum))*bg.prob(b)
    }
  }else{
    liklihood <- function(s,b){
      (1/sqrt(peak.sum))*exp(-(s+b-peak.sum)^2/(2*peak.sum))*bg.prob(b)
    }
  }

  ## posterior
  post <- function(s){

    ## Need to integrate out the background dependence
    lowerlim <- max(0,bg.sum-4*bg.sum.unc)
    upperlim <- bg.sum+4*bg.sum.unc
    integrand <- function(b,signal){
      liklihood(signal,b)*Prior(signal,b)
    }
    ## Integrate over b for each s
    sapply(s,function(s){
      integrate(integrand,lowerlim,upperlim,signal=s)$value
    })
    
  }

  ## Calculate the posterior from 0 to 5*counts
  ## limits to calculate posterior over
  ## can be a peak, or an upper limit
  lowerlim <- max(0,abs(peak.sum-bg.sum)-3*sqrt(abs(peak.sum-bg.sum))-3*bg.sum.unc)
  upperlim <- abs(peak.sum-bg.sum)+3*sqrt(abs(peak.sum-bg.sum))+3*bg.sum.unc
    ##2*peak.sum ##6*sqrt(2*peak.sum)
  grid <- seq(from=lowerlim,to=upperlim,length.out=1000)
  posterior <- post(grid)
  posterior <- posterior/sum(posterior)

  ## Now we need the cumulative posterior
  cumpost <- cumsum(posterior)
  ##plot(grid,cumpost,type="l")
  
  ## Make a spline to interpolate
  interp <- splinefun(cumpost,grid)

  points <- interp(c(0,0.16,0.5,0.84,0.9))

  return(points)
  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Find the positions of a peak
peakpos <- function(E,C){

  Centroid <- sum(E*C)/sum(C)
  Unc <- sum(E*C)/(sum((E-Centroid)^2)*sum(C))
  Centroid <- c(Centroid,sqrt(Unc)*(E[2]-E[1]))

  #*diff(range(E)) 
  return(Centroid)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Go to specified energy
goto <- function(){

  cat("\nEnter Energy: ")
  E <- scan(nlines=1,quiet=TRUE)

  point <- findInterval(E-1:E+1,Energy)

  xrange <<- c(point[1]-100,point[1]+100)
  yrange <<- c(0,max(spec[xrange[1]:xrange[2],2]*1.05))
  doplot(xrange,yrange)

}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Save the calibration to a file
savecal <- function(){

  ## replace .dat or .txt with .cal to save calibration file
  calfile <- sub(".dat",".cal",filename)
  ##calfile <- sub(".txt",".cal",filename)

  ## Write the calibration to that file
  write.table(Cal,file=calfile,sep="\t",row.names=FALSE,col.names=FALSE)

}
load.cal <- function(){

  ## replace .dat or .txt with .cal to load calibration file
  calfile <- sub(".dat",".cal",filename,fixed=TRUE)
  ##calfile <- sub(".txt",".cal",filename)
  
  ## Read calibration from calibration file
  Cal <- read.table(calfile,header=FALSE,skip=1)
  
  ## Change the Energy
##  Energy <<- spec[,1]*Cal[2] + Cal[1]
  Energy <<- (Cal[3,1]*spec[,1]^2 + Cal[2,1]*spec[,1] + Cal[1,1])
  dEnergy <- sqrt(Cal[3,2]^2*spec[,1]^4 + Cal[2,2]^2*spec[,1]^2
                       + Cal[1,2]^2)
  
  replot()
}



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Find peaks in spectrum
peakfind <- function(){

  ## make a pts point moving average of the spectrum
  pts <- 3
  mvavg <- filter(spec[,2], rep(1, pts)/pts)

  ##lines(1:length(mvavg),mvavg,col="red")
  
  ## Arrays of first and second differentials
  ##diff1 <- spec[2:length(spec[,1]),2]-spec[1:(length(spec[,1])-1),2]
  ##diff2 <- diff1[2:length(diff1)]-diff1[1:(length(diff1)-1)]
  diff1 <- mvavg[2:length(mvavg)]-mvavg[1:(length(mvavg)-1)]
  diff2 <- (diff1[2:length(diff1)]-diff1[1:(length(diff1)-1)])

  ##  lines(spec[4:length(spec[,2]),1],diff2,col="red")
  
  ## Easy, the peak position is when the second differential is large and negative
  peakpos <- spec[diff2 < -peakSensitivity,1]+1

  ## For multiple channel peaks, find the beginning and end
  startpeak <- spec[peakpos[diff2[peakpos-1] > -peakSensitivity],1]+1
  endpeak <- spec[peakpos[diff2[peakpos+1] > -peakSensitivity],1]+1
  peakpos <- 2+(startpeak+endpeak)/2
  
  ##cat("Peaks at chn =\n")
  ##cat(c(spec[peakpos,1]),"\n")

  text(peakpos,spec[peakpos,2],
       labels=format(Energy[peakpos],digits=2,nsmall=0),
       pos=3,col="red")
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Select a single region
select.single <- function(draw=TRUE){

  ## Get user to point on plot
  cat("Click on peak to integrate\n")

  positions <- locator(2)
  x <- c(ceiling(positions$x[1]),ceiling(positions$x[2]))-1

  ## Plot the region stelected
  xx <- c(unlist(lapply(spec[x[1]:(x[2]+1),1],rep,2)))
  xx <- c(xx[2:(length(xx)-1)],rev(xx[2:(length(xx)-1)]))
  yy <- c(unlist(lapply(spec[x[1]:(x[2]),2],rep,2)),
          rep(0.01,2*(1+x[2]-x[1])))
  yy[yy==0] <- 0.01
  if(draw)polygon(xx, yy, col="darkred", border=FALSE)
  
  ## Return the index, channel, and counts
  selected.region <- cbind(spec[x[1]:x[2],],Energy[x[1]:x[2]])

  return(selected.region)
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fit background using two regions
select.triple <- function(){

  cat("Click first background region\n")
  bg.xpositions <- locator(2)$x

  ## Plot the region stelected
  x <- c(ceiling(bg.xpositions[1]),ceiling(bg.xpositions[2]))-1
  xx <- c(unlist(lapply(spec[x[1]:(x[2]+1),1],rep,2)))
  xx <- c(xx[2:(length(xx)-1)],rev(xx[2:(length(xx)-1)]))
  yy <- c(unlist(lapply(spec[x[1]:(x[2]),2],rep,2)),
          rep(0.01,2*(1+x[2]-x[1])))
  yy[yy==0] <- 0.01
  polygon(xx, yy, col="darkgreen", border=FALSE)


  ## Get the second peak
  cat("Click second background region\n")
  bg.xpositions <- c(bg.xpositions,locator(2)$x)

  ## Plot the region stelected
  x <- c(ceiling(bg.xpositions[3]),ceiling(bg.xpositions[4]))-1
  xx <- c(unlist(lapply(spec[x[1]:(x[2]+1),1],rep,2)))
  xx <- c(xx[2:(length(xx)-1)],rev(xx[2:(length(xx)-1)]))
  yy <- c(unlist(lapply(spec[x[1]:(x[2]),2],rep,2)),
          rep(0.01,2*(1+x[2]-x[1])))
  yy[yy==0] <- 0.01
  polygon(xx, yy, col="darkgreen", border=FALSE)

  
  ## round
  bg.x <- ceiling(bg.xpositions)-1
  bg.x.sel <- bg.x

  ## Draw them
##  lines(x=spec[bg.x[1]:(bg.x[2]+1),1],y=spec[bg.x[1]:(bg.x[2]+1),2],
##        type="h",col="green",lwd=2)
##  lines(x=spec[bg.x[1]:(bg.x[2]+1),1],y=spec[bg.x[1]:(bg.x[2]+1),2],
##        type="s",col="green",lwd=2)
##  lines(x=spec[bg.x[3]:(bg.x[4]+1),1],y=spec[bg.x[3]:(bg.x[4]+1),2],
##        type="h",col="green",lwd=2)
##  lines(x=spec[bg.x[3]:(bg.x[4]+1),1],y=spec[bg.x[3]:(bg.x[4]+1),2],
##        type="s",col="green",lwd=2)

  ## Get the two ranges, counts
  bg.x<-c(bg.x[1]:bg.x[2],bg.x[3]:bg.x[4])
  bg.counts <- spec[bg.x,2]

  ## Fit with a generalised linear model, with poisson statistics
  ##bg.fit <<- glm(bg.counts~bg.x,family=poisson("identity"))
  y2 <- mean(bg.counts[bg.x>=bg.x.sel[3]])
  y1 <- mean(bg.counts[bg.x<=bg.x.sel[2]])
  x2 <- mean(bg.x.sel[3:4])
  x1 <- mean(bg.x.sel[2:1])
  slope <- (y2-y1)/(x2-x1)
  intercept <- y1 - x1*slope
  bg.fit <<- glm(bg.counts~bg.x,family=poisson("identity"),
                 start=c(intercept,slope))
  
  ## Get the fit parameters and uncertainties
  bg.par <- as.double(bg.fit$coefficients)
  bg.par.unc <- summary(bg.fit)$coefficients[3:4]
  
  ## Now, draw the background line
  linex <- c(bg.x[1]+0.5,bg.x[length(bg.x)]+0.5)
  liney <- predict(bg.fit,data.frame(bg.x=linex))
  lines(x=linex,y=liney,col="green",lwd=2)

  ## ----------
  ## Next, get the peak area
  cat("Click peak region\n")
  peak.xpositions <- locator(2)$x
  peak.x <- ceiling(peak.xpositions)-1

  ## Draw it
##  lines(x=spec[peak.x[1]:(peak.x[2]+1),1],y=spec[peak.x[1]:(peak.x[2]+1),2],
##        type="h",col="red",lwd=2)
##  lines(x=spec[peak.x[1]:(peak.x[2]+1),1],y=spec[peak.x[1]:(peak.x[2]+1),2],
##        type="s",col="red",lwd=2)
  x <- peak.x
  xx <- c(unlist(lapply(spec[x[1]:(x[2]+1),1],rep,2)))
  xx <- c(xx[2:(length(xx)-1)],rev(xx[2:(length(xx)-1)]))
  yy <- c(unlist(lapply(spec[x[1]:(x[2]),2],rep,2)),
          rep(0.01,2*(1+x[2]-x[1])))
  yy[yy==0] <- 0.01
  polygon(xx, yy, col="darkred", border=FALSE)
  lines(x=linex,y=liney,col="green",lwd=2)

  
  ## Get the range and counts
  peak.x <- peak.x[1]:peak.x[2]
  peak.counts <- spec[peak.x,2]

  ## Finally, the summed total(s+b) and bkgd (b)
  bg.counts <- predict(bg.fit,data.frame(bg.x=peak.x))
  bg.counts.unc <- predict(bg.fit,data.frame(bg.x=peak.x),se.fit=TRUE)$se.fit
  
  selected.region <- cbind(spec[peak.x,],bg.counts,bg.counts.unc,Energy[peak.x])

  return(selected.region)
}

