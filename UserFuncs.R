## User Defined Functions
## - use this file to define your own functions
##---------------------------------------------------------------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MCMCGausOriginal <- function(npeaks=2){
    npeaks <<- npeaks
    source("MCMCGaus_original.R")
}

MCMCGaus <- function(npeaks=2){
    npeaks <<- npeaks
    source("MCMCGaus.R")
}

MCMCExpModGaus <- function(npeaks=1,rightTail=FALSE){
    npeaks <<- npeaks
    rightTail <<- rightTail
    source("MCMCExpModGaus.R")
}

MCMCExpModGausConstrained <- function(npeaks=1,rightTail=FALSE){
    npeaks <<- npeaks
    rightTail <<- rightTail
    source("MCMCExpModGausConstrained.R")
}

MCMCLogNorm <- function(npeaks=1){
    npeaks <<- npeaks
    source("MCMCLogNorm.R")
}

MCMCExpModGausAndGaus <- function(npeaks=2,nGaus=1,nEMG=1,rightTail=FALSE){
    npeaks <<- npeaks
    nGaus <<- nGaus
    nEMG <<- nEMG
    rightTail <<- rightTail
    #types <<- types
    source("MCMCExpModGausAndGaus.R")
}

MCMCLogNormAndGaus <- function(npeaks=2,nGaus=1,nLogNorm=1){
    npeaks <<- npeaks
    nGaus <<- nGaus
    nLogNorm <<- nLogNorm
    #types <<- types
    source("MCMCLogNormAndGaus.R")
}

## More complicated calibration that uses fits to Brho from elsewhere
load.cal <- function(){

  ## replace .dat or .txt with .cal to load calibration file
##  calfile <- sub(".dat",".cal",filename,fixed=TRUE)

  ## Read the correction fit
#  load(calfile)
#  newdata <- data.frame(c.chn=21*spec[,1]+1847.5)
#  fitrho <- predict(Fullfit,newdata=newdata,
#                    interval="confidence",se.fit=TRUE)$fit
#  uncert <- apply(fitrho,2,function(x)x/fitrho[,1])
#
#  ## Get Brho matrix that includes calibration uncertainties
#  Brho <- uncert*Energy

  if(file.exists(fieldfile)){
    load(calfile)
    ## remember to scale Brho by B for this spectrum
    newdata <- data.frame(f.rho=Energy/Field)
    Brho.new <- Energy*predict(fitcorr,newdata=newdata,
                               interval="confidence",se.fit=TRUE)$fit
    uncert <- apply(Brho.new,2,function(x)x/Brho.new[,1])
    ##  uncert <- cbind(rep(1,length(Energy)),rep(1,length(Energy)),
    ##                  rep(1,length(Energy)))
    Brho <- uncert*Brho.new
  } else {
    Brho <- cbind(Energy,Energy,Energy)
  }
  
  #fileName <- tclvalue(tkgetOpenFile(filetypes=
  #                                   "{{BrhotoEx file} {.dat}}"))
  ## So far, we've only calibrated the focal plane. Now we need to
  ## load the 22Ne(6Li,d) Brho data to find energies
  BrhotoEx <- function(x){

    ## figure out the BrhotoEx file name
    tmp <- unlist(strsplit(x=filename,split="data/"))
    dir <- tmp[1]
    el <- unlist(strsplit(x=tmp[2],split="-"))[1]
    el <- unlist(strsplit(x=el,split="[1-9]"))
    el <- el[length(el)]
    fileName <- paste(dir,"analysis/MasterCal/BrhotoEx-",el,".dat",sep="")
    
    data <- read.table(fileName,skip=5)
    ## Cut out and use only the angle we care about
    tmp <- unlist(strsplit(x=filename,split="deg"))[1]
    tmp <- substr(tmp,nchar(tmp)-1,nchar(tmp))
    angle <- abs(as.double(tmp))+0.6555
    data <- data[data[,3]==angle,]
    
    data.front <- data[seq(from=1,to=length(data[,1]),by=2),]
    data.back <- data[seq(from=2,to=length(data[,1]),by=2),]
    
    ## make splines to interpolate between Brho values at front and back
    spl.front <- splinefun(data.front[,11]/1000,data.front[,2])
    spl.back <- splinefun(data.back[,11]/1000,data.back[,2])

    ## Finally, from a Brho value, interpolate to find average energy
    res <- (spl.front(x)+spl.back(x))/2
    return(res)
    
  }

  Energy <<- 1000*BrhotoEx(Brho[,1])
  dEnergy <<- sqrt(1000*rowMeans(abs(sapply(2:3,function(x)
                                            BrhotoEx(Brho[,1])-BrhotoEx(Brho[,x])))))

  replot()
  


}

uncal <- function(){
  Energy <<- Energybkup
  replot()
}

## Fit (and plot) a gaussian
net <- function(){

  ## Select the side bands and peak
  region <- select.triple()
  chn <- region[,1]
  peak.counts <- region[,2]
  bg.counts <- region[,3]
  bg.counts.unc <- region[,4]
  E <- region[,5]

  ## The background corrected counts
  C <- peak.counts-bg.counts

  ## Calculate where the peak is
  pos <- peakpos(E,peak.counts-bg.counts)
  chn.pos <- which.min(abs(pos[1]-Energy))
  
  ## Calculate the posterior and get important points
  peak.sum <- sum(peak.counts)
  bg.sum <- sum(bg.counts)
  bg.sum.unc <- sqrt(bg.sum+sum(bg.counts.unc^2))

  points <- upperlim(peak.sum,bg.sum,bg.sum.unc)

  cat("\n--------------------------\n",
      "--- Using binned data  ---\n",
      " Peak at E = ",pos[1]," +- ",sqrt(pos[2]^2+dEnergy[chn.pos]^2),
      " keV\n",sep="")
  
  ##  cat(points,"\n")
  cat(" Net Counts = ",points[3]," +- ",(points[4]-points[2])/2,"\n",
      "--------------------------\n",sep="")

  ## Now calculate a Gausian
  ## Find weighted mean and variance
  mu <- weighted.mean(chn,w=C)
  sigma <- sqrt(sum(C*(chn-mu)^2)/sum(C))
  area <- points[3]

  ## Construct the Gaussian
  x <- seq(from=chn[1]-5,to=chn[length(chn)]+5,length.out=100)
  gaus <- function(x,mu,sigma){
    exp(-0.5*((x-mu)/sigma)^2)/(sigma*sqrt(2*pi))
  }
  y <- gaus(x,mu,sigma)
  y <- y*area
  ## add the background
  y <- y+predict(bg.fit,newdata=data.frame(bg.x=x))
#  lines(x,y,col="blue",lwd=2)

  ## calculate and print chisq
  chifn <- function(par){
    area <- par[1]
    mu <- par[2]
    sigma <- par[3]

    y <- predict(bg.fit,newdata=data.frame(bg.x=chn))+area*gaus(chn,mu,sigma)
    chisq <- sum((peak.counts-y)^2/peak.counts)
    return(chisq)
  }
  chisq <- chifn(c(area,mu,sigma))
  cat("\nStarting chisq = ",chisq,"\n",sep="")
  pars <- c(area,mu,sigma)
  fit <- optim(pars,chifn,hessian=TRUE)
  cat("\nFinishing chisq = ",fit$value,"\n",sep="")
  
  ## Print weighted mean and variance
  pars.fitted <- fit$par
  std.err <- sqrt(diag(solve(fit$hessian)))/pars.fitted

  mu <- approx(x=chn,y=E,xout=pars.fitted[2])$y
  mu.unc <- mu*std.err[2]
  sigma <- pars.fitted[3]*(E[2]-E[1])
  sigma.unc <- sigma*std.err[3]
  area <- pars.fitted[1]
  area.unc <- area*std.err[1]
  cat("\n--------------------------\n",
      "---  Using Gaussian fit  ---\n",
      "Peak at E = ",mu," +- ",mu.unc," with FWHM of ",sigma*2*sqrt(2*log(2)),
      " +- ",sigma.unc," keV \n",
      "Area = ",area," +- ",area.unc,"\n",sep="")

  y <- gaus(x,pars.fitted[2],pars.fitted[3])
  y <- y*area
  ## add the background
  y <- y+predict(bg.fit,newdata=data.frame(bg.x=x))
  lines(x,y,col="blue",lwd=2)

  
}


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Function to fit a n-Gaussian to disentangle multiple peaks
##
fitnGaus <- function(nPeaks){

  BkgdType <- "Linear"

  ##--------------------------------------
  ## Function for background + some peaks
  nGaus <- function(par){
    ## Remember, R calculates things as whole arrays, so ThCouts will
    ##  be calculated for each channel in the range
    ##
    ## first two parameters are background
    ## theoretical bin counts
    ThCounts <- BkgdFn(par[1:2])
    ##
    if(nPeaks > 0){
      ## Next parameter is peak width (same for all peaks)
##      Width <- par[3]
      ## Need to loop over peaks 
      for(PeakNo in 1:nPeaks){
        Width <- par[PeakNo+2]
        amp <- par[nPeaks+PeakNo+2]
        pos <- par[2+PeakNo+2*nPeaks]
        ThCounts <- ThCounts + amp*exp(-0.5*((Chn-pos)/Width)^2)
      }
    }
    ##
    return(ThCounts)
  }


  BkgdFn <- function(Bkgd){
    ##    
    if(BkgdType=="Quadratic"){
      Bxx  <- Bkgd[1]
      BInt <- Bkgd[2]
      BkgdCounts <- Bxx*(Chn-BInt)^2
    } else if(BkgdType=="Linear"){
      BSlope <- Bkgd[1]
      BInt <- Bkgd[2]
      BkgdCounts <- BInt + BSlope*Chn
    }
    ##
    return(BkgdCounts)
  }


  
  ## first get user to input the range to fit over
  region <- select.single()

  ## Only looks at truncated spart of the spectrum
  Chn <- region[,1]
  Counts <- region[,2]
  minChn <- Chn[1]
  maxChn <- Chn[length(Chn)]
  ## The maximum peak amplitude (really rough guess)
  Amax <- max(Counts)


  ## Initial guess for background
  if(BkgdType=="Quadratic")
    Bkgd <- c(max(Counts)/(max(Chn)-min(Chn))^2,max(Chn))
  if(BkgdType=="Linear")
    Bkgd <- c(0,mean(Counts[1:4]))

  ## Initial guesses for peak positions (spread them out evenly)
##  peakAmp <- rep(Amax-Bkgd[2],nPeaks)
##  peakPos <- sapply(1:nPeaks,function(i)Chn[1]+
##                    (i*(Chn[length(Chn)]-Chn[1])/(nPeaks+1)))
  cat("Guess where the peaks are\n")
  guess <- locator(nPeaks)
  peakPos <- guess$x
  peakAmp <- guess$y
  ## Guess at width (in energy - keV)
  Width <- rep(3,nPeaks)

  ##--------------------------------------
  ## Chi-Squared LogLiklihood
  chisqFn <- function(par){
    ThCounts <- nGaus(par) ## get the theoretical counts using fit
                              ## parameters
    ## Gaussian
    den <- Counts
    den[Counts==0] <- 1
    chisq <- sum(((ThCounts-Counts)^2)/den)
    ## Least-squares
    ## chisq <- sum((ThCounts-Counts)^2)
    ## Poisson
    ## if(any(ThCounts < 0) || any(ThCounts == 0))return(1e10)
    ## good <- Counts>0
    ## chisq <- -2*log(prod((((ThCounts^Counts)/factorial(Counts))*
    ##                    exp(-ThCounts))[good]))
    ## print(chisq)
    ## Cash statistic
    ## chisq <- -2*sum(ThCounts-Counts*log(ThCounts))
    ## if(is.nan(chisq))chisq <- 1e10
    ## print(chisq)
    ## G-Test
    ## good <- Counts>0
    ## chisq <- 2*sum(Counts[good]*log(Counts[good]/ThCounts[good]))
    ## if(is.nan(chisq))chisq <- 1e10
    ## print(chisq)
    ##
    ## if anything other than bkgd slope is negative, give a huge chisq
    if(BkgdType=="Quadratic")
      if(min(par[1:length(par)])<0)chisq <- 1e10
    if(BkgdType=="Linear")
      if(min(par[2:length(par)])<0)chisq <- 1e10
    ##
    ## if peak positions are outside the window, make chisq massive
    if(nPeaks>1 ){
      if(min(tail(par,n=nPeaks)) < minChn) chisq <- chisq+1e5
      if(max(tail(par,n=nPeaks)) > maxChn) chisq <- chisq+1e5
    }
    ##
    ## Keep background minimum above maximum
 ##   if(par[2]<maxChn)chisq <- chisq+1e5
    ##
    ## Keep width less than max width
    if(nPeaks>0)
      if(par[3] > 10 || par[3] < 3) chisq <- 1e7
    ##
    return(chisq)
  }

  # initial starting parameters for the fit
  par <- c(Bkgd,Width,peakAmp,peakPos)
##    fit <- optim(par,chisqFn,method="BFGS",hessian=TRUE)
  fit <- optim(par,chisqFn,method="Nelder-Mead",hessian=TRUE)
  fittedPar <<- fit$par
  chisq <- fit$value
  lines(Chn,BkgdFn(fittedPar[1:2]),col="green",lwd=2)
  lines(Chn,nGaus(fittedPar),col="black",lwd=2)
  ## Draw the individual peaks
  gaus <- function(x,mu,sigma){
    exp(-0.5*((x-mu)/sigma)^2)/(sigma*sqrt(2*pi))
  }
  for(i in 1:nPeaks){
    pos <- fittedPar[2+2*nPeaks+i]
    width <- fittedPar[2+i]
    I <- fittedPar[2+nPeaks+i]
    lines(Chn,I*gaus(Chn,pos,width)/max(gaus(Chn,pos,width)),
          col="blue",lwd=2)
  }

  ## Standard errors
  std.err <- sqrt(diag(solve(fit$hessian)))
  for(i in 1:nPeaks){
      
      pos <- fittedPar[2+2*nPeaks+i]
      scale <- diff(Energy[c(floor(pos),ceiling(pos))])
      dpos <- std.err[2+2*nPeaks+i]*scale
      width <- fittedPar[2+i]*scale
      dwidth <- std.err[2+i]*scale
      ## Interpolate channels to get energy
      E <- approx(x=c(floor(pos),ceiling(pos)),
                  y=Energy[c(floor(pos),ceiling(pos))],
                  xout=pos)$y
      dE <- sqrt(dEnergy[floor(pos)]^2 + dpos^2)
      A <- fittedPar[i+2+nPeaks]
      sig <- fittedPar[2+i]
      area <- A*sig*sqrt(2*pi)
      err <- area*sqrt(std.err[2+i]/sig^2 + std.err[i+2+nPeaks]/A^2)
    cat("For peak ",i,":\n",
        " Position = ",E," +- ",dE,"\n",
        " Width = ",width," +- ",dwidth,"\n",
        " Area = ",area," +- ",err,"\n\n",
        sep="")
  }
  cat("Chisq = ",chisq,"\n",sep="")
  
}
fitnGaus_bkup <- function(nPeaks){

  BkgdType <- "Linear"

  ##--------------------------------------
  ## Function for background + some peaks
  nGaus <- function(par){
    ## Remember, R calculates things as whole arrays, so ThCouts will
    ##  be calculated for each channel in the range
    ##
    ## first two parameters are background
    ## theoretical bin counts
    ThCounts <- BkgdFn(par[1:2])
    ##
    if(nPeaks > 0){
      ## Next parameter is peak width (same for all peaks)
      Width <- par[3]
      ## Need to loop over peaks 
      for(PeakNo in 1:nPeaks){
        amp <- par[PeakNo+3]
        pos <- par[3+PeakNo+nPeaks]
        ThCounts <- ThCounts + amp*exp(-0.5*((Chn-pos)/Width)^2)
      }
    }
    ##
    return(ThCounts)
  }


  BkgdFn <- function(Bkgd){
    ##    
    if(BkgdType=="Quadratic"){
      Bxx  <- Bkgd[1]
      BInt <- Bkgd[2]
      BkgdCounts <- Bxx*(Chn-BInt)^2
    } else if(BkgdType=="Linear"){
      BSlope <- Bkgd[1]
      BInt <- Bkgd[2]
      BkgdCounts <- BInt + BSlope*Chn
    }
    ##
    return(BkgdCounts)
  }


  
  ## first get user to input the range to fit over
  region <- select.single()

  ## Only looks at truncated spart of the spectrum
  Chn <- region[,1]
  Counts <- region[,2]
  minChn <- Chn[1]
  maxChn <- Chn[length(Chn)]
  ## The maximum peak amplitude (really rough guess)
  Amax <- max(Counts)

  ## Guess at width (in energy - keV)
  Width <- 10

  ## Initial guess for background
  if(BkgdType=="Quadratic")
    Bkgd <- c(max(Counts)/(max(Chn)-min(Chn))^2,max(Chn))
  if(BkgdType=="Linear")
    Bkgd <- c(0,mean(Counts))

  ## Initial guesses for peak positions (spread them out evenly)
  peakAmp <- rep(Amax-Bkgd[2],nPeaks)
  peakPos <- sapply(1:nPeaks,function(i)Chn[1]+
                    (i*(Chn[length(Chn)]-Chn[1])/(nPeaks+1)))
    
  ##--------------------------------------
  ## Chi-Squared LogLiklihood
  chisqFn <- function(par){
    ThCounts <- nGaus(par) ## get the theoretical counts using fit
                              ## parameters
    ## Gaussian
    den <- Counts
    den[Counts==0] <- 1
    chisq <- sum(((ThCounts-Counts)^2)/den)
    ## Least-squares
    ##chisq <- sum((ThCounts-Counts)^2)
    ## Poisson
    ##if(any(ThCounts < 0) || any(ThCounts == 0))return(1e10)
    ##chisq <- log(prod(((ThCounts^Counts)/factorial(Counts))*
    ##                  exp(-ThCounts)))
    ## Cash statistic
    ##chisq <- -2*sum(ThCounts-Counts*log(ThCounts))
    ##    if(is.nan(chisq))chisq <- 1e10
    ##print(chisq)
    ##
    ## if anything other than bkgd slope is negative, give a huge chisq
    if(BkgdType=="Quadratic")
      if(min(par[1:length(par)])<0)chisq <- 1e10
    if(BkgdType=="Linear")
      if(min(par[2:length(par)])<0)chisq <- 1e10
    ##
    ## if peak positions are outside the window, make chisq massive
    if(nPeaks>1 ){
      if(min(tail(par,n=nPeaks)) < minChn) chisq <- chisq+1e5
      if(max(tail(par,n=nPeaks)) > maxChn) chisq <- chisq+1e5
    }
    ##
    ## Keep background minimum above maximum
 ##   if(par[2]<maxChn)chisq <- chisq+1e5
    ##
    ## Keep width less than max width
    if(nPeaks>0)
      if(par[3] > 20 || par[3] < 5) chisq <- 1e7
    ##
    return(chisq)
  }

  # initial starting parameters for the fit
  par <- c(Bkgd,Width,peakAmp,peakPos)
  ##  fit <- optim(par,chisqFn,method="BFGS")
  fit <- optim(par,chisqFn,method="Nelder-Mead",hessian=TRUE)
  fittedPar <<- fit$par
  chisq <- fit$value
  lines(Chn,BkgdFn(fittedPar[1:2]),col="green",lwd=2)
  lines(Chn,nGaus(fittedPar),col="blue",lwd=2)
  ## Draw the individual peaks
  gaus <- function(x,mu,sigma){
    exp(-0.5*((x-mu)/sigma)^2)/(sigma*sqrt(2*pi))
  }
  for(i in 1:nPeaks){
    pos <- fittedPar[3+nPeaks+i]
    width <- fittedPar[3]
    I <- fittedPar[3+i]
    lines(Chn,I*gaus(Chn,pos,width)/max(gaus(Chn,pos,width)),
          col="blue",lwd=2)
  }

  ## Standard errors
  std.err <- sqrt(diag(solve(fit$hessian)))
  for(i in 1:nPeaks){
    pos <- fittedPar[3+nPeaks+i]
    dpos <- std.err[3+nPeaks+i]
    width <- fittedPar[3]
    dwidth <- std.err[3]
    ## Interpolate channels to get energy
    E <- approx(x=c(floor(pos),ceiling(pos)),
                     y=Energy[c(floor(pos),ceiling(pos))],
                     xout=pos)$y
    dE <- sqrt(dEnergy[floor(pos)]^2 + dpos^2)
    A <- fittedPar[i+3]
    sig <- fittedPar[3]
    area <- A*sig*sqrt(2*pi)
    err <- sqrt(2*pi*(sig^2*std.err[i+3]^2 + A^2*std.err[3]^2))
    cat("For peak ",i,":\n",
        " Position = ",E," +- ",dE,"\n",
        " Width = ",width," +- ",dwidth,"\n",
        " Area = ",area," +- ",err,"\n\n",
        sep="")
  }
  cat("Chisq = ",chisq,"\n",sep="")
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Function to fit a n-Gaussian to disentangle multiple peaks
##
fitnGausKP <- function(nPeaks){

  BkgdType <- "Linear"

  ##--------------------------------------
  ## Function for background + some peaks
  nGaus <- function(par){
    ## Remember, R calculates things as whole arrays, so ThCouts will
    ##  be calculated for each channel in the range
    ##
    ## first two parameters are background
    ## theoretical bin counts
    ThCounts <- BkgdFn(par[1:2])
    ##
    if(nPeaks > 0){
      ## Next parameter is peak width (same for all peaks)
      Width <- par[3]
      ## Need to loop over peaks 
      for(PeakNo in 1:nPeaks){
        amp <- par[PeakNo+3]
        pos <- ppos[PeakNo] ##par[3+PeakNo+nPeaks]
        ThCounts <- ThCounts + amp*exp(-0.5*((Chn-pos)/Width)^2)
      }
    }
    ##
    return(ThCounts)
  }


  BkgdFn <- function(Bkgd){
    ##    
    if(BkgdType=="Quadratic"){
      Bxx  <- Bkgd[1]
      BInt <- Bkgd[2]
      BkgdCounts <- Bxx*(Chn-BInt)^2
    } else if(BkgdType=="Linear"){
      BSlope <- Bkgd[1]
      BInt <- Bkgd[2]
      BkgdCounts <- BInt + BSlope*Chn
    }
    ##
    return(BkgdCounts)
  }


  
  ## first get user to input the range to fit over
  region <- select.single()

  ## Enter Peak Energies
  cat("Enter peak energies for",nPeaks,"peaks in keV\n")
  ##ppos <- scan(nlines=nPeaks)
  ##ppos <- c(10340.7,10398.5,10481.0)

  ## Convert peak position into channel
  ppos <- approx(Energy,spec[,1],ppos)$y
  
  ## Only looks at truncated spart of the spectrum
  Chn <- region[,1]
  Counts <- region[,2]
  minChn <- Chn[1]
  maxChn <- Chn[length(Chn)]
  ## The maximum peak amplitude (really rough guess)
  Amax <- max(Counts)

  ## Guess at width (in energy - keV)
  Width <- 10

  ## Initial guess for background
  if(BkgdType=="Quadratic")
    Bkgd <- c(max(Counts)/(max(Chn)-min(Chn))^2,max(Chn))
  if(BkgdType=="Linear")
    Bkgd <- c(0,mean(Counts))

  ## Initial guesses for peak positions (spread them out evenly)
  peakAmp <- rep(Amax-Bkgd[2],nPeaks)
  peakPos <- sapply(1:nPeaks,function(i)Chn[1]+
                    (i*(Chn[length(Chn)]-Chn[1])/(nPeaks+1)))
    
  ##--------------------------------------
  ## Chi-Squared LogLiklihood
  chisqFn <- function(par){
    ThCounts <- nGaus(par) ## get the theoretical counts using fit
                              ## parameters
    ## Gaussian
    den <- Counts
    den[Counts==0] <- 1
    chisq <- sum(((ThCounts-Counts)^2)/den)
    ## Least-squares
    ##chisq <- sum((ThCounts-Counts)^2)
    ## Poisson
    ##if(any(ThCounts < 0) || any(ThCounts == 0))return(1e10)
    ##chisq <- log(prod(((ThCounts^Counts)/factorial(Counts))*
    ##                  exp(-ThCounts)))
    ## Cash statistic
    ##chisq <- -2*sum(ThCounts-Counts*log(ThCounts))
    ##    if(is.nan(chisq))chisq <- 1e10
    ##print(chisq)
    ##
    ## if anything other than bkgd slope is negative, give a huge chisq
    if(BkgdType=="Quadratic")
      if(min(par[1:length(par)])<0)chisq <- 1e10
    if(BkgdType=="Linear")
      if(min(par[2:length(par)])<0)chisq <- 1e10
    ##
    ## if peak positions are outside the window, make chisq massive
    if(nPeaks>1 ){
      if(min(tail(par,n=nPeaks)) < minChn) chisq <- chisq+1e5
      if(max(tail(par,n=nPeaks)) > maxChn) chisq <- chisq+1e5
    }
    ##
    ## Keep background minimum above maximum
    if(par[2]<maxChn)chisq <- chisq+1e5
    ##
    ## Keep width less than max width
    if(nPeaks>0)
      if(par[3] > 40 || par[3] < 5) chisq <- 1e7
    ##
    return(chisq)
  }

  # initial starting parameters for the fit
  par <- c(Bkgd,Width,peakAmp)
  ##  fit <- optim(par,chisqFn,method="BFGS")
  fit <- optim(par,chisqFn,method="Nelder-Mead",hessian=TRUE)
  fittedPar <<- fit$par
  lines(Chn,BkgdFn(fittedPar[1:2]),col="green",lwd=2)
  lines(Chn,nGaus(fittedPar),col="blue",lwd=2)
  ## Draw the individual peaks
  gaus <- function(x,mu,sigma){
    exp(-0.5*((x-mu)/sigma)^2)/(sigma*sqrt(2*pi))
  }
  for(i in 1:nPeaks){
    pos <- ppos[i] ##fittedPar[3+nPeaks+i]
    width <- fittedPar[3]
    I <- fittedPar[3+i]
    lines(Chn,I*gaus(Chn,pos,width)/max(gaus(Chn,pos,width)),
          col="blue",lwd=2)
  }

  ## Standard errors
  std.err <- sqrt(diag(solve(fit$hessian)))
  for(i in 1:nPeaks){
    pos <- ppos[i] ##fittedPar[3+nPeaks+i]
##    dpos <- std.err[3+nPeaks+i]
    width <- fittedPar[3]
    dwidth <- std.err[3]
    ## Interpolate channels to get energy
##    E <- approx(x=c(floor(pos),ceiling(pos)),
##                     y=Energy[c(floor(pos),ceiling(pos))],
##                     xout=pos)$y
    E <- approx(spec[,1],Energy,pos)$y
##    dE <- sqrt(dEnergy[floor(pos)]^2 + dpos^2)
    A <- fittedPar[i+3]
    sig <- fittedPar[3]
    area <- A*sig*sqrt(2*pi)
    err <- sqrt(2*pi*(sig^2*std.err[i+3]^2 + A^2*std.err[3]^2))
    cat("For peak ",i,":\n",
        " Position = ",E,"\n",
        " Width = ",width," +- ",dwidth,"\n",
        " Area = ",area," +- ",err,"\n\n",
        sep="")
  }
  
}

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Function to fit a n-Gaussian to disentangle multiple peaks
##
fitShoulder <- function(){

  nPeaks <- 1  ## Only one peak for now
  
  ##--------------------------------------
  ## Function for background + some peaks
  nGaus <- function(par){
    ## Remember, R calculates things as whole arrays, so ThCouts will
    ##  be calculated for each channel in the range
    ##
    ## first two parameters are background
    ## theoretical bin counts
    ThCounts <- BkgdFn(par[1:3])
    ##
    if(nPeaks > 0){
      ## Next parameter is peak width (same for all peaks)
      Width <- par[4]
      ## Need to loop over peaks 
      for(PeakNo in 1:nPeaks){
        amp <- par[PeakNo+4]
        pos <- par[4+PeakNo+nPeaks]
        ThCounts <- ThCounts + amp*exp(-0.5*((Chn-pos)/Width)^2)
      }
    }
    ##
    return(ThCounts)
  }


  BkgdFn <- function(Bkgd){
    ##    
##    BkgdCounts <- Bkgd[1]+exp(-Bkgd[2]*(Chn-Bkgd[3])^2)
    BkgdCounts <- bg+Bkgd[3]*exp(-((Chn-Bkgd[2])/Bkgd[1])^2)
    ##
    return(BkgdCounts)
  }


  
  ## first get user to input the range to fit over
  region <- select.single()
  cat("Guess where the peak is\n")
  guess <- locator(1)$x
  cat("Give a range that the peak is allowed\n")
  allowed <- locator(2)$x
  cat("Region for flat background\n")
  bg <- mean(locator(2)$y)
  
  ## Only looks at truncated spart of the spectrum
  Chn <- region[,1]
  Counts <- region[,2]
  minChn <- Chn[1]
  maxChn <- Chn[length(Chn)]
  ## The maximum peak amplitude (really rough guess)
  Amax <- max(Counts)

  ## Guess at width (in energy - keV)
  Width <- 10

  ## Initial guess for background
  Bkgd <- c(Width*2,max(Chn),max(Counts))

  ## Initial guesses for peak positions (spread them out evenly)
  peakAmp <- rep(max(Counts[Chn>allowed[1] & Chn<allowed[2]])-bg,nPeaks)
  ##  peakAmp <- rep(50,nPeaks)
  ##  peakPos <- sapply(1:nPeaks,function(i)Chn[1]+
  ##                    (i*(Chn[length(Chn)]-Chn[1])/(nPeaks+1)))
  peakPos <- rep(guess,nPeaks)
    
  ##--------------------------------------
  ## Chi-Squared LogLiklihood
  chisqFn <- function(par){
    ThCounts <- nGaus(par) ## get the theoretical counts using fit
                              ## parameters
    chisq <- sum(((ThCounts-Counts)/sqrt(Counts))^2)
    ##
    ## if anything other than bkgd slope is negative, give a huge chisq
    if(min(par)<0)chisq <- 1e10
    ##
    ## if peak positions are outside the window, make chisq massive
    if(nPeaks>0 ){
      if(min(tail(par,n=nPeaks)) < allowed[1]) chisq <- chisq+1e10
      if(max(tail(par,n=nPeaks)) > allowed[2]) chisq <- chisq+1e10
    }
    ##
    ## Keep background minimum above maximum
    ##    if(par[2]<maxChn)chisq <- chisq+1e5
    ##
    ## Keep width less than max width
    ##    if(nPeaks>0)
    ##      if(par[3] > 20^2 || par[3] < 8^2) chisq <- 1e7
    ##
    if(par[4] < 5)chisq <- chisq+1e7
    ##
    return(chisq)
  }
  ##--------------------------------------
  ## Chi-Squared LogLiklihood for background fit
  bgchisqFn <- function(par){
    ThCounts <- BkgdFn(par) ## get the theoretical counts using fit
                            ## parameters
    chisq <- sum(((ThCounts-Counts)/sqrt(Counts))^2)
    if(min(par)<0)chisq <- 1e10
    ##
    return(chisq)
  }

  ## initial starting parameters for the fit
  par <- c(Bkgd,Width,peakAmp,peakPos)
  ##  fit <- optim(par,chisqFn,method="BFGS")
  ## Fit background first
  par <- Bkgd
  bgfit <- optim(par,bgchisqFn,method="Nelder-Mead")
  par <- c(bgfit$par,Width,peakAmp,peakPos)
  ##  par[1] <- par[1]-10
  fit <- optim(par,chisqFn,method="Nelder-Mead",hessian=TRUE)
  fittedPar <<- fit$par
  lines(Chn,BkgdFn(fittedPar[1:3]),col="green",lwd=2)
  lines(Chn,nGaus(fittedPar),col="blue",lwd=2)
  ##  lines(Chn,BkgdFn(par[1:4]),col="green")
  ##  lines(Chn,nGaus(par),col="blue")
  ## Draw the individual peaks
  gaus <- function(x,mu,sigma){
    exp(-0.5*((x-mu)/sigma)^2)/(sigma*sqrt(2*pi))
  }
  for(i in 1:nPeaks){
    pos <- fittedPar[4+nPeaks+i]
    width <- fittedPar[4]
    I <- fittedPar[4+i]
    lines(Chn,I*gaus(Chn,pos,width)/max(gaus(Chn,pos,width)),
          col="blue",lwd=2)
  }

  ## Standard errors
  std.err <- sqrt(diag(solve(fit$hessian)))
  for(i in 1:nPeaks){
    pos <- fittedPar[4+nPeaks+i]
    dpos <- std.err[4+nPeaks+i]
    width <- fittedPar[4]
    dwidth <- std.err[4]
    ## Interpolate channels to get energy
    E <- approx(x=c(floor(pos),ceiling(pos)),
                     y=Energy[c(floor(pos),ceiling(pos))],
                     xout=pos)$y
    dE <- sqrt(dEnergy[floor(pos)]^2 + dpos^2)
    A <- fittedPar[i+4]
    sig <- fittedPar[4]
    area <- A*sig*sqrt(2*pi)
    err <- sqrt(2*pi*(sig^2*std.err[i+4]^2 + A^2*std.err[4]^2))
    cat("For peak ",i,":\n",
        " Position = ",E," +- ",dE,"\n",
        " Width = ",width," +- ",dwidth,"\n",
        " Area = ",area," +- ",err,"\n\n",
        sep="")
  }
  
}

## Load some reference lines of where peaks are expected
load.lines <- function(){

  linefile <- sub(".dat",".lines",filename,fixed=TRUE)

  Lines <- read.table(linefile,skip=1,header=FALSE)

  pos <- approx(Energy,spec[,1],xout=Lines[,1])$y
  RefLines <<- pos
    abline(v=pos,lty=2,col="red")
    text(x=pos,y=yrange[2],Lines[,2],col="red",pos=1,xpd=NA)
  print(Lines)
}



cat("\n",
    "To fit multiple peaks: fitnGaus(n=2)","\n",
    "To fit peaks with know energy: fitnGausKP(n=2)","\n",
    "To fit peak and shoulder: fitShoulder()","\n",
    "To fit peaks with MCMC: MCMCGaus(npeaks=2)","\n",
    "\n",sep="")
