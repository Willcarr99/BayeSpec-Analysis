## Use what I've learned to fit Gaussians
## 
library(rethinking)
library(prettyGraphs)
library(emg) # Exponentially modified Gaussian package
library(pracma) # erfc package (error functions)
library(xlsx) # Write data to .xlsx file
#library(tibble) # tibble package or all-encompassing tidyverse package

cat("\n",
    "To MCMC fit multiple types of peaks: MCMCExpModGausAndGaus(npeaks=2,nGaus=1,nEMG=1)","\n",
    "\n",sep="")

## Function to generate n Gaussians
mcmc.nGaus <- function(x,nGaus,intens,cent,sig){
    g <- sapply(1:nGaus, function(i){
        c <- cent[i]
        s <- sig[i]
        I <- abs(intens[i])
        I*exp(-0.5*((x-c)^2)/s^2)
    })
    ##lines(x,rowSums(g),col=add.alpha(rangi2,0.2))
    rowSums(g)
}

## Function to generate n EMG's (further modified to switch the tail from the right side to the left)
## Remember to get the tail on the left, x -> -x and c -> -c
mcmc.nEMG <- function(x,nEMG,intens,cent,sig,lamb){
    g <- sapply(1:nEMG, function(i){
        c <- cent[i]
        s <- sig[i]
        l <- lamb[i]
        I <- abs(intens[i])
        # Exponentially modified Gaussian (tail on the right, x -> x, c -> c)
        #0.5 * l * exp((l/2) * ((2*c) + (l*(s**2)) - (2*x))) * erfc((1 / (sqrt(2) * s)) * (c + (l*(s**2)) - x)) 
        #demg(x,mu=c,sigma=s,lambda=l,log=FALSE)
        # Exponentially modified Gaussian (tail on the left, x -> -x, c -> -c)
        #0.5 * l * exp((l/2) * (-(2*c) + (l*(s**2)) - (2*x))) * erfc((1 / (sqrt(2) * s)) * (-c + (l*(s**2)) - x))
        #demg(-x,mu=-c,sigma=s,lambda=l,log=FALSE)

        # In terms of intensity (tail on right, x -> x, c -> c)
        if (rightTail){
            I*s*l*sqrt(pi/2)*exp((0.5*(s*l)**2) - ((x-c)*l)) * erfc((1/sqrt(2))*((s*l) - ((x-c)/s)))
        } else {
        # In terms of intensity (tail on left, x -> -x, c -> -c)
            I*s*l*sqrt(pi/2)*exp((0.5*(s*l)**2) - ((c-x)*l)) * erfc((1/sqrt(2))*((s*l) - ((c-x)/s)))
        }
    })
    rowSums(g)
}

## Background function
mcmc.BG <- function(x,intens,slope){
    b <- intens + slope*(x-min(x))
    b[b<0] <- 0
    b
}

mcmc.erfcx_eq <- function(x,s,l){
    erfcx(x) - (sqrt(2/pi) / (s*l))
}

mcmc.erfcxinvroot <- function(s,l){
    uniroot(mcmc.erfcx_eq, lower = -8, upper = 10, s=s, l=l)$root
}

## plot the spectrum
mcmc.plotspec <- function(Chn,trial,ylim=NULL,log=''){
    X11(width=10,height=10)
    par(mfrow=c(2,1))
    if(is.null(ylim))ylim <- c(0,1.25*max(trial))
    plot(Chn+0.5,trial,type='n',xlab="Channel",ylab="Counts",yaxs='i',
         ylim=ylim,log=log)
    lines(Chn+0.5,trial,type='S')
}

## Extract a few samples to draw the fit
plotfit <- function(model,Chn,trial,areas_Gaus,areas_EMG,modes,ymodes,npeaks,MCSamples=TRUE,PlotExpModGausAndGaus=TRUE,PlotMode=TRUE,PlotyMode=TRUE){

    if(MCSamples){
        post <- extract.samples(model,n=50)
        
        n_Gaus <- dim(post$cent_Gaus)[2] # number of Gaus peaks or NULL if 1
        n_EMG <- dim(post$cent_EMG)[2] # number of EMG peaks or NULL if 1

        ################################ TESTING ############################################

        ## If there is both 1 Gaussian and 1 EMG peak
        if(is.null(n_Gaus) & is.null(n_EMG)){
            ## Draw those samples
            for(i in 1:length(post$cent_Gaus)){
                lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_Gaus[i],post$cent_Gaus[i],post$sig_Gaus[i]),
                      col=add.alpha("red",0.1))
                lines(Chn,mcmc.nEMG(Chn,n=1,post$intens_EMG[i],post$cent_EMG[i],post$sig_EMG[i],post$lamb[i]),
                      col=add.alpha("orange",0.1))
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=1,post$intens_Gaus[i],post$cent_Gaus[i],post$sig_Gaus[i])+
                          mcmc.nEMG(Chn,n=1,post$intens_EMG[i],post$cent_EMG[i],post$sig_EMG[i],post$lamb[i]),
                      col=add.alpha(rangi2,0.2))
                if(PlotExpModGausAndGaus){
                    lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_EMG[i],post$cent_EMG[i],post$sig_EMG[i]),
                          col=add.alpha("green",0.2))
                }
            }
            if(PlotMode){
                abline(v = modes[1] - modes[2], lty=2)
                abline(v = modes[1])
                abline(v = modes[1] + modes[2], lty=2)                    
            }
            if(PlotyMode){
                abline(h = ymodes[1] - ymodes[2], lty=2)
                abline(h = ymodes[1])
                abline(h = ymodes[1] + ymodes[2], lty=2)    
            }
        ## If there is 1 Gaussian, but more than 1 EMG peak
        } else if (is.null(n_Gaus) & !is.null(n_EMG)){
            ## Draw those samples
            for(i in 1:dim(post$cent_EMG)[1]){
                for(j in 1:n_EMG){
                    lines(Chn,mcmc.nEMG(Chn,n=1,post$intens_EMG[i,j],post$cent_EMG[i,j],post$sig_EMG[i,j],post$lamb[i,j]),
                          col=add.alpha("orange",0.1))
                    if(PlotExpModGausAndGaus){ # Plot the EMG as a Gaussian (parameters are given for the original unmodified Gaussian)
                        lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_EMG[i,j],post$cent_EMG[i,j],post$sig_EMG[i,j]),
                              col=add.alpha("green",0.2))
                    }
                }
                lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_Gaus[i],post$cent_Gaus[i],post$sig_Gaus[i]),
                      col=add.alpha("red",0.1))
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=1,post$intens_Gaus[i],post$cent_Gaus[i],post$sig_Gaus[i])+
                          mcmc.nEMG(Chn,n=n_EMG,post$intens_EMG[i,],post$cent_EMG[i,],post$sig_EMG[i,],post$lamb[i,]),
                      col=add.alpha(rangi2,0.2))
            }
            if(PlotMode){
                for(i in 1:n_EMG){
                    abline(v = modes[i,1] - modes[i,2], lty=2)
                    abline(v = modes[i,1])
                    abline(v = modes[i,1] + modes[i,2], lty=2)
                }            
            }
            if(PlotyMode){
                for(i in 1:n_EMG){
                    abline(h = ymodes[i,1] - ymodes[i,2], lty=2)
                    abline(h = ymodes[i,1])
                    abline(h = ymodes[i,1] + ymodes[i,2], lty=2)
                }
            }

        ## If there is 1 EMG, but more than 1 Gaussian peak
        } else if (!is.null(n_Gaus) & is.null(n_EMG)){
            ## Draw those samples
            for(i in 1:dim(post$cent_Gaus)[1]){
                for(j in 1:n_Gaus){
                    lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_Gaus[i,j],post$cent_Gaus[i,j],post$sig_Gaus[i,j]),
                          col=add.alpha("red",0.1))
                }
                lines(Chn,mcmc.nEMG(Chn,n=1,post$intens_EMG[i],post$cent_EMG[i],post$sig_EMG[i],post$lamb[i]),
                      col=add.alpha("orange",0.1))
                if(PlotExpModGausAndGaus){ # Plot the EMG as a Gaussian (parameters are given for the original unmodified Gaussian)
                    lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_EMG[i],post$cent_EMG[i],post$sig_EMG[i]),
                          col=add.alpha("green",0.2))
                }
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=n_Gaus,post$intens_Gaus[i,],post$cent_Gaus[i,],post$sig_Gaus[i,])+
                          mcmc.nEMG(Chn,n=1,post$intens_EMG[i],post$cent_EMG[i],post$sig_EMG[i],post$lamb[i]),
                      col=add.alpha(rangi2,0.2))
            }
            if(PlotMode){
                abline(v = modes[1] - modes[2], lty=2)
                abline(v = modes[1])
                abline(v = modes[1] + modes[2], lty=2)                    
            }
            if(PlotyMode){
                abline(h = ymodes[1] - ymodes[2], lty=2)
                abline(h = ymodes[1])
                abline(h = ymodes[1] + ymodes[2], lty=2)    
            }
        ## If there is more than 1 Gaussian and more than 1 EMG peak
        } else {
            ## Draw those samples
            for(i in 1:dim(post$cent_Gaus)[1]){ # number of samples
                for(j in 1:n_EMG){
                    lines(Chn,mcmc.nEMG(Chn,n=1,post$intens_EMG[i,j],post$cent_EMG[i,j],post$sig_EMG[i,j],post$lamb[i,j]),
                          col=add.alpha("orange",0.1))
                    if(PlotExpModGausAndGaus){ # Plot the EMG as a Gaussian (parameters are given for the original unmodified Gaussian)
                        lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_EMG[i,j],post$cent_EMG[i,j],post$sig_EMG[i,j]),
                              col=add.alpha("green",0.2))
                    }
                }
                for(j in 1:n_Gaus){
                    lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_Gaus[i,j],post$cent_Gaus[i,j],post$sig_Gaus[i,j]),
                          col=add.alpha("red",0.1))
                }
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=n_Gaus,post$intens_Gaus[i,],post$cent_Gaus[i,],post$sig_Gaus[i,])+
                          mcmc.nEMG(Chn,n=n_EMG,post$intens_EMG[i,],post$cent_EMG[i,],post$sig_EMG[i,],post$lamb[i,]),
                      col=add.alpha(rangi2,0.2))
            }
            if(PlotMode){
                for(i in 1:n_EMG){
                    abline(v = modes[i,1] - modes[i,2], lty=2)
                    abline(v = modes[i,1])
                    abline(v = modes[i,1] + modes[i,2], lty=2)
                }            
            }
            if(PlotyMode){
                for(i in 1:n_EMG){
                    abline(h = ymodes[i,1] - ymodes[i,2], lty=2)
                    abline(h = ymodes[i,1])
                    abline(h = ymodes[i,1] + ymodes[i,2], lty=2)
                }
            }
        }
        #######################################################################################

############# FIX THIS!!!!!!!!!!!!!!!!! #####################

    } else {
        post <- extract.samples(model,n=5000)
        
        n_Gaus <- dim(post$cent_Gaus)[2] # number of Gaus peaks or NULL if 1
        n_EMG <- dim(post$cent_EMG)[2] # number of EMG peaks or NULL if 1

        ## Draw those samples
        ##for(i in 1:dim(post$intens)[1]){
        bg.intens <- median(post$bg.intens) 
        bg.slope <- median(post$bg.slope)
        intens_Gaus <- apply(post$intens_Gaus,2,median)
        cent_Gaus <- apply(post$cent_Gaus,2,median)
        sig_Gaus <- apply(post$sig_Gaus,2,median)
        intens_EMG <- apply(post$intens_EMG,2,median)
        cent_EMG <- apply(post$cent_EMG,2,median)
        sig_EMG <- apply(post$sig_EMG,2,median)
        lamb <- apply(post$lamb,2,median)

        ##print(bg.intens,bg.slope,intens,cent,sig)
        
        for(j in 1:n_EMG){
            lines(Chn,mcmc.nEMG(Chn,n=1,intens_EMG[j],cent_EMG[j],sig_EMG[j],lamb[j]),
                  col=add.alpha("red",1))
            if(PlotExpModGausAndGaus){ # Plot the EMG as a Gaussian (parameters are given for the original unmodified Gaussian)
                lines(Chn,mcmc.nGaus(Chn,n=1,intens_EMG[j],cent_EMG[j],sig_EMG[j]),
                      col=add.alpha("green",0.2))
            }
        }
        for(j in 1:n_Gaus){
            lines(Chn,mcmc.nGaus(Chn,n=1,intens_Gaus[j],cent_Gaus[j],sig_Gaus[j]),
                  col=add.alpha("purple",1))
        }
        lines(Chn,mcmc.BG(Chn,bg.intens,bg.slope)+
                  mcmc.nGaus(Chn,n=n_Gaus,intens_Gaus,cent_Gaus,sig_Gaus),
                  mcmc.nEMG(Chn,n=n_EMG,intens_EMG,cent_EMG,sig_EMG,lamb),
              col=add.alpha(rangi2,1.0))
    }
    ## Plot the spectrum on top
    ##lines(Chn+0.5,trial,type='S')
    box()

    cents <- precis(mod,par=c("cent_Gaus"),depth=2)
    xx <- grconvertX(0.05,from="npc",to="user")
    yy <- grconvertY(0.9,from="npc",to="user")
    ##arrows(x0=cents[,3],x1=cents[,4],y0=yy,col=rangi2,code=3,angle=90,length=0.05)
    Etext <- format(cents[,1],digits=2,nsmall=2)
    dEtext <- format(cents[,2],digits=2,nsmall=2)
    ##text(xx,yy,labels="Energies:",pos=3)
    ##text(cents[,1],yy,labels=paste(Etext,"+-",dEtext),pos=3)

    yy <- grconvertY(0.8,from="npc",to="user")
    Atext <- format(areas_Gaus$mean,nsmall=1,digits=1)
    dAtext <- format(areas_Gaus$sd,nsmall=1,digits=1)
    ##text(xx,yy,labels="Areas:",pos=3)
    ##text(cents[,1],yy,labels=paste(Atext,"+-",dAtext),pos=3)

    ## Now plot the model parameters
    plot(precis(model,depth=2,pars=c("bg.intens","bg.slope","sig_Gaus","intens_Gaus")))
##    plot(precis(model,depth=2,pars=c("cent")))
##    y0 <- 1.75+1*((n:1)-1)
##    y1 <- 2.25+1*((n:1)-1)
##    segments(x0=Sig,y0=y0,y1=y1,lty=2,col=rangi2)
##    y0 <- 1.75+n+1*((n:1)-1)
##    y1 <- 2.25+n+1*((n:1)-1)
##    segments(x0=Cent,y0=y0,y1=y1,lty=2,col=rangi2)
##    y0 <- 1.75+2*n+1*((n:1)-1)
##    y1 <- 2.25+2*n+1*((n:1)-1)
##    segments(x0=Intens,y0=y0,y1=y1,lty=2,col=rangi2)
##    y0 <- 1.75+3*n
##    y1 <- 2.25+3*n
##    segments(x0=bg.Intens,y0=y0,y1=y1,lty=2,col=rangi2)
##    y0 <- 0.75
##    y1 <- 1.25
##    segments(x0=bg.Slope,y0=y0,y1=y1,lty=2,col=rangi2)
    
}

cat("\nSelect plotting region!\n\n")

region <- select.single(draw=FALSE)
Chn <- region[,1]
trial <- region[,2]

cat("\nY-limits\n\n")
yy <- locator(2)$y

mcmc.plotspec(Chn,trial,ylim=c(max(0,min(yy)),max(yy)),log='')

cat("\nSelect fitting region!\n\n")

region <- select.single(draw=FALSE)
Chn <- region[,1]
trial <- region[,2]

## Make a data frame of the portion of the spectrum we're looking at
d <- data.frame(Chn,trial)

## Use input to get starting parameters
cat("\nClick on the Gaussians!\n\n")
guess_Gaus <- locator(nGaus)
g.cent_Gaus <- guess_Gaus$x
g.intens_Gaus <- guess_Gaus$y
g.sig_Gaus <- rep(5,nGaus)

cat("\nClick on the Exponentially-Modified Gaussians (EMG)!\n\n")
guess_EMG <- locator(nEMG)
#g.cent_EMG <- guess$x
#g.intens_EMG <- guess$y

g.bg.intens <- max(1,min(trial))
g.bg.slope <- 0 # 0

g.sig_EMG <- rep(3.5,nEMG)
g.lamb <- rep(0.1,nEMG) # 0.35 # lambda of 1 means nearly gaussian. Closer to 0 means a broader tail

## Peak x is ExpGauss mode and Peak y is ExpGauss ymode | Convert these to Gaussian component centroid and Gaussian component intensity
erfcxinvroot_guess <- mapply(mcmc.erfcxinvroot, g.sig_EMG, g.lamb)
if (rightTail){
    centfun <- function(mode, sig, lamb, z){mode - (lamb * sig^2) + (sqrt(2) * sig * z)}
} else {
    centfun <- function(mode, sig, lamb, z){mode + (lamb * sig^2) - (sqrt(2) * sig * z)}
}
intensfun <- function(ymode, mode, cent, sig){ymode * exp(-0.5*((cent - mode) / sig)^2)}
g.cent_EMG <- centfun(guess_EMG$x, g.sig_EMG, g.lamb, erfcxinvroot_guess)
g.intens_EMG <- intensfun(guess_EMG$y, guess_EMG$x, g.cent_EMG, g.sig_EMG)

## Prior Gaussian standard deviations
g.bg.intens.sd <- 0.01 # 0.01
g.bg.slope.sd <- 0.001 # 0.0001

g.intens_Gaus.sd <- rep(10,nGaus) # 1
g.cent_Gaus.sd <- rep(0.1,nGaus) # 0.1
g.sig_Gaus.sd <- rep(0.01,nGaus) # 0.01

g.intens_EMG.sd <- rep(10,nEMG) # 10
g.cent_EMG.sd <- rep(0.1,nEMG) # 0.1

g.sig_EMG.sd <- rep(0.1,nEMG) # 0.1
g.lamb.sd <- rep(0.01,nEMG) # 0.003

## Fit the bkg + nGaus + nEMG model
mod <- quap(
    alist(
        trial ~ dnorm(bkg_nGaus_nEMG),
        bkg_nGaus_nEMG <- mcmc.BG(Chn,bg.intens,bg.slope)+
                    mcmc.nGaus(Chn,n=nGaus,intens=intens_Gaus,cent=cent_Gaus,sig=sig_Gaus)+
                    mcmc.nEMG(Chn,n=nEMG,intens=intens_EMG,cent=cent_EMG,sig=sig_EMG,lamb=lamb),
        
        ## Automatically generate priors from clicks
        intens_Gaus ~ dnorm(g.intens_Gaus, g.intens_Gaus.sd),
        cent_Gaus ~ dnorm(g.cent_Gaus, g.cent_Gaus.sd),
        sig_Gaus ~ dnorm(g.sig_Gaus, g.sig_Gaus.sd),
        intens_EMG ~ dnorm(g.intens_EMG, g.intens_EMG.sd),
        cent_EMG ~ dnorm(g.cent_EMG, g.cent_EMG.sd),
        sig_EMG ~ dnorm(g.sig_EMG, g.sig_EMG.sd),
        lamb ~ dnorm(g.lamb, g.lamb.sd),
        bg.intens ~ dnorm(g.bg.intens, g.bg.intens.sd),
        bg.slope ~ dnorm(g.bg.slope, g.bg.slope.sd)

        ## Manually generate priors
        #intens_Gaus ~ dnorm(c(25,27,40),0.01),
        #cent_Gaus ~ dnorm(c(2317,2332,2411),0.1),
        #sig_Gaus ~ dnorm(c(4.2,4.3,5),0.1),
        #intens_EMG ~ dnorm(4401,1),
        #cent_EMG ~ dnorm(2365,0.001),
        #sig_EMG ~ dnorm(3.91,0.001),
        #lamb ~ dnorm(0.278,0.001),
        #bg.intens ~ dnorm(g.bg.intens,0.002),
        #bg.slope ~ dnorm(0,0.0001)

    ),
    data=d,
    ## Automatically generate initial guesses from clicks
    start=list(bg.intens=g.bg.intens,bg.slope=g.bg.slope,
               intens_Gaus=g.intens_Gaus,cent_Gaus=g.cent_Gaus,sig_Gaus=g.sig_Gaus,
               intens_EMG=g.intens_EMG,cent_EMG=g.cent_EMG,sig_EMG=g.sig_EMG,lamb=g.lamb)

    ## Manually generate initial guesses
    #start=list(bg.intens=g.bg.intens,bg.slope=g.bg.slope,
    #            intens_Gaus=c(25,27,40),cent_Gaus=c(2317,2332,2411),sig_Gaus=c(4.2,4.3,5),
    #            intens_EMG=9361,cent_EMG=2372,sig_EMG=4,lamb=0.25)
)

post <- extract.samples(mod)

areafun <- function(i,s)sqrt(2*pi)*i*abs(s)

erfcxinvrootsamp <- mapply(mcmc.erfcxinvroot, post$sig_EMG, post$lamb)
if (rightTail){
    modefun <- function(z,mu,sig,lamb){mu + (lamb*sig^2) - (sqrt(2)*sig*z)}
} else {
    modefun <- function(z,mu,sig,lamb){mu - (lamb*sig^2) + (sqrt(2)*sig*z)}
}
ymodefun <- function(mode, i, mu, sig){i * exp(-0.5 * ((mu - mode) / sig)^2)}

asamp_Gaus <- areafun(post$intens_Gaus,post$sig_Gaus)
asamp_EMG <- areafun(post$intens_EMG,post$sig_EMG)

modesamp <- modefun(erfcxinvrootsamp, post$cent_EMG, post$sig_EMG, post$lamb)
ymodesamp <- ymodefun(modesamp, post$intens_EMG, post$cent_EMG, post$sig_EMG)

if(is.null(dim(asamp_Gaus))){
    a.mean_Gaus <- mean(asamp_Gaus)
    a.sd_Gaus <- sd(asamp_Gaus)
    sig.mean_Gaus <- mean(post$sig_Gaus)
    sig.sd_Gaus <- sd(post$sig_Gaus)
    centroid.mean_Gaus <- mean(post$cent_Gaus)
    centroid.sd_Gaus <- sd(post$cent_Gaus)
}else{
    a.mean_Gaus <- apply(asamp_Gaus,2,mean)
    a.sd_Gaus <- apply(asamp_Gaus,2,sd)
    sig.mean_Gaus <- apply(post$sig_Gaus,2,mean)
    sig.sd_Gaus <- apply(post$sig_Gaus,2,sd)
    centroid.mean_Gaus <- apply(post$cent_Gaus,2,mean)
    centroid.sd_Gaus <- apply(post$cent_Gaus,2,sd)
}
if(is.null(dim(asamp_EMG))){
    a.mean_EMG <- mean(asamp_EMG)
    a.sd_EMG <- sd(asamp_EMG)
    sig.mean_EMG <- mean(post$sig_EMG)
    sig.sd_EMG <- sd(post$sig_EMG)
    centroid.mean_EMG <- mean(post$cent_EMG)
    centroid.sd_EMG <- sd(post$cent_EMG)
    lamb.mean <- mean(post$lamb)
    lamb.sd <- sd(post$lamb)
    mode.mean <- mean(modesamp)
    mode.sd <- sd(modesamp)
    ymode.mean <- mean(ymodesamp)
    ymode.sd <- sd(ymodesamp)
}else{
    a.mean_EMG <- apply(asamp_EMG,2,mean)
    a.sd_EMG <- apply(asamp_EMG,2,sd)
    sig.mean_EMG <- apply(post$sig_EMG,2,mean)
    sig.sd_EMG <- apply(post$sig_EMG,2,sd)
    centroid.mean_EMG <- apply(post$cent_EMG,2,mean)
    centroid.sd_EMG <- apply(post$cent_EMG,2,sd)
    lamb.mean <- apply(post$lamb,2,mean)
    lamb.sd <- apply(post$lamb,2,sd)
    mode.mean <- apply(modesamp,2,mean)
    mode.sd <- apply(modesamp,2,sd)
    ymode.mean <- apply(ymodesamp,2,mean)
    ymode.sd <- apply(ymodesamp,2,sd)
}
areas_Gaus <- data.frame(mean=a.mean_Gaus,sd=a.sd_Gaus)
areas_EMG <- data.frame(mean=a.mean_EMG,sd=a.sd_EMG)
modes <- data.frame(modemean=mode.mean,modesd=mode.sd)
ymodes <- data.frame(ymodemean=ymode.mean,ymodesd=ymode.sd)

## Summarise the whole thing!
print_decimals <- 4
mod_df <- precis(mod,depth=2)
mod_df <- mod_df[-c(3,4)] # Delete the 5.5% and 94.5% confidence interval columns
names(mod_df)[1] <- "post_mean" # Change name of column mean to post_mean
names(mod_df)[2] <- "post_sd" # Change name of column sd to post_sd

# Create data frame columns for priors
prior_mean <- c(format(round(g.bg.intens,print_decimals),nsmall=print_decimals), format(round(g.bg.slope,print_decimals),nsmall=print_decimals))
prior_sd <- c(format(round(g.bg.intens.sd,print_decimals),nsmall=print_decimals), format(round(g.bg.slope.sd,print_decimals),nsmall=print_decimals))
for(i in 1:nGaus){
    prior_mean <- append(prior_mean, format(round(g.intens_Gaus[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.intens_Gaus.sd[i],print_decimals),nsmall=print_decimals))
}
for(i in 1:nGaus){
    prior_mean <- append(prior_mean, format(round(g.cent_Gaus[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.cent_Gaus.sd[i],print_decimals),nsmall=print_decimals))
}
for(i in 1:nGaus){
    prior_mean <- append(prior_mean, format(round(g.sig_Gaus[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.sig_Gaus.sd[i],print_decimals),nsmall=print_decimals))
}
for(i in 1:nEMG){
    prior_mean <- append(prior_mean, format(round(g.intens_EMG[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.intens_EMG.sd[i],print_decimals),nsmall=print_decimals))
}
for(i in 1:nEMG){
    prior_mean <- append(prior_mean, format(round(g.cent_EMG[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.cent_EMG.sd[i],print_decimals),nsmall=print_decimals))
}
for(i in 1:nEMG){
    prior_mean <- append(prior_mean, format(round(g.sig_EMG[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.sig_EMG.sd[i],print_decimals),nsmall=print_decimals))
}
for(i in 1:nEMG){
    prior_mean <- append(prior_mean, format(round(g.lamb[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.lamb.sd[i],print_decimals),nsmall=print_decimals))
}
mod_df$prior_mean <- prior_mean
mod_df$prior_sd <- prior_sd

cat("-----------------------------------\n")
cat("Posteriors and Priors [dnorm(mean, sd)]:\n")
print(mod_df)
cat("\nGaussian areas:\n")
print(areas_Gaus)
cat("\nExponentially-Modified Gaussian (EMG) areas/modes:\n")
areas_modes_EMG <- data.frame(areamean=a.mean_EMG,areasd=a.sd_EMG,modemean=mode.mean,modesd=mode.sd)
print(areas_modes_EMG)
cat("-----------------------------------\n\n")

plotfit(mod,Chn,trial,areas_Gaus,areas_EMG,modes,ymodes,npeaks=npeaks,TRUE,FALSE,TRUE,FALSE)

############################################################
## Will's Addition: Append areas and centroids to text file
############################################################

## Ask to append fit data to text file

##filepath <- "./39K-3He-d-40Ca/ElasticScattering/Data/fitdata_FP.txt"
#filepath <- "./39K-3He-d-40Ca/2019_08_26/Transfer/Data/fitdata_summed_transfer.txt"
#
#saveanswer <- readline(paste("Append data to",filepath,"?","[y,n]","\n",sep=" "))
#if(saveanswer == "y" | saveanswer == "yes" | saveanswer == "Y" | saveanswer == "Yes" | saveanswer == "YES" | saveanswer == "yES"){
#    ## Get areas, centroids (and their sd's, i.e. uncertainties) in a dataframe
#    if(is.null(dim(post$cent))){
#        centroid.mean <- mean(post$cent)
#        centroid.sd <- sd(post$cent)
#    }else{
#        centroid.mean <- apply(post$cent,2,mean)
#        centroid.sd <- apply(post$cent,2,sd)
#    }
#    if(is.null(dim(post$sig))){
#        sig.mean <- mean(post$sig)
#        sig.sd <- sd(post$sig)
#    }else{
#        sig.mean <- apply(post$sig,2,mean)
#        sig.sd <- apply(post$sig,2,sd)
#    }
#
#    fitdata <- data.frame(area=a.mean,sdarea=a.sd,sig=sig.mean,sdsig=sig.sd,cent=centroid.mean,sdcent=centroid.sd)
#    
#    ## Ask user if they want to append a header line before the data is appended (e.g. # 15 deg, run 88)
#    headinganswer <- readline("Enter a header line appended to .txt file before data [To not include header line, enter: n]\n")
#    if(headinganswer != "n" & headinganswer != "no" & headinganswer != "N" & headinganswer != "No" & headinganswer != "NO" & headinganswer != "nO"){
#        write(headinganswer,filepath,append=TRUE)
#    }
#
#    ## Write data to text file
#    write.table(fitdata,filepath,append=TRUE,sep=", ",row.names=FALSE,col.names=FALSE)
#    cat(paste("Fit data saved to",filepath,"\n",sep=" "))
#}else{
#    cat("Fit data not saved\n")
#}

##############################################################################
## Will's Addition: Append areas, sigs, cents, lambs, and modes to spreadsheet
##############################################################################

## File to append fit parameter data to
#filepath <- "/home/wcfox/engedaq_TUNL/events/2019-08-26_39K-3He-d/test_sheet.xlsx"
filepath <- "/home/wcfox/engedaq_TUNL/events/2019-08-26_39K-3He-d/fitdata_Transfer.xlsx"
#filepath <- "/home/wcfox/engedaq_TUNL/events/2019-08-26_39K-3He-d/fitdata_ES.xlsx"

## Column in spreadsheet to start appending data to
startCol <- 8 # FP Transfer fits
#startCol <- 27 # Si Transfer fits
#startCol <- 7 # FP ES fits
#startCol <- 25 # Si ES fits

saveanswer <- readline(paste("Add data to",filepath,"?","[y,n]","\n",sep=" "))
if(saveanswer == "y" | saveanswer == "yes" | saveanswer == "Y" | saveanswer == "Yes" | saveanswer == "YES" | saveanswer == "yES"){
    ## Get areas, sigs, cents, lambs (and their sd's, i.e. uncertainties), and EMG_modes in a dataframe

    fitdata_Gaus <- data.frame(area=a.mean_Gaus,sdarea=a.sd_Gaus,sig=sig.mean_Gaus,sdsig=sig.sd_Gaus,cent=centroid.mean_Gaus,sdcent=centroid.sd_Gaus)
    if (rightTail){
        fitdata_EMG <- data.frame(area=a.mean_EMG,sdarea=a.sd_EMG,sig=sig.mean_EMG,sdsig=sig.sd_EMG,cent=centroid.mean_EMG,sdcent=centroid.sd_EMG,lambda=lamb.mean,sdlambda=lamb.sd,EMGmode=mode.mean,sdEMGmode=mode.sd,RightTail="yes")
    } else {
        fitdata_EMG <- data.frame(area=a.mean_EMG,sdarea=a.sd_EMG,sig=sig.mean_EMG,sdsig=sig.sd_EMG,cent=centroid.mean_EMG,sdcent=centroid.sd_EMG,lambda=lamb.mean,sdlambda=lamb.sd,EMGmode=mode.mean,sdEMGmode=mode.sd)
    }

    wb <- loadWorkbook(filepath)
    sheets <- getSheets(wb)
    sheet <- sheets[[1]]

    ## Ask user for cell input. Which cell should the data start from (each successive number will be pasted in the next column, same row)
    if(is.null(dim(asamp_Gaus))){ # 1 peak
        rowanswer_Gaus <- readline("Enter row number to append Gaussian data. [Type \"n\" to not record this peak]\n")
        #colanswer <- readline("Enter column number.\n")
        if (rowanswer_Gaus != "n"){
            addDataFrame(fitdata_Gaus, sheet, col.names = FALSE, row.names = FALSE, startRow = rowanswer_Gaus, startColumn = startCol)
        }
    }else{ # Number of peaks = dim(asamp_Gaus)[2]
        for(i in 1:dim(asamp_Gaus)[2]){
            rowanswer_Gaus <- readline(paste("Enter row number to append Gaussian data for peak ",i,"/",dim(asamp_Gaus)[2],". [Type \"n\" to not record this peak]","\n",sep=""))
            #colanswer <- readline(paste("Enter column number to append data for peak ",i,"/",dim(asamp_Gaus)[2],".","\n",sep=""))
            if (rowanswer_Gaus != "n"){
                addDataFrame(fitdata_Gaus[i,], sheet, col.names = FALSE, row.names = FALSE, startRow = rowanswer_Gaus, startColumn = startCol)
            }
        }
    }

    if(is.null(dim(asamp_EMG))){ # 1 peak
        rowanswer_EMG <- readline("Enter row number to append EMG data. [Type \"n\" to not record this peak]\n")
        #colanswer <- readline("Enter column number.\n")
        if (rowanswer_EMG != "n"){
            addDataFrame(fitdata_EMG, sheet, col.names = FALSE, row.names = FALSE, startRow = rowanswer_EMG, startColumn = startCol)
        }
    }else{ # Number of peaks = dim(asamp_EMG)[2]
        for(i in 1:dim(asamp_EMG)[2]){
            rowanswer_EMG <- readline(paste("Enter row number to append EMG data for peak ",i,"/",dim(asamp_EMG)[2],". [Type \"n\" to not record this peak]","\n",sep=""))
            #colanswer <- readline(paste("Enter column number to append data for peak ",i,"/",dim(asamp_EMG)[2],".","\n",sep=""))
            if (rowanswer_EMG != "n"){
                addDataFrame(fitdata_EMG[i,], sheet, col.names = FALSE, row.names = FALSE, startRow = rowanswer_EMG, startColumn = startCol)
            }
        }
    }
    saveWorkbook(wb, filepath)

    ## Print whether data is saved
    cat(paste("Fit data saved to",filepath,"\n",sep=" "))
}else{
    cat("Fit data not saved\n")
}
