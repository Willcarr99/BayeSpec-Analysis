## Use what I've learned to fit Gaussians
## 
library(rethinking)
library(prettyGraphs)
library(emg) # Exponentially modified Gaussian package
library(pracma) # erfc package (error functions)
library(xlsx) # Write data to .xlsx file

cat("\n",
    "To MCMC fit multiple peaks: MCMCExpModGaus(npeaks=2)","\n",
    "\n",sep="")

## Function to generate n Gaussians
mcmc.nGaus <- function(x,npeaks,intens,cent,sig){
    g <- sapply(1:npeaks, function(i){
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
mcmc.nEMG <- function(x,npeaks,intens,cent,sig,lamb){
    g <- sapply(1:npeaks, function(i){
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
plotfit <- function(model,Chn,trial,areas,modes,ymodes,npeaks,MCSamples=TRUE,PlotExpModGausAndGaus=TRUE,PlotMode=TRUE,PlotyMode=TRUE){

    if(MCSamples){
        post <- extract.samples(model,n=50)
        
        n <- dim(post$cent)[2]
        ## For one peak
        if(is.null(n)){
            n <- 1
            ## Draw those samples
            for(i in 1:length(post$cent)){
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nEMG(Chn,n=n,post$intens[i],post$cent[i],post$sig[i],post$lamb[i]),
                      col=add.alpha(rangi2,0.2))
                if(PlotExpModGausAndGaus){
                    lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                              mcmc.nGaus(Chn,n=n,post$intens[i],post$cent[i],post$sig[i]),
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
        } else {
            ## Draw those samples
            for(i in 1:dim(post$cent)[1]){
                for(j in 1:n){
                    lines(Chn,mcmc.nEMG(Chn,n=1,post$intens[i,j],post$cent[i,j],post$sig[i,j],post$lamb[i,j]),
                          col=add.alpha("red",0.1))
                    if(PlotExpModGausAndGaus){
                        lines(Chn,mcmc.nGaus(Chn,n=1,post$intens[i,j],post$cent[i,j],post$sig[i,j]),
                              col=add.alpha("green",0.2))
                    }
                }
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nEMG(Chn,n=n,post$intens[i,],post$cent[i,],post$sig[i,],post$lamb[i,]),
                      col=add.alpha(rangi2,0.2))
            }
            if(PlotMode){
                for(i in 1:n){
                    abline(v = modes[i,1] - modes[i,2], lty=2)
                    abline(v = modes[i,1])
                    abline(v = modes[i,1] + modes[i,2], lty=2)
                }            
            }
            if(PlotyMode){
                for(i in 1:n){
                    abline(h = ymodes[i,1] - ymodes[i,2], lty=2)
                    abline(h = ymodes[i,1])
                    abline(h = ymodes[i,1] + ymodes[i,2], lty=2)
                }
            }
        }
    } else {
        post <- extract.samples(model,n=5000)
        
        n <- dim(post$cent)[2]
        ## For one peak
        if(is.null(n)){
            n <- 1
            ## Get the median value
            bg.intens <- median(post$bg.intens) 
            bg.slope <- median(post$bg.slope)
            intens <- median(post$intens)
            cent <- median(post$cent)
            sig <- median(post$sig)
            lamb <- median(post$lamb)

            
            ## Draw those samples
            ##for(i in 1:length(post$intens)){
            lines(Chn,mcmc.BG(Chn,bg.intens,bg.slope)+
                      mcmc.nEMG(Chn,n=n,intens,cent,sig,lamb),
                  col=add.alpha(rangi2,1))
            if(PlotExpModGausAndGaus){
                lines(Chn,mcmc.BG(Chn,post$bg.intens,post$bg.slope)+
                          mcmc.nGaus(Chn,n=n,intens,cent,sig),
                      col=add.alpha("green",0.2))
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
            ##           }
        } else {
            ## Draw those samples
            ##for(i in 1:dim(post$intens)[1]){
            bg.intens <- median(post$bg.intens) 
            bg.slope <- median(post$bg.slope)
            intens <- apply(post$intens,2,median)
            cent <- apply(post$cent,2,median)
            sig <- apply(post$sig,2,median)
            lamb <- apply(post$lamb,2,median)

            ##print(bg.intens,bg.slope,intens,cent,sig)
            
            for(j in 1:n){
                lines(Chn,mcmc.nEMG(Chn,n=1,intens[j],cent[j],sig[j],lamb[j]),
                      col=add.alpha("red",1))
                if(PlotExpModGausAndGaus){
                    lines(Chn,mcmc.nGaus(Chn,n=1,intens[j],cent[j],sig[j]),
                          col=add.alpha("green",0.2))
                }
            }
            lines(Chn,mcmc.BG(Chn,bg.intens,bg.slope)+
                      mcmc.nEMG(Chn,n=n,intens,cent,sig,lamb),
                  col=add.alpha(rangi2,1.0))
            if(PlotMode){
                for(i in 1:n){
                    abline(v = modes[i,1] - modes[i,2], lty=2)
                    abline(v = modes[i,1])
                    abline(v = modes[i,1] + modes[i,2], lty=2)
                }            
            }
            if(PlotyMode){
                for(i in 1:n){
                    abline(h = ymodes[i,1] - ymodes[i,2], lty=2)
                    abline(h = ymodes[i,1])
                    abline(h = ymodes[i,1] + ymodes[i,2], lty=2)
                }
            }
        }
    }

    ## Plot the spectrum on top
    #box()

    #cents <- precis(mod,par=c("cent"),depth=2)
    #xx <- grconvertX(0.05,from="npc",to="user")
    #yy <- grconvertY(0.9,from="npc",to="user")

    #Etext <- format(cents[,1],digits=2,nsmall=2)
    #dEtext <- format(cents[,2],digits=2,nsmall=2)


    #yy <- grconvertY(0.8,from="npc",to="user")
    #Atext <- format(areas$mean,nsmall=1,digits=1)
    #dAtext <- format(areas$sd,nsmall=1,digits=1)
    
    ## Now plot the model parameters
    #plot(precis(model,depth=2,pars=c("bg.intens","bg.slope","sig","intens")))

    
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
cat("\nClick on the peaks!\n\n")
guess <- locator(npeaks)

#g.cent <- guess$x
#g.intens <- guess$y
g.bg.intens <- max(1,min(trial))
g.bg.slope <- 0

g.sig <- rep(5,npeaks) # 5
g.lamb <- rep(0.09,npeaks) # 0.1 # lambda of 1 means nearly gaussian. Closer to 0 means a broader tail

## Peak x is ExpGauss mode and Peak y is ExpGauss ymode | Convert these to Gaussian component centroid and Gaussian component intensity
erfcxinvroot_guess <- mapply(mcmc.erfcxinvroot, g.sig, g.lamb)
if (rightTail){
    centfun <- function(mode, sig, lamb, z){mode - (lamb * sig^2) + (sqrt(2) * sig * z)}
} else {
    centfun <- function(mode, sig, lamb, z){mode + (lamb * sig^2) - (sqrt(2) * sig * z)}
}
intensfun <- function(ymode, mode, cent, sig){ymode * exp(-0.5*((cent - mode) / sig)^2)}
g.cent <- centfun(guess$x, g.sig, g.lamb, erfcxinvroot_guess)
g.intens <- intensfun(guess$y, guess$x, g.cent, g.sig)

## Prior Gaussian standard deviations:
# Default - g.bg.intens.sd = 0.1, g.bg.slope.sd = 0.1, g.intens.sd = 10, g.cent.sd = 0.1, g.sig.sd = 1, g.lamb.sd = 0.01
g.bg.intens.sd <- 0.1 # 0.1
g.bg.slope.sd <- 0.001 # 0.1

g.intens.sd <- rep(10,npeaks) # 10
g.cent.sd <- rep(0.1,npeaks) # 0.1

g.sig.sd <- rep(0.2,npeaks) # 1
g.lamb.sd <- rep(0.1,npeaks) # 0.01

## Fit the bkg + nEMG model
mod <- quap(
    alist(
        trial ~ dnorm(bkg_nEMG),
        bkg_nEMG <- mcmc.BG(Chn,bg.intens,bg.slope)+
                    mcmc.nEMG(Chn,n=npeaks,intens=intens,cent=cent,sig=sig,lamb=lamb),
        intens ~ dnorm(g.intens, g.intens.sd),
        cent ~ dnorm(g.cent, g.cent.sd),
        sig ~ dlnorm(g.sig, g.sig.sd),
        lamb ~ dnorm(g.lamb, g.lamb.sd),
        bg.intens ~ dnorm(g.bg.intens, g.bg.intens.sd),
        bg.slope ~ dnorm(g.bg.slope, g.bg.slope.sd)
    ),
    data=d,
    start=list(bg.intens=g.bg.intens,bg.slope=g.bg.slope,
               intens=g.intens,cent=g.cent,sig=g.sig,lamb=g.lamb)
)

post <- extract.samples(mod)

areafun <- function(i,s)sqrt(2*pi)*i*abs(s)

erfcxinvrootsamp <- mapply(mcmc.erfcxinvroot, post$sig, post$lamb)
if (rightTail){
    modefun <- function(z,mu,sig,lamb){mu + (lamb*sig^2) - (sqrt(2)*sig*z)}
} else {
    modefun <- function(z,mu,sig,lamb){mu - (lamb*sig^2) + (sqrt(2)*sig*z)}
}
ymodefun <- function(mode, i, mu, sig){i * exp(-0.5 * ((mu - mode) / sig)^2)}

asamp <- areafun(post$intens,post$sig)
modesamp <- modefun(erfcxinvrootsamp, post$cent, post$sig, post$lamb)
ymodesamp <- ymodefun(modesamp, post$intens, post$cent, post$sig)

if(is.null(dim(asamp))){
    a.mean <- mean(asamp)
    a.sd <- sd(asamp)
    sig.mean <- mean(post$sig)
    sig.sd <- sd(post$sig)
    centroid.mean <- mean(post$cent)
    centroid.sd <- sd(post$cent)
    lamb.mean <- mean(post$lamb)
    lamb.sd <- sd(post$lamb)
    mode.mean <- mean(modesamp)
    mode.sd <- sd(modesamp)
    ymode.mean <- mean(ymodesamp)
    ymode.sd <- sd(ymodesamp)
}else{
    a.mean <- apply(asamp,2,mean)
    a.sd <- apply(asamp,2,sd)
    sig.mean <- apply(post$sig,2,mean)
    sig.sd <- apply(post$sig,2,sd)
    centroid.mean <- apply(post$cent,2,mean)
    centroid.sd <- apply(post$cent,2,sd)
    lamb.mean <- apply(post$lamb,2,mean)
    lamb.sd <- apply(post$lamb,2,sd)
    mode.mean <- apply(modesamp,2,mean)
    mode.sd <- apply(modesamp,2,sd)
    ymode.mean <- apply(ymodesamp,2,mean)
    ymode.sd <- apply(ymodesamp,2,sd)
}
areas <- data.frame(areamean=a.mean,areasd=a.sd)
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
for(i in 1:npeaks){
    prior_mean <- append(prior_mean, format(round(g.intens[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.intens.sd[i],print_decimals),nsmall=print_decimals))
}
for(i in 1:npeaks){
    prior_mean <- append(prior_mean, format(round(g.cent[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.cent.sd[i],print_decimals),nsmall=print_decimals))
}
for(i in 1:npeaks){
    prior_mean <- append(prior_mean, format(round(g.sig[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.sig.sd[i],print_decimals),nsmall=print_decimals))
}
for(i in 1:npeaks){
    prior_mean <- append(prior_mean, format(round(g.lamb[i],print_decimals),nsmall=print_decimals))
    prior_sd <- append(prior_sd, format(round(g.lamb.sd[i],print_decimals),nsmall=print_decimals))
}
mod_df$prior_mean <- prior_mean
mod_df$prior_sd <- prior_sd

cat("-----------------------------------\n")
cat("Posteriors and Priors [dnorm(mean, sd)]:\n")
print(mod_df)
cat("\n\nExponentially-Modified Gaussian (EMG) areas and modes:\n")
areas_modes <- data.frame(areamean=a.mean,areasd=a.sd,modemean=mode.mean,modesd=mode.sd)
print(areas_modes)
cat("-----------------------------------\n\n")

plotfit(mod,Chn,trial,areas,modes,ymodes,npeaks=npeaks,TRUE,TRUE,TRUE,FALSE)

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

    if (rightTail){
        fitdata <- data.frame(area=a.mean,sdarea=a.sd,sig=sig.mean,sdsig=sig.sd,cent=centroid.mean,sdcent=centroid.sd,lambda=lamb.mean,sdlambda=lamb.sd,EMGmode=mode.mean,sdEMGmode=mode.sd,RightTail="yes")
    } else {
        fitdata <- data.frame(area=a.mean,sdarea=a.sd,sig=sig.mean,sdsig=sig.sd,cent=centroid.mean,sdcent=centroid.sd,lambda=lamb.mean,sdlambda=lamb.sd,EMGmode=mode.mean,sdEMGmode=mode.sd)
    }

    ## Ask user for cell input. Which cell should the data start from (each successive number will be pasted in the next column, same row)
    if(is.null(dim(asamp))){ # 1 peak
        wb <- loadWorkbook(filepath)
        sheets <- getSheets(wb)
        sheet <- sheets[[1]]
        rowanswer <- readline("Enter row number to append data.\n")
        #colanswer <- readline("Enter column number.\n")
        addDataFrame(fitdata, sheet, col.names = FALSE, row.names = FALSE, startRow = rowanswer, startColumn = startCol)
        saveWorkbook(wb, filepath)
    }else{ # Number of peaks = dim(asamp)[2]
        wb <- loadWorkbook(filepath)
        sheets <- getSheets(wb)
        sheet <- sheets[[1]]
        for(i in 1:dim(asamp)[2]){
            rowanswer <- readline(paste("Enter row number to append data for peak ",i,"/",dim(asamp)[2],". [Type \"n\" to not record this peak]","\n",sep=""))
            #colanswer <- readline(paste("Enter column number to append data for peak ",i,"/",dim(asamp)[2],".","\n",sep=""))
            if (rowanswer != "n"){
                addDataFrame(fitdata[i,], sheet, col.names = FALSE, row.names = FALSE, startRow = rowanswer, startColumn = startCol)
            }
        }
        saveWorkbook(wb, filepath)
    }

    ## Print whether data is saved
    cat(paste("Fit data saved to",filepath,"\n",sep=" "))
}else{
    cat("Fit data not saved\n")
}
