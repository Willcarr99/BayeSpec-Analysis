## Use what I've learned to fit Gaussians
## 
library(rethinking)
library(prettyGraphs)

cat("\n",
    "To MCMC fit multiple types of peaks: MCMCLogNormAndGaus(npeaks=2,nGaus=1,nLogNorm=1)","\n",
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

## Function to generate n Log-Norms (modified to switch the tail from the right side to the left)
## Remember to get the tail on the left, x -> (2*mode) - x
mcmc.nLogNorm <- function(x,nLogNorm,intens,mode,sig){
    g <- sapply(1:nLogNorm, function(i){
        mu <- log(mode[i]) + (sig[i])**2
        A <- abs(intens[i]) * exp(mu - ((sig[i])**2 / 2))
        s <- sig[i]
        # Tail on the right side
        #(A/x) * exp(-(log(x) - mu)**2 / (2 * (s)**2))
        # Tail on the left side: x -> (2*mode) - x. Argument of ln(2*mode - x) becomes negative if x > 2*mode, and division by zero if x = 2*mode.
        #if (x < 2*mode[i]){
        (A/((2 * mode[i]) - x)) * exp(-(log((2 * mode[i]) - x) - mu)**2 / (2 * (s)**2))
        #}
    })
    rowSums(g)
}

## Background function
mcmc.BG <- function(x,intens,slope){
    b <- intens + slope*(x-min(x))
    b[b<0] <- 0
    b
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
plotfit <- function(model,Chn,trial,areas_Gaus,areas_LogNorm,npeaks,MCSamples=TRUE){

    if(MCSamples){
        post <- extract.samples(model,n=50)
        
        n_Gaus <- dim(post$cent_Gaus)[2] # number of Gaus peaks or NULL if 1
        n_LogNorm <- dim(post$mode_LogNorm)[2] # number of LogNorm peaks or NULL if 1

        ################################ TESTING ############################################

        ## If there is both 1 Gaussian and 1 LogNorm peak
        if(is.null(n_Gaus) & is.null(n_LogNorm)){
            ## Draw those samples
            for(i in 1:length(post$cent_Gaus)){
                lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_Gaus[i],post$cent_Gaus[i],post$sig_Gaus[i]),
                      col=add.alpha("red",0.1))
                lines(Chn,mcmc.nLogNorm(Chn,n=1,post$intens_LogNorm[i],post$mode_LogNorm[i],post$sig_LogNorm[i]),
                      col=add.alpha("orange",0.1))
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=1,post$intens_Gaus[i],post$cent_Gaus[i],post$sig_Gaus[i])+
                          mcmc.nLogNorm(Chn,n=1,post$intens_LogNorm[i],post$mode_LogNorm[i],post$sig_LogNorm[i]),
                      col=add.alpha(rangi2,0.2))
            }
        ## If there is 1 Gaussian, but more than 1 LogNorm peak
        } else if (is.null(n_Gaus) & !is.null(n_LogNorm)){
            ## Draw those samples
            for(i in 1:dim(post$mode_LogNorm)[1]){
                for(j in 1:n_LogNorm){
                    lines(Chn,mcmc.nLogNorm(Chn,n=1,post$intens_LogNorm[i,j],post$mode_LogNorm[i,j],post$sig_LogNorm[i,j]),
                          col=add.alpha("orange",0.1))
                }
                lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_Gaus[i],post$cent_Gaus[i],post$sig_Gaus[i]),
                      col=add.alpha("red",0.1))
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=1,post$intens_Gaus[i],post$cent_Gaus[i],post$sig_Gaus[i])+
                          mcmc.nLogNorm(Chn,n=n_LogNorm,post$intens_LogNorm[i,],post$mode_LogNorm[i,],post$sig_LogNorm[i,]),
                      col=add.alpha(rangi2,0.2))
            }

        ## If there is 1 LogNorm, but more than 1 Gaussian peak
        } else if (!is.null(n_Gaus) & is.null(n_LogNorm)){
            ## Draw those samples
            for(i in 1:dim(post$cent_Gaus)[1]){
                for(j in 1:n_Gaus){
                    lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_Gaus[i,j],post$cent_Gaus[i,j],post$sig_Gaus[i,j]),
                          col=add.alpha("red",0.1))
                }
                lines(Chn,mcmc.nLogNorm(Chn,n=1,post$intens_LogNorm[i],post$mode_LogNorm[i],post$sig_LogNorm[i]),
                      col=add.alpha("orange",0.1))
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=n_Gaus,post$intens_Gaus[i,],post$cent_Gaus[i,],post$sig_Gaus[i,])+
                          mcmc.nLogNorm(Chn,n=1,post$intens_LogNorm[i],post$mode_LogNorm[i],post$sig_LogNorm[i]),
                      col=add.alpha(rangi2,0.2))
            }
        ## If there is more than 1 Gaussian and more than 1 LogNorm peak
        } else {
            ## Draw those samples
            for(i in 1:dim(post$cent_Gaus)[1]){ # number of samples
                for(j in 1:n_LogNorm){
                    lines(Chn,mcmc.nLogNorm(Chn,n=1,post$intens_LogNorm[i,j],post$mode_LogNorm[i,j],post$sig_LogNorm[i,j]),
                          col=add.alpha("orange",0.1))
                }
                for(j in 1:n_Gaus){
                    lines(Chn,mcmc.nGaus(Chn,n=1,post$intens_Gaus[i,j],post$cent_Gaus[i,j],post$sig_Gaus[i,j]),
                          col=add.alpha("red",0.1))
                }
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=n_Gaus,post$intens_Gaus[i,],post$cent_Gaus[i,],post$sig_Gaus[i,])+
                          mcmc.nLogNorm(Chn,n=n_LogNorm,post$intens_LogNorm[i,],post$mode_LogNorm[i,],post$sig_LogNorm[i,]),
                      col=add.alpha(rangi2,0.2))
            }
        }
        #######################################################################################

############# FIX THIS!!!!!!!!!!!!!!!!! #####################

    } else {
        post <- extract.samples(model,n=5000)
        
        n_Gaus <- dim(post$cent_Gaus)[2] # number of Gaus peaks or NULL if 1
        n_LogNorm <- dim(post$mode_LogNorm)[2] # number of LogNorm peaks or NULL if 1

        ## Draw those samples
        ##for(i in 1:dim(post$intens)[1]){
        bg.intens <- median(post$bg.intens) 
        bg.slope <- median(post$bg.slope)
        intens_Gaus <- apply(post$intens_Gaus,2,median)
        cent_Gaus <- apply(post$cent_Gaus,2,median)
        sig_Gaus <- apply(post$sig_Gaus,2,median)
        intens_LogNorm <- apply(post$intens_LogNorm,2,median)
        mode_LogNorm <- apply(post$mode_LogNorm,2,median)
        sig_LogNorm <- apply(post$sig_LogNorm,2,median)

        ##print(bg.intens,bg.slope,intens,cent,sig)
        
        for(j in 1:n_LogNorm){
            lines(Chn,mcmc.nLogNorm(Chn,n=1,intens_LogNorm[j],mode_LogNorm[j],sig_LogNorm[j]),
                  col=add.alpha("red",1))
        }
        for(j in 1:n_Gaus){
            lines(Chn,mcmc.nGaus(Chn,n=1,intens_Gaus[j],cent_Gaus[j],sig_Gaus[j]),
                  col=add.alpha("purple",1))
        }
        lines(Chn,mcmc.BG(Chn,bg.intens,bg.slope)+
                  mcmc.nGaus(Chn,n=n_Gaus,intens_Gaus,cent_Gaus,sig_Gaus),
                  mcmc.nLogNorm(Chn,n=n_LogNorm,intens_LogNorm,mode_LogNorm,sig_LogNorm),
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
guess <- locator(nGaus)
g.cent_Gaus <- guess$x
g.intens_Gaus <- guess$y
g.sig_Gaus <- rep(5,nGaus)

cat("\nClick on the Log-Normal Distributions!\n\n")
guess <- locator(nLogNorm)
g.mode_LogNorm <- guess$x
g.intens_LogNorm <- guess$y
g.sig_LogNorm <- rep(0.1,nLogNorm)

g.bg.intens <- max(1,min(trial))
g.bg.slope <- 0

## Fit the bkg + nGaus model
#mod <- quap(
#    alist(
#        trial ~ dnorm(lambda),
#        lambda <- mcmc.BG(Chn,bg.intens,bg.slope)+
#                       mcmc.nGaus(Chn,n=npeaks,intens=intens,cent=cent,sig=sig),
#        intens ~ dnorm(g.intens,20),
#        cent ~ dnorm(g.cent,5),
#        sig ~ dnorm(g.sig,2),
#        bg.intens ~ dnorm(g.bg.intens,10),
#        bg.slope ~ dnorm(0,0.1)
#    ),
#    data=d,
#    start=list(bg.intens=g.bg.intens,bg.slope=g.bg.slope,
#               intens=g.intens,cent=g.cent,sig=g.sig)
#)

## Fit the bkg + nGaus + nLogNorm model
mod <- quap(
    alist(
        trial ~ dnorm(bkg_nGaus_nLogNorm),
        bkg_nGaus_nLogNorm <- mcmc.BG(Chn,bg.intens,bg.slope)+
                    mcmc.nGaus(Chn,n=nGaus,intens=intens_Gaus,cent=cent_Gaus,sig=sig_Gaus)+
                    mcmc.nLogNorm(Chn,n=nLogNorm,intens=intens_LogNorm,mode=mode_LogNorm,sig=sig_LogNorm),
        
        ## Automatically generate priors from clicks
        intens_Gaus ~ dnorm(g.intens_Gaus,10), #20
        cent_Gaus ~ dnorm(g.cent_Gaus,0.1), #5
        sig_Gaus ~ dnorm(g.sig_Gaus,1), # 2
        intens_LogNorm ~ dnorm(g.intens_LogNorm,10), #20
        mode_LogNorm ~ dnorm(g.mode_LogNorm,0.1), #5
        sig_LogNorm ~ dnorm(g.sig_LogNorm,0.01), #2
        bg.intens ~ dnorm(g.bg.intens,2), #10
        bg.slope ~ dnorm(0,0.01) # 0.1

        ## Manually generate priors
        #intens_Gaus ~ dnorm(c(35,60,45),2),
        #cent_Gaus ~ dnorm(c(2320,2336,2413),0.01),
        #sig_Gaus ~ dnorm(c(5,5,5),0.01),
        #intens_LogNorm ~ dnorm(7200,0.01),
        #mode_LogNorm ~ dnorm(2370,0.1)
        #sig_LogNorm ~ dnorm(4,0.0001),
        #bg.intens ~ dnorm(g.bg.intens,0.002),
        #bg.slope ~ dnorm(0,0.000001)

    ),
    data=d,
    ## Automatically generate initial guesses from clicks
    start=list(bg.intens=g.bg.intens,bg.slope=g.bg.slope,
               intens_Gaus=g.intens_Gaus,cent_Gaus=g.cent_Gaus,sig_Gaus=g.sig_Gaus,
               intens_LogNorm=g.intens_LogNorm,mode_LogNorm=g.mode_LogNorm,sig_LogNorm=g.sig_LogNorm)

    ## Manually generate initial guesses
    #start=list(bg.intens=g.bg.intens,bg.slope=g.bg.slope,
    #            intens_Gaus=c(35,60,45),cent_Gaus=c(2320,2336,2413),sig_Gaus=c(5,5,5),
    #            intens_LogNorm=9361,mode_LogNorm=2372,sig_LogNorm=4)
)

post <- extract.samples(mod)

areafun <- function(i,s)sqrt(2*pi)*i*abs(s)
logNorm_A_fun <- function(i,mode,sig) i*exp(log(mode) + (sig)**2 - ((sig)**2 / 2)) # i is intens here
asamp_Gaus <- areafun(post$intens_Gaus,post$sig_Gaus)
asamp_LogNorm <- areafun(logNorm_A_fun(post$intens_LogNorm,post$mode_LogNorm,post$sig_LogNorm),post$sig_LogNorm) # i param is LogNorm A here, not intens

if(is.null(dim(asamp_Gaus))){
    a.mean_Gaus <- mean(asamp_Gaus)
    a.sd_Gaus <- sd(asamp_Gaus)
}else{
    a.mean_Gaus <- apply(asamp_Gaus,2,mean)
    a.sd_Gaus <- apply(asamp_Gaus,2,sd)
}
if(is.null(dim(asamp_LogNorm))){
    a.mean_LogNorm <- mean(asamp_LogNorm)
    a.sd_LogNorm <- sd(asamp_LogNorm)
}else{
    a.mean_LogNorm <- apply(asamp_LogNorm,2,mean)
    a.sd_LogNorm <- apply(asamp_LogNorm,2,sd)
}
areas_Gaus <- data.frame(mean=a.mean_Gaus,sd=a.sd_Gaus)
areas_LogNorm <- data.frame(mean=a.mean_LogNorm,sd=a.sd_LogNorm)

## Summarise the whole thing!
cat("-----------------------------------\n")
print(precis(mod,depth=2))
cat("Gaussian areas:\n")
print(areas_Gaus)
cat("Log-Normal areas:\n")
print(areas_LogNorm)
cat("-----------------------------------\n")

plotfit(mod,Chn,trial,areas_Gaus,areas_LogNorm,npeaks=npeaks,TRUE)

############################################################
## Will's Addition: Append areas and centroids to text file
############################################################

## Ask to append fit data to text file

#filepath <- "./39K-3He-d-40Ca/ElasticScattering/Data/fitdata_FP.txt"
filepath <- "./39K-3He-d-40Ca/2019_08_26/Transfer/Data/fitdata_summed_transfer.txt"

saveanswer <- readline(paste("Append data to",filepath,"?","[y,n]","\n",sep=" "))
if(saveanswer == "y" | saveanswer == "yes" | saveanswer == "Y" | saveanswer == "Yes" | saveanswer == "YES" | saveanswer == "yES"){
    ## Get areas, centroids (and their sd's, i.e. uncertainties) in a dataframe
    if(is.null(dim(post$cent))){
        centroid.mean <- mean(post$cent)
        centroid.sd <- sd(post$cent)
    }else{
        centroid.mean <- apply(post$cent,2,mean)
        centroid.sd <- apply(post$cent,2,sd)
    }
    if(is.null(dim(post$sig))){
        sig.mean <- mean(post$sig)
        sig.sd <- sd(post$sig)
    }else{
        sig.mean <- apply(post$sig,2,mean)
        sig.sd <- apply(post$sig,2,sd)
    }

    fitdata <- data.frame(area=a.mean,sdarea=a.sd,sig=sig.mean,sdsig=sig.sd,cent=centroid.mean,sdcent=centroid.sd)
    
    ## Ask user if they want to append a header line before the data is appended (e.g. # 15 deg, run 88)
    headinganswer <- readline("Enter a header line appended to .txt file before data [To not include header line, enter: n]\n")
    if(headinganswer != "n" & headinganswer != "no" & headinganswer != "N" & headinganswer != "No" & headinganswer != "NO" & headinganswer != "nO"){
        write(headinganswer,filepath,append=TRUE)
    }

    ## Write data to text file
    write.table(fitdata,filepath,append=TRUE,sep=", ",row.names=FALSE,col.names=FALSE)
    cat(paste("Fit data saved to",filepath,"\n",sep=" "))
}else{
    cat("Fit data not saved\n")
}
