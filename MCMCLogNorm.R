## Use what I've learned to fit Gaussians
## 
library(rethinking)
library(prettyGraphs)

cat("\n",
    "To MCMC fit multiple peaks: MCMCLogNorm(npeaks=2)","\n",
    "\n",sep="")

## Function to generate n Log-Norms (modified to switch the tail from the right side to the left)
## Remember to get the tail on the left, x -> (2*mode) - x
mcmc.nLogNorm <- function(x,nLogNorm,intens,mode,sig){
    g <- sapply(1:nLogNorm, function(i){
        mu <- log(mode[i]) + (sig[i])**2
        A <- abs(intens[i]) * exp(mu - ((sig[i])**2 / 2))
        s <- sig[i]
        # Tail on the right side
        (A/x) * exp(-(log(x) - mu)**2 / (2 * (s)**2))
        # Tail on the left side: x -> (2*mode) - x. Argument of ln(2*mode - x) becomes negative if x > 2*mode, and division by zero if x = 2*mode.
        #if (x < 2*mode[i]){
        #(A/((2 * mode[i]) - x)) * exp(-(log((2 * mode[i]) - x) - mu)**2 / (2 * (s)**2))
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
plotfit <- function(model,Chn,trial,areas,npeaks,MCSamples=TRUE){

    if(MCSamples){
        post <- extract.samples(model,n=50)
        
        n <- dim(post$mode)[2]
        ## For one peak
        if(is.null(n)){
            n <- 1
            ## Draw those samples
            for(i in 1:length(post$mode)){
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nLogNorm(Chn,n=n,post$intens[i],post$mode[i],post$sig[i]),
                      col=add.alpha(rangi2,0.2))
            }
        } else {
            ## Draw those samples
            for(i in 1:dim(post$mode)[1]){
                for(j in 1:n)
                    lines(Chn,mcmc.nLogNorm(Chn,n=1,post$intens[i,j],post$mode[i,j],post$sig[i,j]),
                          col=add.alpha("red",0.1))
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nLogNorm(Chn,n=n,post$intens[i,],post$mode[i,],post$sig[i,]),
                      col=add.alpha(rangi2,0.2))
            }
        }
    } else {
        post <- extract.samples(model,n=5000)
        
        n <- dim(post$mode)[2]
        ## For one peak
        if(is.null(n)){
            n <- 1
            ## Get the median value
            bg.intens <- median(post$bg.intens) 
            bg.slope <- median(post$bg.slope)
            intens <- median(post$intens)
            mode <- median(post$mode)
            sig <- median(post$sig)
            
            ## Draw those samples
            ##for(i in 1:length(post$intens)){
            lines(Chn,mcmc.BG(Chn,bg.intens,bg.slope)+
                      mcmc.nLogNorm(Chn,n=n,intens,mode,sig),
                  col=add.alpha(rangi2,1))    
            ##           }
        } else {
            ## Draw those samples
            ##for(i in 1:dim(post$intens)[1]){
            bg.intens <- median(post$bg.intens) 
            bg.slope <- median(post$bg.slope)
            intens <- apply(post$intens,2,median)
            mode <- apply(post$mode,2,median)
            sig <- apply(post$sig,2,median)

            ##print(bg.intens,bg.slope,intens,cent,sig)
            
            for(j in 1:n)
                lines(Chn,mcmc.nLogNorm(Chn,n=1,intens[j],mode[j],sig[j]),
                      col=add.alpha("red",1))
            lines(Chn,mcmc.BG(Chn,bg.intens,bg.slope)+
                      mcmc.nLogNorm(Chn,n=n,intens,mode,sig),
                  col=add.alpha(rangi2,1.0))
        
        }

    }
    
    ## Plot the spectrum on top
    ##lines(Chn+0.5,trial,type='S')
    box()

    modes <- precis(mod,par=c("mode"),depth=2)
    xx <- grconvertX(0.05,from="npc",to="user")
    yy <- grconvertY(0.9,from="npc",to="user")
    ##arrows(x0=cents[,3],x1=cents[,4],y0=yy,col=rangi2,code=3,angle=90,length=0.05)
    Etext <- format(modes[,1],digits=2,nsmall=2)
    dEtext <- format(modes[,2],digits=2,nsmall=2)
    ##text(xx,yy,labels="Energies:",pos=3)
    ##text(cents[,1],yy,labels=paste(Etext,"+-",dEtext),pos=3)

    yy <- grconvertY(0.8,from="npc",to="user")
    Atext <- format(areas$mean,nsmall=1,digits=1)
    dAtext <- format(areas$sd,nsmall=1,digits=1)
    ##text(xx,yy,labels="Areas:",pos=3)
    ##text(cents[,1],yy,labels=paste(Atext,"+-",dAtext),pos=3)
    
    ## Now plot the model parameters
    plot(precis(model,depth=2,pars=c("bg.intens","bg.slope","sig","intens")))
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
cat("\nClick on the peaks!\n\n")
guess <- locator(npeaks)
g.mode <- guess$x
g.intens <- guess$y
g.bg.intens <- max(1,min(trial))
g.bg.slope <- 0
g.sig <- rep(0.1,npeaks)

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

## Fit the bkg + nLogNorm model
mod <- quap(
    alist(
        trial ~ dnorm(bkg_nLogNorm),
        bkg_nLogNorm <- mcmc.BG(Chn,bg.intens,bg.slope)+
                    mcmc.nLogNorm(Chn,n=npeaks,intens=intens,mode=mode,sig=sig),
        intens ~ dnorm(g.intens,20),
        mode ~ dnorm(g.mode,0.5),
        sig ~ dnorm(g.sig,0.1), # 2
        bg.intens ~ dnorm(g.bg.intens,0.001),
        bg.slope ~ dnorm(0,0.001)
    ),
    data=d,
    start=list(bg.intens=g.bg.intens,bg.slope=g.bg.slope,
               intens=g.intens,mode=g.mode,sig=g.sig)
)

post <- extract.samples(mod)

areafun <- function(i,s)sqrt(2*pi)*i*abs(s)
logNorm_A_fun <- function(i,mode,sig) i*exp(log(mode) + (sig)**2 - ((sig)**2 / 2)) # i is intens here
asamp <- areafun(logNorm_A_fun(post$intens,post$mode,post$sig),post$sig) # i param is LogNorm A here, not intens

if(is.null(dim(asamp))){
    a.mean <- mean(asamp)
    a.sd <- sd(asamp)
}else{
    a.mean <- apply(asamp,2,mean)
    a.sd <- apply(asamp,2,sd)
}
areas <- data.frame(mean=a.mean,sd=a.sd)

## Summarise the whole thing!
cat("-----------------------------------\n")
print(precis(mod,depth=2))
print(areas)
cat("-----------------------------------\n")

plotfit(mod,Chn,trial,areas,npeaks=npeaks,TRUE)

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
