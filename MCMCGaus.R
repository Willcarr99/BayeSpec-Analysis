## Use what I've learned to fit Gaussians
## 
library(rethinking)
library(prettyGraphs)
library(xlsx) # Write data to .xlsx file

cat("\n",
    "To MCMC fit multiple peaks: MCMCGaus(n=2)","\n",
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
        
        n <- dim(post$cent)[2]
        ## For one peak
        if(is.null(n)){
            n <- 1
            ## Draw those samples
            for(i in 1:length(post$intens)){
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=n,post$intens[i],post$cent[i],post$sig[i]),
                      col=add.alpha(rangi2,0.2))
            }
        } else {
            ## Draw those samples
            for(i in 1:dim(post$intens)[1]){
                for(j in 1:n)
                    lines(Chn,mcmc.nGaus(Chn,n=1,post$intens[i,j],post$cent[i,j],post$sig[i,j]),
                          col=add.alpha("red",0.1))
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=n,post$intens[i,],post$cent[i,],post$sig[i,]),
                      col=add.alpha(rangi2,0.2))
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

            
            ## Draw those samples
            ##for(i in 1:length(post$intens)){
            lines(Chn,mcmc.BG(Chn,bg.intens,bg.slope)+
                      mcmc.nGaus(Chn,n=n,intens,cent,sig),
                  col=add.alpha(rangi2,1))
            ##           }
        } else {
            ## Draw those samples
            ##for(i in 1:dim(post$intens)[1]){
            bg.intens <- median(post$bg.intens) 
            bg.slope <- median(post$bg.slope)
            intens <- apply(post$intens,2,median)
            cent <- apply(post$cent,2,median)
            sig <- apply(post$sig,2,median)

            ##print(bg.intens,bg.slope,intens,cent,sig)
            
            for(j in 1:n)
                lines(Chn,mcmc.nGaus(Chn,n=1,intens[j],cent[j],sig[j]),
                      col=add.alpha("red",1))
            lines(Chn,mcmc.BG(Chn,bg.intens,bg.slope)+
                      mcmc.nGaus(Chn,n=n,intens,cent,sig),
                  col=add.alpha(rangi2,1.0))
        
        }

    }
    
    ## Plot the spectrum on top
    ##lines(Chn+0.5,trial,type='S')
    box()

    cents <- precis(mod,par=c("cent"),depth=2)
    xx <- grconvertX(0.05,from="npc",to="user")
    yy <- grconvertY(0.9,from="npc",to="user")
    ##arrows(x0=cents[,3],x1=cents[,4],y0=yy,col=rangi2,code=3,angle=90,length=0.05)
    Etext <- format(cents[,1],digits=2,nsmall=2)
    dEtext <- format(cents[,2],digits=2,nsmall=2)
    ##text(xx,yy,labels="Energies:",pos=3)
    ##text(cents[,1],yy,labels=paste(Etext,"+-",dEtext),pos=3)

    yy <- grconvertY(0.8,from="npc",to="user")
    Atext <- format(areas$mean,nsmall=1,digits=1)
    dAtext <- format(areas$sd,nsmall=1,digits=1)
    ##text(xx,yy,labels="Areas:",pos=3)
    ##text(cents[,1],yy,labels=paste(Atext,"+-",dAtext),pos=3)
    
    ## Now plot the model parameters
    plot(precis(model,depth=2,pars=c("bg.intens","bg.slope","intens","sig")))
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
##trial <- trial+10

## Make a data frame of the portion of the spectrum we're looking at
d <- data.frame(Chn,trial)


## Use input to get starting parameters
cat("\nClick on the peaks!\n\n")
guess <- locator(npeaks)
g.cent <- guess$x
g.intens <- guess$y

g.bg.intens <- max(1,min(trial))
g.bg.slope <- 0 # 0
g.sig <- rep(5,npeaks) # 5

g.bg.intens.sd <- 1 # 0.01
g.bg.slope.sd <- 0.1 # 0.1

g.intens.sd <- rep(1,npeaks) # 0.01
g.cent.sd <- rep(0.5,npeaks) # 0.5
g.sig.sd <- rep(1,npeaks) # 0.1

## Fit the model using MAP
mod <- quap(
    alist(
        trial ~ dnorm( lambda),
        lambda <- mcmc.BG(Chn,bg.intens,bg.slope)+
                       mcmc.nGaus(Chn,n=npeaks,intens=intens,cent=cent,sig=sig),
        intens ~ dnorm(g.intens, g.intens.sd),
        cent ~ dnorm(g.cent, g.cent.sd),
        sig ~ dnorm(g.sig, g.sig.sd),
        bg.intens ~ dnorm(g.bg.intens, g.bg.intens.sd),
        bg.slope ~ dnorm(g.bg.slope, g.bg.slope.sd)
    ),
    data=d,
    start=list(bg.intens=g.bg.intens,bg.slope=g.bg.slope,
               intens=g.intens,cent=g.cent,sig=g.sig)
)
##g.sig=1.95
##mod <- map(
##    alist(
##        trial ~ dpois( mcmc.BG(Chn,bg.intens,bg.slope)+
##                       mcmc.nGaus(Chn,n=npeaks,intens=intens,cent=cent,sig=sig)),
##        intens ~ dnorm(g.intens,20),
##        cent ~ dnorm(g.cent,10),
##        sig ~ dnorm(g.sig,0.06),
##        bg.intens ~ dnorm(g.bg.intens,10),
##        bg.slope ~ dnorm(0,10)
##    ),
##    data=d,start=list(bg.intens=g.bg.intens,bg.slope=g.bg.slope,
##                      intens=g.intens,cent=g.cent,sig=g.sig)
##)


post <- extract.samples(mod)

areafun <- function(i,s)sqrt(2*pi)*i*abs(s)
asamp <- areafun(post$intens,post$sig)
if(is.null(dim(asamp))){
    a.mean <- mean(asamp)
    a.sd <- sd(asamp)
    sig.mean <- mean(post$sig)
    sig.sd <- sd(post$sig)
    centroid.mean <- mean(post$cent)
    centroid.sd <- sd(post$cent)
}else{
    a.mean <- apply(asamp,2,mean)
    a.sd <- apply(asamp,2,sd)
    sig.mean <- apply(post$sig,2,mean)
    sig.sd <- apply(post$sig,2,sd)
    centroid.mean <- apply(post$cent,2,mean)
    centroid.sd <- apply(post$cent,2,sd)
}
areas <- data.frame(mean=a.mean,sd=a.sd)

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
mod_df$prior_mean <- prior_mean
mod_df$prior_sd <- prior_sd

cat("-----------------------------------\n")
cat("Posteriors and Priors [dnorm(mean, sd)]:\n")
print(mod_df)
cat("\nGaussian areas:\n")
print(areas)
cat("-----------------------------------\n")

plotfit(mod,Chn,trial,areas,npeaks=npeaks,TRUE)

############################################################
## Will's Addition: Append areas and centroids to text file
############################################################

## Ask to append fit data to text file
#
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
#   ## Ask user if they want to append a header line before the data is appended (e.g. # 15 deg, run 88)
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

    fitdata <- data.frame(area=a.mean,sdarea=a.sd,sig=sig.mean,sdsig=sig.sd,cent=centroid.mean,sdcent=centroid.sd)

    ## Ask user for cell input. Which cell should the data start from (each successive fit parameter will be pasted in the next column, same row)
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
