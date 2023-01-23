## Use what I've learned to fit Gaussians
## 
library(rethinking)
library(prettyGraphs)

cat("\n",
    "To MCMC fit multiple peaks: MCMCGaus(n=2)","\n",
    "\n",sep="")

## Function to generate n Gaussians
mcmc.nGaus <- function(x,npeaks,intens,cent,sig){
    g <- sapply(1:npeaks, function(i){
        c <- cent[i]
        s <- sig  ##[i]
        I <- abs(intens[i])
        I*exp(-0.5*((x-c)^2)/s^2)
    })
##    lines(x,rowSums(g),col=add.alpha(rangi2,0.2))
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
                          mcmc.nGaus(Chn,n=n,post$intens[i],post$cent[i],post$sig),#[i]),
                      col=add.alpha(rangi2,0.2))
            }
        } else {
            ## Draw those samples
            for(i in 1:dim(post$intens)[1]){
                for(j in 1:n)
                    lines(Chn,mcmc.nGaus(Chn,n=1,post$intens[i,j],post$cent[i,j],post$sig[i]),##,j]),
                          col=add.alpha("red",0.1))
                lines(Chn,mcmc.BG(Chn,post$bg.intens[i],post$bg.slope[i])+
                          mcmc.nGaus(Chn,n=n,post$intens[i,],post$cent[i,],post$sig[i]), ##,]),
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
                      mcmc.nGaus(Chn,n=n,intens,cent,sig),#[i]),
                  col=add.alpha(rangi2,1))
            ##           }
        } else {
            ## Draw those samples
            ##for(i in 1:dim(post$intens)[1]){
            bg.intens <- median(post$bg.intens) 
            bg.slope <- median(post$bg.slope)
            intens <- apply(post$intens,2,median)
            cent <- apply(post$cent,2,median)
            sig <- median(post$sig)

            ##print(bg.intens,bg.slope,intens,cent,sig)
            
            for(j in 1:n)
                lines(Chn,mcmc.nGaus(Chn,n=1,intens[j],cent[j],sig),##,j]),
                      col=add.alpha("red",1))
            lines(Chn,mcmc.BG(Chn,bg.intens,bg.slope)+
                      mcmc.nGaus(Chn,n=n,intens,cent,sig), ##,]),
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
g.bg.slope <- 0
g.sig <- 5  ##rep(10,npeaks)
##cat("\nYour guess was\n")
##print(guess)
                                        #npeaks <- npeaks


## Fit the model using MAP
mod <- quap(
    alist(
        trial ~ dnorm( lambda),
        lambda <- mcmc.BG(Chn,bg.intens,bg.slope)+
                       mcmc.nGaus(Chn,n=npeaks,intens=intens,cent=cent,sig=sig),
        intens ~ dnorm(g.intens,20),
        cent ~ dnorm(g.cent,5),
        sig ~ dnorm(g.sig,2), #2
        bg.intens ~ dnorm(g.bg.intens,10),
        bg.slope ~ dnorm(0,0.1)
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
