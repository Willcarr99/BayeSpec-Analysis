#### Spectrum analysis program

##----------------------------------------------------------------
## Definition of functions to be used
source("Funcs.R")

source("UserFuncs.R")

##source("MCMCGaus.R")

source("GUI.R")

##----------------------------------------------------------------
## Definition of global control parameters

## The correction fit to account for quadratic focal plane
calfile <- "../MasterCal/FullCalCorr.cal"
#calfile <- "../MasterCal/CarbonCal.cal"

## The uncertainty to use for the background
BackgroundUncert="Gaussian"
## Bayes prior normalisation power
BayesNorm=0.5

## Confidence level for upper limits
ConfUL <- 0.9

## Plot scale
logscale <- ""

## Peak finder sensitivity
peakSensitivity=300


##----------------------------------------------------------------
## Graphics parameters
newpar <- list(bg="black",    # black background
               fg="white",   # White plot colour
               col="white",
               col.axis="white",
               col.lab="white",
               col.main="white",
               col.sub="white",
               cex.axis=1.2,  # axis annotation scale
               cex.lab=1.3
               )

##while( TRUE ) Sys.sleep(1)
