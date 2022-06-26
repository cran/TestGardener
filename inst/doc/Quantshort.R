## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#", out.width='//textwidth')
system.file(package="TestGardener")

## ----setup--------------------------------------------------------------------
library(TestGardener)

## -----------------------------------------------------------------------------
titlestr  <- "SweSAT-Q: 80 math analysis items"

## -----------------------------------------------------------------------------
U <- scan("U_5000.txt","o")
N <- length(U) # Number of examinees

## -----------------------------------------------------------------------------
Uvector <- as.integer(unlist(stringr::str_split(U,"")))
n       <- length(Uvector)/N # Number of items
U       <- matrix(Uvector,N,n,byrow=TRUE)

## -----------------------------------------------------------------------------
key   <- scan("key_5000.txt","o")
key <- as.integer(unlist(stringr::str_split(key,"")))

## -----------------------------------------------------------------------------
noption <- rep(4,n)
noption[23:28] <- 5
noption[63:68] <- 5
for (i in 1:n) {
  if (any(U[,i] > noption[i])) {
    noption[i]  <- noption[i] + 1 # Add one option for invalid responses
    U[U[,i] > noption[i],i] <- noption[i]
  }
}

## -----------------------------------------------------------------------------
grbg <- noption 

## -----------------------------------------------------------------------------
ScoreList <- list() # option scores
for (item in 1:n){
  scorei <- rep(0,noption[item])
  scorei[key[item]] <- 1
  ScoreList[[item]] <- scorei
}

## -----------------------------------------------------------------------------
optList <- list(itemLab=NULL, optLab=NULL, optScr=ScoreList)

## -----------------------------------------------------------------------------
Math_dataList <- TestGardener::make.dataList(U, key, optList, grbg)
names(Math_dataList)

## -----------------------------------------------------------------------------
hist(Math_dataList$scrvec, Math_dataList$scrrng[2], xlab="Sum Score",
     main=titlestr)

## -----------------------------------------------------------------------------
theta     <- Math_dataList$percntrnk
thetaQnt  <- Math_dataList$thetaQnt
WfdResult <- TestGardener::Wbinsmth(theta, Math_dataList)

## -----------------------------------------------------------------------------
WfdList <- WfdResult$WfdList
binctr  <- WfdResult$aves
Qvec    <- c(5,25,50,75,95)
indfine <- seq(0,100,len=101)
TestGardener::ICC.plot(indfine, WfdList, Math_dataList, Qvec, binctr,  
                       Wrng=c(0,3), plotindex=1)


## -----------------------------------------------------------------------------
ncycle <- 10

## -----------------------------------------------------------------------------
AnalyzeResult <- TestGardener::Analyze(theta, thetaQnt, Math_dataList, 
                                       ncycle, itdisp=FALSE, verbose=FALSE) 

## -----------------------------------------------------------------------------
parList      <- AnalyzeResult$parList
meanHsave    <- AnalyzeResult$meanHsave
arclengthvec <- AnalyzeResult$arclengthsave

## -----------------------------------------------------------------------------
cycleno <- 1:ncycle
plot(cycleno,    meanHsave[cycleno], type="b", lwd=2, xlab="Cycle Number", ylab="Mean Fit")

## -----------------------------------------------------------------------------
cycleno <- 1:ncycle
plot(cycleno, arclengthvec[cycleno], type="b", lwd=2, xlab="Cycle Number", ylab="Arclength")

## -----------------------------------------------------------------------------
icycle <- ncycle
Math_parListi  <- parList[[icycle]]

## -----------------------------------------------------------------------------
WfdList    <- Math_parListi$WfdList
theta      <- Math_parListi$theta
Qvec       <- Math_parListi$Qvec
binctr     <- Math_parListi$binctr

## -----------------------------------------------------------------------------
Math_infoList <- TestGardener::theta2arclen(theta, Qvec, WfdList, binctr)

## -----------------------------------------------------------------------------
arclength     <- Math_infoList$arclength
arclengthvec  <- Math_infoList$arclengthvec
arclengthfd   <- Math_infoList$arclengthfd
theta_al      <- Math_infoList$theta_al
thetafine_al  <- Math_infoList$thetafine.rng
Qvec_al       <- Math_infoList$Qvec_al
binctr_al     <- Math_infoList$binctr_al
Wfd_info      <- Math_infoList$Wfd_info
Wdim          <- Math_infoList$Wdim

## -----------------------------------------------------------------------------
indfine <- seq(0,100,len=101)
TestGardener::ICC.plot(indfine, WfdList, Math_dataList, Qvec, binctr,  
                       data_point=TRUE, plotType=c("P", "W"), Wrng=c(0,3), plotindex=1)

## -----------------------------------------------------------------------------
indfine <- seq(0,100,len=101)
TestGardener::ICC.plot(arclengthvec, WfdList, Math_dataList, Qvec_al, binctr_al,  
                       data_point=TRUE, plotType=c("P", "W"), Wrng=c(0,3), plotindex=1)

## -----------------------------------------------------------------------------
TestGardener::density_plot(theta, c(0,100), Qvec, xlabstr="Score index", 
                           titlestr="Theta Density", scrnbasis=15)

## -----------------------------------------------------------------------------
TestGardener::density_plot(theta_al, c(0,arclength), Qvec_al, xlabstr="Arclength", 
                           titlestr="Arc length Density",  scrnbasis=15)

## -----------------------------------------------------------------------------
indfine <- seq(0,100,len=101)
mufine <- TestGardener::testscore(indfine, WfdList, optList)
TestGardener::mu.plot(mufine, Math_dataList$scrrng, titlestr)

## -----------------------------------------------------------------------------
print(paste("Arc length =", round(arclength,2)))
TestGardener::ArcLength.plot(arclength, arclengthvec, titlestr)

## -----------------------------------------------------------------------------
Result <- TestGardener::Wpca.plot(arclength, WfdList, Math_dataList$Wdim, titlestr=titlestr)

## -----------------------------------------------------------------------------
TestGardener::Hfuns.plot(theta, WfdList, U, plotindex=1)

