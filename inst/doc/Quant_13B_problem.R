## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#", out.width='//textwidth')
system.file(package="TestGardener")

## -----------------------------------------------------------------------------
titlestr  <- "SweSAT-Q: 24 math analysis items and 1000 examinees"

## -----------------------------------------------------------------------------
U <- as.matrix(read.table("Quant_13B_problem_U.txt"))
N <- nrow(U) # Number of examinees 
n <- ncol(U) # Number of items

## -----------------------------------------------------------------------------
key   <- scan("Quant_13B_problem_key.txt","o")
key <- as.numeric(unlist(stringr::str_split(key,"")))

## -----------------------------------------------------------------------------
noption <- rep(5,n)

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
Math_dataList <- TestGardener::make.dataList(U, key, optList)
names(Math_dataList)

## ----fig.width = 7------------------------------------------------------------
hist(Math_dataList$scrvec, Math_dataList$scrrng[2], xlab="Sum Score",
     main=titlestr)

## -----------------------------------------------------------------------------
theta     <- Math_dataList$percntrnk
thetaQnt  <- Math_dataList$thetaQnt
WfdResult <- TestGardener::Wbinsmth(theta, Math_dataList)

## ----fig.width = 7------------------------------------------------------------
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
HCycle <- matrix(0,ncycle,2)
for (icycle in 1:ncycle) {
    HCycle[icycle,1] <- parList[[icycle]]$meanH
    HCycle[icycle,2] <- parList[[icycle]]$arclength
}

## ----fig.width = 7------------------------------------------------------------
plot(1:ncycle, HCycle[,1], type="b", lwd=2,
     xlab="Cycle Number", ylab="Mean H")

## ----fig.width = 7------------------------------------------------------------
plot(1:ncycle, HCycle[,2], type="b", lwd=2, 
     xlab="Cycle Number", ylab="Arc Length")

## -----------------------------------------------------------------------------
Math_parListi  <- parList[[ncycle]]

## -----------------------------------------------------------------------------
WfdList    <- Math_parListi$WfdList
theta      <- Math_parListi$theta
Qvec       <- Math_parListi$Qvec
Hval       <- Math_parListi$Hval
DHval      <- Math_parListi$DHval
D2Hval     <- Math_parListi$D2Hval
binctr     <- Math_parListi$binctr
indfine    <- seq(0,100,len=101)

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

## ----fig.width = 7------------------------------------------------------------
indfine <- seq(0,100,len=101)
TestGardener::density_plot(theta, c(0,100), Qvec, xlabstr="Score index", 
                           titlestr="Theta Density", scrnbasis=15)

## ----fig.width = 7------------------------------------------------------------
TestGardener::density_plot(theta_al, c(0,arclength), Qvec_al, xlabstr="Arclength", 
                           titlestr="Arc length Density",  scrnbasis=15)

## ----fig.width = 7,fig.height = 5---------------------------------------------
indfine <- seq(0,100,len=101)
TestGardener::ICC.plot(indfine, WfdList, Math_dataList, Qvec, binctr,  
                       data_point=TRUE, plotType=c("P", "W"), Wrng=c(0,3), plotindex=1)

## ----fig.width = 7,fig.height = 5---------------------------------------------
TestGardener::ICC.plot(arclengthvec, WfdList, Math_dataList, Qvec_al, binctr_al,  
                       data_point=TRUE, plotType=c("P", "W"), Wrng=c(0,3), plotindex=1)

## -----------------------------------------------------------------------------
indfine <- seq(0,100,len=101)
mufine <- TestGardener::testscore(indfine, WfdList, optList)
TestGardener::mu.plot(mufine, Math_dataList$scrrng, titlestr)

## ----fig.width = 7------------------------------------------------------------
print(paste("Arc length =", round(arclength,2)))
TestGardener::ArcLength.plot(arclength, arclengthvec, titlestr)

## ---- setup-------------------------------------------------------------------
library(rgl)
options(rgl.useNULL = TRUE) # Suppress the separate window.

## ----fig.width = 7,webgl=TRUE-------------------------------------------------
Result <- TestGardener::Wpca.plot(arclength, WfdList, Math_dataList$Wdim, 3,
                                  rotate=FALSE, titlestr=titlestr)
rglwidget()

## ----fig.width = 7------------------------------------------------------------
TestGardener::Hfuns.plot(indfine, theta, WfdList, U, plotindex=1:5)

