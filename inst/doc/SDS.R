## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#", out.width='//textwidth')
system.file(package="TestGardener")

## ----setup--------------------------------------------------------------------
library(TestGardener)

## -----------------------------------------------------------------------------
titlestr  <- "Symptom Distress"

## -----------------------------------------------------------------------------
U         <- scan("SDS.txt","o")
# U         <- scan(paste(getwd(),"/data/SDS.txt",sep=""),"o")
U         <- matrix(U,473,2,byrow=TRUE)
U         <- U[,2]
N         <- length(U) # Number of examinees
Umat      <- as.integer(unlist(stringr::str_split(U,"")))
n         <- length(Umat)/N # Number of items
U         <- matrix(Umat,N,n,byrow=TRUE)

## -----------------------------------------------------------------------------
key <- NULL

## -----------------------------------------------------------------------------
noption <- matrix(5,n,1)
for (i in 1:n)
{
  if (any(U[,i] > noption[i]))
  {
    noption[i]  <- noption[i] + 1 # Add one option for invalid responses
    U[U[,i] >= noption[i],i] <- noption[i]
  }
}

## -----------------------------------------------------------------------------
optScore <- list() # option scores
for (item in 1:n){
  scorei <- c(0:4,0)
  optScore[[item]] <- scorei
}

## -----------------------------------------------------------------------------
itemVec <- c("Inability to sleep", "Fatigue", "Bowel symptoms", "Breathing symptoms",
             "Coughing", "Inability to concentrate", "Intensity of nausea",
             "Frequency of nausea", "Intensity of pain", "Frequency of pain",
             "Bad outlook on life", "Loss of appetite", "Poor appearance")

## -----------------------------------------------------------------------------
optVec <-c("0", "1", "2", "3", "4", " ")
optLab <- list()
for (i in 1:n)
{
  optLab[[i]] = optVec
}

## -----------------------------------------------------------------------------
optList <- list(itemLab=itemVec, optLab=optLab, optScr=optScore)

## -----------------------------------------------------------------------------
scrrng = c(0,37)

## -----------------------------------------------------------------------------
SDS_dataList <- TestGardener::make.dataList(U, key, optList, scrrng=scrrng)

## -----------------------------------------------------------------------------
hist(SDS_dataList$scrvec, SDS_dataList$scrrng[2], xlab="Sum Score",
     main=titlestr)

## -----------------------------------------------------------------------------
theta     <- SDS_dataList$percntrnk
thetaQnt  <- SDS_dataList$thetaQnt
WfdResult <- TestGardener::Wbinsmth(theta, SDS_dataList)

## -----------------------------------------------------------------------------
WfdList <- WfdResult$WfdList
binctr  <- WfdResult$aves
Qvec    <- c(5,25,50,75,95)
indfine <- seq(0,100,len=101)
TestGardener::ICC.plot(indfine, WfdList, SDS_dataList, Qvec, binctr,  Wrng=c(0,3), plotindex=1)

## -----------------------------------------------------------------------------
ncycle <- 10

## -----------------------------------------------------------------------------
AnalyzeResult <- TestGardener::Analyze(theta, thetaQnt, SDS_dataList, ncycle, itdisp=FALSE) 

## -----------------------------------------------------------------------------
parList  <- AnalyzeResult$parList
meanHsave <- AnalyzeResult$meanHsave

## -----------------------------------------------------------------------------
cycleno <- 1:ncycle
plot(cycleno,meanHsave[cycleno], type="b", lwd=2, xlab="Cycle Number")

## -----------------------------------------------------------------------------
icycle <- 10
SDS_parList  <- parList[[icycle]]

## -----------------------------------------------------------------------------
WfdList    <- SDS_parList$WfdList
theta      <- SDS_parList$theta
Qvec       <- SDS_parList$Qvec
binctr     <- SDS_parList$binctr

## -----------------------------------------------------------------------------
SMS_infoList <- TestGardener::theta2arclen(theta, Qvec, WfdList, binctr)

## -----------------------------------------------------------------------------
arclength     <- SMS_infoList$arclength
arclengthvec  <- SMS_infoList$arclengthvec
arclengthfd   <- SMS_infoList$arclengthfd
theta_al      <- SMS_infoList$theta_al
thetafine_al  <- SMS_infoList$thetafine.rng
Qvec_al       <- SMS_infoList$Qvec_al
binctr_al     <- SMS_infoList$binctr_al
Wfd_info      <- SMS_infoList$Wfd_info
Wdim          <- SMS_infoList$Wdim

## -----------------------------------------------------------------------------
TestGardener::ICC.plot(indfine, WfdList, SDS_dataList, Qvec, binctr,  Wrng=c(0,3), plotindex=1)

## -----------------------------------------------------------------------------
theta_in <- theta[theta > 0 & theta < 100]
dens.plot <- TestGardener::scoreDensity(theta_in, ttlstr=titlestr)
print(dens.plot)

## -----------------------------------------------------------------------------
indfine <- seq(0,100,len=101)
mufine  <- TestGardener::testscore(indfine, WfdList, optList)
TestGardener::mu.plot(mufine, SDS_dataList$scrrng, titlestr)

## -----------------------------------------------------------------------------
print(paste("Arc length =", round(arclength,2)))
TestGardener::ArcLength.plot(arclength, arclengthvec, titlestr)

## -----------------------------------------------------------------------------
Result <- TestGardener::Wpca.plot(arclength, WfdList, SDS_dataList$Wdim, titlestr=titlestr)

## -----------------------------------------------------------------------------
TestGardener::Sensitivity.plot(arclengthvec, WfdList, Qvec, SDS_dataList, 
                               titlestr=titlestr, plotindex=8)

## -----------------------------------------------------------------------------
Result <- TestGardener::Power.plot(arclengthvec, WfdList, Qvec_al, SDS_dataList, 
                                   plotindex=9, height=0.3)

## -----------------------------------------------------------------------------
Result <- TestGardener::Power.plot(arclengthvec, WfdList, Qvec_al, SDS_dataList, 
                                   plotindex=8, height=0.3)

## -----------------------------------------------------------------------------
TestGardener::Hfuns.plot(theta, WfdList, U, plotindex=1:5)

