## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#", out.width='//textwidth')
system.file(package="TestGardener")

## -----------------------------------------------------------------------------
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
chartList <- SDS_dataList$chartList
WfdResult <- TestGardener::Wbinsmth(theta, SDS_dataList, thetaQnt, chartList)

## -----------------------------------------------------------------------------
WfdList <- WfdResult$WfdList
binctr  <- WfdResult$aves
Qvec    <- c(5,25,50,75,95)
TestGardener::Wbinsmth.plot(binctr, Qvec, WfdList, SDS_dataList, Wrng=c(0,3), plotindex=1)

## -----------------------------------------------------------------------------
ncycle <- 10

## -----------------------------------------------------------------------------
AnalyzeResult <- TestGardener::Analyze(theta, thetaQnt, SDS_dataList, ncycle=ncycle, itdisp=TRUE) 

## -----------------------------------------------------------------------------
parList  <- AnalyzeResult$parList
meanHvec <- AnalyzeResult$meanHvec

## -----------------------------------------------------------------------------
cycleno <- 1:ncycle
plot(cycleno,meanHvec[cycleno], type="b", lwd=2, xlab="Cycle Number")

## -----------------------------------------------------------------------------
icycle <- 10
SDS_parListi  <- parList[[icycle]]

## -----------------------------------------------------------------------------
WfdList    <- SDS_parListi$WfdList
theta      <- SDS_parListi$theta
Qvec       <- SDS_parListi$Qvec
binctr     <- SDS_parListi$binctr
arclength  <- SDS_parListi$arclength
alfine     <- SDS_parListi$alfine

## -----------------------------------------------------------------------------
TestGardener::Wbinsmth.plot(binctr, Qvec, WfdList, SDS_dataList, Wrng=c(0,3), plotindex=8)

## -----------------------------------------------------------------------------
ttllab     <- paste(titlestr,": percent rank", sep="")
scrrng     <- c(0,100)
theta_in   <- theta[theta > 0 & theta < 100]
indden10   <- TestGardener::scoreDensity(theta_in, scrrng, ttlstr=ttllab)

## -----------------------------------------------------------------------------
mu <- testscore(theta, WfdList, optList)
ttllab <- paste(titlestr,": expected score", sep="")
muden  <- TestGardener::scoreDensity(mu, SDS_dataList$scrrng, ttlstr=ttllab) 

## -----------------------------------------------------------------------------
print(paste("Arc length =", round(arclength,2)))
TestGardener::ArcLength.plot(arclength, alfine, titlestr)

## -----------------------------------------------------------------------------
Result <- TestGardener::Wpca.plot(arclength, WfdList, SDS_dataList$Wdim, titlestr=titlestr)

## -----------------------------------------------------------------------------
TestGardener::Sensitivity.plot(WfdList, Qvec, SDS_dataList, titlestr=titlestr, plotindex=8)

## -----------------------------------------------------------------------------
Result <- TestGardener::Power.plot(WfdList, Qvec, SDS_dataList, plotindex=9, height=0.3)

## -----------------------------------------------------------------------------
Result <- TestGardener::Power.plot(WfdList, Qvec, SDS_dataList, plotindex=8, height=0.3)

## -----------------------------------------------------------------------------
TestGardener::Hfuns.plot(theta, WfdList, U, plotindex=1:5)

