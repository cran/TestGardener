## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#", out.width='//textwidth')
system.file(package="TestGardener")

## -----------------------------------------------------------------------------
titlestr  <- "SweSAT-Q: 24 math analysis items"

## -----------------------------------------------------------------------------
U     <- scan("Ushort.txt","o")
# U     <- scan(paste(getwd(),"/data/Ushort.txt",sep=""),"o")
N <- length(U) # Number of examinees

## -----------------------------------------------------------------------------
Uvector <- as.integer(unlist(stringr::str_split(U,"")))
n       <- length(Uvector)/N # Number of items
U       <- matrix(Uvector,N,n,byrow=TRUE)

## -----------------------------------------------------------------------------
key   <- scan("keyshort.txt","o")
# key   <- scan(paste(getwd(),"/data/keyshort.txt",sep=""),"o")
key <- as.integer(unlist(stringr::str_split(key,"")))

## -----------------------------------------------------------------------------
noption <- rep(4,n)
for (j in 1:n) {
  if (any(U[,j] > noption[j])) {
    noption[j]  <- noption[j] + 1 # Add one option for invalid responses
    U[U[,j] > noption[j],j] <- noption[j]
  }
}

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
Quant_dataList <- TestGardener::make.dataList(U, key, optList)
names(Quant_dataList)

## -----------------------------------------------------------------------------
hist(Quant_dataList$scrvec, Quant_dataList$scrrng[2], xlab="Sum Score",
     main=titlestr)

## -----------------------------------------------------------------------------
theta     <- Quant_dataList$percntrnk
thetaQnt  <- Quant_dataList$thetaQnt
chartList <- Quant_dataList$chartList
WfdResult <- TestGardener::Wbinsmth(theta, Quant_dataList, thetaQnt, chartList)

## -----------------------------------------------------------------------------
WfdList <- WfdResult$WfdList
binctr  <- WfdResult$aves
Qvec    <- c(5,25,50,75,95)
TestGardener::Wbinsmth.plot(binctr, Qvec, WfdList, Quant_dataList, Wrng=c(0,3), plotindex=1)

## -----------------------------------------------------------------------------
ncycle <- 10

## -----------------------------------------------------------------------------
AnalyzeResult <- TestGardener::Analyze(theta, thetaQnt, Quant_dataList, ncycle=ncycle, itdisp=FALSE) 

## -----------------------------------------------------------------------------
parList  <- AnalyzeResult$parList
meanHvec <- AnalyzeResult$meanHvec

## -----------------------------------------------------------------------------
cycleno <- 1:ncycle
plot(cycleno,meanHvec[cycleno], type="b", lwd=2, xlab="Cycle Number")

## -----------------------------------------------------------------------------
icycle <- 10
Quant_parListi  <- parList[[icycle]]

## -----------------------------------------------------------------------------
WfdList    <- Quant_parListi$WfdList
theta      <- Quant_parListi$theta
Qvec       <- Quant_parListi$Qvec
binctr     <- Quant_parListi$binctr
arclength  <- Quant_parListi$arclength
alfine     <- Quant_parListi$alfine

## -----------------------------------------------------------------------------
TestGardener::Wbinsmth.plot(binctr, Qvec, WfdList, Quant_dataList, Wrng=c(0,3), plotindex=1)

## -----------------------------------------------------------------------------
ttllab     <- paste(titlestr,": percent rank", sep="")
scrrng     <- c(0,100)
theta_in   <- theta[theta > 0 & theta < 100]
indden10   <- TestGardener::scoreDensity(theta_in, scrrng, ttlstr=ttllab)

## -----------------------------------------------------------------------------
mu <- TestGardener::testscore(theta, WfdList, optList)
ttllab <- paste(titlestr,": expected score", sep="")
muden  <- TestGardener::scoreDensity(mu, Quant_dataList$scrrng, ttlstr=ttllab) 
print(muden)

## -----------------------------------------------------------------------------
indfine <- seq(0,100,len=101)
mufine <- TestGardener::testscore(indfine, WfdList, optList)
TestGardener::mu.plot(mufine, Quant_dataList$scrrng, titlestr)

## -----------------------------------------------------------------------------
print(paste("Arc length =", round(arclength,2)))
TestGardener::ArcLength.plot(arclength, alfine, titlestr)

## -----------------------------------------------------------------------------
Result <- TestGardener::Wpca.plot(arclength, WfdList, Quant_dataList$Wdim, titlestr=titlestr)

## -----------------------------------------------------------------------------
TestGardener::Sensitivity.plot(WfdList, Qvec, Quant_dataList, plotindex=1)

## -----------------------------------------------------------------------------
Result <- TestGardener::Power.plot(WfdList, Qvec, Quant_dataList, plotindex=1, height=0.25)

## -----------------------------------------------------------------------------
Result <- TestGardener::Power.plot(WfdList, Qvec, Quant_dataList, plotindex=21, height=0.25)

## -----------------------------------------------------------------------------
TestGardener::Wbinsmth.plot(binctr, Qvec, WfdList, Quant_dataList, Wrng=c(0,3), plotindex=21)

## -----------------------------------------------------------------------------
TestGardener::Sensitivity.plot(WfdList, Qvec, Quant_dataList, titlestr=titlestr, plotindex=21)

## -----------------------------------------------------------------------------
TestGardener::Hfuns.plot(theta, WfdList, U, plotindex=1)

