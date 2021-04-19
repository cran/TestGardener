
#  -------------  preliminaries  -----------------

library(TestGardener)
#  set working directory to TestGardener

# ----------- read in data and key ----------- 

titlestr  <- "SweSAT-Q: 24 math analysis items"

U         <- scan("data/Ushort.txt", "o")
N         <- length(U) # Number of examinees
Umat      <- as.integer(unlist(stringr::str_split(U,"")))
n         <- length(Umat)/N # Number of items
U         <- matrix(Umat,N,n,byrow=TRUE)

Quant_U <- U

knostring <- scan("data/keyshort.txt", "o")
key       <- as.integer(unlist(stringr::str_split(knostring,"")))

Quant_key <- key

# Define the coding of invalid responses.
# Any value in U greater than number of legitimate options is treated as missing.

noption <- rep(0,n)
for (i in 1:n)
{
  noption[i]  <- length(unique(U[,i]))
}

# summary count for each option

for (i in 1:n)
{
  print(paste("Item: ",i,sep=""))
  for (m in 1:noption[i])
  {
    print(paste("Option: ",m," N= ",sum(U[,i]==m),sep=""))
  }
}

# --------- Define the option score values for each item ---------

ScoreList <- list() # option scores
for (item in 1:n){
  scorei <- rep(0,noption[item])
  scorei[key[item]] <- 1
  ScoreList[[item]] <- scorei
}

optList <- list(itemLab=NULL, optLab=NULL, optScr=ScoreList)

# ----------------  Initialization Steps  ------------------------

Quant_dataList <- make.dataList(U, key, optList, titlestr=titlestr)

#  plot the sum scores as a histogram 

hist(Quant_dataList$scrvec, Quant_dataList$scrrng[2], xlab="Sum Score",
     main=titlestr)

#  compute the initial option surprisal curves using the 
#  percentage ranks as initial estimates of theta

theta     <- Quant_dataList$percntrnk
thetaQnt  <- Quant_dataList$thetaQnt
chartList <- Quant_dataList$chartList

WfdResult <- Wbinsmth(theta, Quant_dataList, thetaQnt, chartList)

#  Plot the initial option proability and surprisal curves

WfdList <- WfdResult$WfdList
binctr  <- WfdResult$aves
Qvec    <- c(5,25,50,75,95)

Wbinsmth.plot(binctr, Qvec, WfdList, Quant_dataList, Wrng=c(0,3))

# ---------------  Optimal scoring: cycle of smoothing/theta estimation  ------------

#  Set number of cycles and the cell array to containing the parameter

ncycle=10

#  ----------------------------------------------------------------------------
#                      Proceed through the cycles
#  ----------------------------------------------------------------------------

AnalyzeResult <- Analyze(theta, thetaQnt, Quant_dataList, ncycle=ncycle, itdisp=FALSE) 

parList  <- AnalyzeResult$parList
meanHvec <- AnalyzeResult$meanHvec
  
#  ----------------------------------------------------------------------------
#              Plot the average H value, meanHsave, over cycles
#  ----------------------------------------------------------------------------

cycleno <- 1:ncycle
plot(cycleno,meanHvec[cycleno], type="b", lwd=2, xlab="Cycle Number")

#  select cycle for plotting

icycle <- 10

Quant_parListi  <- parList[[icycle]]

WfdList    <- Quant_parListi$WfdList
Qvec       <- Quant_parListi$Qvec
binctr     <- Quant_parListi$binctr
theta      <- Quant_parListi$theta
theta_al   <- Quant_parListi$theta_al
arclength  <- Quant_parListi$arclength
alfine     <- Quant_parListi$alfine

#  ----------------------------------------------------------------------------
#                   Plot surprisal curves for each test question
#  ----------------------------------------------------------------------------

#  plot both the probability and surprisal curves along with data points

Wbinsmth.plot(binctr, Qvec, WfdList, Quant_dataList, Wrng=c(0,3))

#  ----------------------------------------------------------------------------
#                         Plot density of theta
#  ----------------------------------------------------------------------------

ttllab     <- paste(titlestr,": percent rank", sep="")
scrrng     <- c(0,100)
indden10   <- scoreDensity(theta, scrrng, ttlstr=ttllab)

#  ----------------------------------------------------------------------------
#      Compute expected test scores for all examinees
#      Plot expected test scores and expected test score over mesh
#  ----------------------------------------------------------------------------

mu <- testscore(theta, WfdList, optList)
ttllab <- paste(titlestr,": expected score", sep="")
muden  <- scoreDensity(mu, Quant_dataList$scrrng, ttlstr=ttllab) 

#  compute expected score for each value in the fine mesh of theta values

indfine <- seq(0,100,len=101)
mufine <- testscore(indfine, WfdList, optList)

mu.plot(mufine, Quant_dataList$scrrng, titlestr)

#  ----------------------------------------------------------------------------
#         Compute arc length over a fine mesh of theta values and plot
#  ----------------------------------------------------------------------------

#  print length of the test effort curve

print(paste("Arc length =", round(arclength,2)))

#  plot arc length over fine mesh

ArcLength.plot(arclength, alfine, titlestr)

#  ----------------------------------------------------------------------------
#  Display test effort curve projected into its first two principal components
#  ----------------------------------------------------------------------------

# nharm=2
Wpca.plot(arclength, WfdList, Quant_dataList$Wdim, titlestr=titlestr)

# nharm=3
Wpca.plot(arclength, WfdList, Quant_dataList$Wdim, 3, dodge = 1.005, titlestr=titlestr)

#  ----------------------------------------------------------------------------
#                          Display sensitivity curves
#  ----------------------------------------------------------------------------

#  This code needs to put in a legend and indication of right answer if 
#  scoring is multiple choice

Sensitivity.plot(WfdList, Qvec, Quant_dataList)
  
#  ----------------------------------------------------------------------------
#                          Display power curves
#  ----------------------------------------------------------------------------

Power.plot(WfdList, Qvec, Quant_dataList, height=0.25)

#  ----------------------------------------------------------------------------
#                          Display entropy curves
#  ----------------------------------------------------------------------------

Entropy.plot(WfdList, Qvec, Quant_dataList, height=1)

#  ----------------------------------------------------------------------------
#                   Display H, DH and D2H curves for selected examinees
#  ----------------------------------------------------------------------------

index <- 1:5

Hfuns.plot(theta, WfdList, Quant_dataList$U, plotindex=1:5)

#  ----------------------------------------------------------------------------
#             simulate data samples and analyze simulated samples
#  ----------------------------------------------------------------------------

simList <- dataSimulation(Quant_dataList, Quant_parListi, nsample=500)
dataSimulation.plot(simList, Qvec)
