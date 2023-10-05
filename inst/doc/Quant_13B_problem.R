## ---- chunk=1,, echo = FALSE--------------------------------------------------
knitr::opts_chunk$set(collapse=TRUE, comment="#", out.width='//textwidth')
system.file(package="TestGardener")

## ----define title, chunk=2----------------------------------------------------
titlestr  <- "SweSAT-Q: 24 math analysis items and 1000 examinees"

## ----read U, chunk=3----------------------------------------------------------
U <- as.matrix(read.table("Quant_13B_problem_U.txt"))
N <- nrow(U) # Number of examinees 
n <- ncol(U) # Number of items

## ----read key, chunk=4--------------------------------------------------------
key   <- scan("Quant_13B_problem_key.txt","o")
key <- as.numeric(unlist(stringr::str_split(key,"")))

## ----define numbers of options, chunk=5---------------------------------------
noption <- rep(5,n)

## ----define option weights, chunk=6-------------------------------------------
ScoreList <- list() # option scores
for (item in 1:n){
  scorei <- rep(0,noption[item])
  scorei[key[item]] <- 1
  ScoreList[[item]] <- scorei
}

## ----define option labels, chunk=7--------------------------------------------
optList <- list(itemLab=NULL, optLab=NULL, optScr=ScoreList)

## ----make the dataList object, chunk=8----------------------------------------
Math_dataList <- TestGardener::make.dataList(U, key, optList)
names(Math_dataList)

## ----plot a histogram of sum scores, chunk=9, fig.width = 7-------------------
hist(Math_dataList$scrvec, Math_dataList$scrrng[2], xlab="Sum Score",
     main=titlestr)

## ----initialize the analysis, chunk=10----------------------------------------
theta     <- Math_dataList$percntrnk
thetaQnt  <- Math_dataList$thetaQnt
WfdResult <- TestGardener::Wbinsmth(theta, Math_dataList)

## ----plot initial curves, chunk=11, fig.width = 7-----------------------------
WfdList <- WfdResult$WfdList
binctr  <- WfdResult$aves
Qvec    <- c(5,25,50,75,95)
indfine <- seq(0,100,len=101)
TestGardener::ICC.plot(indfine, WfdList, Math_dataList, Qvec, binctr,  
                       Wrng=c(0,5), plotindex=1, plotType="P")

## ----set number of cycles, chunk=12-------------------------------------------
ncycle <- 10

## ----chunk=13-----------------------------------------------------------------
AnalyzeResult <- TestGardener::Analyze(theta, thetaQnt, Math_dataList, 
                                       ncycle, itdisp=FALSE, verbose=FALSE) 

## ----set up cycle progress for mean fit and arc length, chunk=14--------------
parList      <- AnalyzeResult$parList
HCycle <- matrix(0,ncycle,2)
for (icycle in 1:ncycle) {
    HCycle[icycle,1] <- parList[[icycle]]$meanH
    HCycle[icycle,2] <- parList[[icycle]]$arclength
}

## ----plot the cycle progress of the mean data fit, chunk=15, fig.width = 7----
plot(1:ncycle, HCycle[,1], type="b", lwd=2,
     xlab="Cycle Number", ylab="Mean H")

## ----plot the cycle progress of arclength, chunk=16, fig.width = 7------------
plot(1:ncycle, HCycle[,2], type="b", lwd=2, 
     xlab="Cycle Number", ylab="Arc Length")

## ----define parList, chunk=17-------------------------------------------------
Math_parListi  <- parList[[ncycle]]

## ----output parList objects, chunk=18-----------------------------------------
WfdList    <- Math_parListi$WfdList
theta      <- Math_parListi$theta
Qvec       <- Math_parListi$Qvec
binctr     <- Math_parListi$binctr
indfine    <- seq(0,100,len=101)

## ----define infoList, chunk=19------------------------------------------------
Math_infoList <- TestGardener::theta2arclen(theta, Qvec, WfdList, binctr)

## ----output objects in infoList, chunk=20-------------------------------------
arclength     <- Math_infoList$arclength
arclengthvec  <- Math_infoList$arclengthvec
WfdList_al    <- Math_infoList$arclength
theta_al      <- Math_infoList$theta_al
Qvec_al       <- Math_infoList$Qvec_al
binctr_al     <- Math_infoList$binctr_al

## ----plot score index density, chunk=21, fig.width = 7------------------------
TestGardener::density_plot(theta, c(0,100), Qvec, xlabstr="Score index", 
                           titlestr="Theta Density", scrnbasis=15)

## ----plot arclength density,chunk=22,fig.width = 7----------------------------
TestGardener::density_plot(theta_al, c(0,arclength), Qvec_al, xlabstr="Arclength", 
                           titlestr="Arc length Density",  scrnbasis=15)

## ----plot curves over score index, chunk=23,  fig.width = 7, fig.height = 7----
TestGardener::ICC.plot(indfine, WfdList, Math_dataList, Qvec, binctr,  
                       data_point=TRUE, plotType=c("P", "W"), Wrng=c(0,3), plotindex=1)

## ----plot curves as functions of arclength, chunk=24, fig.width = 7, fig.height = 5----
TestGardener::ICC.plot(arclengthvec, WfdList, Math_dataList, Qvec_al, binctr_al,  
                       data_point=TRUE, plotType=c("P", "W"), Wrng=c(0,5),            
                       plotindex=1)

## ----plot expected score index against arc length, chunk=25-------------------
indfine <- seq(0,100,len=101)
mufine <- TestGardener::testscore(indfine, WfdList, optList)
TestGardener::mu.plot(mufine, Math_dataList$scrrng, titlestr)

## ----plot arc length against score index, chunk=26, fig.width = 7-------------
print(paste("Arc length =", round(arclength,2)))
TestGardener::ArcLength.plot(arclength, arclengthvec, titlestr)

## ----compute and display test information curve, chunk=27, eval=TRUE, fig.width = 7, webgl=TRUE----
Result <- TestGardener::Wpca.plot(WfdList, nharm=3,
                                  rotate=FALSE, titlestr = titlestr)
print(Result$pcaplt)
print("Percentagess of variance for principal components:")
print(round(100*Result$varpropvarmx,1))
print("Total percentage of variation:")
print(round(sum(100*Result$varpropvarmx,1)))


## ----plot the first five examinee data fits, chunk=28, eval = TRUE, fig.width = 7----
TestGardener::Hfuns.plot(indfine, theta, WfdList, U, plotindex=1:5)

