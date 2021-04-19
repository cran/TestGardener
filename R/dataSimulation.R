dataSimulation <- function(dataList, parList, nsample=1000) {
  #  Simulate data, analyze each simulated sample, and display results
  
  #  Last modified on 8 February 2021
  
  #  info from dataList
  
  optList <- dataList$optList
  Wdim    <- dataList$Wdim
  
  #  info from parList
  
  WfdList   <- parList$WfdList
  logdensfd <- parList$logdensfd
  Qvec      <- parList$Qvec
  n         <- length(WfdList)
  
  ##  Generate simulated data
  
  scrfine <- dataList$scrfine
  
  #  define population theta values from scrfine
  
  theta.pop   <- seq(0,100,len=51)  #  percent rank seed
  ntheta.pop  <- length(theta.pop)
  
  #  compute population values for mu and arc length
  
  mu.pop        <- testscore(theta.pop, WfdList, optList)
  Result        <- theta2arclen(theta.pop, WfdList, Wdim)
  al.pop        <- Result$theta_al 
  arclength     <- Result$arclength 
  arclengthfine <- Result$arclengthfine
  
  al.pop <- al.pop*100/arclength
  
  #  arrays to store simulated data
  
  Usave <- array(0,c(ntheta.pop,n,nsample))
  
  for (isample in 1:nsample) {
    if (round(isample/100)*100 == isample) {
      print(paste("Simulation Sample",isample))
    }
    Usave[,,isample] <- Usimulate(n, theta.pop, WfdList)
  }
  
  #  Temporarily removing the parallel feature until a later stage in debugging for CRAN check
  
  # #----------------------------------------------------------
  # #        Run this in parallel
  # #---------------------------------------------------------
  # 
  # library(foreach)
  # library(doParallel)
  # library(parallel)
  # 
  # #setup parallel backend to use many processors
  # cores=detectCores()
  # cl <- makeCluster(cores[1]-1) #not to overload your computer
  # registerDoParallel(cl)
  # 
  # start_time <- Sys.time()
  # finalMatrix <- foreach(isample=1:nsample, .combine=rbind) %dopar% {
  #   tempMatrix <- Usimulate(n, indfine, theta.pop, WfdList)
  #   
  #   tempMatrix #Equivalent to finalMatrix <- cbind(finalMatrix, tempMatrix)
  # }
  # end_time <- Sys.time()
  # end_time - start_time
  # 
  # #stop cluster
  # stopCluster(cl)
  
  ##  Analyze each sample for each cycle
  
  #  set up matrices to save(results
  
  sumscrsave <- matrix(0,ntheta.pop, nsample)
  thetasave  <- matrix(0,ntheta.pop, nsample)
  musave     <- matrix(0,ntheta.pop, nsample)
  alsave     <- matrix(0,ntheta.pop, nsample)
  
  Hvalsave   <- matrix(0,ntheta.pop, nsample)
  DHvalsave  <- matrix(0,ntheta.pop, nsample)
  D2Hvalsave <- matrix(0,ntheta.pop, nsample)
  
  alphascrsave <- matrix(0,nsample)  #  Cronbach's alpha for sum score
  alphamusave  <- matrix(0,nsample)  #  Cronbach's alpha for sum score
  
  #  analyze the samples
  
  for (isample in 1:nsample) {
    if (round(isample/100)*100 == isample) {
      print(paste("Analysis Sample", isample))
    }
    Umati <- Usave[,,isample]
    dataList$U <- Umati
    #  sum scores for this sample
    scrmat    <- matrix(0,ntheta.pop,n)
    for (j in 1:ntheta.pop) {
      for (i in 1:n) {
        scorei <- optList$optScr[[i]]
        scrmat[j,i]    <- scorei[Umati[j,i]]
      }
    }
    scrveci <- apply(scrmat,1,sum)   
    vartot  <- var(scrveci)
    scrmean <- apply(scrmat,2,mean)
    scrvar  <- apply((scrmat-scrmean)^2,2,mean)
    alphascr <- n*(1-sum(scrvar)/vartot)/(n-1)
    #  optimal scores for multi-option analysis
    Result <- thetascorefn(theta.pop, dataList, WfdList, optList)
    thetaveci  <- Result$theta 
    muveci     <- Result$mu 
    Hvalveci   <- Result$Hval 
    DHvalveci  <- Result$DHval 
    D2Hvalveci <- Result$D2Hval 
    alveci     <- Result$arclen
    alphamu    <- Result$alphamu
    
    sumscrsave[,isample]  <- scrveci
    thetasave[,isample]   <- thetaveci     
    musave[,isample]      <- muveci
    Hvalsave[,isample]    <-   Hvalveci
    DHvalsave[,isample]   <-  DHvalveci
    D2Hvalsave[,isample]  <- D2Hvalveci
    alsave[,isample]      <- alveci
    alphascrsave[isample] <- alphascr
    alphamusave[isample]  <- alphamu
  }
  
  #  List  object simList to save the results
  
  simList <- list(
    sumscr   = sumscrsave,
    theta    = thetasave,
    mu       = musave,
    al       = alsave,
    alphascr = alphascrsave,
    alphamu  = alphamusave,
    thepop   = theta.pop,
    mupop    = mu.pop,
    alpop    = al.pop,
    n        = n,
    Qvec     = Qvec
  )
  
  ##  Set up comparisons between score types using tables
  
  simList = scorePerformance(dataList, simList)
  
  return(simList)
  
}

#  ---------------------------------------------------------------------------

thetascorefn <- function(theta.pop, dataList, WfdList, optList) {
  simList <- thetafun(theta.pop, WfdList, dataList)
  theta   <- simList$theta
  Hval    <- simList$Hval
  DHval   <- simList$DHval
  D2Hval  <- simList$D2Hval
  N <- length(theta.pop)
  n <- length(WfdList)
  mumat <- matrix(0,N,n)
  for (item in 1:n)
  {
    SListi <- WfdList[[item]]
    Sfdi  <- SListi$Wfd
    Mi    <- SListi$M
    if (Mi == 1)
    {
      stop("Mi = 1.  Binary data should use Mi = 2.")
    } else {
      Smati <- eval.surp(theta, Sfdi)
      Pmati <- exp(-Smati*log(Mi))
      scri  <- matrix(optList$optScr[[item]], N, Mi, byrow=TRUE)
      mumat[,item] <- rowSums(scri * Pmati)
    }
  }
  # mu <- testscore(theta, WfdList, optList)
  
  Wdim      <- sum(dataList$noption)
  alList    <- theta2arclen(theta.pop, WfdList, Wdim)
  theta_al  <- alList$theta_al
  arclength <- alList$arclength
  arclen    <- theta_al*100/arclength
  
  mu      <- apply(mumat,1,sum)   
  vartot  <- var(mu)
  mumean  <- apply(mumat,2,mean)
  muvar   <- apply((mumat-mumean)^2,2,mean)
  alphamu <- n*(1-sum(muvar)/vartot)/(n-1)
  
  
  return(list(theta = theta, mu = mu, theta_al = theta_al, 
              Hval = Hval, DHval = DHval, D2Hval = D2Hval, 
              arclen = arclen, alphamu = alphamu))
}


