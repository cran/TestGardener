dataSimulation.plot <- function(simList, Qvec,
                                ttlsz=NULL,axisttl=NULL,axistxt=NULL,lgdlab=NULL) {
  
  RMSErng = c(0,ceiling(max(simList$thetaRMSE)))
  biasrng = c(floor(min(simList$thetabias)),ceiling(max(simList$thetabias)))
  
  thetaplot <- RMSEbiasplot.theta(simList$thetaRMSE, simList$thetabias, 
                     RMSErng, biasrng, Qvec,
                     ttlsz,axisttl,axistxt,lgdlab)
  print(thetaplot)
  readline(prompt = ". Press [enter] to see next plot")
  
  RMSErng = c(0,ceiling(max(simList$sumscrRMSE)))
  biasrng = c(floor(min(simList$sumscrbias)),ceiling(max(simList$sumscrbias)))
  
  muplot <- RMSEbiasplot.mu(simList$sumscrRMSE, simList$muRMSE, 
                  simList$sumscrbias, simList$mubias, 
                  RMSErng, biasrng, Qvec,
                  ttlsz,axisttl,axistxt,lgdlab)
  print(muplot)
  return(list(thetaplot=thetaplot, muplot=muplot))
}

#  ---------------------------------------------------------------------------

RMSEbiasplot.theta = function(thetaRMSE, thetabias, RMSErng, biasrng, Qvec,
                              ttlsz=NULL,axisttl=NULL,axistxt=NULL,lgdlab=NULL) {
  
  indfine = seq(0,100,len=101)
  
  default_size=16
  default_size1=12
  default_size2=10
  My_Theme = ggplot2::theme(
    plot.title = ggplot2::element_text(size = ifelse(is.null(ttlsz),default_size,ttlsz)),
    axis.title = ggplot2::element_text(size = ifelse(is.null(axisttl),default_size1,axisttl)),
    axis.text = ggplot2::element_text(size = ifelse(is.null(axistxt),default_size2,axistxt)))
  
  #  ------------  RMSE results  ------------
  df <- data.frame(x=indfine,y=thetaRMSE)
  vec <- sqrt(thetaRMSE^2 - thetabias^2)
  df2 <- data.frame(x=indfine,y=vec)
  pt1 <- ggplot2::ggplot(df, ggplot2::aes(indfine, thetaRMSE, color='Model'))+
    ggplot2::geom_line(size=2)+
    ggplot2::geom_line(data=df2, ggplot2::aes(indfine, vec, color='Variance'), 
                       size=2, linetype = "dashed")+
    ggplot2::scale_color_manual(name="",
                                labels = c("Model", 
                                           "Variance"), 
                                values=c("blue","red"))+
    ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed")+
    ggplot2::xlim(0,100)+
    ggplot2::ylim(RMSErng)+
    ggplot2::ylab("Percentile Index RMSE")+
    ggplot2::xlab("") + 
    ggplot2::theme(legend.position="top")
  pt1 <- pt1 + My_Theme
    
  #  ------------  bias results  ------------
  df <- data.frame(x=indfine,y=thetabias)
  pt2 <- ggplot2::ggplot(df, ggplot2::aes(indfine, thetabias))+
    ggplot2::geom_line(color='blue', size=2)+
    ggplot2::geom_hline(yintercept = 0,  color="black", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed")+
    ggplot2::xlim(0,100)+
    ggplot2::ylim(biasrng)+
    ggplot2::ylab("Percentile Index Bias") +
    ggplot2::xlab("True Percentile Index")
    
  pt2 <- pt2 + My_Theme
  
  p <- ggpubr::ggarrange(pt1, pt2, ncol = 1, nrow = 2)
  return(p)
}

#  ---------------------------------------------------------------------------

RMSEbiasplot.mu <- function(sumscrRMSE, muRMSE, sumscrbias, mubias, 
                            RMSErng, biasrng, Qvec,
                            ttlsz=NULL,axisttl=NULL,axistxt=NULL,lgdlab=NULL) {
  
  indfine <- seq(0,100,len=101)
  indrng  <- c(0,100)
  
  default_size=16
  default_size1=12
  default_size2=10
  My_Theme = ggplot2::theme(
    plot.title = ggplot2::element_text(size = ifelse(is.null(ttlsz),default_size,ttlsz)),
    axis.title = ggplot2::element_text(size = ifelse(is.null(axisttl),default_size1,axisttl)),
    axis.text = ggplot2::element_text(size = ifelse(is.null(axistxt),default_size2,axistxt)))
  
  #  ------------  RMSE results  ------------
  df <- data.frame(x=indfine,y=muRMSE)
  df2 <- data.frame(x=indfine,y=sumscrRMSE)
  pt1 <- ggplot2::ggplot(df, ggplot2::aes(indfine, muRMSE, linetype = "Model"))+
    ggplot2::geom_line(color='blue', size=2)+
    ggplot2::geom_line(data=df2, ggplot2::aes(indfine, sumscrRMSE, linetype = "Sum score"), 
                       color='blue', size=2)+
    ggplot2::scale_linetype_manual(name="",
                                   labels = c("Model", 
                                              "Sum score"), 
                                   values=c("solid","dashed"))+
    ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed")+
    ggplot2::xlim(indrng)+
    ggplot2::ylim(RMSErng)+
    ggplot2::ylab("Test Score RMSE")+
    ggplot2::xlab("")+ 
    ggplot2::theme(legend.position="top")+
    ggplot2::guides(linetype=ggplot2::guide_legend(keywidth = 5, keyheight = 1))
  pt1 <- pt1 + My_Theme 
  # plot(indfine, muRMSE, type="l", lty=1, col=1, lwd=2, 
  #      xlim=indrng, ylim=RMSErng, 
  #      xlab="True Percentile Index", ylab="Test Score RMSE")
  # lines(indfine, sumscrRMSE, col=1, lty=2,  lwd=2)
  # for (k in 1:length(Qvec)) {
  #   lines(c(Qvec[k],Qvec[k]), RMSErng, col=1, lty=2, lwd=1)
  # }
  #  ------------  bias results  ------------
  df <- data.frame(x=indfine,y=mubias)
  df2 <- data.frame(x=indfine,y=sumscrbias)
  pt2 <- ggplot2::ggplot(df, ggplot2::aes(indfine, mubias))+
    ggplot2::geom_line(color='blue', size=2)+
    ggplot2::geom_line(data=df2, ggplot2::aes(indfine, sumscrbias), 
                       color='blue', size=2, linetype = "dashed")+
    ggplot2::scale_linetype_manual(name="",
                                   labels = c("Model", 
                                              "Sum score"), 
                                   values=c("solid","dashed"))+
    ggplot2::geom_hline(yintercept = 0,  color="black", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = Qvec, color="black", linetype = "dashed")+
    ggplot2::xlim(0,100)+
    ggplot2::ylim(biasrng)+
    ggplot2::ylab("Test Score Bias")+
    ggplot2::xlab("True Percentile Index")
  pt2 <- pt2 + My_Theme 
  
  # plot(indfine, mubias,  type="l", lty=1, col=1, lwd=2, 
  #      xlim=c(0,100), ylim=biasrng, 
  #      xlab="True Percentile Index", ylab="Test Score Bias")  
  # lines(indfine, sumscrbias,  col=1, lty=2,  lwd=2)
  # lines(indrng, c(0,0), col=1, lty=2, lwd=1);
  # for (k in 1:length(Qvec)) {
  #   lines(c(Qvec[k],Qvec[k]), biasrng, col=1, lty=2, lwd=1)
  # }
  # text(Qvec[1],biasrng[1]+0.1,' 5%')
  # text(Qvec[2],biasrng[1]+0.1,'25%')
  # text(Qvec[3],biasrng[1]+0.1,'50%')
  # text(Qvec[4],biasrng[1]+0.1,'75%')
  # text(Qvec[5],biasrng[1]+0.1,'95%')
  
  p <- ggpubr::ggarrange(pt1, pt2, ncol = 1, nrow = 2)
  return(p)
}
