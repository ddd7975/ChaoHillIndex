ChaoHill <-
function(dat, datatype=c("abundance", "incidence"), from=0, to=2, interval=0.1,
         B=200, conf=0.95, detail=c(TRUE, FALSE)){ # for real data estimation
  if (is.matrix(dat) == T || is.data.frame(dat) == T){
    if (ncol(dat) != 1 & nrow(dat) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(dat) == 1){
      data <- data[, 1]
    } else {
      data <- data[1, ]
    }
  }
  q <- seq(from, to, by=interval)
  if (datatype == "abundance"){
    y <- dat; n <- sum(dat)
  }else{
    y <- dat[-1]; n <- dat[1]
  }
  ###------------------------------------
  conf.reg=function(x,LCL,UCL,...) {
    x.sort <- order(x)
    x <- x[x.sort]
    LCL <- LCL[x.sort]
    UCL <- UCL[x.sort]
    polygon(c(x,rev(x)),c(LCL,rev(UCL)), ...)
  }
  ###------------------------------------ 
  z <- -qnorm((1 - conf)/2)
  D <- length(which(y != 0)); U <- sum(y)
  out.pro <- qDPRO(dat, q, B, datatype)
  ci.pro.l <- round(out.pro[1,] - z*out.pro[2,], 3)
  ci.pro.u <- round(out.pro[1,] + z*out.pro[2,], 3)
  out.mle <- qDMLE(dat, q, B, datatype)
  ci.mle.l <- round(out.mle[1,] - z*out.mle[2,], 3)
  ci.mle.u <- round(out.mle[1,] + z*out.mle[2,], 3)
  plot(q, out.pro[1,], type="l", lty=5, lwd=2, col=2, 
       xlab="index q",
       ylab=expression(paste(""^q, "D")), 
       main = "Hill's Number",
       ylim=c(min(ci.pro.l, ci.mle.l), max(ci.pro.u, ci.mle.u)))
  points(q, out.mle[1,], type="l", lty=1, lwd=2)
  
  #  points(q, qdmle.s, type = "l", lty = 3, lwd = 2, col = 4)
  conf.reg(q, ci.mle.l, ci.mle.u, col=adjustcolor("red", 0.25), border=NA)
  #  points(q, qdpro.s, type = "l", lty = 5, lwd = 2, col = 2)
  conf.reg(q, ci.pro.l, ci.pro.u, col=adjustcolor("green", 0.25), border=NA)
  legend("topright", legend=c("MLE", "Proposed"),
         lty=c(1, 5), lwd=c(1, 2), col=c(1, "red"), cex=1, bty="n")
  if (detail == TRUE){
    table.est <- data.frame(matrix(c(out.mle[1,], out.pro[1,]),
                                   nrow=2, byrow=T))
  }else{
    table.est <- data.frame(matrix(c(out.mle[1,][which(round(q)==q)], out.pro[1,][which(round(q)==q)]),
                                   nrow=2, byrow=T))
  }
  rownames(table.est)=c("Obersved", "Chao_2013")
  if (detail==T){
    colnames(table.est)=paste("q=", q, sep="")    
  }else{
    colnames(table.est)=c("q = 0", "q = 1", "q = 2")
  }
  
  if(detail==T){
    table.sd <- data.frame(matrix(c(out.mle[2,], out.pro[2,]),
                                  nrow=2, byrow=T))
  }else{
    table.sd <- data.frame(matrix(c(out.mle[2,][which(round(q)==q)], out.pro[2,][which(round(q)==q)]),
                                nrow=2, byrow=T))
  }
  rownames(table.sd)=c("Obsersved", "Chao_2013")
  if(detail==T){
    colnames(table.sd)=paste("q=", q, sep="")  
    }else{
    colnames(table.sd)=c("q = 0", "q = 1", "q = 2")
  }
  return(list(EST = table.est, SD = table.sd))
}
