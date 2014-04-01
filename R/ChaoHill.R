ChaoHill <-
function(dat, datatype=c("abundance", "incidence"), from=0, to=2, interval=0.1,
         B=200, conf=0.95, detail=c(TRUE, FALSE)){ # for real data estimation
  if (is.matrix(dat) == T || is.data.frame(dat) == T){
    if (ncol(dat) != 1 & nrow(dat) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(dat) == 1){
      dat <- dat[, 1]
    } else {
      dat <- dat[1, ]
    }
  }
  dat <- as.numeric(dat)
  q <- seq(from, to, by=interval)
  if (datatype == "abundance"){
    y <- dat; n <- sum(dat)
  }else{
    y <- dat[-1]; n <- dat[1]
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
#   plot(q, out.pro[1,], type="l", lty=5, lwd=2, col=2, 
#        xlab="index q",
#        ylab=expression(paste(""^q, "D")), 
#        main = "Hill's Number",
#        ylim=c(min(ci.pro.l, ci.mle.l), max(ci.pro.u, ci.mle.u)))
#   points(q, out.mle[1,], type="l", lty=1, lwd=2)
#   
#   #  points(q, qdmle.s, type = "l", lty = 3, lwd = 2, col = 4)
#   conf.reg(q, ci.mle.l, ci.mle.u, col=adjustcolor("red", 0.25), border=NA)
#   #  points(q, qdpro.s, type = "l", lty = 5, lwd = 2, col = 2)
#   conf.reg(q, ci.pro.l, ci.pro.u, col=adjustcolor("green", 0.25), border=NA)
#   legend("topright", legend=c("MLE", "Proposed"),
#          lty=c(1, 5), lwd=c(1, 2), col=c(1, "red"), cex=1, bty="n")
  if (detail == TRUE){
    table.est <- data.frame(matrix(c(out.mle[1,], out.pro[1,]),
                                   nrow=2, byrow=T))
    table.sd <- data.frame(matrix(c(out.mle[2,], out.pro[2,]),
                                  nrow=2, byrow=T))
    table.lci <- table.est - z*table.sd
    table.uci <- table.est + z*table.sd
    colnames(table.est) <- colnames(table.sd) <- colnames(table.lci) <- colnames(table.uci) <- paste("q=", q, sep="")    
  }else{
    table.est <- data.frame(matrix(c(out.mle[1,][which((q %in% c(0, 1, 2))==T)],
                                     out.pro[1,][which((q %in% c(0, 1, 2))==T)]),
                                   nrow=2, byrow=T))
    table.sd <- data.frame(matrix(c(out.mle[2,][which((q %in% c(0, 1, 2))==T)], 
                                    out.pro[2,][which((q %in% c(0, 1, 2))==T)]),
                                  nrow=2, byrow=T))
    table.lci <- table.est - z*table.sd
    table.uci <- table.est + z*table.sd
    colnames(table.est) <- colnames(table.sd) <- colnames(table.lci) <- colnames(table.uci) <- c("q=0", "q=1", "q=2")
  }
  rownames(table.est) <- rownames(table.sd) <- rownames(table.lci) <- rownames(table.uci) <- c("Observed", "Chao_2013")
  return(list(EST = round(table.est, 3), 
              SD = round(table.sd, 3),
              LCI = round(table.lci, 3),
              UCI = round(table.uci, 3)))
}
