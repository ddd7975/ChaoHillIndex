qDMLE <-
function(dat, q, B, datatype=c("abundance", "incidence")){ # dat=(n, y1, y2,..., ys), q can be a vector.
  dat <- as.numeric(dat)
  qDmle <- function(y, q, n, datatype){
    if(datatype == "abundance"){
     p <- y/n; p <- p[which(p != 0)] 
     if (q == 1){
       D_mle <- exp(-sum(p/sum(p)*log(p/sum(p))))  
     } else {
       D_mle <- sum((p/sum(p))^q)^(1/(1 - q))
     }
    }else{
      p <- y/n; p <- p[which(p != 0)]
      if (q == 1){
        D_mle <- exp(-sum(p/sum(p)*log(p/sum(p))))  
      } else {
        D_mle <- sum((p/sum(p))^q)^(1/(1 - q))
      }
    }
    return(D_mle)
  }
  ibootstrap <- function(y, fun, B, q, n, datatype){      # improved bootstrap
    if (datatype == "abundance"){
      f <- function(r, data){sum(data == r)}
      y <- y[which(y != 0)]
      f1 <- f(1, y); f2 <- f(2, y)
      
      if (f2 > 0){
        alfa1 <- 2*f2/((n - 1)*f1 + 2*f2)
      } else if (f2 == 0 & f1 != 0){
        alfa1 <- 2/((n - 1)*(f1 - 1) + 1)
      } else {
        alfa1 <- 1  
      }
      Chat <- 1 - f1/n*(1 - alfa1)
      
      W <- (1 - Chat)/sum(y/n*(1 - y/n)^n)
      pi.hat <- y/n*(1 - W*(1 - y/n)^n)
      
      if (f2 > 0){
        f0_hat <- (n - 1)/n*f1^2/(2*f2)
      }else{
        f0_hat <- (n - 1)/n*f1*(f1 - 1)/(2*(f2 + 1))
      }
      f0_hat <- ceiling(f0_hat)
      pi.hat.r <- rep((1 - Chat)/f0_hat, f0_hat)         
      pi.star <- c(pi.hat, pi.hat.r)     
      
      set.seed(456)
      Y <- matrix(rmultinom(B, n, pi.star), ncol=B)
      
      b <- 1:B
      sd(sapply(b, function(b){fun(Y[, b], q, n, datatype)}))
      
    }else{
      Q <- function(r, data){sum(data == r)} 
      y <- y[which(y != 0)]
      q1 <- Q(1, y); q2 <- Q(2, y)
      U <- sum(y)
      
      if (q2 > 0){
        alfa1 <- 2*q2/((n - 1)*q1 + 2*q2)
      } else if (q2 == 0 & q1 != 0){
        alfa1 <- 2/((n - 1)*(q1 - 1) + 1)
      } else {
        alfa1 <- 1  
      }
      Chat <- 1 - q1/U*(1 - alfa1)
      
      W <- U/n*(1 - Chat)/sum(y/n*(1 - y/n)^n)
      pi.hat <- y/n*(1 - W*(1 - y/n)^n)
      
      if (q2 > 0){
        Q0_hat <- (n - 1)/n*q1^2/(2*q2)
      }else{
        Q0_hat <- (n - 1)/n*q1*(q1 - 1)/(2*(q2 + 1))
      }
      Q0_hat <- ceiling(Q0_hat)
      pi.hat.r <- rep(U/n*(1 - Chat)/Q0_hat, Q0_hat)         
      pi.star <- c(pi.hat, pi.hat.r)     
      
      set.seed(456)
      Y <- matrix(rbinom(length(pi.star) * B, n, pi.star), ncol = B)
      
      b <- 1:B
      sd(sapply(b, function(b){fun(Y[, b], q, n, datatype)}))
      
    }
  }
  if (datatype == "abundance"){
    y <- dat; n <- sum(dat)
  }else{
    y <- dat[-1]; n <- dat[1]
  }
  #estimate
  i <- 1:length(q)
  est.value <- sapply(i, function(i) qDmle(y, q[i], n, datatype))
  est.sd <- sapply(i, function(i) ibootstrap(y, qDmle, B, q[i], n, datatype))
  return(matrix(c(est.value, est.sd), nrow=2, byrow=T))
}
