\name{HillPlot}
\alias{HillPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HillPlot(dat, from, to, interval, B, conf, datatype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{from}{
%%     ~~Describe \code{from} here~~
}
  \item{to}{
%%     ~~Describe \code{to} here~~
}
  \item{interval}{
%%     ~~Describe \code{interval} here~~
}
  \item{B}{
%%     ~~Describe \code{B} here~~
}
  \item{conf}{
%%     ~~Describe \code{conf} here~~
}
  \item{datatype}{
%%     ~~Describe \code{datatype} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dat, from, to, interval, B, conf, datatype) 
{
    q <- seq(from, to, by = interval)
    if (datatype == "abundance") {
        y <- dat
        n <- sum(dat)
    }
    else {
        y <- dat[-1]
        n <- dat[1]
    }
    conf.reg = function(x, LCL, UCL, ...) {
        x.sort <- order(x)
        x <- x[x.sort]
        LCL <- LCL[x.sort]
        UCL <- UCL[x.sort]
        polygon(c(x, rev(x)), c(LCL, rev(UCL)), ...)
    }
    z <- -qnorm((1 - conf)/2)
    D <- length(which(y != 0))
    U <- sum(y)
    out.pro <- qDPRO(dat, q, B, datatype)
    ci.pro.l <- out.pro[1, ] - z * out.pro[2, ]
    ci.pro.u <- out.pro[1, ] + z * out.pro[2, ]
    out.mle <- qDMLE(dat, q, B, datatype)
    ci.mle.l <- out.mle[1, ] - z * out.mle[2, ]
    ci.mle.u <- out.mle[1, ] + z * out.mle[2, ]
    plot(q, out.pro[1, ], type = "l", lty = 5, lwd = 2, col = 2, 
        xlab = "index q", ylab = expression(paste(""^q, "D")), 
        main = "Hill's Number", ylim = c(min(ci.pro.l, ci.mle.l), 
            max(ci.pro.u, ci.mle.u)))
    points(q, out.mle[1, ], type = "l", lty = 1, lwd = 2)
    conf.reg(q, ci.mle.l, ci.mle.u, col = adjustcolor("red", 
        0.25), border = NA)
    conf.reg(q, ci.pro.l, ci.pro.u, col = adjustcolor("green", 
        0.25), border = NA)
    legend("topright", legend = c("MLE", "Proposed"), lty = c(1, 
        5), lwd = c(1, 2), col = c(1, "red"), cex = 1, bty = "n")
    table.est <- matrix(c(out.mle[1, ][which(round(q) == q)], 
        out.pro[1, ][which(round(q) == q)]), nrow = 2, byrow = T)
    rownames(table.est) = c("MLE", "Pro")
    colnames(table.est) = c("q=0", "q=1", "q=2")
    table.sd <- matrix(c(out.mle[2, ][which(round(q) == q)], 
        out.pro[2, ][which(round(q) == q)]), nrow = 2, byrow = T)
    rownames(table.sd) = c("MLE", "Pro")
    colnames(table.sd) = c("q=0", "q=1", "q=2")
    return(list(EST = data.frame(table.est), SD = data.frame(table.sd)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
