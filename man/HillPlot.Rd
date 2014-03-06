\name{HillPlot}
\alias{HillPlot}
\title{
Estimation of Hill's number
}
\description{
\kbd{Hillplot} is a function to provide a Species richness estimator and curve.
}
\usage{
HillPlot(dat, from, to, interval, B, conf, datatype)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
 a vector of species abundance or incidence frequency. If \code{datatype = "incidence"}, then the input format of first entry should be total number of sampling units, and followed by species incidence frequency.}
}
  \item{from}{
the start of index q
}
  \item{to}{
the end of index q
}
  \item{interval}{
decide the interval of sequence q
}
  \item{B}{
number of bootstrap
}
  \item{conf}{
positive number \eqn{\le} 1. "\code{conf}" specifies the confidence level for confidence   interval. The default is \code{0.95}.
}
  \item{datatype}{
the data type of input data. That is individual-based abundance data (\code{datatype = "abundance"}) or sample-based incidence data (\code{datatype = "incidence"}).
}
\value{
\kbd{HillPlot} returns a table of Hill's number estimators, their standard error and 95 \% confidence interval which the \code{method} you choose.
}
\references{

}
\author{
K.S. Tseng ; Y.H. Lee : Anne Chao\cr
Maintainer : Anne Chao \email{chao@stat.nthu.edu.tw}
}
\examples{
}