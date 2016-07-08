\name{sim.partially.matched}
\alias{sim.partially.matched}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Partially Matched Two-Sample Data
}
\description{
Generates partially matched two-sample data for Monte Carlo studies.
}
\usage{
sim.partially.matched(m, n.x, n.y, 
distr = c("normal", "logistic", "student", "mixnormal", "gamma"), params, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
Number of pairs.
}
  \item{n.x}{
Number of extra Xs.
}
  \item{n.y}{
Number of extra Ys.
}
  \item{distr}{
Distributions.
}
  \item{params}{
Parameters. See details.
}
  \item{seed}{
Seed for random number generator.
}
}
\details{
The parameters are interpreted in the following block of code:
    if (distr %in% c("normal","student","logistic")) 
        params=c(loc.2=paramss[1],rho=paramss[2],scale.2=paramss[3])
    else if (distr=="mixnormal") 
        params=c(p.1=paramss[1],p.2=paramss[2],sd.n=paramss[3])
    else if (distr=="gamma") 
        params=c(loc.2=paramss[1],shape.1=paramss[2],shape.2=paramss[3],rate.1=paramss[4],rate.2=paramss[5],rho=paramss[6])
    
}
\value{
A list with the following components:
  \item{X}{m sample 1 that pair with Y}
  \item{Y}{m sample 2 that pair with X}
  \item{Xprime}{n.x sample 1}
  \item{Yprime}{n.y sample 2}
}
\seealso{
 \code{\link{r2sample}}
}
\examples{
\dontrun{
    dat=sim.partially.matched(m=m,n.x=0,n.y=0,distr=distr,params=params,seed=seed)
    X=dat$X; Y=dat$Y; Yprime=dat$Yprime    
}

}
