\name{sim.partially.matched}
\alias{r2sample}
\alias{sim.paired.with.replicates}
\alias{sim.partially.matched}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Paired, Independent, or Partially Matched Two-Sample Data
}
\description{
sim.partially.matched generates partially matched two-sample data.  for Monte Carlo studies. 
r2sample is a wrapper for sim.partially.matched and generates indepenent two-sample data.
}
\usage{
sim.partially.matched(m, n.x, n.y, 
distr = c("normal","logistic","student","mixnormal","gamma","lognormal","beta",
    "uniform","hybrid1","hybrid2","doublexp"), params, seed)

r2sample(m, n, 
distr = c("normal", "logistic", "student", "mixnormal"), params, seed)

sim.paired.with.replicates(m, meanRatio, sdRatio, within.sd, type, hyp, distr, seed)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
Number of pairs.
}
  \item{n}{
Number of Ys.
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
Named vector. See details.
}
  \item{seed}{
Seed for random number generator.
}
  \item{meanRatio}{
%%     ~~Describe \code{meanRatio} here~~
}
  \item{sdRatio}{
%%     ~~Describe \code{sdRatio} here~~
}
  \item{within.sd}{
%%     ~~Describe \code{within.sd} here~~
}
  \item{type}{
%%     ~~Describe \code{type} here~~
}
  \item{hyp}{
%%     ~~Describe \code{hyp} here~~
}
}
\details{
If the distribution is in \code{c("normal","student","logistic")}, \code{params} should have three fields: loc.2, rho and scale.2; loc.1 is set to 0 and scale.1 is set to 1.

If the distribution is mixnormal, \code{params} should have three fields: p.1, p.2 and sd.n. 

If the distribution is gamma, \code{params} should have fix fields: loc.2, shape.1, shape.2, rate.1, rate.2 and rho.

For details on bivariate logistic distribution, see \code{rbilogistic}
    
}
\value{
sim.partially.matched return a list with the following components:
  \item{X}{m sample 1 that pair with Y}
  \item{Y}{m sample 2 that pair with X}
  \item{Xprime}{n.x sample 1}
  \item{Yprime}{n.y sample 2}
r2sample returns a list with the following components:
  \item{X}{m sample 1 that are independent of Y}
  \item{Y}{n sample 2 that are independent of X}
}
\examples{
dat=sim.partially.matched(m=10,n.x=5,n.y=4,distr="normal",
    params=c("loc.2"=0,"rho"=0,"scale.2"=1),seed=1)
X=dat$X; Y=dat$Y; Yprime=dat$Yprime    

dat=sim.partially.matched(m=10,n.x=5,n.y=4,distr="logistic",
    params=c("loc.2"=0,"rho"=0,"scale.2"=1),seed=1)
X=dat$X; Y=dat$Y; Yprime=dat$Yprime    


}
