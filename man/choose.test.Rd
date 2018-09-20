\name{choose.test}
\alias{choose.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make Recommendations on the Most Powerful Test to Use
}
\description{
Performs simulations to compare the power of different tests
}
\usage{
choose.test(Xpaired, Ypaired, Xextra = NULL, Yextra = NULL, mc.rep = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xpaired}{
%%     ~~Describe \code{Xpaired} here~~
}
  \item{Ypaired}{
%%     ~~Describe \code{Ypaired} here~~
}
  \item{Xextra}{
%%     ~~Describe \code{Xextra} here~~
}
  \item{Yextra}{
%%     ~~Describe \code{Yextra} here~~
}
  \item{mc.rep}{
%%     ~~Describe \code{mc.rep} here~~
}
}
\examples{


# There are unpaired observations from both samples
dat=sim.partially.matched(m=20,n.x=40,n.y=5,distr="normal",
    params=c(loc.2=.8,rho=.1,scale.2=1),seed=1)
choose.test(dat$X, dat$Y, dat$Xprime, dat$Yprime)

## There are unpaired observations from only one sample
#dat=sim.partially.matched(m=20,n.x=0,n.y=10,distr="normal",
#    params=c(loc.2=.5,rho=.8,scale.2=1),seed=1)
#choose.test(dat$X, dat$Y, dat$Xprime, dat$Yprime)


}
