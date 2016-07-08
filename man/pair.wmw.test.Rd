\name{pair.wmw.test}
\alias{pair.wmw.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
WMW test for paired data
}
\description{
Performs a WMW-type test of the strong null for paired data.
}
\usage{
pair.wmw.test(X,Y
, alternative = c("two.sided", "less", "greater"), correct = TRUE, perm=NULL, mc.rep=1e4
, method=c("exact.2","large.0","large","exact","exact.0","exact.1","exact.3"), trace=0
, mode=c("test","var"), useC=TRUE)
}
\arguments{
  \item{X}{
Sample 1.
}
  \item{Y}{
Sample 2. 
}
  \item{alternative}{
Alternative hypothesis.
}
  \item{correct}{
Whether to apply continuity correction.
}
  \item{perm}{
Whether to use permutation distribution or normal approximation to find p-value. See details.
}
  \item{mc.rep}{
Number of Monte Carlo replicates for permutation test.
}
  \item{method}{
Choices of test statistics.
}
  \item{trace}{
Print debug message when positive.
}
  \item{mode}{
For development used.
}
  \item{useC}{
For development used.
}
}
\details{
When perm is NULL, if (min(m,n)>=20) normal approximatino is used to find p value, otherwise permutation test is used.
When permutation test is used, if the number of possible permutations is less than mc.rep, a test statistic is computed for all permutations; otherwise, Monte Carlo is done.

}
\value{
P value for now.
}
\references{
Under prep.
}
\examples{

dat=sim.partially.matched(m=15,n.x=0,n.y=20,distr="mixnormal",params=c(p.1=0.3,p.2=0.3),seed=1)
X=dat$X; Y=dat$Y
pair.wmw.test(X, Y, perm=TRUE,  method="large.0", trace=1)
pair.wmw.test(X, Y, perm=FALSE, method="large.0", trace=1)

}
