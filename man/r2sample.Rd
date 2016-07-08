\name{r2sample}
\alias{r2sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate Two-Sample Data
}
\description{
Generates independent two-sample data for Monte Carlo studies.
}
\usage{
r2sample(m, n, distr = c("normal", "logistic", "student", "mixnormal"), params, seed)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
Number of Xs.
}
  \item{n}{
Number of Ys.
}
  \item{distr}{
Distributions.
}
  \item{params}{
See params and details in sim.partially.matched man page.
}
  \item{seed}{
See for random number generator.
}
}
\details{
r2sample is a wrapper for sim.partially.matched.
}
\value{
A list with the following components:
  \item{X}{m sample 1 that are independent of Y}
  \item{Y}{n sample 2 that are independent of X}
}
\seealso{
 \code{\link{sim.partially.matched}}
}
