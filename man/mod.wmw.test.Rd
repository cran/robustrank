\name{mod.wmw.test}
\alias{mod.wmw.test}
\title{
Modified Wilcoxon-Mann-Whitney Test
}
\description{
Also known as the Fligner-Policello test.
}
\usage{
mod.wmw.test(X, Y, alternative = c("two.sided", "less", "greater"),
         correct = TRUE, perm = NULL, mc.rep = 10000, method =
         c("combine", "comb2", "fp", "wmw", "fplarge", "nsm3"),
         verbose = FALSE, mode = c("test", "var"), useC = TRUE)
}
\arguments{
  \item{X}{
Samples from population 1.
}
  \item{Y}{
Samples from population 2.
}
  \item{alternative}{
Directon of the alternative hypothesis.
}
  \item{correct}{
Whether to do continutiy correction.
}
  \item{perm}{
Boolean, whether to do permutation to get p-value or use normal approximation. See details.
}
  \item{mc.rep}{
Default number of replicates when doing permutation. See details.
}
  \item{method}{
For development.
}
  \item{verbose}{
For development. Print some debug info.
}
  \item{mode}{
For development. 
}
  \item{useC}{
For development. Run C or R implementation.
}
}
\details{
When perm is null, we will compute permutation-based p values if either sample size is less than 20 and compute normal approximation-based p values otherwise.
When doing permuation, if the possible number of combinations is less than mc.rep, every possible configuration is done.
}
\value{
A p value for now.
}
\references{
manuscript in preperation
}
\examples{

# Example 4.1, Hollander, Wolfe and Chicken (2014) Nonparameteric Statistics
X <- c(0.80, 0.83, 1.89, 1.04, 1.45, 1.38, 1.91, 1.64, 0.73, 1.46)
Y <- c(1.15, 0.88, 0.90, 0.74, 1.21)
mod.wmw.test(X, Y, method="wmw", alternative="greater")
mod.wmw.test(X, Y, method="combine", alternative="greater", verbose=1)

# Section 4.1 Problem 1, Hollander et al. 
X=c(1651,1112,102.4,100,67.6,65.9,64.7,39.6,31.0)
Y=c(48.1,48.0,45.5,41.7,35.4,34.3,32.4,29.1,27.3,18.9,6.6,5.2,4.7)
mod.wmw.test(X, Y, method="wmw")
mod.wmw.test(X, Y, method="combine", verbose=1)

# Section 4.1 Problem 5, Hollander et al. 
X=c(12 ,44 ,34 ,14 ,9  ,19 ,156,23 ,13 ,11 ,47 ,26 ,14 ,33 ,15 ,62 ,5  ,8  ,0  ,154,146)
Y=c(37,39,30,7,13,139, 45,25,16,146,94,16,23,1,290,169,62,145,36, 20, 13)
mod.wmw.test(X, Y, method="wmw", alternative="less")
mod.wmw.test(X, Y, method="combine", alternative="less", verbose=1)

# Section 4.1 Problem 15, Hollander et al. 
X=c(0.19,0.14,0.02,0.44,0.37)
Y=c(0.89,0.76,0.63,0.69,0.58,0.79,0.02,0.79)
mod.wmw.test(X, Y, method="wmw")
mod.wmw.test(X, Y, method="combine", verbose=1)

# Table 4.7, Hollander et al. 
X=c(297,340,325,227,277,337,250,290)
Y=c(293,291,289,430,510,353,318)
mod.wmw.test(X, Y, method="wmw", alternative="less")
mod.wmw.test(X, Y, method="combine", alternative="less", verbose=1)


}
\keyword{ Mann-Whitney U test }
\keyword{ Wilcoxon rank sum test }
