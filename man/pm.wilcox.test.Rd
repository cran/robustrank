\name{pm.wilcox.test}
\alias{pm.wilcox.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wilcoxon test for Partially Matched Two Sample Data
}
\description{
Performs rank-based two sample test for partially matched two sample data by combining information from matched and unmatched data
}
\usage{
pm.wilcox.test(Xpaired, Ypaired, Xextra = NULL, Yextra = NULL,
 alternative = c("two.sided", "less", "greater"),
 method = c("SR-MW", "MW-MW", "all"), mode = c("test",
 "var", "power.study"), useC = FALSE, correct = NULL,
 verbose=FALSE)
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
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
}
  \item{method}{
String. SR-MW is recommended, all is for development only.
}
  \item{mode}{
String. Do not change it to var, for developpment only.
}
  \item{useC}{
Boolean. Do not set it to TRUE, for development only.
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{correct}{
Continuity correction.
}
}
\details{
If Xpaired and Ypaired have NAs, the corresponding unpaired data in Ypaired and Xpaired will be combined with Yextra and Xextra.

}
\value{
An htest object.
}
\examples{

set.seed(1)
z=rnorm(20, sd=0.5) # induces correlation between X and Y
X=rnorm(20)+z
Y=rnorm(20,mean=0.8)+z
X[1:10]=NA
boxplot(X,Y,names=c("X","Y"))

pm.wilcox.test(X,Y)
# for comparison
wilcox.test(X,Y,paired=TRUE)
wilcox.test(X,Y,paired=FALSE)# often a conservative test due to the correlation

# no paired data
Y1=Y
Y1[11:20]=NA
pm.wilcox.test(X,Y1)
# should match the following
wilcox.test(X,Y1,paired=FALSE)

# only 1 pair of matched data
Y1=Y
Y1[12:20]=NA
pm.wilcox.test(X,Y1)


}
