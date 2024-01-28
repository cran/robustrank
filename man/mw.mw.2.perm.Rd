\name{mw.mw.2.perm}
\alias{mw.mw.2.perm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A Test that Combines WMW for Paired Data and WMW for Unpaired Data
}
\description{
Use permutation-based reference distribution to obtain p values for a test that combines WMW for paired data and WMW for unpaired data
}
\usage{
mw.mw.2.perm(X, Y, Xprime, Yprime, .corr, mc.rep = 10000, 
    alternative = c("two.sided", "less", "greater"), verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
X
}
  \item{Y}{
Y
}
  \item{Xprime}{
Xprime
}
  \item{Yprime}{
Yprime
}
  \item{.corr}{
.corr
}
  \item{mc.rep}{
mc.rep
}
  \item{alternative}{
alternative
}
  \item{verbose}{
verbose
}
}
