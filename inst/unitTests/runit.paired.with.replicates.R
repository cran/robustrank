library("RUnit")
library("robustrank")


test.paired.with.replicates <- function() {

  suppressWarnings(RNGversion("3.5.0"))
RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-3
if(file.exists("C:/_checkReproducibility")) tolerance=1e-6
verbose=FALSE

dat=sim.paired.with.replicates(m=15, meanRatio=10, sdRatio=0, within.sd=.3, type=1, hyp=0, distr="logistic", seed=1)
checkEqualsNumeric(abs(sum(unlist(dat))), 38.86667, tolerance=tolerance)
checkEqualsNumeric(multinom.test(dat$X, dat$Y, useC=TRUE,  alternative = "two.sided", correct = FALSE, perm=TRUE, mc.rep=101), 0.8910891, tolerance=tolerance)
checkEqualsNumeric(multinom.test(dat$X, dat$Y, useC=FALSE, alternative = "two.sided", correct = FALSE, perm=TRUE, mc.rep=101), 0.8910891, tolerance=tolerance)

dat=sim.paired.with.replicates(m=15, meanRatio=1, sdRatio=0, within.sd=.3, type=1, hyp=0, distr="logistic", seed=1)
checkEqualsNumeric(multinom.test(dat$X, dat$Y, useC=TRUE,  alternative = "two.sided", correct = FALSE, perm=TRUE, mc.rep=101), 0.5841584, tolerance=tolerance)
checkEqualsNumeric(multinom.test(dat$X, dat$Y, useC=FALSE, alternative = "two.sided", correct = FALSE, perm=TRUE, mc.rep=101), 0.5841584, tolerance=tolerance)

dat=sim.paired.with.replicates(m=15, meanRatio=11, sdRatio=0, within.sd=.3, type=1, hyp=0, distr="logistic", seed=1)
checkEqualsNumeric( wmw.paired.replicates.test(dat$X, dat$Y, useC=FALSE, alternative = "two.sided", correct = FALSE, perm=TRUE, mc.rep=101), 
                    wmw.paired.replicates.test(dat$X, dat$Y, useC=TRUE,  alternative = "two.sided", correct = FALSE, perm=TRUE, mc.rep=101), tolerance=tolerance)



}
