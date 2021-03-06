library("RUnit")
library("robustrank")


test.pair.wmw <- function() {

  suppressWarnings(RNGversion("3.5.0"))
RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-3
if(file.exists("C:/_checkReproducibility")) tolerance=1e-6

verbose=FALSE
dat=sim.partially.matched(m=15,n.x=0,n.y=20,distr="mixnormal",params=c(p.1=0.3,p.2=0.3),seed=1)
X=dat$X; Y=dat$Y
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, mode="var", useC=FALSE, verbose=verbose), c(0.59555556,0.08880423,0.08344974,0.06187795,0.23651217,0.07689877,0.08519244,0.08614423), tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, verbose=verbose, perm=FALSE, method="large.0")$p.value, pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, verbose=verbose, perm=FALSE, method="large.0")$p.value, tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, verbose=verbose, perm=FALSE, method="large")$p.value,   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, verbose=verbose, perm=FALSE, method="large")$p.value, tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, verbose=verbose, perm=FALSE, method="exact")$p.value,   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, verbose=verbose, perm=FALSE, method="exact")$p.value, tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, verbose=verbose, perm=FALSE, method="exact.0")$p.value,   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, verbose=verbose, perm=FALSE, method="exact.0")$p.value, tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, verbose=verbose, perm=FALSE, method="exact.1")$p.value,   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, verbose=verbose, perm=FALSE, method="exact.1")$p.value, tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, verbose=verbose, perm=FALSE, method="exact.2")$p.value,   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, verbose=verbose, perm=FALSE, method="exact.2")$p.value, tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, verbose=verbose, perm=FALSE, method="exact.3")$p.value,   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, verbose=verbose, perm=FALSE, method="exact.3")$p.value, tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, verbose=verbose, perm=TRUE,  method="large", mc.rep=10)$p.value, pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, verbose=verbose, perm=TRUE, method="large", mc.rep=10)$p.value, tolerance=tolerance)

dat=sim.partially.matched(m=5,n.x=0,n.y=20,distr="mixnormal",params=c(p.1=0.3,p.2=0.3),seed=1)
X=dat$X; Y=dat$Y
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, verbose=verbose, perm=TRUE, method="large", mc.rep=4e4)$p.value, pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, verbose=verbose, perm=TRUE, method="large", mc.rep=4e4)$p.value, tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=TRUE, useC=TRUE, verbose=verbose, perm=TRUE, method="large", mc.rep=4e4)$p.value, pair.wmw.test(X, Y, correct=TRUE, useC=FALSE, verbose=verbose, perm=TRUE, method="large", mc.rep=4e4)$p.value, tolerance=tolerance)


}
