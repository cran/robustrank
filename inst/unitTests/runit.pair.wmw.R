library("RUnit")
library("robustrank")


test.pair.wmw <- function() {

RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-3
if(file.exists("C:/_checkReproducibility")) tolerance=1e-6

trace=0
dat=sim.partially.matched(m=15,n.x=0,n.y=20,distr="mixnormal",params=c(p.1=0.3,p.2=0.3),seed=1)
X=dat$X; Y=dat$Y
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, mode="var", useC=FALSE, trace=trace), c(0.59555556,0.08880423,0.08344974,0.06187795,0.23651217,0.07689877,0.08519244,0.08614423), tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, trace=trace, perm=FALSE, method="large.0"), pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, trace=trace, perm=FALSE, method="large.0"), tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, trace=trace, perm=FALSE, method="large"),   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, trace=trace, perm=FALSE, method="large"), tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, trace=trace, perm=FALSE, method="exact"),   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, trace=trace, perm=FALSE, method="exact"), tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, trace=trace, perm=FALSE, method="exact.0"),   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, trace=trace, perm=FALSE, method="exact.0"), tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, trace=trace, perm=FALSE, method="exact.1"),   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, trace=trace, perm=FALSE, method="exact.1"), tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, trace=trace, perm=FALSE, method="exact.2"),   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, trace=trace, perm=FALSE, method="exact.2"), tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, trace=trace, perm=FALSE, method="exact.3"),   pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, trace=trace, perm=FALSE, method="exact.3"), tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, trace=trace, perm=TRUE,  method="large", mc.rep=10), pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, trace=trace, perm=TRUE, method="large", mc.rep=10), tolerance=tolerance)

dat=sim.partially.matched(m=5,n.x=0,n.y=20,distr="mixnormal",params=c(p.1=0.3,p.2=0.3),seed=1)
X=dat$X; Y=dat$Y
checkEqualsNumeric(pair.wmw.test(X, Y, correct=FALSE, useC=TRUE, trace=trace, perm=TRUE, method="large", mc.rep=4e4), pair.wmw.test(X, Y, correct=FALSE, useC=FALSE, trace=trace, perm=TRUE, method="large", mc.rep=4e4), tolerance=tolerance)
checkEqualsNumeric(pair.wmw.test(X, Y, correct=TRUE, useC=TRUE, trace=trace, perm=TRUE, method="large", mc.rep=4e4), pair.wmw.test(X, Y, correct=TRUE, useC=FALSE, trace=trace, perm=TRUE, method="large", mc.rep=4e4), tolerance=tolerance)


}
