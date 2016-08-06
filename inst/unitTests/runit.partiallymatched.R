library("RUnit")
library("robustrank")


test.pm <- function() {

RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-3
if(file.exists("C:/_checkReproducibility")) tolerance=1e-6

dat=sim.partially.matched(m=10,n.x=10,n.y=10,distr="mixnormal",params=c(p.1=0.3,p.2=0.3),seed=1)
X=dat$X; Y=dat$Y; Xprime=dat$Xprime; Yprime=dat$Yprime        
checkEqualsNumeric(pm.wilcox.test(Y, X, Yprime, Xprime, trace=1)[,2], c(0.6870638, 0.6195517, 0.4839197, 0.4715029, 0.6084627, 0.6596283, 0.3746707, 0.4446945, 0.4015319, 0.4671714, 0.8810059, 0.5866820), tolerance=tolerance)
#checkEqualsNumeric(pm.wilcox.test(Y, X, NULL, Xprime, useC=FALSE, trace=1)["mw.mw.3",1], pm.wilcox.test(Y, X, NULL, Xprime, correct=FALSE, useC=TRUE , trace=1)^2, tolerance=tolerance)



#checkEqualsNumeric(sign.mw.test(X, Y, Xprime), 0.1422232, tolerance=tolerance)
#checkEqualsNumeric(mw.mw.test(Y, X, Xprime), 0.6065135, tolerance=tolerance)

#
#dat=sim.partially.matched(m=10,n.x=10,n.y=12,distr="mixnormal",params=c(p.1=0.3,p.2=0.7),seed=1)
#X=dat$X; Y=dat$Y; Xprime=dat$Xprime; Yprime=dat$Yprime
#checkEqualsNumeric(sign.mw.test(X, Y, Xprime), 0.3007986, tolerance=tolerance)
#checkEqualsNumeric(mw.mw.test(Y, X, Xprime), 0.01713436, tolerance=tolerance)
#checkEqualsNumeric(mw.mw.test(X, Y, Yprime), 0.002023006, tolerance=tolerance)



}
