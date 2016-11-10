library("RUnit")
library("robustrank")

test.pm <- function() {

RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-3
if(file.exists("C:/_checkReproducibility")) tolerance=1e-6

dat=sim.partially.matched(m=10,n.x=10,n.y=10,distr="mixnormal",params=c(p.1=0.3,p.2=0.3),seed=1)
X=dat$X; Y=dat$Y; Xprime=dat$Xprime; Yprime=dat$Yprime       

# call wilcox.test
checkEqualsNumeric(pm.wilcox.test(X,Y, trace=1,method="all")$p.value, 0.4316406, tolerance=tolerance) 
# switch two samples
checkEqualsNumeric(pm.wilcox.test(X,Y,Xprime,NULL,trace=1,method="all")[,2], pm.wilcox.test(Y,X,NULL,Xprime,trace=1,method="all")[,2], tolerance=tolerance)  
checkEqualsNumeric(pm.wilcox.test(X,Y,Xprime,NULL,trace=1,method="all")[,1], -pm.wilcox.test(Y,X,NULL,Xprime,trace=1,method="all")[,1], tolerance=tolerance)  
# create extra samples
X1=X; X1[1]=NA
checkEqualsNumeric(pm.wilcox.test(X1,Y,trace=1,method="all")[,2], c(0.23337562,NA,0.09152589,NA,0.20020911,NA,0.12883094,NA,0.16885113,NA), tolerance=tolerance)

# useC. Note that it does not implement var.U.p.0.high
checkEqualsNumeric(0.1133, pm.wilcox.test(X,Y, NULL, Yprime, useC=TRUE , trace=1,method="all")^2, tolerance=tolerance) 

# SR-MW, default
checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, Yprime)$statistic, 0.5270253,tolerance=tolerance) 
# MW-MW
checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, Yprime, method="MW-MW")$statistic, 0.5324883,tolerance=tolerance) 

# two different input methods
checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, method="MW-MW")$statistic, pm.wilcox.test(c(X,Xprime), c(Y,rep(NA,length(Xprime))), method="MW-MW")$statistic,tolerance=tolerance)

# two-sided vs one-sided
checkEqualsNumeric(pm.wilcox.test(Y,X, NULL, Xprime, trace=1,method="all",alternative="two.sided")[,2], c(0.6870638, 0.6195517, 0.4839197, 0.4715029, 0.6022422, 0.6596283, 0.3709689, 0.4446945, 0.4000556, 0.4671714), tolerance=tolerance)
checkEqualsNumeric(pm.wilcox.test(Y,X, NULL, Xprime, trace=1,method="all",alternative="greater")[,2],   c(0.6870638, 0.6195517, 0.7580401, 0.4715029, 0.6022422, 0.6596283, 0.8145156, 0.4446945, 0.4000556, 0.4671714), tolerance=tolerance)
checkEqualsNumeric(pm.wilcox.test(Y,X, NULL, Xprime, trace=1,method="all",alternative="less")[,2],      c(0.6870638, 0.6195517, 0.2419599, 0.4715029, 0.6022422, 0.6596283, 0.1854844, 0.4446945, 0.4000556, 0.4671714), tolerance=tolerance)

checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, Yprime, useC=FALSE, trace=1,method="all",alternative="two.sided")[,2], c(0.5981760, 0.5943879, 0.4795192),tolerance=tolerance)
checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, Yprime, useC=FALSE, trace=1,method="all",alternative="greater")[,2],   c(0.7009120, 0.7028061, 0.4795192),tolerance=tolerance)
checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, Yprime, useC=FALSE, trace=1,method="all",alternative="less")[,2],      c(0.2990880, 0.2971939, 0.4795192),tolerance=tolerance)
# is.switched and alternative
checkEqualsNumeric(pm.wilcox.test(X,Y, NULL, Yprime, useC=FALSE, trace=1,method="all",alternative="less")[,2], pm.wilcox.test(Y,X, Yprime, NULL, useC=FALSE, trace=1,method="all",alternative="greater")[,2], tolerance=tolerance)


## additional functions
#checkEqualsNumeric(sign.mw.test(X, Y, Xprime), 0.1422232, tolerance=tolerance)
#checkEqualsNumeric(mw.mw.test(Y, X, Xprime), 0.6065135, tolerance=tolerance)

#dat=sim.partially.matched(m=10,n.x=10,n.y=12,distr="mixnormal",params=c(p.1=0.3,p.2=0.7),seed=1)
#X=dat$X; Y=dat$Y; Xprime=dat$Xprime; Yprime=dat$Yprime
#checkEqualsNumeric(sign.mw.test(X, Y, Xprime), 0.3007986, tolerance=tolerance)
#checkEqualsNumeric(mw.mw.test(Y, X, Xprime), 0.01713436, tolerance=tolerance)
#checkEqualsNumeric(mw.mw.test(X, Y, Yprime), 0.002023006, tolerance=tolerance)

}
