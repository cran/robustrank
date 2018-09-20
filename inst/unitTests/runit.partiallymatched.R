library("RUnit")
library("robustrank")

test.pm <- function() {

RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-3
if(file.exists("C:/_checkReproducibility")) tolerance=1e-6

dat=sim.partially.matched(m=10,n.x=10,n.y=10,distr="mixnormal",params=c(p.1=0.3,p.2=0.3),seed=1)
X=dat$X; Y=dat$Y; Xprime=dat$Xprime; Yprime=dat$Yprime       

# call wilcox.test
checkEqualsNumeric(pm.wilcox.test(X,Y, verbose=1,method="all")$p.value, 0.4316406, tolerance=tolerance) 
# switch two samples
checkEqualsNumeric(pm.wilcox.test(X,Y,Xprime,NULL,verbose=1,method="all")[,2], pm.wilcox.test(Y,X,NULL,Xprime,verbose=1,method="all")[,2], tolerance=tolerance)  
checkEqualsNumeric(pm.wilcox.test(X,Y,Xprime,NULL,verbose=1,method="all")[,1], -pm.wilcox.test(Y,X,NULL,Xprime,verbose=1,method="all")[,1], tolerance=tolerance)  
# create extra samples
X1=X; X1[1]=NA
checkEqualsNumeric(pm.wilcox.test(X1,Y,verbose=1,method="all")[,2], c(0.3089915,NA,NA,0.1310076,0.3319130,NA,NA,0.2124292,NA,0.1787910), tolerance=tolerance)

# useC. Note that it does not implement var.U.p.0.high
checkEqualsNumeric(0.1133, pm.wilcox.test(X,Y, NULL, Yprime, useC=TRUE , verbose=1,method="all")^2, tolerance=tolerance) 

# SR-MW, default
checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, Yprime)$statistic, 0.523538,tolerance=tolerance) 
# MW-MW
checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, Yprime, method="MW-MW")$statistic, 0.5242973,tolerance=tolerance) 

# two different input methods
checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, method="MW-MW")$statistic, pm.wilcox.test(c(X,Xprime), c(Y,rep(NA,length(Xprime))), method="MW-MW")$statistic,tolerance=tolerance)

# two-sided vs one-sided
checkEqualsNumeric(pm.wilcox.test(Y,X, NULL, Xprime, verbose=1,method="all",alternative="two.sided")[,2], c(0.7171121,0.6525842,0.5019282,0.5149631,0.6299976,0.6776026,0.4823612,0.4429210,0.5073727,0.4102135), tolerance=tolerance)
checkEqualsNumeric(pm.wilcox.test(Y,X, NULL, Xprime, verbose=1,method="all",alternative="greater")[,2],   c(0.7171121,0.6525842,0.5019282,0.7425185,0.6299976,0.6776026,0.4823612,0.4429210,0.5073727,0.7948933), tolerance=tolerance)
checkEqualsNumeric(pm.wilcox.test(Y,X, NULL, Xprime, verbose=1,method="all",alternative="less")[,2],      c(0.6564530,0.5863457,0.4421231,0.2269602,0.5721599,0.6388692,0.4087645,0.3597855,0.4288016,0.1670716), tolerance=tolerance)

checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, Yprime, useC=FALSE, verbose=1,method="all",alternative="two.sided")[c("sr.mw.l0","mw.mw.l0","mw.mw.00"),2], c(0.6005999,0.6000718,0.4887205),tolerance=tolerance)
checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, Yprime, useC=FALSE, verbose=1,method="all",alternative="greater")[c("sr.mw.l0","mw.mw.l0","mw.mw.00"),2],   c(0.6745656,0.6770202,0.5264579),tolerance=tolerance)
checkEqualsNumeric(pm.wilcox.test(X,Y, Xprime, Yprime, useC=FALSE, verbose=1,method="all",alternative="less")[c("sr.mw.l0","mw.mw.l0","mw.mw.00"),2],      c(0.2737552,0.2723921,0.4349615),tolerance=tolerance)
# is.switched and alternative
checkEqualsNumeric(pm.wilcox.test(X,Y, NULL, Yprime, useC=FALSE, verbose=1,method="all",alternative="less")[,2], pm.wilcox.test(Y,X, Yprime, NULL, useC=FALSE, verbose=1,method="all",alternative="greater")[,2], tolerance=tolerance)


## additional functions
#checkEqualsNumeric(sign.mw.test(X, Y, Xprime), 0.1422232, tolerance=tolerance)
#checkEqualsNumeric(mw.mw.test(Y, X, Xprime), 0.6065135, tolerance=tolerance)

#dat=sim.partially.matched(m=10,n.x=10,n.y=12,distr="mixnormal",params=c(p.1=0.3,p.2=0.7),seed=1)
#X=dat$X; Y=dat$Y; Xprime=dat$Xprime; Yprime=dat$Yprime
#checkEqualsNumeric(sign.mw.test(X, Y, Xprime), 0.3007986, tolerance=tolerance)
#checkEqualsNumeric(mw.mw.test(Y, X, Xprime), 0.01713436, tolerance=tolerance)
#checkEqualsNumeric(mw.mw.test(X, Y, Yprime), 0.002023006, tolerance=tolerance)

}
