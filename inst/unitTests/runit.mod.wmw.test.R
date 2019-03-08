test.mod.wmw.test <- function() {

  suppressWarnings(RNGversion("3.5.0"))
library("RUnit")
library("robustrank")
RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-2
if(file.exists("D:/gdrive/3software/_checkReproducibility")) tolerance=1e-6 # under i386, the U computed by R and by C are not identical and that causes somep problem
verbose=FALSE


# WMW
set.seed(1); X=rnorm(50); Y=rnorm(75,mean=1,sd=1.4)
checkEqualsNumeric(mod.wmw.test(Y, X, correct=TRUE,  method="wmw", perm=FALSE, verbose=verbose), wilcox.test(Y, X, correct=TRUE)$p.value,  tolerance=tolerance)
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="wmw", perm=FALSE, verbose=verbose), wilcox.test(Y, X, correct=FALSE)$p.value, tolerance=tolerance)
checkEqualsNumeric(mod.wmw.test(Y, X, correct=TRUE,  method="wmw", perm=FALSE, alternative="greater", verbose=verbose), wilcox.test(Y, X, correct=TRUE, alternative="greater")$p.value,  tolerance=tolerance)
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE,  method="wmw", perm=FALSE, alternative="greater", verbose=verbose), wilcox.test(Y, X, correct=FALSE, alternative="greater")$p.value,  tolerance=tolerance)
Y=c(Y,X[1:10]) # create ties
checkEqualsNumeric(mod.wmw.test(Y, X, correct=TRUE,  method="wmw", perm=FALSE, verbose=verbose), wilcox.test(Y, X, correct=TRUE)$p.value,  tolerance=tolerance)

# FP
set.seed(1); X=rnorm(5); Y=rnorm(10,mean=1,sd=1.4)
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="fp", perm=FALSE, verbose=verbose), 0.08068789, tolerance=tolerance) # exact
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="fp", perm=TRUE, mc.rep=1e4, verbose=verbose), 0.1165501, tolerance=tolerance) # exact
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="fp", perm=TRUE, mc.rep=1e2, verbose=verbose), 0.06, tolerance=tolerance) # MC
Y=c(Y,X[1:3]) # create ties NOTE that the R implementation of the following is not correct
checkEqualsNumeric(mod.wmw.test(Y, X, correct=TRUE,  method="fp", perm=FALSE, verbose=verbose), 0.35205,  tolerance=tolerance)
checkEqualsNumeric(mod.wmw.test(Y, X, correct=TRUE,  method="fp", perm=TRUE, mc.rep=1e4, verbose=verbose),  0.3311158,  tolerance=tolerance)
mod.wmw.test(Y, X, correct=TRUE,  method="fp", perm=TRUE, mc.rep=1e1, verbose=verbose, useC=FALSE)


# same as the next commented block but with hard-coded values
set.seed(1); X=rnorm(5); Y=rnorm(10,mean=1,sd=1.4)
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=FALSE, alternative="less", verbose=verbose), 0.9568903, tolerance=tolerance) 
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=FALSE, verbose=verbose), 0.08621937, tolerance=tolerance) 
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=TRUE, alternative="less", verbose=verbose), 0.9407259, tolerance=tolerance) 
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=TRUE, verbose=verbose), 0.1252081, tolerance=tolerance) 
# ties
set.seed(1); X=rnorm(5); Y=rnorm(10,mean=1,sd=1.4); Y=c(Y,X[1:3]) 
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=FALSE, alternative="less", verbose=verbose), 0.836123, tolerance=tolerance) 
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=FALSE, verbose=verbose), 0.327754, tolerance=tolerance) 
checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=TRUE, alternative="less", verbose=verbose), 0.8288982, tolerance=tolerance) 

## the following, except the 4th should hold true, but it is slow and requires NSM3
#library(NSM3)
#set.seed(1); X=rnorm(5); Y=rnorm(10,mean=1,sd=1.4)
#checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=FALSE, alternative="less"), pFligPoli(Y,X,method="Asymptotic")$p.val, tolerance=tolerance) 
#checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=FALSE), pFligPoli(Y,X,method="Asymptotic")$two.sided, tolerance=tolerance) 
#checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=TRUE, alternative="less"), pFligPoli(Y,X,method="Exact")$p.val, tolerance=tolerance) 
#checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=TRUE), pFligPoli(Y,X,method="Exact")$two.sided, tolerance=tolerance) # this does not hold b/c the def of two-sided p value are different between the two functions
## ties
#set.seed(1); X=rnorm(5); Y=rnorm(10,mean=1,sd=1.4); Y=c(Y,X[1:3]) 
#checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=FALSE, alternative="less"), pFligPoli(Y,X,method="Asymptotic")$p.val, tolerance=tolerance) 
#checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=FALSE), pFligPoli(Y,X,method="Asymptotic")$two.sided, tolerance=tolerance) 
#checkEqualsNumeric(mod.wmw.test(Y, X, correct=FALSE, method="nsm3", perm=TRUE, alternative="less"), pFligPoli(Y,X,method="Exact")$p.val, tolerance=tolerance) 



#library("coin"); set.seed(1); coin::wilcox_test(v~g, data.frame(v=c(X,Y),g=c(rep("1",5),rep("2",10))), distribution = approximate(B = 10000)) # stat is the same, but p value different, which is ok

}


#library("robustrank")
#set.seed(1)
## no ties
#X=c(1,2,3); Y=c(2.1,4,5)
#mod.wmw.test(X, Y, mc.rep=2, perm=F, method="wmw", verbose=verbose)
#wilcox.test(Y, X, correct=TRUE, exact=F)
## ties
#X=c(1,2,3); Y=c(2,4,5)
#mod.wmw.test(X, Y, mc.rep=2, perm=F, method="wmw")
#wilcox.test(Y, X, correct=TRUE, exact=F)
#
