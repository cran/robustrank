test.mw.mw.2.perm <- function() {


library("RUnit")
library("robustrank")
RNGkind("Mersenne-Twister", "Inversion")
tolerance=1e-3
if(file.exists("C:/_checkReproducibility")) tolerance=1e-6


dat=sim.partially.matched(m=6,n.x=3,n.y=2,distr="normal",params=c(loc.2=1,rho=.5,scale.2=1),seed=1)
X=dat$X; Y=dat$Y; Xprime=dat$Xprime; Yprime=dat$Yprime       
out=mw.mw.2.perm(X, Y, Xprime, Yprime, 2, verbose=T)
checkEqualsNumeric(c(out$stat, out$p.val), c(2.819392, 0.003125), tolerance=tolerance)



# pm.wilcox.test result is potentially different b/c it uses exact variance estimate of U.p. if we comment out var.U.p.0=var.U.p.0.high, we could compare the two
dat=sim.partially.matched(m=10,n.x=7,n.y=5,distr="mixnormal",params=c(p.1=0.3,p.2=0.3),seed=1)
X=dat$X; Y=dat$Y; Xprime=dat$Xprime; Yprime=dat$Yprime       

z=.Call("mw_mw_2_perm", X, Y, Xprime, Yprime, .corr=as.integer(2), 0, as.integer(1), as.integer(1)); z
out=pm.wilcox.test(X,Y,Xprime,Yprime, method="all",correct=T)["mw.mw.2.perm",1] 
checkEqualsNumeric(z,out, tolerance=tolerance)

z=.Call("mw_mw_2_perm", X, Y, Xprime, Yprime, .corr=0, 0, as.integer(1), as.integer(1)); z
checkEqualsNumeric(z,0.04116935, tolerance=tolerance)


dat=sim.partially.matched(m=10,n.x=12,n.y=5,distr="mixnormal",params=c(p.1=0.3,p.2=0.3),seed=1)
X=dat$X; Y=dat$Y; Xprime=dat$Xprime; Yprime=dat$Yprime       

z=.Call("mw_mw_2_perm", X, Y, Xprime, Yprime, .corr=as.integer(2), 0, as.integer(1), as.integer(1)); z
checkEqualsNumeric(z,-2.647543, tolerance=tolerance)

z=.Call("mw_mw_2_perm", X, Y, Xprime, Yprime, .corr=0, 0, as.integer(1), as.integer(1)); z
checkEqualsNumeric(z,-2.714718, tolerance=tolerance)


}
