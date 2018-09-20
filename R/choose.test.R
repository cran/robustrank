# perform power study to decide which method to use
choose.test=function(Xpaired, Ypaired, Xextra=NULL, Yextra=NULL, mc.rep=1000) {
    m=length(Xpaired)
    n=length(Yextra)
    l=length(Xextra)
    if(n==0 & l==0 ) stop("There are no unpaired observations.")
    
    scale.2=sd(Ypaired)/sd(Xpaired)
    loc.2=(mean(Ypaired)-mean(Xpaired))/sd(Xpaired)
    rho=cor(Ypaired,Xpaired, method="pearson")
    
    res=sapply (1:mc.rep, function(seed) {
        dat=sim.partially.matched(m=m,n.x=l,n.y=n,distr="normal",params=c(loc.2=loc.2,rho=rho,scale.2=scale.2),seed=seed)
        test.sr=wilcox.test(dat$X, dat$Y, paired=T,conf.int=F) # sign rank test using only paired data
        test.pm=pm.wilcox.test (dat$X, dat$Y, Xextra=dat$Xprime, Yextra=dat$Yprime, method="all", mode="power.study", correct=F)
        if(n==0 | l==0) c(SR=test.sr$p.val, test.pm) else test.pm
    })
    apply(res, 1, function(x) round(mean(x<0.05)*100))    
}
