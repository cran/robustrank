# generates m pairs of X and Y and independent m X samples (called Xprime) and n Y samples (called Yprime)
sim.partially.matched=function(m, n.x, n.y, distr=c("normal","logistic","student","mixnormal","gamma"), params, seed){
    set.seed(seed)    
    N=m+n.x+n.y    
    distr=match.arg(distr)
    if (distr %in% c("normal","student","logistic")) {
        loc.2=params["loc.2"]
        scale.2=params["scale.2"]
        rho=params["rho"]; if (is.na(rho)) rho=0        
        if (distr=="normal") {
            dat.all=mvtnorm::rmvnorm(n = N, mean=c(0,loc.2),  sigma=matrix(c(1,scale.2*rho,scale.2*rho,scale.2^2),nrow=2))        
        } else if (distr=="student") {        
            dat.all=mvtnorm::rmvt   (n = N, delta=c(0,loc.2), sigma=matrix(c(1,scale.2*rho,scale.2*rho,scale.2^2),nrow=2), df=4)        
        } else if (distr=="logistic") {
            dat.all=rbilogistic(N, 0, loc.2, scale.1=1, scale.2=scale.2, rho) 
        }
    } else if (distr=="mixnormal") {        
        p.1=params["p.1"]
        p.2=params["p.2"]
        sd.n=params["sd.n"]; if (is.na(sd.n)) sd.n=1
        sd2=params["sd2"];   if (is.na(sd2)) sd2=0.5
        dat.all=rnorm(N, sd=sd.n) + cbind(rmixnorm (N, mix.p=p.1, mu1=0, mu2=2, sd1=.5, sd2=sd2), rmixnorm (N, mix.p=p.2, mu1=0, mu2=2, sd1=.5, sd2=sd2))            
    } else if (distr=="gamma") {        
        loc.2=params["loc.2"]
        shape.1=params["shape.1"]
        shape.2=params["shape.2"]
        rate.1=params["rate.1"]
        rate.2=params["rate.2"]
        rho=params["rho"]
        dat.all=rbigamma(n=N, shape.1, shape.2, rate.1, rate.2, rho=rho) 
        dat.all[,2]=dat.all[,2]+loc.2
    }    
    index.1=if(n.x==0) NULL else m+1:n.x
    index.2=if(n.y==0) NULL else m+n.x+1:n.y
    list(X=dat.all[1:m,1],            Y=dat.all[1:m,2], 
         Xprime=dat.all[index.1,1],   Ymissing=dat.all[index.1,2], 
         Xmissing=dat.all[index.2,1], Yprime=dat.all[index.2,2])         
}


r2sample=function(m, n, distr=c("normal","logistic","student","mixnormal"), params, seed){
    distr=match.arg(distr)
    dat=sim.partially.matched(m=m, n.x=0, n.y=n, distr, params, seed)[c("X","Yprime")]
    names(dat)=c("X","Y")
    dat
}
