# generates m pairs of X and Y and independent m X samples (called Xprime) and n Y samples (called Yprime)
sim.partially.matched=function(m, n.x, n.y, distr=c("normal","logistic","student","mixnormal","gamma","lognormal","beta","uniform","hybrid1","hybrid2","doublexp"), params, seed){
    set.seed(seed)    
    N=m+n.x+n.y    
    distr=match.arg(distr)
    if (distr %in% c("normal","student","logistic","lognormal","doublexp")) {
        loc.2=params["loc.2"]
        scale.2=params["scale.2"]
        rho=params["rho"]; if (is.na(rho)) rho=0        
        if (distr=="normal") {
            dat.all=mvtnorm::rmvnorm(n = N, mean=c(0,loc.2),  sigma=matrix(c(1,scale.2*rho,scale.2*rho,scale.2^2),nrow=2))        
        } else if (distr=="lognormal") {
            dat.all=exp(mvtnorm::rmvnorm(n = N, mean=c(0,0),  sigma=matrix(c(1,scale.2*rho,scale.2*rho,scale.2^2),nrow=2)) )
            dat.all[,2]=dat.all[,2]+loc.2
        } else if (distr=="student") {        
            dat.all=mvtnorm::rmvt   (n = N, delta=c(0,loc.2), sigma=matrix(c(1,scale.2*rho,scale.2*rho,scale.2^2),nrow=2), df=4)        
        } else if (distr=="logistic") {
            dat.all=rbilogistic(N, 0, loc.2, scale.1=1, scale.2=scale.2, rho) 
        } else if (distr=="doublexp") { # double exponential also known as Laplace
            dat.all=rbidoublexp(N, 0, loc.2, scale.1=1, scale.2=scale.2, rho) 
        }
    } else if (distr=="mixnormal") {        
        p.1=params["p.1"]
        p.2=params["p.2"]
        sd.n=params["sd.n"]; if (is.na(sd.n)) sd.n=1
        sd1=params["sd1"];   if (is.na(sd1)) sd1=0.5
        sd2=params["sd2"];   if (is.na(sd2)) sd2=0.5
        dat.all=rnorm(N, sd=sd.n) + cbind(rmixnorm (N, mix.p=p.1, mu1=0, mu2=2, sd1=sd1, sd2=sd2), rmixnorm (N, mix.p=p.2, mu1=0, mu2=2, sd1=sd1, sd2=sd2))            
    } else if (distr=="gamma") {        
        loc.2=params["loc.2"]
        shape.1=params["shape.1"]
        shape.2=params["shape.2"]
        rate.1=params["rate.1"]
        rate.2=params["rate.2"]
        rho=params["rho"]
        dat.all=rbigamma(n=N, shape.1, shape.2, rate.1, rate.2, rho=rho) 
        dat.all[,2]=dat.all[,2]+loc.2
    } else if (distr=="beta") {        
        shape1.1=params["shape1.1"]
        shape2.1=params["shape2.1"]
        shape1.2=params["shape1.2"]
        shape2.2=params["shape2.2"]
        dat.all=cbind(rbeta(N,shape1.1,shape2.1),rbeta(N,shape1.2,shape2.2))
    } else if (distr=="uniform") {        
        loc.2=params["loc.2"]
        dat.all=cbind(runif(N),runif(N))
        dat.all[,2]=dat.all[,2]+loc.2
    } else if (distr=="hybrid1") {
        # X gamma Y logistic
        loc.2=params["loc.2"]
        toggle=params["toggle"]
        if (toggle==0) {
            dat.all=cbind(rgamma(N,shape=2,rate=0.5), rlogis(N)+loc.2)
        } else if (toggle==1) {
            dat.all=cbind(rlogis(N)+loc.2, rgamma(N,shape=2,rate=0.5))
        }
    } else if (distr=="hybrid2") {
        # X lognormal Y logistic
        loc.2=params["loc.2"]
        toggle=params["toggle"]
        if (toggle==0) {
            dat.all=cbind(exp(rnorm(N)), rlogis(N)+loc.2)
        } else if (toggle==1) {
            dat.all=cbind(rlogis(N)+loc.2, exp(rnorm(N)))
        }
        
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


sim.paired.with.replicates=function(m, meanRatio, sdRatio, within.sd, type, hyp, distr, seed){    
    set.seed(seed)
    ratios=round(rnorm(n=m,mean=meanRatio,sd=sdRatio))    
    tmp=sim.partially.matched(m=m,n.x=0,n.y=0,distr=distr,params=c(loc.2=0,rho=0.5,scale.2=1),seed=seed)
    mother=list()
    
    if(type==1) {
        # infant seq and mother sequences are exchangeable
        infant=numeric(m)
        for (i in 1:m){
            tmp.1=tmp$Y[i]+rnorm(n=ratios[i]+1, sd=within.sd)
            if (hyp==0) {   
                # null: choose a random one to be infant: 1 is a random one
                sel=1
            } else if (hyp==1) {
                # filtering model: choose smaller ones to be infant with higher prob, to be consistent with the direction of 1-sided test
                tmp.1=sort(tmp.1)
                sel=rbern(1, 1/(1:length(tmp.1))^1.05, generalized=TRUE)                
            } else if (hyp==2) {
                # location shift model
                sel=1
                tmp.1[sel]=tmp.1[sel]-.5
            } else stop("wrong hyp")
            infant[i]=tmp.1[sel]
            mother=c(mother, list(tmp.1[-sel]))            
        }
    } else if (type==2) {
        # null only
        # infant seq and mother sequences are not exchangeable
        infant=tmp$X
        for (i in 1:m){
            tmp.1=tmp$Y[i]+rnorm(n=ratios[i], sd=within.sd)
            if (hyp==0) {
                infant[i]=infant[i]+rnorm(n=1, sd=within.sd)
                mother=c(mother, list(tmp.1))
            } else stop("for type 2, only hyp 0 is implemented")
        }        
    } else stop ("wrong type")
    list(X=infant, Y=mother)
}
