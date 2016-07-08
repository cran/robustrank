# mode: var is used to examine variance estimators, test is used to perform testing
# comb2: use combine to compute stat, use wmw to compute permutation distribution
mod.wmw.test=function(X,Y, alternative = c("two.sided", "less", "greater"), correct = TRUE, perm=NULL, mc.rep=1e4, method=c("combine","comb2","fp","wmw","fplarge","nsm3"), trace=0 
    , mode=c("test","var"), useC=TRUE # develeopment parameters
) {
    
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    mode <- match.arg(mode)
    
    m=length(X)
    n=length(Y)
    N=m+n
    
    # if not converted to double, .Call may throw errors for integer X/Y
    X=as.double(X[!is.na(X)])
    Y=as.double(Y[!is.na(Y)])
    
    if(!correct) .corr=0 else .corr=switch(alternative, "two.sided" = 2, "greater" = -1, "less" = 1)         
    if(useC) {
        z=.Call("mod_wmw_test", X, Y, .corr, switch(method,"wmw"=0,"fp"=1,"fplarge"=2,"combine"=3,"comb2"=3,"nsm3"=4), as.integer(1), as.integer(1))
    } else {
        #the following R implementation may give a different result under i386 software on 64-bit machine, b/c C is always 64 bit, but R is 32 bit
        z=compute.Z(X,Y, alternative, correct, ifelse(method=="comb2","combine",method)); 
    }    
    
    # for developing manuscript
    if(mode=="var") return(compute.Z(X,Y, alternative, correct, ifelse(method=="comb2","combine",method), return.all=TRUE)) 
        
    p.method=NULL
    if(is.null(perm)) {
        if (min(m,n)>=20) p.method="asymptotic" 
    } else if (!perm) {
        p.method="asymptotic" 
    }
    if(is.null(p.method)) {
        if (choose(N,m)<=mc.rep) p.method="exact" else p.method="Monte Carlo"
    }
    if(trace==1) print(p.method)
    
    if (p.method=="asymptotic") {
       # one sided p values are opposite to wilcox.test because the way U is computed here as the rank of Y's
        pval=switch(alternative,
           "two.sided" = 2 * pnorm(abs(z), lower.tail=FALSE),
           "less" = pnorm(z, lower.tail=FALSE), # X<Y
           "greater" = pnorm(z, lower.tail=TRUE) #X>Y
        )
        
    } else {
        xy=c(X,Y)
        if(p.method=="exact"){        
            indexes=combn(1:N, m) 
            if (useC) {
                z.mc=.Call("mod_wmw_test", X, Y, .corr, switch(method,"wmw"=0,"fp"=1,"fplarge"=2,"combine"=3,"comb2"=0,"nsm3"=4), as.integer(0), indexes)
            } else {
                z.mc=apply(indexes, 2, function(ind.x) {
                    ind.y=setdiff(1:N, ind.x)
                    compute.Z(xy[ind.x], xy[ind.y], alternative=alternative, correct=correct, method=ifelse(method=="comb2","wmw",method))
                })
            }
## numerical instability
#t=2075
#print((z-z.mc[t]))
#print(which(z==z.mc))
            
        } else if (p.method=="Monte Carlo") {
            #### save rng state before set.seed in order to restore before exiting this function
            save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
            if (class(save.seed)=="try-error") { set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
    
            set.seed(1)
            if(useC) {
                z.mc=.Call("mod_wmw_test", X, Y, .corr, switch(method,"wmw"=0,"fp"=1,"fplarge"=2,"combine"=3,"comb2"=0,"nsm3"=4), as.integer(mc.rep), as.integer(1))
            } else {
                z.mc=sapply(1:mc.rep, function(b) {
                    ind.x=sample(1:N, m)
                    ind.y=setdiff(1:N,ind.x)
                    compute.Z(xy[ind.x], xy[ind.y], alternative=alternative, correct=correct, method=ifelse(method=="comb2","wmw",method))
                })            
            }
                
            #### restore rng state 
            assign(".Random.seed", save.seed, .GlobalEnv)     
        }        
        
        numerical.stability.margin=1e-10 # e.g. in C, ((6.6)/(100))+((0.45)/(25)) != ((6.84)/(100))+((0.39)/(25)), some permutations should have identical stat as the obs., but end up not
        pval=switch(alternative,
           "two.sided" = mean(z.mc>=abs(z)*(1-numerical.stability.margin)) + mean(z.mc<= -abs(z)*(1-numerical.stability.margin)),
           "less" = mean(z.mc>=z*(1-sign(z)*numerical.stability.margin)), # X<Y
           "greater" = mean(z.mc<=z*(1+sign(z)*numerical.stability.margin)) #X>Y
        )
    }
    
    pval
}


# R implementation, for debugging use
# ties are properly handled for WMW variance, but not properly handled for estimated Var
compute.Z=function(X,Y, alternative, correct, method, return.all=FALSE) {
    m=length(X)
    n=length(Y)
    N=m+n
    
    r <- rank(c(X, Y))
    U=(sum(r[(m+1):N]) - n*(n+1)/2)/m/n    
    correction <- switch(alternative,
       "two.sided" = sign(U-0.5) * 0.5/m/n,
       "greater" = 0.5/m/n,
       "less" = -0.5/m/n
    )
    if(correct) U=U-correction
    
    # V_mn
    var.wmw=1/12/m + 1/12/n + 1/12/m/n    
    # ties
    NTIES <- table(r)
    var.wmw <- var.wmw - sum(NTIES^3-NTIES)/(12*m*n*N*(N - 1))
    
    # estimated var
    Fx=ecdf(X); Fy=ecdf(Y)
    a=var(Fy(X)); b=var(Fx(Y))
    var.exact=a/m*(1-1/n) + b/n*(1-1/m) + mean(Fy(X))/m * mean(Fx(Y))/n
    var.large=a/m + b/n # large sample approx
    var.nsm3=a/m*(1-1/m) + b/n*(1-1/n) + mean(Fy(X))/m * mean(Fx(Y))/n
    
    # combined test
    var.3=min(var.wmw, var.exact)
    
    z=(U-0.5)/sqrt(switch(method, "fp"=var.exact, "fplarge"=var.large, "combine"=var.3, "wmw"=var.wmw, "nsm3"=var.nsm3))
    #print(z)
    
    if(return.all) return(c(U=U, var.wmw=var.wmw, var.large=var.large, var.exact=var.exact, var.nsm3=var.nsm3, z=z)) else return(z)
}


mod.wsr.test=function(X,Y, debug=FALSE) {
    
    m=length(X)
    n=length(Y)
    N=m+n
    
    delta=Y-X    
    abs.delta=abs(delta)
    sgn.delta=sign(delta)
    pos.delta=ifelse(delta>0,1,0)
    
    Fdelta=ecdf(delta)
    h.1=1-Fdelta(-delta)
#    # alternatively
#    h.1.f=function(i) (sum(delta+delta[i]>0) - sum(delta[i]>0))/(length(delta)-1)
#    h.1=numeric(m)
#    for (i in 1:m) h.1[i]=h.1.f(i)    
    
    r.1 <- rank(abs.delta)    
    W.plus = mean(pos.delta*r.1)
    W.plus.2=mean(sgn.delta*r.1) # definition in van der Vaart
    var.W.plus = m * var(h.1) # var(W.plus.2) is 4 x var(W.plus)
    var.W.plus.0 = (m+1)*(2*m+1)/(24*m)
    
    if(debug) {
        return(c(W.plus=W.plus, var.W.plus=var.W.plus, var.W.plus.0=var.W.plus.0))
    }    
    
    # p value
    pchisq( (W.plus-m/4)^2/var.W.plus, df=1, lower.tail=FALSE)
        
}
