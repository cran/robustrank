stat=function(a,b,m,n) {
    #lr.stat(a,b,m,n)
    lik.null.stat(a,b,m,n)
}


lr.stat = function(a,b,m,n) {
    c=m-a
    d=n-b
    r=a+b
    s=c+d
    ifelse(a>0,a*log(a/r),0) + ifelse(b>0,b*log(b/r),0) + ifelse(c>0,c*log(c/s),0) + ifelse(d>0,d*log(d/s),0)
}

lik.null.stat = function(a,b,m,n) {
    c=m-a
    d=n-b
    r=a+b
    s=c+d
    # add neg sign so that large is good for alternative
    -(-lgamma(a+1) - lgamma(b+1) - lgamma(c+1) - lgamma(d+1) + ifelse(r>0,r*log(r),0) + ifelse(s>0,s*log(s),0))
}

binom.lr.test=function(a,b,m,n, useC=FALSE, verbose=FALSE, mc.rep=1e4) {
    
    N=m+n; c=m-a; d=n-b; r=a+b; s=c+d
    p.method=ifelse(choose(N,m)<=mc.rep, "exact", "Monte Carlo")
    if(verbose) print(p.method)
    
    z=stat(a,b,m,n)
    
    if(p.method=="exact"){ 
        indexes=combn(1:N, m)
        if (useC) {
            z.mc=.Call("binom_lr_test", a, b, m, n, as.integer(0), as.integer(indexes))
        } else {
            yy = c(rep(1,a), rep(0,c), rep(1,b), rep(0,d))
            z.mc=apply(indexes, 2, function(ind.x) {  
                ind.y=setdiff(1:N, ind.x)
                a.1=sum(yy[ind.x])
                b.1=sum(yy[ind.y])
                stat(a.1,b.1,m,n)
            })
        }
        
    } else if (p.method=="Monte Carlo") {
        #### save rng state before set.seed in order to restore before exiting this function
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") { set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
    
        set.seed(1)
        if(useC) {
            #z.mc=.Call("pair_wmw_test", X, Y, .corr, switch(method, "large.0"=-1, "large"=-2, "exact"=-3, "exact.0"=0, "exact.1"=1, "exact.2"=2, "exact.3"=3), as.integer(mc.rep), as.integer(1))
        } else {
            yy = c(rep(1,a), rep(0,c), rep(1,b), rep(0,d))
            z.mc=sapply(1:mc.rep, function(b) {
                ind.x=sample(1:N, m)
                ind.y=setdiff(1:N, ind.x)
                a.1=sum(yy[ind.x])
                b.1=sum(yy[ind.y])
                stat(a.1,b.1,m,n)
            })            
        }
            
        #### restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)     
    }        
    
    # z.mc may be NaN b/c var estimate may be less than 0
    # we treat NaN here as inf b/c a small value of var.est means large z
    z.mc[is.nan(z.mc)] = Inf 
    
    numerical.stability.margin=1e-10 # e.g. in C, ((6.6)/(100))+((0.45)/(25)) != ((6.84)/(100))+((0.39)/(25)), some permutations should have identical stat as the obs., but end up not
    pval = mean(z.mc>=z*(1+numerical.stability.margin)) # since z is lr, it is always negative
}


#library(kyotil)
#set.seed(1)
#n.2=10
#B=1e4
##pp=rbind(c(.5,.5), c(.5,.7), c(.5,.9))
#pp=rbind(c(.5,.5))
#for (n.1 in c(5)) { # 5, 10, 20
#    tab=NULL
#    for (j in 1:nrow(pp)) {
#        myprint(n.1, j)
#        p.1=pp[j,1]; p.2=pp[j,2]
#        X=cbind(rbinom(B, n.1, p.1), rbinom(B, n.2, p.2))    
#        pvals=sapply(1:B, function(i) binom.lr.test(X[i,1], X[i,2],n.1,n.2, mc.rep=1e4))
#        pow=mean(pvals < 0.05); pow
#        
#        tab=rbind(tab, c(p.1, p.2, pow))
#    }
#    colnames(tab)=c("p.1","p.2","new")
#}
#mytex(tab)
#
#
# n.1=5: size 0.0286 using null model likelihood
# n.1=5: size 0.028 using lr (corrected version )
