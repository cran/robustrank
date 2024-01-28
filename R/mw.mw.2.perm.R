#.corr=2; mc.rep=1e4; alternative = "two.sided"
mw.mw.2.perm=function(X, Y, Xprime, Yprime, .corr, mc.rep=1e4, alternative = c("two.sided", "less", "greater"), verbose=FALSE ) {
    
    alternative <- match.arg(alternative)
    m=length(X); l=length(Xprime); n=length(Yprime)
        
    n.perm = 2^m * choose(l+n, l)
    p.method=ifelse(n.perm<=mc.rep, "exact", "Monte Carlo")
    if(verbose) myprint(n.perm, p.method)
    
    if(p.method=="exact"){ 
        # indexes to resample paired data
        # binarys is adapted from Thomas Lumley https://stat.ethz.ch/pipermail/r-help/2000-February/010141.html
        binarys<-function(i,m) {a<-2^((m-1):0); b<-2*a; sapply(i,simplify="array", function(x) as.integer((x %% b)>=a))}
        xy.idx=binarys(0:(2^m-1), m)
        
        # indexes to resample unpaired data
        # each col is a vector of l+n, where 1 indicates X, 0 indicates Y
        xyprime.idx=combn(l+n,l,function(x) {out=rep(0,l+n); out[x]=1; out })        
        
        z.mc=.Call("mw_mw_2_perm", X, Y, Xprime, Yprime, .corr, as.integer(-1), as.integer(xy.idx), as.integer(xyprime.idx))
        
    } else if (p.method=="Monte Carlo") {
        #### save rng state before set.seed in order to restore before exiting this function
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (inherits(save.seed,"try-error")) { set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
        
        set.seed(1)
        z.mc=.Call("mw_mw_2_perm", X, Y, Xprime, Yprime, .corr, as.integer(mc.rep), as.integer(1), as.integer(1))
     }           
    
    # z.mc may be NaN b/c var estimate may be less than 0
    # we treat NaN here as inf b/c a small value of var.est means large z
    z=z.mc[1]
    z.mc=z.mc[-1]
    z.mc[is.nan(z.mc)] = Inf 
    
    numerical.stability.margin=1e-10 # e.g. in C, ((6.6)/(100))+((0.45)/(25)) != ((6.84)/(100))+((0.39)/(25)), some permutations should have identical stat as the obs., but end up not
    pval=switch(alternative,
       "two.sided" = mean(z.mc>=abs(z)*(1-numerical.stability.margin)) + mean(z.mc<= -abs(z)*(1-numerical.stability.margin)),
       "less" = mean(z.mc>=z*(1-sign(z)*numerical.stability.margin)), # X<Y
       "greater" = mean(z.mc<=z*(1+sign(z)*numerical.stability.margin)) #X>Y
    )
    
    res=list()
    class(res)=c("htest",class(res))
    res$statistic=z
    res$p.value=pval
    res$alternative=alternative
    return(res)
}
