# X=dat$X; Y= dat$Y; useC=TRUE; alternative = "two.sided"; correct = FALSE; perm=TRUE; mc.rep=101; trace=1
multinom.test=function(X,Y
    , alternative = c("two.sided", "less", "greater"), correct = FALSE, perm=NULL, mc.rep=1e4
    , method=c("exact.2","large.0","large","exact","exact.0","exact.1","exact.3"), trace=0
    , mode=c("test","var")
    , useC=TRUE)
{
    
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    mode <- match.arg(mode)
    
    if(!correct) .corr=0 else .corr=switch(alternative, "two.sided" = 2, "greater" = -1, "less" = 1)         
    # remove pairs that have missing values
    # as.double converts int to num (double), calls to C function may throw errors for integer X/Y
    stopifnot(length(X)==length(Y))    
    X=as.double(X[!is.na(X) & !is.na(Y)])
    Y=Y[!is.na(X) & !is.na(Y)]
    for (i in 1:length(X)) Y[[i]]=as.double(Y[[i]][!is.na(Y[[i]])])
    len.Y=sapply(Y,length)
    if (any(0==len.Y)) {
        X=X[!0==len.Y]
        Y=Y[!0==len.Y]
    }
    m=length(X)
    n=m
    len.Y=sapply(Y,length)
    
    # for developing manuscript
    #if(mode=="var") return(compute.pair.wmw.Z(X,Y, alternative, correct, method, return.all=TRUE)) 
    
    W=compute.multinom.stat(X,Y, alternative, correct, FALSE)
                
    p.method=NULL # asymptotic, exact, Monte Carlo
    if(is.null(perm)) {
        if (min(m,n)>=50) p.method="asymptotic" 
    } else if (!perm) {
        p.method="asymptotic" 
    }
    if(is.null(p.method)) p.method="Monte Carlo"
    if(trace==1) print(p.method)
    
    if (p.method=="asymptotic") {
       # one sided p values are opposite to wilcox.test because the way U is computed here as the rank of Y's
       stop("not implemented yet")
        pval=switch(alternative,
           "two.sided" = 2 * pnorm(abs(W), lower.tail=FALSE),
           "less" = pnorm(W, lower.tail=FALSE), # X<Y
           "greater" = pnorm(W, lower.tail=TRUE) #X>Y
        )
        
    } else {
        #### save rng state before set.seed in order to restore before exiting this function
        save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
        if (class(save.seed)=="try-error") { set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
    
        # permutation method like that in Elena except two-sided
        if(useC) {
            set.seed(1) # so that it can be repeated
            W.mc=.Call("multinom_test", X, unlist(Y), len.Y, .corr, 1, as.integer(mc.rep))
        } else {    
            set.seed(1) # so that it can be repeated
            W.mc=sapply (1:mc.rep, function (b) compute.multinom.stat(X,Y, alternative, correct, TRUE) )
        }        
            
        #### restore rng state 
        assign(".Random.seed", save.seed, .GlobalEnv)     
            
        
        # W.mc may be NaN b/c var estimate may be less than 0
        # we treat NaN here as inf b/c a small value of var.est means large W
        W.mc[is.nan(W.mc)] = Inf 
        
        numerical.stability.margin=1e-10 # e.g. in C, ((6.6)/(100))+((0.45)/(25)) != ((6.84)/(100))+((0.39)/(25)), some permutations should have identical stat as the obs., but end up not
        pval=switch(alternative,
           "two.sided" = mean(W.mc>=abs(W)*(1-numerical.stability.margin)) + mean(W.mc<= -abs(W)*(1-numerical.stability.margin)),
           "less" = mean(W.mc>=W*(1-sign(W)*numerical.stability.margin)), # X<Y
           "greater" = mean(W.mc<=W*(1+sign(W)*numerical.stability.margin)) #X>Y
        )
    }
    
    pval
}

# R implementation, for debugging use
# X is a vector of length m, Y a list of length m, sel is a vector of length m
compute.multinom.stat=function(X,Y, alternative, correct, to.sample) {
    
    m=length(X)
    stat=0
    mu=0
    for (i in 1:m) {
        all=c(X[i], Y[[i]])
        if (to.sample) {
            sel=sample(1:length(all), 1)
            #stat=stat+ sum(rank(all)[sel])-1 # this stat matches wilcox.test(all[sel], all[-sel], correct=F)$statistic             
            # there does not seem to be a need to rank since we are just picking a random number out
            stat=stat+ sum(sel)-1 
        } else {
            sel=1
            stat=stat+ sum(rank(all)[sel])-1 # this stat matches wilcox.test(all[sel], all[-sel], correct=F)$statistic             
        }
        mu=mu+1*(length(all)-1)/2
    }
    stat=stat-mu
    stat
    
}    


#        # one-sided test is not of main interest
#        # permutation method as in Elena, one-sided: rank within pair, then sum
#        W=0; for (i in 1:nInfant) W=W+sum(rank(c(infant[i], mother[[i]]))[-1]) 
#        W.B=sapply (1:B, function (b) {
#            set.seed(b)
#            # for each infant/mother pair, randomly choose one to be infant
#            stat=0
#            for (i in 1:nInfant) {
#                all=c(infant[i], mother[[i]])
#                sel=sample(1:length(all),1)
#                stat=stat+ sum(rank(all)[-sel]) 
#            }
#            stat
#        })    
#        p.elena=mean(W.B>=W); p.elena
#        cumsum(rev(table(W.B)))/1e4
        
