# X=dat$X; Y= dat$Y; useC=TRUE; alternative = "two.sided"; correct = FALSE; perm=TRUE; mc.rep=3; verbose=1
wmw.paired.replicates.test=function(X,Y
    , alternative = c("two.sided", "less", "greater"), correct = FALSE, perm=NULL, mc.rep=1e4
    , method=c("exact.2","large.0","large","exact","exact.0","exact.1","exact.3"), verbose=FALSE
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
    
    W=wilcox.test(X, unlist(Y), correct=F)$statistic
    mu=m*sum(len.Y)/2 # this is the mean of stat returned by wilcox.test 
    W=W-mu
                
    p.method=NULL # asymptotic, exact, Monte Carlo
    if(is.null(perm)) {
        if (min(m,n)>=50) p.method="asymptotic" 
    } else if (!perm) {
        p.method="asymptotic" 
    }
    if(is.null(p.method)) p.method="Monte Carlo"
    if(verbose) print(p.method)
    
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
            W.mc=.Call("wmw_paired_replicates", X, unlist(Y), len.Y, .corr, 1, as.integer(mc.rep))
        } else {    
            set.seed(1) # so that it can be repeated
            W.mc=sapply (1:mc.rep, function (b) compute.wmw.paired.replicates.stat(X,Y, alternative, correct) )
        }
        W.mc=W.mc-mu        
            
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
compute.wmw.paired.replicates.stat=function(X,Y, alternative, correct) {    
    m=length(X)
    for (i in 1:m) {
        all=c(Y[[i]], X[i])
        sel=sample(1:length(all), 1)
        X[i]=all[sel]
        Y[[i]]=all[-sel]
    }
    stat=wilcox.test(X, unlist(Y), correct=F)$statistic    
}    


#        set.seed(1)
#        W=wilcox.test(infant, unlist(mother), correct=F)$statistic
#        W.B=sapply (1:B, function (b) {
#            # for each infant/mother pair, randomly choose one to be infant
#            for (i in 1:nInfant) {
#                all=c(infant[i], mother[[i]])
#                sel=sample(1:length(all),1)
#                infant[i]=all[sel]
#                mother[[i]]=all[-sel]
#            }
#            wilcox.test(infant, unlist(mother), correct=F)$statistic
#        })    
#        mu=nInfant*sum(ratios)/2 # this is the mean of stat returned by wilcox.test 
#        W=W-mu
#        W.B=W.B-mu
#        p.W=mean(W.B>=abs(W))+mean(W.B<=-abs(W))
