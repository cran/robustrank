# method="large.0"; perm=T; alternative = "two.sided"; correct=T; mc.rep=1e4; trace=1; mode="test"; useC=T; 
pair.wmw.test=function(X,Y
    , alternative = c("two.sided", "less", "greater"), correct = TRUE, perm=NULL, mc.rep=1e4
    , method=c("exact.2","large.0","large","exact","exact.0","exact.1","exact.3"), trace=0
    , mode=c("test","var"), useC=TRUE)
{
    
    alternative <- match.arg(alternative)
    method <- match.arg(method)
    mode <- match.arg(mode)
    
    if(!correct) .corr=0 else .corr=switch(alternative, "two.sided" = 2, "greater" = -1, "less" = 1)     
    
    # remove pairs that have missing values
    # if not converted to double, calls to C function may throw errors for integer X/Y
    stopifnot(length(X)==length(Y))    
    X=as.double(X[!is.na(X) & !is.na(Y)])
    Y=as.double(Y[!is.na(X) & !is.na(Y)])
    m=length(X)
    n=m
    
    # for developing manuscript
    if(mode=="var") return(compute.pair.wmw.Z(X,Y, alternative, correct, method, return.all=TRUE)) 
    
    if(useC) {
        z=.Call("pair_wmw_test", X, Y, .corr, switch(method, "large.0"=-1, "large"=-2, "exact"=-3, "exact.0"=0, "exact.1"=1, "exact.2"=2, "exact.3"=3), as.integer(1), as.integer(1))
    } else {
        #the following R implementation may give a different result under i386 software on 64-bit machine, b/c C is always 64 bit, but R is 32 bit
        z=compute.pair.wmw.Z(X,Y, alternative, correct, method); 
    }        
        
        
    p.method=NULL # asymptotic, exact, Monte Carlo
    if(is.null(perm)) {
        if (min(m,n)>=50) p.method="asymptotic" 
    } else if (!perm) {
        p.method="asymptotic" 
    }
    if(is.null(p.method)) p.method=ifelse(2^m<=mc.rep, "exact", "Monte Carlo")
    if(trace==1) print(p.method)
    
    if (p.method=="asymptotic") {
       # one sided p values are opposite to wilcox.test because the way U is computed here as the rank of Y's
        pval=switch(alternative,
           "two.sided" = 2 * pnorm(abs(z), lower.tail=FALSE),
           "less" = pnorm(z, lower.tail=FALSE), # X<Y
           "greater" = pnorm(z, lower.tail=TRUE) #X>Y
        )
        
    } else {
        if(p.method=="exact"){ 
            # binarys is adapted from Thomas Lumley https://stat.ethz.ch/pipermail/r-help/2000-February/010141.html
            binarys<-function(i,m) {
                a<-2^((m-1):0)
                b<-2*a
                sapply(i,simplify="array", function(x) as.integer((x %% b)>=a))
            }
            indexes=binarys(0:(2^m-1), m)
            if (useC) {
                z.mc=.Call("pair_wmw_test", X, Y, .corr, switch(method, "large.0"=-1, "large"=-2, "exact"=-3, "exact.0"=0, "exact.1"=1, "exact.2"=2, "exact.3"=3), as.integer(0), indexes)
            } else {
                z.mc=apply(indexes, 2, function(ind) {
                    compute.pair.wmw.Z(ifelse(ind==0,X,Y), ifelse(ind==1,X,Y), alternative=alternative, correct=correct, method=method)
                })
            }
            
        } else if (p.method=="Monte Carlo") {
            #### save rng state before set.seed in order to restore before exiting this function
            save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
            if (class(save.seed)=="try-error") { set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
    
            set.seed(1)
            if(useC) {
                z.mc=.Call("pair_wmw_test", X, Y, .corr, switch(method, "large.0"=-1, "large"=-2, "exact"=-3, "exact.0"=0, "exact.1"=1, "exact.2"=2, "exact.3"=3), as.integer(mc.rep), as.integer(1))
            } else {
                z.mc=sapply(1:mc.rep, function(b) {
                    ind=runif(m)
                    compute.pair.wmw.Z(ifelse(ind<0.5,X,Y), ifelse(ind>=0.5,X,Y), alternative=alternative, correct=correct, method=method)
                })            
            }
                
            #### restore rng state 
            assign(".Random.seed", save.seed, .GlobalEnv)     
        }        
        
        # z.mc may be NaN b/c var estimate may be less than 0
        # we treat NaN here as inf b/c a small value of var.est means large z
        z.mc[is.nan(z.mc)] = Inf 
        
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
compute.pair.wmw.Z=function(X,Y, alternative, correct, method, return.all=FALSE) {
      
    m=length(X)
    n=m
    
    # empirical distribution function estimated from observed
    Fx=ecdf(X)
    Fy=ecdf(Y)
    FyX=Fy(X); SyX=1-FyX
    FxY=Fx(Y); SxY=1-FxY
    
    tao=mean(Y>X)
    theta=mean(outer(1:m,1:n, function(i,j) ifelse(i!=j,Y[j]>X[i],NA)), na.rm=TRUE) # note that theta!=FxY here b/c i==j not counted
    theta.sq=theta*theta
    
    r <- rank(c(X, Y))
    U=(sum(r[m+1:n]) - n*(n+1)/2) /m/n    
    correction <- switch(alternative,
       "two.sided" = sign(U-0.5) * 0.5/m/n,
       "greater" = 0.5/m/n,
       "less" = -0.5/m/n
    )
    if(correct) U=U-correction
    #U.lin=m*N/(m+N)*(-mean(Fyyprime(X)) + mean(Fx(YYprime))) # linear expansion for debugging purpose, but this is ridiculous b/c mean(Fyyprime(X)) = 1 - mean(Fx(YYprime)), so this has variance 0
    
    #### variance
    
    # to get more unbiased estimates, as in variance. a caveat is that when the covariance is high, the true sample size may be smaller
    # in addition, it is not quite clear what the proper finite sample adjustment should be here b/c there are three indexes: i, j, k
    # it is not clear how var.large does it
    #finite.adj=m/(m-1) 
    finite.adj=1
    
    # lower order terms
    # tmp variables are needed b/c they are used in both var.exact and var.exact.0
    tmp.1=mean(outer(1:m,1:n, function(i,j) ifelse(i!=j,Y[i]>X[i] & Y[j]>X[i],NA)), na.rm=TRUE) 
    tmp.2=mean(outer(1:m,1:n, function(i,j) ifelse(i!=j,Y[i]>X[i] & Y[i]>X[j],NA)), na.rm=TRUE)
    tmp.3=mean(outer(1:m,1:n, function(i,j) ifelse(i!=j,Y[j]>X[i] & Y[i]>X[j],NA)), na.rm=TRUE)
    v.f=function(theta.1, tao.1) {
        finite.adj.1 = ifelse(theta.1==1/2 & tao.1==1/2, 1, finite.adj) # if under the null, we are not estimating some terms, so there is no need to adjust
        A  =m*        tao.1*(1-tao.1) *finite.adj.1
        B  =m*(m-1)*  (tmp.1 - tao*theta + tmp.2 - tao*theta) *finite.adj
        C.1=m*(m-1)*  {theta.1*(1-theta.1) *finite.adj.1  + (tmp.3 - theta*theta)*finite.adj}
        A + 2*B + C.1
    }
        
    # VarF_{X}(Y_{j}) and VarF_{Y}(X_{i})
    tmp.5=outer.3.mean(1:m,1:m,1:m, function(i,j,k) ifelse(i==j|i==k|j==k,NA,(Y[j]>X[i])*(Y[k]>X[i])) )
    tmp.6=outer.3.mean(1:m,1:m,1:m, function(i,j,k) ifelse(i==j|i==k|j==k,NA,(Y[j]>X[i])*(Y[j]>X[k])) )
    C.2.var=    m*(m-1)*(m-2)* (tmp.5 - theta.sq + tmp.6 - theta.sq) *finite.adj 
    C.2.var.0=  m*(m-1)*(m-2)* (1/12 + 1/12)
    
    # -Cov{F_{Y}(X_{i}),F_{X}(Y_{i})} * 2
    tmp.4=outer.3.mean(1:m,1:m,1:m, function(i,j,k) ifelse(i==j|i==k|j==k,NA,(Y[j]>X[i])*(Y[i]>X[k])) )    
    C.2.cov  =  m*(m-1)*(m-2)* (tmp.4 - theta.sq) * 2 *finite.adj 
    C.2.cov.0=  m*(m-1)*(m-2)* (tmp.4 - 1/4     ) * 2 *finite.adj # reduces bias by a little, but inflates variability by a great deal
    
    #var.exact  =(C.2.var   + C.2.cov + v.f(1/2,1/2))/m^3 # difference between this and next line is fairly small
    var.exact  =(C.2.var   + C.2.cov + v.f(theta,tao))/m^3 # only difference from var.exact.0 is in the first term, like the relationship between var.large and var.large.0
    # the following only works under null
    var.exact.0=(C.2.var.0 + C.2.cov.0+v.f(theta,tao))/m^3 
    var.exact.1=(C.2.var.0 + C.2.cov + v.f(theta,tao))/m^3 # use est theta in C.2.cov leads to a good bias-variance tradeoff
#    var.exact.x=(C.2.var.0 + C.2.cov)/m^3 # no lower order term, does not do as well, dropped from laster MC studies
    # correct for bias in theta^2 est
    # one-step
    C.2.cov.1=  m*(m-1)*(m-2)* (tmp.4 - (theta.sq-var.exact.1/m)) * 2 
    var.exact.2=(C.2.var.0 + C.2.cov.1+v.f(theta,tao))/m^3 # the best so far
#    # two-step, fairly small improvement over one-step, dropped from further studies
#    C.2.cov.2=  m*(m-1)*(m-2)* (tmp.4 - (theta.sq-var.exact.2/m)) * 2 
#    var.exact.x=(C.2.var.0 + C.2.cov.2+v.f(theta,tao))/m^3 
    # also correct for bias in lower order terms, not much of an improvement either
    var.exact.3=(C.2.var.0 + C.2.cov.1+v.f(1/2,1/2) + m*(m-1)*var.exact.1/m )/m^3
    
    # note that it is not just the higher order terms, it also includes some lower order terms
    var.large=  var(Fy(X)) + var(Fx(Y)) - 2*cov(Fy(X), Fx(Y))
    var.large.0=1/12       + 1/12       - 2*cov(Fy(X), Fx(Y))    
    
#    # part of var.exact, performance not so good
#    var.large=  (C.2.var   + C.2.cov)/m^3 
#    var.large.0=(C.2.var.0 + C.2.cov)/m^3 # change C.2.cov.0 to C.2.cov reduces variability of the estimated covariance, but intro bias, so we use C.2.cov instead
    
    z=sqrt(m)*(U-0.5)/sqrt(get(paste("var.",method,sep="")))
    
    if(!return.all) return(z) else return(c(U=U, var.large=var.large, var.large.0=var.large.0, 
        var.exact=var.exact, var.exact.0=var.exact.0, var.exact.1=var.exact.1, var.exact.2=var.exact.2, var.exact.3=var.exact.3)) 
}

outer.3.mean=function(x,y,z,func) {
    allidx = expand.grid(x,y,z)
    mean(func(allidx[,1],allidx[,2],allidx[,3]),na.rm=TRUE)
}
