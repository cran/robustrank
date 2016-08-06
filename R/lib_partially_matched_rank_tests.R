# useC only implements mw.mw.3 for now
pm.wilcox.test = function(Xpaired, Ypaired, Xextra, Yextra
#    , alternative = c("two.sided", "less", "greater"), correct = TRUE
    , method=c("sr.mw.20", "mw.mw.20", "mw.mw.00", "all")
    , mode=c("test","var"), trace=0)
{
    
    method=match.arg(method)
    mode <- match.arg(mode)
    #alternative <- match.arg(alternative)
    alternative="two.sided"
    
    correct=TRUE
    if(!correct) .corr=0 else .corr=switch(alternative, "two.sided" = 2, "greater" = -1, "less" = 1)     
    
    stopifnot(length(Xpaired)==length(Ypaired))    
    m=length(Xpaired)
    n=length(Yextra)
    lx=length(Xextra) # so that it does not get mixed up with number 1
    N=m+n
    alpha=n/N
    beta=N/(m+N)
    
    # for developing manuscript
    if(mode=="var") return(compute.pm.wilcox.z(Xpaired,Ypaired,Xextra,Yextra, alternative, correct, method, return.all=TRUE)) 
    
    useC=FALSE
    mc.rep=1
    
#    if(useC) {
#        z=.Call("pm_wmw_test", Xpaired, Ypaired, Yextra, .corr, switch(method, "mw.mw.3"=3), as.integer(1))
#    } else {
#        #the following R implementation may give a different result under i386 software on 64-bit machine, b/c C is always 64 bit, but R is 32 bit
        z=compute.pm.wilcox.z(Xpaired,Ypaired,Xextra,Yextra, alternative, correct, method)
#    }        
    
    # no permutation/Monte Carlo method yet
    p.method="asymptotic" 
    
    if (p.method=="asymptotic") {
       # one sided p values are opposite to wilcox.test because the way U is computed here as the rank of Ypaired's
        pval=switch(alternative,
           "two.sided" = z, #2 * pnorm(abs(z), lower.tail=FALSE),
           "less" = pnorm(z, lower.tail=FALSE), # Xpaired<Ypaired
           "greater" = pnorm(z, lower.tail=TRUE) #Xpaired>Ypaired
        )
        
    } else {
        if (p.method=="Monte Carlo") {
            #### save rng state before set.seed in order to restore before exiting this function
            save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
            if (class(save.seed)=="try-error") { set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
    
            set.seed(1)
            if(useC) {
                z.mc=.Call("pm_wmw_test", Xpaired, Ypaired, .corr, switch(method, "mw.mw.3"=3), as.integer(mc.rep))
            } else {
                z.mc=sapply(1:mc.rep, function(b) {
                    ind=runif(m)
                    compute.pair.wmw.Z(ifelse(ind<0.5,Xpaired,Ypaired), ifelse(ind>=0.5,Xpaired,Ypaired), alternative=alternative, correct=correct, method=method)
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
           "less" = mean(z.mc>=z*(1-sign(z)*numerical.stability.margin)), # Xpaired<Ypaired
           "greater" = mean(z.mc<=z*(1+sign(z)*numerical.stability.margin)) #Xpaired>Ypaired
        )
    }
    
    pval
}        

compute.pm.wilcox.z=function(Xpaired,Ypaired,Xextra,Yextra, alternative, correct, method, return.all=FALSE) {
    m=length(Xpaired)
    n=length(Yextra)
    lx=length(Xextra)
    N=m+n
    alpha=n/N
    beta=N/(m+N)
    
    YYprime=c(Ypaired,Yextra)
    
    delta=Ypaired-Xpaired    
    abs.delta=abs(delta)
    sgn.delta=sign(delta)
    pos.delta=ifelse(delta>0,1,0)
    
    # empirical distribution function estimated from observed
    Fx=ecdf(Xpaired)
    Fy=ecdf(Ypaired)
    Fyprime=ecdf(Yextra)
    Fyyprime=ecdf(YYprime)
    Fplus=ecdf(abs.delta)
    Fxyyprime=ecdf(c(Xpaired,Ypaired,Yextra))
    cov.G.F=cov(Fy(Xpaired), Fx(Ypaired)) 
    cov.G.F.1=cov(Fyyprime(Xpaired), Fx(Ypaired)) 
    
    # compute h.1 for W.plus
    # Vectorize the following does not help speed b/c the loop is still implemented in R
    h.1.f=function(i) (sum(delta+delta[i]>0) - sum(delta[i]>0))/(length(delta)-1)
    h.1=numeric(m)
    for (i in 1:m) h.1[i]=h.1.f(i)
#    #using outer is actually slower probably b/c memory operations
#    delta.pair.plus=outer(delta, delta, FUN="+")
#    diag(delta.pair.plus)=NA # we don't want to count delta_i + delta_j
#    h.1.0=apply(delta.pair.plus>0, 1, mean, na.rm=TRUE)
#    all(h.1==h.1.0)
    
    # Wilcoxon signed rank stat and var
    r.1 <- rank(abs.delta)    
    W.plus = mean(pos.delta*r.1)
    W.plus.2=mean(sgn.delta*r.1)
    var.W.plus = m * var(h.1) # var(W.plus.2) is 4 x var(W.plus)
    var.W.plus.0 = (m+1)*(2*m+1)/(24*m)
    
    # Wilcoxon-Mann-Whitney stat and var
    r.2 <- rank(c(Xpaired, Yextra))
    W.mw=sum(r.2[(m+1):N])/(N+1)
    var.W.mw=m * (alpha^2*var(Fyprime(Xpaired)) + alpha*(1-alpha)*var(Fx(Yextra)))
    var.W.mw.0=m*n/(12*(N+1))
        
    # covariance between W.plus and W.mw    
    cov. = -m*alpha * cov(h.1, Fy(Xpaired))
    rho=cov./sqrt(var.W.plus*var.W.mw)
    rho.0 = - 6*sqrt(alpha) * cov(Fplus(abs.delta) * sgn.delta, Fx(Xpaired))
    cov.0 =  rho.0  * sqrt(var.W.plus.0*var.W.mw.0) # good
    # we may also put in the population mean in the estimate, it reduces bias but increases var greatly, not good
    rho.0a= - 6*sqrt(alpha) * mean(Fplus(abs.delta) * sgn.delta * Fx(Xpaired)) 
    cov.0a = rho.0a * sqrt(var.W.plus.0*var.W.mw.0) 
    # other choices tried that do not make big differences
    # cov.0d = (2*rho-rho.0) * sqrt(var.W.plus.0*var.W.mw.0) # 
    #rho.0b = - 6*sqrt(alpha) * cov(Fplus(abs.delta) * sgn.delta, Fyprime(Xpaired))
    #cov.0  = rho    * sqrt(var.W.plus.0*var.W.mw.0) # 
    
    # MW-MW1 and MW-MW2
    U.p=(sum(rank(c(Xpaired, Ypaired))[m+1:m]) - m*(m+1)/2)/(m*m)
    U.mw=(W.mw*(N+1)-n*(n+1)/2)/(m*n)
    var.U.p.0=1/6-2*cov.G.F
    var.U.mw.0=1/(12*alpha)
    cov.U.mw.p.0=1/12-cov.G.F
    var.U.p =  var(Fy(Xpaired)) + var(Fx(Ypaired)) - 2*cov.G.F
    var.U.mw=  var(Fyprime(Xpaired)) + (1/alpha-1)*var(Fx(Yextra))
    cov.U.mw.p=var(Fy(Xpaired)) - cov.G.F
    
    # MW-MW0, another way to get this is to combine MW-MW1 and MW-MW2
    T = sum(rank(c(Xpaired, Ypaired, Yextra))[(m+1):(m+N)])/(m+N+1)
    var.T.0=(m*N/(m+N))^2 * (1/12/m + 1/12/N - 2*cov.G.F.1/N) 
    var.T=  (m*N/(m+N))^2 * (var(Fyyprime(Xpaired))/m + var(Fx(YYprime))/N - 2*cov.G.F.1/N)
    #var.T=(m*N/(m+N))^2 * (var(Fyyprime(Xpaired))/m + var(Fx(YYprime))/N - 2*cov(Fyyprime(Xpaired), Fx(Ypaired))/N)
    #T.lin=m*N/(m+N)*(-mean(Fyyprime(Xpaired)) + mean(Fx(YYprime))) # linear expansion, note that mean(Fyyprime(Xpaired)) and mean(Fx(YYprime)) are perfectly correlated in data
    
    tests=NULL
    
    ####################################################
    # SR-MW tests, several versions of rho are used
    
    vec=c(W.plus-(m+1)/4, W.mw-n/2)
    
    stat.1 = try(c(vec %*% solve(matrix(c(var.W.plus.0, cov.0, cov.0, var.W.mw.0),2) , vec)), silent=TRUE)  # solve may fail
    if (inherits(stat.1, "try-error")) stat.1<-pval.1<-NA else pval.1=pchisq(stat.1, df=2, lower.tail=FALSE)
    tests=rbind(tests, "sr.mw.10"=c(stat.1, pval.1)) #SR-MW$_{1}$
    
#    stat.1 = try(c(vec %*% solve(matrix(c(var.W.plus.0, cov.0a, cov.0a, var.W.mw.0),2) , vec)), silent=TRUE)  # solve may fail
#    if (inherits(stat.1, "try-error")) stat.1<-pval.1<-NA else pval.1=pchisq(stat.1, df=2, lower.tail=FALSE)
#    tests=rbind(tests, "SR-MW.1a"=c(stat.1, pval.1))
    
    stat.1 = try(c(vec %*% solve(matrix(c(var.W.plus, cov., cov., var.W.mw),2) , vec)), silent=TRUE)  # solve may fail
    if (inherits(stat.1, "try-error")) stat.1<-pval.1<-NA else pval.1=pchisq(stat.1, df=2, lower.tail=FALSE)
    tests=rbind(tests, "sr.mw.11"=c(stat.1, pval.1))
    
    comb=c(1, var.W.plus.0/var.W.mw.0)
    stat.1 = (comb %*% vec)^2 / (comb %*% matrix(c(var.W.plus.0, cov.0, cov.0, var.W.mw.0),2) %*% comb)      
    tests=rbind(tests, "sr.mw.20"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
    
    comb=c(1, var.W.plus/var.W.mw)
    stat.1 = (comb %*% vec)^2 / (comb %*% matrix(c(var.W.plus, cov., cov., var.W.mw),2) %*% comb)      
    tests=rbind(tests, "sr.mw.21"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
    
    
    ####################################################
    # MW-MW tests
        
#    stat.1 = (T-N/2)^2/var.T.0
#    tests=rbind(tests, "mw.mw.00"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
#        
#    stat.1 = (T-N/2)^2/var.T
#    tests=rbind(tests, "mw.mw.01"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
        
    vec=sqrt(m)*c(U.p-1/2, U.mw-1/2)    
    
    stat.1 = try(c(vec %*% solve(matrix(c(var.U.p.0, cov.U.mw.p.0, cov.U.mw.p.0, var.U.mw.0),2) , vec)), silent=TRUE)  # solve may fail
    if (inherits(stat.1, "try-error")) stat.1<-pval.1<-NA else pval.1=pchisq(stat.1, df=2, lower.tail=FALSE)
    tests=rbind(tests, "mw.mw.10"=c(stat.1, pval.1))
    
    stat.1 = try(c(vec %*% solve(matrix(c(var.U.p, cov.U.mw.p, cov.U.mw.p, var.U.mw),2) , vec)), silent=TRUE)  # solve may fail
    if (inherits(stat.1, "try-error")) stat.1<-pval.1<-NA else pval.1=pchisq(stat.1, df=2, lower.tail=FALSE)
    tests=rbind(tests, "mw.mw.11"=c(stat.1, pval.1))
    
    comb=c(1, var.U.p.0/var.U.mw.0)
    stat.1 = (comb %*% vec)^2 / (comb %*% matrix(c(var.U.p.0, cov.U.mw.p.0, cov.U.mw.p.0, var.U.mw.0),2) %*% comb)      
    tests=rbind(tests, "mw.mw.20"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
    
    comb=c(1, var.U.p/var.U.mw)
    stat.1 = (comb %*% vec)^2 / (comb %*% matrix(c(var.U.p, cov.U.mw.p, cov.U.mw.p, var.U.mw),2) %*% comb)      
    tests=rbind(tests, "mw.mw.21"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
    
    comb=c(1, n/m)
    stat.1 = (comb %*% vec)^2 / (comb %*% matrix(c(var.U.p.0, cov.U.mw.p.0, cov.U.mw.p.0, var.U.mw.0),2) %*% comb)      
    tests=rbind(tests, "mw.mw.00"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
    
    comb=c(1, n/m)
    stat.1 = (comb %*% vec)^2 / (comb %*% matrix(c(var.U.p, cov.U.mw.p, cov.U.mw.p, var.U.mw),2) %*% comb)      
    tests=rbind(tests, "mw.mw.01"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
    
    # using both Yextra and Xextra
    if (!is.null(Xextra)) {
        U.3=(sum(rank(c(Xextra, Ypaired))[lx+1:m]) - m*(m+1)/2)/(lx*m)
        U.4=(sum(rank(c(Xextra, Yextra))[lx+1:n]) - n*(n+1)/2)/(lx*n)
    
        vec=sqrt(m)*c(U.p-1/2, U.mw-1/2, U.3-1/2, U.4-1/2)    
        v.tmp=diag(c(var.U.p.0, var.U.mw.0, ((m/lx+1)/12), ((m/lx+m/n)/12)))
        v.tmp[2,1]<-v.tmp[1,2]<- cov.U.mw.p.0
        v.tmp[3,1]<-v.tmp[1,3]<- cov.U.mw.p.0
        v.tmp[4,1]<-v.tmp[1,4]<- 0
        v.tmp[3,2]<-v.tmp[2,3]<- -cov.G.F
        v.tmp[4,2]<-v.tmp[2,4]<- m/(12*n)
        v.tmp[4,3]<-v.tmp[3,4]<- m/(12*lx)
    
        comb=c(1, 0, 0, var.U.p.0/((m/lx+m/n)/12))
        stat.1 = (comb %*% vec)^2 / (comb %*% v.tmp %*% comb)      
        tests=rbind(tests, "mw.mw.20a"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
    
        comb=c(1, var.U.p.0/var.U.mw.0, var.U.p.0/((m/lx+1)/12), var.U.p.0/((m/lx+m/n)/12))
        stat.1 = (comb %*% vec)^2 / (comb %*% v.tmp %*% comb)      
        tests=rbind(tests, "mw.mw.20b"=c(stat.1, pchisq(stat.1, df=1, lower.tail=FALSE)))
    }    
    
    if(!return.all) 
        return(tests)
    else 
        return(c(W.plus=W.plus, W.mw=W.mw, var.W.plus=var.W.plus, var.W.mw=var.W.mw, cov.=cov., rho=rho, rho.0=rho.0, rho.0a=rho.0a, T=T, var.T=var.T, 
                 U.p=sqrt(m)*U.p, U.mw=sqrt(m)*U.mw, var.U.mw.0=var.U.mw.0, var.U.p.0=var.U.p.0, cov.U.mw.p.0=cov.U.mw.p.0 ))
    
}


# naively or simply combining sign test with MW, treat the two as independent
part.naive=function(Xpaired,Ypaired,Xprime){
    stopifnot(length(Xpaired)==length(Ypaired))
    test.1=wilcox.test(Xpaired, Ypaired, paired=T)# wilcox test on the paired data
    test.2=wilcox.test(Xprime, Ypaired) # unpaired test
    test.1
    test.2
    partial.p.1 = pchisq(qnorm(test.1$p.val)**2 + qnorm(test.2$p.val)**2, df=2, lower.tail=FALSE)                
}

sign.mw.test = function(Xpaired, Ypaired, Z) {    
    
    stopifnot(length(Xpaired)==length(Ypaired))
    m=length(Xpaired)
    n=length(Z)
    N=m+n
    
    # a signed test for the paired data
    t.1=mean(Xpaired>Ypaired)    
    if (t.1==0) t.1=t.1+.5/m else if (t.1==1) t.1=t.1-.5/m # continuity correction
    var.t.1=t.1 * (1-t.1) / m
    var.t.1
    
    r <- rank(c(Ypaired, Z))
    t.2=(sum(r[(n+1):N]) - m*(m+1)/2)/n/m # U or AUC
    if (t.2==0) t.2=t.2+.5/n/m else if (t.2==1) t.2=t.2-.5/n/m
    Fy=ecdf(Ypaired)
    Fz=ecdf(Z)
    var.t.2=var(Fz(Ypaired))/n + var(Fy(Z))/m
    
    # covariance
    cov.1.2=1/m* {mean (outer(1:n, 1:m, function (i,j) as.integer(Z[i]>Ypaired[j]) * as.integer(Xpaired[j]>Ypaired[j]) )) - mean (outer(1:n, 1:m, function (i,j) as.integer(Z[i]>Ypaired[j]) )) * t.1}
    
    # var matrix
    Sigma=matrix(c(var.t.1,cov.1.2,cov.1.2,var.t.2),2,2)
    # p-value
    pval=try(pchisq(c(c(t.1-0.5, t.2-0.5) %*% solve(Sigma, c(t.1-0.5, t.2-0.5))), df=2, lower.tail=FALSE))
    ifelse (inherits(pval, "try-error"), NA, pval)
    
}

mw.mw.test=function(Xpaired,Ypaired,Yextra) {
    
    stopifnot(length(Ypaired)==length(Xpaired))
    m=length(Xpaired)
    n=length(Yextra)
    N=m+n
    
    t.1=mean(Ypaired>Xpaired)    
    if (t.1==0) t.11=t.1+.5/m else if (t.1==1) t.11=t.1-.5/m else t.11=t.1  # continuity correction
    var.t.1=t.11 * (1-t.11)/m
    var.t.1
    
    r <- rank(c(Xpaired, Yextra))
    t.2=(sum(r[(m+1):N]) - n*(n+1)/2)/m/n # U or AUC
    if (t.2==0) t.2=t.2+.5/m/n else if (t.2==1) t.2=t.2-.5/m/n
    Fx=ecdf(Xpaired)
    Fyprime=ecdf(Yextra)
    var.t.2=var(Fyprime(Xpaired))/m + var(Fx(Yextra))/n # same result as pROC::var.auc.delong
    
    mean.x.gr.y = mean( outer(1:m,1:m, function(i,j) ifelse(i==j,NA,Ypaired[i]>Xpaired[j])), na.rm=TRUE )
    mean.z.gr.y=mean( outer(1:m,1:n, function(j,k) Yextra[k]>Xpaired[j]) )
    
    #tmp = mean( multi.outer(function(i,j,k) ifelse(i==j,NA,(Ypaired[i]>Xpaired[j])*(Yextra[k]>Xpaired[j])), 1:m,1:m,1:n), na.rm=TRUE )
    tmp   =outer.3.mean(1:m,1:m,1:n, function(i,j,k) ifelse(i==j,NA,(Ypaired[i]>Xpaired[j])*(Yextra[k]>Xpaired[j])) )
    tmp.1 =outer.3.mean(1:m,1:m,1:m, function(i,j,k) ifelse(i==j|i==k|j==k,NA,(Ypaired[i]>Xpaired[j])*(Ypaired[i]>Xpaired[k])) )
    tmp.11=outer.3.mean(1:m,1:m,1:m, function(i,j,k) ifelse(i==j|i==k|j==k,NA,(Ypaired[i]>Xpaired[j])*(Ypaired[k]>Xpaired[i])) )
    tmp.2= outer.3.mean(1:m,1:m,1:m, function(i,j,k) ifelse(i==j|i==k|j==k,NA,(Ypaired[i]>Xpaired[j])*(Ypaired[k]>Xpaired[j])) )
    tmp.21=tmp.11 # it can be shown that tmp.21 equals tmp.11. it helps to draw pictures
    #tmp.21=outer.3.mean(1:m,1:m,1:m, function(i,j,k) ifelse(i==j|i==k|j==k,NA,(Ypaired[i]>Xpaired[j])*(Ypaired[j]>Xpaired[k])) ) 
    
    # variance
    Sigma=(
        #second term
          m^2*n^2*var.t.2
        #third term
        + 2*m*n* {mean( outer(1:m, 1:n, function(j,k) (Ypaired[j]>Xpaired[j])*(Yextra[k]>Xpaired[j])) ) - t.1 * mean.z.gr.y}    
        + 2*m*n*(m-1)* (tmp - mean.x.gr.y * mean.z.gr.y)
        #first term
        + m* (m*var.t.1) # m*var.t.1 is the variance of (Ypaired[i]>Xpaired[i])
        + m*(m-1)* (mean.x.gr.y*(1-mean.x.gr.y))
        + m*(m-1)*(m-2)* (tmp.1 - mean.x.gr.y^2)
        + m*(m-1)*(m-2)* (tmp.11- mean.x.gr.y^2)
        + m*(m-1)*(m-2)* (tmp.2 - mean.x.gr.y^2)
        + m*(m-1)*(m-2)* (tmp.21- mean.x.gr.y^2)
        + 2*m*(m-1)* {mean( outer(1:m,1:m, function(i,j) ifelse(i==j,NA,(Ypaired[i]>Xpaired[i])*(Ypaired[i]>Xpaired[j]))), na.rm=TRUE ) - t.1 * mean.x.gr.y}  
        + 2*m*(m-1)* {mean( outer(1:m,1:m, function(i,j) ifelse(i==j,NA,(Ypaired[i]>Xpaired[i])*(Ypaired[j]>Xpaired[i]))), na.rm=TRUE ) - t.1 * mean.x.gr.y}  
    )/(m*(m+n))^2
    
    # test stat
    test.stat = unname(wilcox.test(c(Ypaired,Yextra),Xpaired)$statistic/m/(m+n) - 0.5)
    # p value
    pchisq( test.stat^2/Sigma, df=1, lower.tail=FALSE)
    
}


#    ####################################################
#    # maximum weighted combination
#    a=seq(1,5,by=.5)
#    ww=cbind(1, c(rev(1/a),a[-1]))
#    # normalize ww
#    ww=ww/sqrt(ww[,1]**2+ww[,2]**2)
#    stat.2=max(abs(ww%*%vec))
#    
#    # save rng state before set.seed in order to restore before exiting this function
#    save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
#    if (class(save.seed)=="try-error") { set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }                        
#    set.seed(1)
#    
#    sam=mvrnorm (mc.n, rep(0,2), matrix(c(1,rho.0,rho.0,1),2)) # mvtnorm::rmvnorm can fail 
#    sam.max=apply(sam %*% t(ww), 1, function(x) max(abs(x)))
#    pval.2 = mean(sam.max>stat.2)                    
#    tests=rbind(tests, "max weighted"=c(stat.2, pval.2))

#    # restore rng state 
#    assign(".Random.seed", save.seed, .GlobalEnv)     
#        } else if (method=="perm"){
#            # save rng state before set.seed in order to restore before exiting this function
#            save.seed <- try(get(".Random.seed", .GlobalEnv), silent=TRUE) 
#            if (class(save.seed)=="try-error") {set.seed(1); save.seed <- get(".Random.seed", .GlobalEnv) }    
#            set.seed(1)    
#            
#            stats=numeric(n.perm)
#            x.yprime=c(Xpaired,Yextra)
#            for (b in 1:n.perm) {            
#                s=sample(N,m,replace=FALSE)
#                # suppose m=n=4 and s=c(3,4,5,6)
#                #
#                # Yprime.b=x.yprime[c(1,2,7,8)]
#                # X.b=x.yprime[c(3,4,5,6)]
#                # Y.b=c(Ypaired[3,4,1,2])
#                # s.1=sample(c(3,4),1)
#                # tmp=X.b[s.1]; X.b[s.1]=Y.b[s.1]; Y.b[s.1]=tmp
#                
#                # more generally
#                #
#                Yprime.b=x.yprime[setdiff(1:N,s)]
#                X.b=x.yprime[s]
#                s.2=intersect(1:m, s) # 3,4
#                s.1=sample(s.2,length(s.2)/2)
#                Y.b=c(Ypaired[c(s.1, setdiff(1:m,s.1))])
#                tmp=X.b[s.1]; X.b[s.1]=Y.b[s.1]; Y.b[s.1]=tmp
#                
#                stats[b]=pm.wilcox.statistic(Y.b-X.b, x.b, Yprime.b)[1]              
#            }
#            #myprint(stat)
#            #print(summary(stats))        
#            pval=mean(stats>stat, na.rm=T)
#            
#            # restore rng state 
#            assign(".Random.seed", save.seed, .GlobalEnv)             
#        }
