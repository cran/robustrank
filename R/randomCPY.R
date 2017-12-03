#generating rejective sample
rejective.sampling = function (N, n, pik) {
    s = sample(N, n, replace = TRUE, prob = pik)
    while (length(unique(s)) != n) {
        s = sample(N, n, replace = TRUE, prob = pik)
    } 
    s   
}

# generate random number from (generalized) Bernoulli dist
# if generalized, output is 1,2,3,...; otherwise, output is 1 or 0
# n is number of random numbers to generate, prob can be one number or a vector
rbern=function(n, prob, generalized=FALSE) {
    if (!generalized){
        # binary outcome
        if (length(prob)==1) {
            rbinom(n,1,prob)
        } else {
            prob=cbind(1:n,prob)[,2]
            rbinom(n,1,prob)
#            sapply(prob, function(p) rbinom(1,1,p))
        }
    } else {
        x=rmultinom(n, 1, prob)
        out = apply(x, 2, function(y) which(y==1))
        out
    }
}

dbern=function(x,prob,log=FALSE){
    out=ifelse(x==1, prob, 1-prob)
    ifelse(log, log(out), out)
}

# correlated Bernoulli
dcorbern=function(x, p, a, log=FALSE){
    out=dbern(x[1],p)
    if (length(x)>1) {
        for (i in 2:length(x)) {
            out = out * dbern(x[i], p+a*(x[i-1]-p))
        }
    }
    ifelse(log, log(out), out)
}


# the integratd density of a normal random variable, whose mean and precision follow a normal gamma distribution. It is a three-parameter t distribution. 
# keywords: normgamma
# when x has length greater than 1 and same.distr is true, the data are considered to be from the same mean and variance
dnorm.norm.gamma = function(x, p, same.distr=FALSE, log=FALSE) {
    mu.0=p[1]; lambda=p[2]; a=p[3]; beta=p[4]
    n=length(x)
    if (!same.distr) {
        if (n>1) {
            mu.0=rep(mu.0, n)
            lambda=rep(lambda, n)
            a=rep(a, n)
            beta=rep(beta, n)
        }
        ll= 1/2*log(lambda/(2*pi*(1+lambda))) + log(gamma(a+1/2)/gamma(a)) + a*log(beta) - (a+1/2)*log(beta+lambda*(x-mu.0)**2/(2*(1+lambda)))
    } else{
        s2=var(x)
        ll= 1/2*log(lambda/((2*pi)**n*(n+lambda))) + log(gamma(a+n/2)/gamma(a)) + a*log(beta) - (a+n/2)*log(beta + n*s2/2 +n*lambda*(mean(x)-mu.0)**2/(2*(1+lambda)))
    }
    ll=unname(ll)
    ifelse(log, ll, exp(ll))
}


# simulation samples from a normal random variable, whose mean and precision follow a normal gamma distribution
rnorm.norm.gamma = function(n, mu.0, lambda, alpha, beta) {
    tao=rgamma(n, shape=alpha, rate=beta)
    mu=rnorm(n, mean=mu.0, sd=sqrt(1/(lambda*tao)))
    rnorm(n, mean=mu, sd=tao**-.5)
}

# simulate correlated normal random variables, correlation is sqrt(alpha)
rnorm.cor = function (n, mu, sd, alpha) {
    out=numeric(n)
    out[1]=rnorm(1, mu, sd)
    if (n>1) 
        for (i in 2:n) {
            out[i]=sqrt(alpha)*out[i-1] + sqrt(1-alpha)*rnorm(1, mu, sd)
        }
    out
}

# mixture normal density funtion 
# mix.p: proportion of first component
dmixnorm=function (x, mix.p, sd1, sd2, log=FALSE){
    out=mix.p*dnorm (x,0,sd1) + (1-mix.p)*dnorm (x,0,sd2)
    if(log) out=log(out)
    out
}

#mix.p=.5; mu1=0; mu2=2; sd1=1; sd2=1
rmixnorm=function (n, mix.p, mu1, mu2, sd1, sd2){
    r=rbern(n, mix.p)
    x.1=rnorm(n, mu1, sd1)
    x.2=rnorm(n, mu2, sd2)
    ifelse(r, x.1, x.2)
}


#n=1000; loc.1<-loc.2<-0; scale.1<-scale.2<-1
rbilogistic=function(n, loc.1, loc.2, scale.1, scale.2, rho) {
    if(rho==0.5) {
        dat=VGAM::rbilogis(n=n, loc1 = loc.1, scale1 = scale.1, loc2 = loc.2, scale2 = scale.2)
    } else if (rho<0.478 & rho> -0.271) {
        if(rho==0) {
            theta=0
        } else{
            rho.f=function(t) (12*(1+t)*dilog(1-t) -24*(1-t)*log(1-t) )/t^2 - 3*(t+12)/t - rho
            theta = uniroot(rho.f, interval=c(-1+1e-3,1-1e-3))$root
        }
        suppressMessages({ myCop <- copula::amhCopula(param=theta) }) # suppress msg: parameter at boundary ==> returning indepCopula()
        myMvd <- copula::mvdc(copula=myCop, margins=c("logis", "logis"), paramMargins=list(list(location=loc.1, scale=scale.1), list(location=loc.2, scale=scale.2)))
        dat <- (copula::rMvdc(n, myMvd))        
    } else {
        # use normal copula outside those windows of rho
        myCop <- copula::normalCopula(param=rho)
        myMvd <- copula::mvdc(copula=myCop, margins=c("logis", "logis"), paramMargins=list(list(location=loc.1, scale=scale.1), list(location=loc.2, scale=scale.2)))
        dat <- (copula::rMvdc(n, myMvd))        
    }
}

rbigamma=function(n, shape.1, shape.2, rate.1, rate.2, rho) {
    myCop <- copula::normalCopula(param=rho)
    myMvd <- copula::mvdc(copula=myCop, margins=c("gamma", "gamma"), paramMargins=list(list(shape=shape.1, rate=rate.1), list(shape=shape.2, rate=rate.2)))
    dat <- (copula::rMvdc(n, myMvd))        
}

dilog <- function(x) {
    flog  <- function(t) log(t) / (1-t)  # singularity at t=1, almost at t=0
    pracma::simpadpt(flog, 1, x, tol = 1e-12)
}

# log.p and lower.tail are not implemented to reduce complexity
rdoublexp=function(n, location=0, scale=1) rexp(n,rate=1/scale)*(rbern(n,1/2)*2-1)+location
ddoublexp=function(x, location=0, scale=1) dexp(abs(x-location),rate=1/scale)/2
qdoublexp=function(p, location=0, scale=1) ifelse(p>0.5,1,-1)*qexp(1-(1-ifelse(p>0.5,p,1-p))*2,rate=1/scale) + location
pdoublexp=function(q, location=0, scale=1) ifelse(q>location, 1-pexp(abs(q-location),rate=1/scale,lower.tail=FALSE)/2, pexp(abs(location-q),rate=1/scale,lower.tail=FALSE)/2)
## check against rmulti::dlaplace etc
#pdoublexp(c(-1,0,1,2,3), location=-1, scale=2);          plaplace(c(-1,0,1,2,3), m=-1, 2)
#ddoublexp(c(-1,0,1,2,3), location=-1, scale=2);          dlaplace(c(-1,0,1,2,3), m=-1, 2)
#qdoublexp(seq(0,1,by=0.1), location=-1, scale=2);      qlaplace(seq(0,1,by=0.1), m=-1, 2)
#plot(density(rdoublexp(1e3, location=-1, scale=2))); lines(density(rlaplace(1e3, m=-1, 2)), col=2)


rbidoublexp=function(n, loc.1, loc.2, scale.1, scale.2, rho) {
    
    # this seems to be necessary; otherwise linux installation complains about no def found for these functions
    rdoublexp=function(n, location=0, scale=1) rexp(n,rate=1/scale)*(rbern(n,1/2)*2-1)+location
    ddoublexp=function(x, location=0, scale=1) dexp(abs(x-location),rate=1/scale)/2
    qdoublexp=function(p, location=0, scale=1) ifelse(p>0.5,1,-1)*qexp(1-(1-ifelse(p>0.5,p,1-p))*2,rate=1/scale) + location
    pdoublexp=function(q, location=0, scale=1) ifelse(q>location, 1-pexp(abs(q-location),rate=1/scale,lower.tail=FALSE)/2, pexp(abs(location-q),rate=1/scale,lower.tail=FALSE)/2)
    
    myCop <- copula::normalCopula(param=rho)
    myMvd <- copula::mvdc(copula=myCop, margins=c("doublexp", "doublexp"), paramMargins=list(list(location=loc.1, scale=scale.1), list(location=loc.2, scale=scale.2)))
    dat <- (copula::rMvdc(n, myMvd))        
}
#x=rbidoublexp(100,0,1,1,1,.5)
#plot(x[,1], x[,2])
