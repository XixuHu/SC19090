#' @title Histogram of mixture of normal distribution 
#' @description Generate a random sample of size 1000 from a normal location mixture. 
#' $$X_1\sim N(0,1),X_2\sim N(3,1),Z=p_1X_1+(1-p_1)X_2$$
#' Graph the histogram of the sample with density superimposed
#' @param p1 The proportion of mixture
#'
#' @return A histogram of the sample with density superimposed
#' @export
#'
#' @examples
#' \dontrun{
#'  gen_Mixnorm(0.25)
#' }
gen_Mixnorm = function(p1){
 n = 1e4
 X1 = rnorm(n,0,1)
 X2 = rnorm(n,3,1)
 r = sample(c(0,1),n,prob=c(1-p1,p1),replace=TRUE)
 Z = r*X1+(1-r)*X2
 hist(Z,prob=TRUE,main=paste(p1,"X1+",1-p1,"X2"),ylim=c(0,0.8))
 y = seq(-4,6,.01)
 lines(y,p1*dnorm(y,0,1)+(1-p1)*dnorm(y,3,1),col=cbPalette[sample(seq(1:8),1,replace=F)],lwd=2)
}

#' @title Generate random samples from Wishart distribution $W_d(\Sigma,n)$
#'
#' @param d the dimension of the square matrix
#' @param n degrees of freedom
#' @param Sigma scale matrix
#' @description This function generate random samples from Wishart distribution based on Bartlett’s decomposition.
#' @details Let T=($T_{ij}$) be a lower triangular $d\times d$ random matrix with independent entries satisfying
#' -$T_{ij}\sim N(0,1),i>j,i.i.d$
#' -$T_{ii}\sim \sqrt{\chi^2(n-i+1)},i=1,...,d$
#' Then the matrix $A=TT^{T}$ has a $W_d(I_d,n)$ distribution
#' Obtain the Choleski factorization $\Sigma=LL^{T}$,where L is a lower triangular
#' Then $LAL^{T}\sim W_d(\Sigma,n)$
#' @return The generated samples from $W_d(\Sigma,n)$
#' @export
#'
#' @examples
#' \dontrun{
#'  Sigma=matrix(c(3,2,2,4),nrow=2)
#'  gen_Wisha(2,5,Sigma)
#' }
gen_Wisha = function (d,n,Sigma){
 T=matrix(0,nrow=d,ncol=d)
 for (i in 1:d){
  for (j in 1:d){
   if (i>j) T[i,j]=rnorm(1)
   else if (i==j) T[i,i]=sqrt(rchisq(1,n-i+1))
   else next
  }
 }
 L=t(chol(Sigma)) #Choleski decomposition
 return(L%*%T%*%t(T)%*%t(L))
}


#' @title Monte Carlo integration with antithetic variables
#' @description This function is defined to estimate $\int_{0}^{1}\frac{e^{-x}}{1+x^2}dx$
#' and test the effect of antithetic variables in variance reduction.
#' @param R number of generated samples to estimate the integration
#' @param anti use antithetic variables or not
#'
#' @return The numerical estimation of $\int_{0}^{1}\frac{e^{-x}}{1+x^2}dx$
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(2222)
#' m=1000;
#' MC.de=MC=numeric(m)
#' for (i in 1:m){
#'  MC.de[i]=MC.Phi(R=1000)
#'  MC[i]=MC.Phi(R=1000,anti=FALSE)
#'  }
#' print("Approximate reduction in variance: ");(var(MC)-var(MC.de))/var(MC)
#' }
MC.Phi = function(R = 1e4, anti = TRUE) {
 u = runif(R/2)
 if (!anti) v = runif(R/2) else v = (1-u)
 u = c(u, v)
 MC = mean(exp(-u)/(1+u^2))
 return(MC)
}

#' @title Computes the sample skewness coefficient
#'
#' @param x samples
#'
#' @return the sample skewness coefficient
#' @export
#'
#' @examples
#' \dontrun{
#'  x=rnorm(100);skewness=sk(x)
#' }
sk = function(x) {
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

#' @title Calculate the standard error of q sample quantile using exact variance formula of skewness
#' @details  The variance of q sample quantile of skewness is $Var(\hat{x_q})=\frac{q(1-q)}{nf^2(x)}$
#' where f is the density of the sampled distribution. 
#' Using the normal approximation for the density, the skewness is asymptotically normal with mean 0 and variance $\frac{6(n-2)}{(n+1)(n+3)}$
#' @param x the quantiles of skewness coefficient
#' @param q quantiles
#' @param n number of samples
#'
#' @return the standard error of q sample quantile of skewness coefficient
#' @export
#'
#' @examples
#' \dontrun{
#'  m=1000;n=100
#'  stat = replicate(m, expr={
#'  x = rnorm(n)
#'  y = sk(x)
#'  })
#'  q = c(.025,.05,.95,.975)
#'  x = quantile(stat,q)
#'  print('The standard error of the estimates:');sd_sk(x,q,n)
#' }
sd_sk = function(x,q,n){
  sv = 6*(n-2)/(n+1)/(n+3)
  var_q=q*(1-q)/n/dnorm(x,mean=0,sd=sqrt(sv))
  return(sqrt(var_q))
}


#' @title The graph of the power of the skewness test of symmetric Beta$(\alpha, \alpha)$ distribution
#'
#' @param a the parameter for symmetric Beta distribution
#'
#' @return a graph of the power of the skewness test
#' @export
#'
#' @examples \dontrun{
#'  a1 = c(seq(0,1,0.1))
#'  pwr_beta(a1)
#'  a2 = c(seq(1,20,2),seq(20,100,5))
#'  pwr_beta(a2)
#' }
pwr_beta = function(a){
 alpha = .1;n = 30;m = 10000;N = length(a)
 pwr = numeric(N)
 cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
 for (j in 1:N) { #for each a
  sktests = numeric(m)
  for (i in 1:m) { 
   x = rbeta(n, a[j], a[j])
   sktests[i] = as.integer(abs(sk(x))>= cv)
  }
  pwr[j] = mean(sktests)
 }
 se = sqrt(pwr * (1-pwr) / m) #add standard errors
 plot(a, pwr, type = "b",xlab = "a",ylab="pwr", ylim = c(0,0.12),pch=20,
 main=paste('a from',min(a),'to',max(a)))
 abline(h = .1, lty = 2)
 lines(a, pwr+se, lty = 3)
 lines(a, pwr-se, lty = 3)
 text(a[ceiling(0.3*N)],0.105,paste("alpha =",alpha))
}

#' @title The graph of the power of the skewness test of t distribution
#'
#' @param v the parameter of t distribution
#'
#' @return a graph of the power of the skewness test
#' @export
#'
#' @examples \dontrun{
#' v=seq(1,20);pwr_t(v)
#' }
pwr_t = function(v){
 alpha = .1;n = 30;m = 2500;N = length(v)
 pwr = numeric(N)
 cv = qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
 for (j in 1:N) { #for each v
  sktests = numeric(m)
  for (i in 1:m) { 
   x = rt(n,v[j])
   sktests[i] = as.integer(abs(sk(x))>= cv)
   }
   pwr[j] = mean(sktests)
 }
 se = sqrt(pwr*(1-pwr) / m) #add standard errors
 plot(v, pwr, type = "b",xlab = "v",ylab="pwr",pch=20,ylim=c(0,1),
 main="Power of skewness test against t(v)")
 abline(h = .1, lty = 2)
 lines(v, pwr+se, lty = 3)
 lines(v, pwr-se, lty = 3)
 text(v[ceiling(0.8*N)],0.8,paste("alpha =",alpha))
}

  
#' @title Compute the maximum number of extreme points
#'
#' @param z the samples
#' @param ix the index after permutation
#' @param n1 the length of the first sample
#'
#' @return the number of extreme points
#' @export
#'

count5num = function(z,ix,n1) {
 z = z[ix]
 x = z[1:n1]
 y = z[-(1:n1)]
 X = x - mean(x)
 Y = y - mean(y)
 outx = sum(X > max(Y)) + sum(X < min(Y))
 outy = sum(Y > max(X)) + sum(Y < min(X))
 return(max(c(outx, outy)))
}

#' @title A permutation test for equal variance based on the maximum number of extreme points
#'
#' @param y the samples
#' @param n1 the length of the first sample
#'
#' @return 1 for reject null hypothesis (equal variance), 0 for accept
#' @export
#'
#' @examples \dontrun{
#' compute the typr-I error for the permutation test
#' n1=20;n2=30;m=500
#' a2 = mean(replicate(m, expr={ 
#' x1 = rnorm(n1);x2 = rnorm(n2)
#' y=c(x1,x2)
#' count5_test(y,n1)
#' }))
#' }
count5_test =function(y,n1){ 
 boot.obj = boot(data = y, statistic = count5num,R = 2000,sim = "permutation",n1=n1)
 tb = c(boot.obj$t0, boot.obj$t)
 p.value = mean(tb>=tb[1]) #larger t favors alternative
 return(as.integer(p.value<0.05))
}

#' @title The original "Count Five” test for equal variance
#'
#' @param x sample 1
#' @param y sample 2
#'
#' @return 1 for reject null hypothesis (equal variance), 0 for accept
#' @export
#'
#' @examples \dontrun{
#' compute the type-I error for the original test
#' n1=20;n2=30;m=500
#' a1 = mean(replicate(m, expr={ 
#'  x = rnorm(n1);y <- rnorm(n2)
#'  x = x - mean(x) #centered by sample mean
#'  y = y - mean(y)
#'  count5test(x, y)
#' }))
#' }
count5test = function(x, y) {
 X = x - mean(x)
 Y = y - mean(y)
 outx = sum(X > max(Y)) + sum(X < min(Y))
 outy = sum(Y > max(X)) + sum(Y < min(X))
 return(as.integer(max(c(outx, outy)) > 5))
}

f = function(x) {
 return(exp(-abs(x)))
}

#' @title A random walk Metropolis sampler for generating the standard Laplace distribution
#'
#' @param sigma the standard deviation of the sample-generating normal distribution
#' @param x0 Initial values for the random walk
#' @param N number of the generated samples
#'
#' @return alist of a Markov Chain sampled from standard Laplace distribution and the times of rejection
#' @export
#'
#' @examples \dontrun{
#' N=2000
#' sigma = c(.05, .5, 2, 16)
#' x0=20
#' R1=random.M(sigma[1],x0,N)
#' R2=random.M(sigma[2],x0,N)
#' R3=random.M(sigma[3],x0,N)
#' R4=random.M(sigma[4],x0,N)
#' }
random.M = function(sigma, x0, N){
 x = numeric(N)
 x[1] = x0
 u = runif(N)
 rej = 0
 for (i in 2:N) {
  y = rnorm(1, x[i-1], sigma)
  if (u[i] <= (f(y) / f(x[i-1]))) x[i] = y 
  else {x[i] = x[i-1];rej=rej+1}
 }
 return(list(x=x, rej=rej))
}



#' @title Compute $\frac{2\Gamma(\frac{k+1}{2})}{\sqrt{\pi k} \Gamma(\frac{k}{2})}\int_0^{c_{k}}(1+\frac{u^2}{k})^{-(k+1)/2}du$
#' @description $c_k=\sqrt{\frac{a^2k}{k+1-a^2}}$
#' @param k integer
#' @param a real number
#'
#' @return the numerical integration
#' @export

ff = function(k,a){
 ck=sqrt(a^2*k/(k+1-a^2))
 int=integrate(function(u){(1+u^2/k)^(-(k+1)/2)},0,ck)$value
 mul=2/sqrt(pi*k)*exp(lgamma((k+1)/2)-lgamma(k/2))
 mul*int
}

#' @title Solve the equation of $\frac{2\Gamma(\frac{k}{2})}{\sqrt{\pi (k-1)} \Gamma(\frac{k-1}{2})}\int_0^{c_{k-1}}(1+\frac{u^2}{k-1})^{-k/2}du$
#' $=\frac{2\Gamma(\frac{k+1}{2})}{\sqrt{\pi k} \Gamma(\frac{k}{2})}\int_0^{c_{k}}(1+\frac{u^2}{k})^{-(k+1)/2}du$
#'
#' @param k integer
#'
#' @return the root of the equation
#' @export
#'
#' @examples \dontrun{
#' Root=matrix(0,2,5)
#' k=c(4,25,100,500,1000)
#' for (i in 1:length(k)){
#'  Root[2,i]=round(solve_eq(k[i]),4)
#' }
#' Root[1,]=k;rownames(Root)=c('k','root')
#' Root
#' }
solve_eq = function(k){
 out=uniroot(function(a){ff(k,a)-ff(k-1,a)},lower=1,upper=2)
 out$root
}



#' @title A multicore version of sapply() (For Windows only)
#'
#' @param cluster 
#' @param X object of the function
#' @param FUN function
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples \dontrun{
#' cores = detectCores()
#' cluster = makePSOCKcluster(cores)
#' boot_i = function(i){
#' r_adj = function(x,id) {
#'   x = x[id,]     #bootstrap sample
#'   res = lm(mpg~wt+disp,data=x)
#'   summary(res)$adj.r.squared
#' }
#' as.numeric(boot::boot(mtcars,r_adj,1)$t)
#' }
#' system.time(sapply(1:500,boot_i))
#' system.time(mcsapply(cluster,1:500,boot_i))
#' }
mcsapply = function(cluster,X,FUN,...){
 res=parLapply(cluster,X,FUN,...) #Use parLapply in Windows
 simplify2array(res)
}

