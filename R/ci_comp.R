#' @title A function used to calculate the log-likelihood of sample mean
#'
#' @param x samples
#' @param mu theoretical mean
#' @details When n is large, np[i] %~~% frac(1,1+frac((bar(y)-mu)(y[i]-mu),S(mu))$$
#' where S(mu)=frac(1,n) sum((y[i]-mu)^2,i==1,n)
#' @return the value of sum(log(np[i](mu)),i==1,n)
#' @export
log_npi = function(x,mu){
 n=length(x)
 Smu=mean((x-mu)^2)
 del=mean(x)-mu
 S=0
 for (i in 1:n){
  d=(1+del*(x[i]-mu)/Smu)^(-1)
  S=S+log(d)
 }
 S
}

#' @title A function that uses uniroot to solve for the (1-alpha) confidence interval of sample mean
#' by the method of empirical likelihood
#'
#' @param y samples
#' @param alpha significance level
#' @param ub upper bound for solving for the root 
#'
#' @return a (1-alpha) confidence interval by the method of empirical likelihood
#' @details For empirical likelihood method, the (1-alpha) CI for mu is
#' {mu|sum(log(np[i](mu)),i==1,n)>-0.5 chi[1-alpha]^2*(1)}
#' Noticing that mu doesn't have explicit formula to calculate, 
#' I solved the equation numerically.
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(7524)
#' n=c(500,5000,50000)
#' lik=matrix(0,ncol=2,nrow=3)
#' rownames(lik)=paste("n=",n);colnames(lik)=c("left","right")
#' up=c(1.23,1.07,1.07)
#' for (i in 1:3){
#'  x=rnorm(n[i])
#'  y=x^2 #turn var into mean
#'  lik[i,]=CI_ELR(y,0.05,up[i])}
#' lik
#' }
CI_ELR = function(y,alpha,ub){
 b= -qchisq(1-alpha,df=1)/2
 out1=uniroot(function(mu){log_npi(y,mu)-b},lower=0.5,upper=1)
 root1=out1$root
 out2=uniroot(function(mu){log_npi(y,mu)-b},lower=1,upper=ub)
 root2=out2$root
 return(round(c(root1,root2),4))
}

#' @title A function that computes the the (1-alpha) confidence interval of sample mean
#' by the method of nonparametric delta method
#'
#' @param y samples
#' @param alpha significance level
#'
#' @return a (1-alpha) confidence interval by the method of nonparametric delta method
#' @details For nonparametric delta method of mean value, hat(tau)^2=frac(1,n)sum((Y[i]-bar{Y})^2,i==1,n)
#'  group('[',list(bar(Y)-frac(z[1-frac(alpha,2)]%*%tau,sqrt(n)),bar(Y)+frac(z[1-frac(alpha,2)]%*%tau,sqrt(n))
#'  ),']')
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(7524)
#' n=c(500,5000,50000)
#' int=matrix(0,ncol=2,nrow=3)
#' rownames(int)=paste("n=",n);colnames(int)=c("left","right")
#' for (i in 1:3){
#'  x=rnorm(n[i])
#'  y=x^2 #turn var into mean
#'  int[i,]=CI_NDM(y,0.05)}
#' int
#' }
CI_NDM = function(y,alpha){
  TFn=mean(y);n=length(y)
  tau=sqrt(mean((y-TFn)^2))
  d=qnorm(1-alpha/2)*tau/sqrt(n)
  return(round(c(TFn-d,TFn+d),4))
  
}






