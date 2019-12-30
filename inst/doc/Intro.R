## ---- echo=FALSE, message=FALSE, warning=TRUE---------------------------------
knitr::opts_chunk$set(comment = "#", warning = TRUE, eval = TRUE, message = FALSE)
library(SC19090)

## -----------------------------------------------------------------------------
library(clinfun,quietly = TRUE)
medical=c(125,136,116,101,105,109,
122,114,132,120,119,127,
128,142,128,134,135,132,140,129)
gr.medical=c(rep(1,6),rep(2,6),rep(3,8))
jonckheere.test(medical,gr.medical)

## -----------------------------------------------------------------------------
data(m)
attach(m)
jh_extend(medical,gr)

## -----------------------------------------------------------------------------
set.seed(7524)
n=c(500,5000,50000)
int=matrix(0,ncol=2,nrow=3)
lik=matrix(0,ncol=2,nrow=3)
rownames(int)=rownames(lik)=paste("n=",n)
colnames(int)=colnames(lik)=c("left","right")
up=c(1.23,1.07,1.07)
for (i in 1:3){
 x=rnorm(n[i])
 y=x^2 #turn var into mean
 int[i,]=CI_NDM(y,0.05)
 lik[i,]=CI_ELR(y,0.05,up[i])
}
print("95% Confidence Interval for Nonparametric Delta Method");int
print("95% Confidence Interval for Empirical Likelihood Ratio");lik

## ----results='hide'-----------------------------------------------------------
library(lattice)
library(dplyr)
library(knitr)

## ----pressure, echo=TRUE------------------------------------------------------
n <- seq(5, 45, 5)
x <- rnorm(sum(n))
y <- factor(rep(n, n), labels=paste("n =", n))
densityplot(~ x | y,
panel = function(x, ...) {
panel.densityplot(x, col="DarkOliveGreen", ...)
panel.mathdensity(dmath=dnorm,
args=list(mean=mean(x), sd=sd(x)),
col="darkblue")
})

## -----------------------------------------------------------------------------
dt <- mtcars[1:5, 1:6]
kable(dt)

## -----------------------------------------------------------------------------
cbPalette = c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
              "#D55E00", "#CC79A7")
gen_Raylei= function(sigma){
 n = 1000;
 u = runif(n);
 x = sqrt(-2*sigma^2*log(1-u));
 hist(x, prob = TRUE,ylim=c(0,0.65),main=paste("Sigma = ",sigma))
 y = seq(0,ceiling(max(x)),.01)
 lines(y,y/sigma^2*exp(-y^2/2/sigma^2),col=cbPalette[sigma],lwd=2)
}
par(mfrow=c(1,1))
for (sigma in 1:4){
 gen_Raylei(sigma)
}

## -----------------------------------------------------------------------------
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
par(mfrow=(c(1,1)))
gen_Mixnorm(0.25)
gen_Mixnorm(0.5)
gen_Mixnorm(0.75)
gen_Mixnorm(1)

## -----------------------------------------------------------------------------
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
Sigma=matrix(c(3,2,2,4),nrow=2)
gen_Wisha(2,5,Sigma)


## -----------------------------------------------------------------------------
set.seed(1234)
m = 1e4;
x = runif(m,min=0,max=pi/3)
theta.hat = mean(sin(x)*pi/3)
re=matrix(c(theta.hat,cos(0)-cos(pi/3)),nrow=1)
re=round(re,4)
colnames(re)=c("Estimate","Exact Value")
re

## -----------------------------------------------------------------------------
set.seed(2222)
MC.Phi <- function(R = 1e4, anti = TRUE) {
  u <- runif(R/2)
  if (!anti) v <- runif(R/2) else v <- (1-u)
  u <- c(u, v)
  MC = mean(exp(-u)/(1+u^2))
  return(MC)
  }
m=1000;
MC.de=MC=numeric(m)
for (i in 1:m){
 MC.de[i]=MC.Phi(R=1000)
 MC[i]=MC.Phi(R=1000,anti=FALSE)
}
print("Approximate reduction in variance: ")
(var(MC)-var(MC.de))/var(MC)


## -----------------------------------------------------------------------------
set.seed(5678)
M = 1e4; k = 5 
r = M/k             #replicates per stratum
N = 50              #number of times to repeat the estimation
T2 = numeric(k)     #Store the results of each stratum
est = matrix(0, N, 2)
g = function(x)exp(-x)/(1+x^2)*(x>0)*(x<1) 
f1 = function(x)exp(-x)/(1-exp(-1))
f2 = function(x)exp(-x)/(exp(-(j-1)/k)-exp(-j/k))

for (i in 1:N) {
  u = runif(M);x = -log(1-u*(1-exp(-1))) #u~U[0,1]; inverse trans; 
  est[i, 1] <- mean(g(x)/f1(x))            
  for(j in 1:k){
   u = runif(r)  
   x = -log(exp(-(j-1)/k)-(exp(-(j-1)/k)-exp(-j/k))*u) 
   T2[j]<-mean(5*g(x)/f2(x))
  }
  est[i,2] <- mean(T2)
}
re = cbind(apply(est,2,mean),apply(est,2,sd))
rownames(re)=c("Important Sampling","Im Sa + Str Sa")
colnames(re)=c("theta_hat","sd")
knitr::kable(re)


## -----------------------------------------------------------------------------
n=20;m=1000
set.seed(5557)


CP = replicate(m,expr={
 x = rchisq(n,df=2)
 d = qt(.975,n-1)*sd(x)/sqrt(n)
 c(mean(x)-d,mean(x)+d)
 })
mean(2> CP[1,]& 2< CP[2,])

UCL = replicate(m, expr ={
  x = rchisq(n,df=2)
 (n-1)*var(x)/qchisq(.05,df=n-1)
})

mean(UCL>4)


## -----------------------------------------------------------------------------
set.seed(1234)
sk = function(x) {
 #computes the sample skewness coeff.
 xbar <- mean(x)
 m3 <- mean((x - xbar)^3)
 m2 <- mean((x - xbar)^2)
 return( m3 / m2^1.5 )
}
m=1000;n=100
stat = replicate(m, expr={
 x = rnorm(n)
 y = sk(x)
})
q = c(.025,.05,.95,.975)
x = quantile(stat,q)

sd_sk = function(x,q,n){
 sv = 6*(n-2)/(n+1)/(n+3)
 var_q=q*(1-q)/n/dnorm(x,mean=0,sd=sqrt(sv))
 return(sqrt(var_q))
}
print('The standard error of the estimates:')
sd_sk(x,q,n)
cv = qnorm(q,0,sqrt(6/n))
quat=rbind(x,cv,x-cv)
rownames(quat)=c('estimated quantiles','quat of large sample approx','bias (est-app)')
knitr::kable(quat)


## -----------------------------------------------------------------------------
set.seed(1234)
sk = function(x) {
  xbar = mean(x)
  m3 = mean((x - xbar)^3)
  m2 = mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

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

a1 = c(seq(0,1,0.1))
a2 = c(seq(1,20,2),seq(20,100,5))
par(mfrow=c(1,1))
pwr_beta(a1);pwr_beta(a2)

## -----------------------------------------------------------------------------
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
v=seq(1,20);pwr_t(v)


## -----------------------------------------------------------------------------
m = 1e4; n = 10; mu0 = 1 
mu.hat = mu.se = p.val1 = numeric(m)
 for(i in 1:m){
 x = rchisq(n,1)
 mu.hat[i] <- mean(x);mu.se[i] <- sd(x)
 p.val1[i] <- t.test(x,mu=1)$p.value
}
T=sqrt(n)*(mu.hat-mu0)/mu.se
p.val2=2*pt(abs(T),(n-1),lower.tail=FALSE)
print(c(mean(p.val1<=0.05),mean(p.val2<=0.05)))

## -----------------------------------------------------------------------------
for(i in 1:m){
 x=runif(n,0,2)
 mu.hat[i] <- mean(x);mu.se[i] <- sd(x)
 p.val1[i] <- t.test(x,mu=1)$p.value
}
T=sqrt(n)*(mu.hat-mu0)/mu.se
p.val2=2*pt(abs(T),(n-1),lower.tail=FALSE)
print(c(mean(p.val1<=0.05),mean(p.val2<=0.05)))

## -----------------------------------------------------------------------------
for(i in 1:m){
 x=rexp(n,1)
 mu.hat[i] <- mean(x);mu.se[i] <- sd(x)
 p.val1[i] <- t.test(x,mu=1)$p.value
}
T=sqrt(n)*(mu.hat-mu0)/mu.se
p.val2=2*pt(abs(T),(n-1),lower.tail=FALSE)
print(c(mean(p.val1<=0.05),mean(p.val2<=0.05)))

## -----------------------------------------------------------------------------
set.seed(8888)
data(scor, package = "bootstrap")
scor=data.matrix(scor)
pairs(scor,main='Scatter plots for each pair of test scores',pch=20)
cor_mat = function(x){
 n=dim(x)[1]
 one=rep(1,n);dim(one)=c(n,1);ONE=diag(n)
 Po=one%*%t(one)/n
 S=t(x)%*%(ONE-Po)%*%x/(n-1)    
 D=c();for (i in 1:dim(x)[2]) {D=c(D,S[i,i])};D=diag(D)
 R=sqrt(solve(D))%*%S%*%sqrt(solve(D))
 colnames(R)=rownames(R)=colnames(x)  
 return(R=round(R,4))
}
cor_mat(scor)

## -----------------------------------------------------------------------------
rho_12 = function(x,i) cor(x[i,1],x[i,2])
rho_34 = function(x,i) cor(x[i,3],x[i,4])
rho_35 = function(x,i) cor(x[i,3],x[i,5])
rho_45 = function(x,i) cor(x[i,4],x[i,5])

r12_b = boot::boot(data=scor,statistic=rho_12,R=2000)
r34_b = boot::boot(data=scor,statistic=rho_34,R=2000)
r35_b = boot::boot(data=scor,statistic=rho_35,R=2000)
r45_b = boot::boot(data=scor,statistic=rho_45,R=2000)

OUT=matrix(c(sd(r12_b$t),sd(r34_b$t),sd(r35_b$t),sd(r45_b$t)),nrow=1,byrow=T)
OUT=round(OUT,4)
colnames(OUT)=c("p_12","p_34","p_35","p_45");rownames(OUT)="standard error"
knitr::kable(OUT)

## -----------------------------------------------------------------------------
SK=0
boot_sk = function(x,i) {
  xbar = mean(x[i])
  m3 = mean((x[i] - xbar)^3)
  m2 = mean((x[i] - xbar)^2)
  return( m3 / m2^1.5 )
}
n=10;m=1000;
CI.n=CI.b=CI.p=matrix(0,m,2)
for(i in 1:m){
 x=rnorm(n)
 b = boot::boot(data=x,statistic=boot_sk, R = 1000)
 CI = boot::boot.ci(b,type=c("norm","basic","perc"))
 CI.n[i,]=CI$norm[2:3]
 CI.b[i,]=CI$basic[4:5]
 CI.p[i,]=CI$percent[4:5]
}
out=out2=out3=matrix(0,nrow=2,ncol=3)
out[1,]=c(mean(CI.n[,1]<=SK & CI.n[,2]>=SK),mean(CI.b[,1]<=SK & CI.b[,2]>=SK),
mean(CI.p[,1]<=SK & CI.p[,2]>=SK))

out2[1,]=c(mean(CI.n[,1]>SK),mean(CI.b[,1]>SK),mean(CI.p[,1]>SK)) 
#miss on the left

out3[1,]=c(mean(CI.n[,2]<SK),mean(CI.b[,2]<SK),mean(CI.p[,2]<SK)) 
#miss on the right

## -----------------------------------------------------------------------------
sk=sqrt(8/5)
Ci.n=Ci.b=Ci.p=matrix(0,m,2)
for(i in 1:m){
 x=rchisq(n,5)
 b = boot::boot(data=x,statistic=boot_sk, R = 1000)
 Ci = boot::boot.ci(b,type=c("norm","basic","perc"))
 Ci.n[i,]=Ci$norm[2:3]
 Ci.b[i,]=Ci$basic[4:5]
 Ci.p[i,]=Ci$percent[4:5]
}
out[2,]=c(mean(Ci.n[,1]<=sk & Ci.n[,2]>=sk),mean(Ci.b[,1]<=sk & Ci.b[,2]>=sk),
mean(Ci.p[,1]<=sk & Ci.p[,2]>=sk))

out2[2,]=c(mean(Ci.n[,1]>sk),mean(Ci.b[,1]>sk),mean(Ci.p[,1]>sk))
#miss on the left

out3[2,]=c(mean(Ci.n[,2]<sk),mean(Ci.b[,2]<sk),mean(Ci.p[,2]<sk))
#miss on the right

rownames(out)=rownames(out2)=rownames(out3)=c("N(0,1)","X^2(5)")
colnames(out)=colnames(out2)=colnames(out3)=c("norm","basic","perc")
print("Empirical Coverage Rates");out
print("The proportion miss on the left");out2
print("The proportion miss on the right");out3

## -----------------------------------------------------------------------------
library(bootstrap,quietly=TRUE)
data(scor)
scor=data.matrix(scor)
Sigma=cov(scor)
lambda=eigen(Sigma)$values

theta_compute = function(lambda){
 return(lambda[1]/sum(lambda))
}

theta_hat = theta_compute(lambda)
n=dim(scor)[1]
theta_jack = numeric(n)
for(i in 1:n){
 Sigma_i=cov(scor[-i,])
 lambda_i=eigen(Sigma_i)$values
 theta_jack[i] = theta_compute(lambda_i)
}

bias_jack = (n-1)*(mean(theta_jack)-theta_hat)
var_jack = (n-1)^2*var(theta_jack)/n
se_jack = sqrt(var_jack)
round(c(original=theta_hat,bias_jack=bias_jack,
se_jack=se_jack),3)


## -----------------------------------------------------------------------------
data(ironslag,package="DAAG")
attach(ironslag)
n = length(magnetic) 
e1 = e2 = e3 = e4 = numeric(n)
for (k in 1:n) {
 y = magnetic[-k]
 x = chemical[-k]

 J1 = lm(y ~ x) 
 yhat1 = J1$coef[1] + J1$coef[2] * chemical[k]
 e1[k] = magnetic[k] - yhat1

 J2 = lm(y ~ x + I(x^2))
 yhat2 = J2$coef[1] + J2$coef[2] * chemical[k] +
 J2$coef[3] * chemical[k]^2
 e2[k] = magnetic[k] - yhat2
 

 J3 = lm(log(y) ~ x)
 logyhat3 = J3$coef[1] + J3$coef[2] * chemical[k]
 yhat3 = exp(logyhat3)
 e3[k] = magnetic[k] - yhat3
 

 J4 = lm(y ~ x + I(x^2)+ I(x^3))
 yhat4 = J4$coef[1] + J4$coef[2] * chemical[k] +
 J4$coef[3] * chemical[k]^2 + J4$coef[4] * chemical[k]^3
 e4[k] = magnetic[k] - yhat4
 

}

round(c(spe_1=mean(e1^2), spe_2=mean(e2^2), 
spe_3=mean(e3^2), spe_4=mean(e4^2)),4)

## -----------------------------------------------------------------------------
y = magnetic
x = chemical
L1 = lm(y ~ x) 
L2 = lm(y ~ x + I(x^2))
L3 = lm(log(y) ~ x)
L4 = lm(y ~ x + I(x^2)+ I(x^3))
round(c(adj_R2_1=summary(L1)$adj.r.squared,
        adj_R2_2=summary(L2)$adj.r.squared,
        adj_R2_3=summary(L3)$adj.r.squared,
        adj_R3_3=summary(L4)$adj.r.squared),4)

## -----------------------------------------------------------------------------
L2 = lm(magnetic ~ chemical + I(chemical^2))
coef=round(L2$coef,3)
paste("Magnetic =",coef[1],coef[2],"*Chemical +" ,coef[3],"*Chemical^2")


## -----------------------------------------------------------------------------
library(boot,quietly = TRUE)
count5num <- function(z,ix,n1) {
 z = z[ix]
 x = z[1:n1]
 y = z[-(1:n1)]
 X <- x - mean(x)
 Y <- y - mean(y)
 outx <- sum(X > max(Y)) + sum(X < min(Y))
 outy <- sum(Y > max(X)) + sum(Y < min(X))
 # return the number of extreme points
 return(max(c(outx, outy)))
}
count5_test =function(y,n1){ #Permutation version
 boot.obj = boot(data = y, statistic = count5num
 ,R = 2000,sim = "permutation",n1=n1)
 tb = c(boot.obj$t0, boot.obj$t)
 p.value = mean(tb>=tb[1]) #larger t favors alternative
 return(as.integer(p.value<0.05))
}
count5test = function(x, y) {
 X = x - mean(x)
 Y = y - mean(y)
 outx = sum(X > max(Y)) + sum(X < min(Y))
 outy = sum(Y > max(X)) + sum(Y < min(X))
 # return 1 (reject) or 0 (do not reject H0)
 return(as.integer(max(c(outx, outy)) > 5))
}
n1=20;n2=30;m=500
a1 = mean(replicate(m, expr={ #Count 5 criterion
 x = rnorm(n1);y <- rnorm(n2)
 x = x - mean(x) #centered by sample mean
 y = y - mean(y)
 count5test(x, y)
}))
a2 = mean(replicate(m, expr={ #Permutation test
 x1 = rnorm(n1);x2 = rnorm(n2)
 y=c(x1,x2)
 count5_test(y,n1)
}))
round(c(c5_t1e=a1,pt_t1e=a2),3)

## -----------------------------------------------------------------------------
n2=50;m=500
a1 = mean(replicate(m, expr={ #Count 5 criterion
 x = rnorm(n1);y <- rnorm(n2)
 x = x - mean(x) #centered by sample mean
 y = y - mean(y)
 count5test(x, y)
}))
a2 = mean(replicate(m, expr={ #Permutation test
 x1 = rnorm(n1);x2 = rnorm(n2)
 y=c(x1,x2)
 count5_test(y,n1)
}))
round(c(c5_t1e=a1,pt_t1e=a2),3)

## ----eval=FALSE---------------------------------------------------------------
#  library(snowfall,quietly=TRUE)
#  Akl = function(x) {
#   d = as.matrix(dist(x))
#   m = rowMeans(d); M = mean(d)
#   a = sweep(d, 1, m); b = sweep(a, 2, m)
#   return(b+M)
#  }
#  dCov = function(x, y) {
#   x = as.matrix(x); y = as.matrix(y)
#   n = nrow(x); m = nrow(y)
#   if (n != m || n < 2) stop("Sample sizes must agree")
#   if (! (all(is.finite(c(x, y)))))
#   stop("Data contains missing or infinite values")
#   A = Akl(x); B = Akl(y)
#   sqrt(mean(A * B))
#  }
#  ndCov2 = function(z, ix, dims) {
#   #dims contains dimensions of x and y
#   p = dims[1]
#   q = dims[2]
#   d = p + q
#   x = z[ , 1:p] #leave x as is
#   y = z[ix, -(1:p)] #permute rows of y
#   return(nrow(z) * dCov(x, y)^2)
#  }
#  dcov.test= function(z){
#   boot.obj = boot(data = z, statistic = ndCov2, R = 999,
#   sim = "permutation", dims = c(2, 2))
#   tb = c(boot.obj$t0, boot.obj$t)
#   return(p.cor = mean(tb>=tb[1]))
#  }
#  run_cop1 = function(n){
#  for (i in 1:10){
#   x = rmvnorm(n,mu=c(0,0),sigma=diag(1,2))
#   e = rmvnorm(n,mu=c(0,0),sigma=diag(1,2))
#   y = x/4+e
#   z = cbind(x,y)
#   p1.values[i,1] = dcov.test(z)
#   p1.values[i,2] = bcov.test(z[,1:2],z[,3:4],R=999,seed=i+1)$p.value
#  }
#   return(colMeans(p1.values<0.1))
#  }
#  p1.values=matrix(0,ncol=2,nrow=10)
#  sfInit(parallel = TRUE, cpus = 10)
#  sfLibrary(mixtools)
#  sfLibrary(Ball)
#  sfLibrary(boot)
#  sfExport("Akl","dCov","ndCov2","dcov.test")
#  sfExport("p1.values")
#  n = seq(50,200,20)
#  result1 = sfLapply(n, run_cop1)  #Model 1
#  result1=matrix(unlist(result1),ncol=2,byrow=T)
#  plot(n,result1[,1],type="b",pch=20,ylab="Power",main="Model 1",ylim=c(0,1.1))
#  lines(n,result1[,2],type="b",lty=2)
#  legend(140,0.3,c("Distance","Ball"),lty=1:2)
#  t1=t(result1);colnames(t1)=as.character(n);rownames(t1)=c("Distance","Ball")
#  knitr::kable(t1)

## ----eval=FALSE---------------------------------------------------------------
#  run_cop2 = function(n){
#  for (i in 1:10){
#   x = rmvnorm(n,mu=c(0,0),sigma=diag(1,2))
#   e = rmvnorm(n,mu=c(0,0),sigma=diag(1,2))
#   y = x/4*e
#   z = cbind(x,y)
#   p2.values[i,1] = dcov.test(z)
#   p2.values[i,2] = bcov.test(z[,1:2],z[,3:4],R=999,seed=i+1)$p.value
#  }
#   return(colMeans(p2.values<0.1))
#  }
#  p2.values=matrix(0,ncol=2,nrow=10)
#  sfExport("p2.values")
#  n = seq(50,200,20)
#  result2 = sfLapply(n, run_cop2)  #Model 2
#  result2=matrix(unlist(result2),ncol=2,byrow=T)
#  plot(n,result2[,1],type="b",pch=20,ylab="Power",main="Model 2",ylim=c(0,1.1))
#  lines(n,result2[,2],type="b",lty=2)
#  legend(140,0.3,c("Distance","Ball"),lty=1:2)
#  t2=t(result2);colnames(t2)=as.character(n);rownames(t2)=c("Distance","Ball")
#  knitr::kable(t2)

## -----------------------------------------------------------------------------
set.seed(8888)
f = function(x) {
  return(exp(-abs(x)))
}
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
N=2000
sigma = c(.05, .5, 2, 16)
x0=20
R1=random.M(sigma[1],x0,N)
R2=random.M(sigma[2],x0,N)
R3=random.M(sigma[3],x0,N)
R4=random.M(sigma[4],x0,N)
#number of candidate points rejected
Rej=cbind(R1$rej, R2$rej, R3$rej, R4$rej)
Rej= round((N-Rej)/N,3)
rownames(Rej)="Accept rates";colnames(Rej)=paste("sigma",sigma)
knitr::kable(Rej)
par(mfrow=c(1,1))
plot(1:N,R1$x,type='l',ylab='N',main=paste("sigma=",sigma[1]))
plot(1:N,R2$x,type='l',ylab='N',main=paste("sigma=",sigma[2]))
plot(1:N,R3$x,type='l',ylab='N',main=paste("sigma=",sigma[3]))
plot(1:N,R4$x,type='l',ylab='N',main=paste("sigma=",sigma[4]))
library(rmutil,quietly=TRUE)
p = c(.05, seq(.1, .9, .1), .95)
Q = qlaplace(p)
R_all = cbind(R1$x, R2$x, R3$x, R4$x)
MC = R_all[501:N, ]
QR = apply(MC, 2, function(x) quantile(x, p))
colnames(QR)=paste("sigma",sigma)
knitr::kable(round(cbind(real_q=Q, QR),3))
#xtable::xtable(round(cbind(real_q=Q, QR),3))
real=rlaplace(N)
par(mfrow=c(1,1))
qqplot(MC[,1],real,pch=20,col='red',main=paste("sigma=",sigma[1]));
qqline(MC[,1],distribution=qlaplace,lwd=1.9)
qqplot(MC[,2],real,pch=20,col='red',main=paste("sigma=",sigma[2]));
qqline(MC[,2],distribution=qlaplace,lwd=1.9)
qqplot(MC[,3],real,pch=20,col='red',main=paste("sigma=",sigma[3]));
qqline(MC[,3],distribution=qlaplace,lwd=1.9)
qqplot(MC[,4],real,pch=20,col='red',main=paste("sigma=",sigma[4]));
qqline(MC[,4],distribution=qlaplace,lwd=1.9)


## -----------------------------------------------------------------------------
a=3
log(exp(a))==exp(log(a));isTRUE(all.equal(log(exp(a)),exp(log(a))))

a=5
log(exp(a))==exp(log(a));isTRUE(all.equal(log(exp(a)),exp(log(a))))

## -----------------------------------------------------------------------------
f = function(k,a){
 ck=sqrt(a^2*k/(k+1-a^2))
 int=integrate(function(u){(1+u^2/k)^(-(k+1)/2)},0,ck)$value
 mul=2/sqrt(pi*k)*exp(lgamma((k+1)/2)-lgamma(k/2))
 mul*int
}
# f = function(k,a){
#   ck=sqrt(a^2*k/(k+1-a^2))
#   pt(ck,df=k)
# }
solve_eq = function(k){
  out=uniroot(function(a){f(k,a)-f(k-1,a)},lower=1,upper=2)
  out$root
}
Root=matrix(0,2,5)
k=c(4,25,100,500,1000)
for (i in 1:length(k)){
  Root[2,i]=round(solve_eq(k[i]),4)
}
Root[1,]=k;rownames(Root)=c('k','root')
Root

## -----------------------------------------------------------------------------
Sk = function(k,a){
 ck=sqrt(a^2*k/(k+1-a^2))
 pt(ck,df=k,lower.tail=FALSE)
}

solve_ist = function(k){
  o=uniroot(function(a){Sk(k,a)-Sk(k-1,a)},lower=1,upper=2)
  o$root
}
Ist=matrix(0,2,5)
k=c(4,25,100,500,1000)
for (i in 1:length(k)){
  Ist[2,i]=round(solve_ist(k[i]),4)
}
Ist[1,]=k;rownames(Ist)=c('k','ist')
Ist

## -----------------------------------------------------------------------------
set.seed(5678)
log_like = function(p,q,r){
 p_t=p^2/(p^2+2*p*r);q_t=q^2/(q^2+2*q*r)
 l=n_A*log(2*p*r)+n_B*log(2*q*r)+2*n_OO*log(r)+n_AB*log(2*p*q)+
   p_t*n_A*log(p/2/r)+q_t*n_B*log(q/2/r)
 l
}
N = 10000 
B = c(.5,.4,.1) 
tol = .Machine$double.eps^0.5
n_A=28;n_B=24;n_OO=41;n_AB=70
n=n_A+n_B+n_OO+n_AB
B.new=B;Log_l=log_like(B[1],B[2],B[3])
for (j in 1:N) {
 p_tmp=B[1]^2/(B[1]^2+2*B[1]*B[3])
 B[1]=(n_AB+n_A+n_A*p_tmp)/n/2
 q_tmp=B[2]^2/(B[2]^2+2*B[2]*B[3])
 B[2]=(n_AB+n_B+n_B*q_tmp)/n/2
 B[3]=1-B[1]-B[2]
 if (sum(abs(B-B.new)/B.new)<tol) break
 B.new = B;Log_l=c(Log_l,log_like(B[1],B[2],B[3]))
}
print(list(pqr = B, iter = j, tol = tol))

## -----------------------------------------------------------------------------
plot(1:j,Log_l,type='l',xlab='iteration',ylab='log-likelihood')

## -----------------------------------------------------------------------------
attach(mtcars)
formulas=list(
 mpg ~ disp,
 mpg ~ I(1 / disp),
 mpg ~ disp + wt,
 mpg ~ I(1 / disp) + wt
)
#Loop
n=length(formulas)
out=vector("list",n)
for (i in 1:n){
 out[[i]]=lm(formulas[[i]],data=mtcars)
}
out
#lapply
lapply(formulas,lm)

## -----------------------------------------------------------------------------
bootstraps = lapply(1:10, function(i) {
 rows = sample(1:nrow(mtcars), rep = TRUE)
 mtcars[rows, ]
})
oUt=vector("list",10)
for (i in 1:10) oUt[[i]]=lm(mpg ~ disp,data=bootstraps[[i]])
oUt[1:3]#Only print results 1-3  
Out=lapply(seq_along(bootstraps), function(i){
 lm(mpg ~ disp,data=bootstraps[[i]])})
Out[1:3]#Only print resuls 1-3
OUT=lapply(bootstraps,lm,formula = mpg ~ disp)#No anonymous function
OUT[1:3]

## -----------------------------------------------------------------------------
rsq = function(mod) summary(mod)$r.squared
unlist(lapply(out,rsq)) #Model in Ex3
unlist(lapply(Out,rsq)) #Model in Ex4


## -----------------------------------------------------------------------------
trials=replicate(100,t.test(rpois(10, 10),rpois(7, 10)),simplify = FALSE)
sapply(trials,function(mod){mod$p.value})
sapply(trials, "[[", "p.value") #No anonymous function

## ----eval=FALSE---------------------------------------------------------------
#  library(parallel)
#  cores = detectCores()
#  cluster = makePSOCKcluster(cores)
#  mcsapply = function(cluster,X,FUN,...){
#   res=parLapply(cluster,X,FUN,...) #Use parLapply in Windows
#   simplify2array(res)
#  }
#  #Example of Parallelisation
#  boot_i = function(i){
#   r_adj = function(x,id) {
#    x = x[id,]     #bootstrap sample
#    res = lm(mpg~wt+disp,data=x)
#    summary(res)$adj.r.squared
#  }
#   as.numeric(boot::boot(mtcars,r_adj,1)$t)
#  }
#  system.time(sapply(1:500,boot_i))
#  system.time(mcsapply(cluster,1:500,boot_i))
#  #In fact, in Windows, we have function *parSapply* directly.

## -----------------------------------------------------------------------------
library(Rcpp)
set.seed(1278)
sourceCpp('C:/Users/H.X.X/Desktop/Statistical Computing/A-19090-2019-12-20/randomMHC.cpp')
N=2000;sigma=2;x0=20
C=randomMHC(sigma,x0,N) #Samples generated by Cpp function
paste('Rej=',C$rej)
plot(1:N,C$x,type='l',xlab='N',ylab='Samples',main='RW Metropolis sampler for generating Laplace dist')

## -----------------------------------------------------------------------------
set.seed(1278)
f = function(x) {
  return(exp(-abs(x)))
}
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
R=random.M(sigma,x0,N) #Samples generated by R function
library(rmutil,quietly=TRUE)
p = c(.05, seq(.1, .9, .1), .95)
real=rlaplace(N) #Real examples of Laplace distribution 
Ex = cbind(C$x, R$x)
Mc = Ex[501:N, ] #burn in the first 500 samples
QR = apply(Mc, 2, function(x) quantile(x, p))
# par(mfrow=c(1,3))
# qqplot(Mc[,1],real,pch=20,col='red',main='Samples generated by Cpp');
# qqline(Mc[,1],distribution=qlaplace,lwd=1.9)
# qqplot(Mc[,2],real,pch=20,col='red',main='Samples generated by R');
# qqline(Mc[,2],distribution=qlaplace,lwd=1.9)
qqplot(Mc[,1],Mc[,2],xlab='Samples generated by Cpp',ylab='Samples generated by R',col='red',pch=20)
abline(a=0,b=1,col='black')

## -----------------------------------------------------------------------------
library(microbenchmark)
time = microbenchmark(MHC=randomMHC(sigma,x0,N),MHR=random.M(sigma,x0,N))
summary(time)[,c(1,3,5,6)]

