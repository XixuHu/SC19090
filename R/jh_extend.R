#' 
#' @title Extended jonckheere function that can handle data with ties.
#' @param x the multi-group data (vector)
#' @param gr the group index for the data (numeric/factor)
#' 
#' @return a list that contains the value of J statistics and p.value
#' @details For Jonkheere-Terpstra test with ties
#' W[ij]=\#(X[ik]<X[jl],k=1,2,...,n[i],l=1,2,...,n[j])+frac(1,2)\#(X[ik]=X[jl],k=1,2,...,n[i],l=1,2,...,n[j])
#' J=sum(W[ij],i<j)
#' frac(J-E[H[0]](J),sqrt(var[H[0]](J)))
#' E[H[0]](J)=frac(N^2-sum(n[i]^2,i==1,k),4)
#' 
#'
#' @examples
#' \dontrun{
#'  data(m)
#'  attach(m)
#'  jh.tie(medical,gr)
#' }
#' @export
jh.tie = function(x,gr){
 n=nlevels(factor(gr))
 nn=as.numeric(table(gr))
 J=0
 for (i in 1:(n-1)){
  for (j in (i+1):n){
   for (k in 1:nn[i]){
    for (l in 1:nn[j]) {
     q=sum(nn[1:(j-1)])+l
     if(i==1) J=J+as.integer(x[k]<x[q])+0.5*as.integer(x[k]==x[q])
     else {
      p=sum(nn[1:(i-1)])+k
      J=J+as.integer(x[p]<x[q])+0.5*as.integer(x[p]==x[q])
     }}}}}
 N=sum(nn);tau=as.numeric(table(x));tau=tau[tau>1] #Tie Statistics
 EJ=(N^2-sum(nn^2))/4
 VarJ=(N*(N-1)*(2*N+5)-sum(nn*(nn-1)*(2*nn+5))-sum(tau*(tau-1)*(2*tau+5)))/72+
      sum(nn*(nn-1)*(nn-2))*sum(tau*(tau-1)*(tau-2))/36/N/(N-1)/(N-2)+
      sum(nn*(nn-1))*sum(tau*(tau-1))/8/N/(N-1)
 p_value=round(pnorm((J-EJ)/sqrt(VarJ),lower.tail=FALSE),6)
 cat("\n")
 cat("        Jonckheere-Terpstra test (data with ties)        \n")
 cat("\n")
 cat("data:\n")
 cat("JT = ",J,", p-value = ",p_value,'\n')
 return(list(JT=J,p.value=p_value))
}

#' 
#' @title Jonckheere function that can handle data both with ties and without ties.
#' @param x the multi-group data (vector)
#' @param gr the group index for the data (numeric/factor)
#' @import clinfun
#' @return a list that contains the value of J statistics and p.value
#' @export
#'
#' @examples
#' \dontrun{
#'  data(m)
#'  attach(m)
#'  jh_extend(medical,gr)
#' }
jh_extend = function(x,gr){
 if(length(unique(x))==length(x)) jonckheere.test(x,gr)
 else jh.tie(x,gr)
}