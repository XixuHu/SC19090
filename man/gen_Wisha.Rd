% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ex.R
\name{gen_Wisha}
\alias{gen_Wisha}
\title{Generate random samples from Wishart distribution W[d](Sigma,n)}
\usage{
gen_Wisha(d, n, Sigma)
}
\arguments{
\item{d}{the dimension of the square matrix}

\item{n}{degrees of freedom}

\item{Sigma}{scale matrix}
}
\value{
The generated samples from W[d](Sigma,n)
}
\description{
This function generate random samples from Wishart distribution based on Bartlett’s decomposition.
}
\details{
Let T=(T[ij]) be a lower triangular d%*%d random matrix with independent entries satisfying
T[ij]%~% N(0,1),i>j,i.i.d
T[ii]%~% sqrt(chi^2*(n-i+1)),i=1,...,d
Then the matrix A=TT^T has a W[d](I[d],n) distribution
Obtain the Choleski factorization Sigma=LL^T,where L is a lower triangular
Then LAL^T %~% W[d](Sigma,n)
}
\examples{
\dontrun{
 Sigma=matrix(c(3,2,2,4),nrow=2)
 gen_Wisha(2,5,Sigma)
}
}
