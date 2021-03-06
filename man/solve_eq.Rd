% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ex.R
\name{solve_eq}
\alias{solve_eq}
\title{Solve the equation ff(k,a)=ff(k-1,a)}
\usage{
solve_eq(k)
}
\arguments{
\item{k}{integer}
}
\value{
the root of the equation
}
\description{
Solve the equation ff(k,a)=ff(k-1,a)
}
\examples{
\dontrun{
Root=matrix(0,2,5)
k=c(4,25,100,500,1000)
for (i in 1:length(k)){
 Root[2,i]=round(solve_eq(k[i]),4)
}
Root[1,]=k;rownames(Root)=c('k','root')
Root
}
}
