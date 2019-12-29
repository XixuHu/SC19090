#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double f(double x) {
	return exp(-fabs(x));
}

//[[Rcpp::export]]
List randomMHC (double sigma, double x0, int N) {
	NumericVector x(N);
	x[0] = x0;  //index start from 0
	NumericVector u = runif(N);
	int rej = 0;
		for (int i = 1; i < N;i++ ) {
			NumericVector y = rnorm(1, x[i-1], sigma);
			if (u[i] <= (f(y[0]) / f(x[i - 1])))
				{
					x[i] = y[0];
				}
			else 
				{ 
					x[i] = x[i-1]; 
					rej ++;
				}
		}
		return List::create(Named("x") = x, Named("rej") = rej);
}