#include <cmath>
#include <Rcpp.h>
using namespace Rcpp;

//[[Rcpp::export]]
double f(double x) {
	return exp(-fabs(x));
}

//' @title A Metropolis-Hastings sampler using Rcpp
//' @description A Metropolis-Hastings sampler using Rcpp
//' @param sigma the standard deviation of the sample-generating normal distribution
//' @param x0 initial values
//' @param N length of Markov Chain
//' @return A Markov China of generating samples
//' @useDynLib SC19090
//' @examples
//' \dontrun{
//' R= randomMHC(2,10,2000)
//' }
//' @export
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