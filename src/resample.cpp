// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;


double ess(NumericVector w)
{
  NumericVector w2 = w*w;
  return 1/std::accumulate(w2.begin(), w2.end(), 0.0);
}



// [[Rcpp::export]]
IntegerVector multinomial_resample(NumericVector w, int n)
{
  IntegerVector ii = seq_len(w.length());
  return RcppArmadillo::sample(ii, n, TRUE, w);
}
