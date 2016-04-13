#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;


double ess(NumericVector w)
{
  NumericVector w2 = w*w;
  return 1/std::accumulate(w2.begin(), w2.end(), 0.0);
}



// [[Rcpp::export]]
IntegerVector multinomial_resample(NumericVector w, int nI)
{
  return RcppArmadillo::sample(seq_len(nI), nI, TRUE, w);
}
