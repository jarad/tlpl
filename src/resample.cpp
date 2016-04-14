// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;



#define STRATIFIED_RESAMPLING 1
#define MULTINOMIAL_RESAMPLING 2
#define SYSTEMATIC_RESAMPLING 3
#define RESIDUAL_RESAMPLING_THEN_STRATIFIED 4
#define RESIDUAL_RESAMPLING_THEN_MULTINOMIAL 5
#define RESIDUAL_RESAMPLING_THEN_SYSTEMATIC 6

double ess(NumericVector w)
{
  NumericVector w2 = w*w;
  return 1/std::accumulate(w2.begin(), w2.end(), 0.0);
}


IntegerVector multinomial_resample(NumericVector w, int n)
{
  IntegerVector ii = seq_len(w.length());
  return RcppArmadillo::sample(ii, n, TRUE, w);
}


// [[Rcpp::export]]
IntegerVector resample(NumericVector w, int n, int method)
{
  switch(method)
  {
//    case STRATIFIED_RESAMPLING:
//      stratified_resample(nW, adWeights, nI, anIndices);
//      break;
    case MULTINOMIAL_RESAMPLING:
      return multinomial_resample(w, n);
      break;
//    case SYSTEMATIC_RESAMPLING:
//      systematic_resample(nW, adWeights, nI, anIndices);
//      break;
//    case RESIDUAL_RESAMPLING_THEN_STRATIFIED:
//      residual_resample(nW, adWeights, nI, anIndices, STRATIFIED_RESAMPLING);
//      break;
//    case RESIDUAL_RESAMPLING_THEN_MULTINOMIAL: 
//      residual_resample(nW, adWeights, nI, anIndices, MULTINOMIAL_RESAMPLING);
//      break;
//    case RESIDUAL_RESAMPLING_THEN_SYSTEMATIC : 
//      residual_resample(nW, adWeights, nI, anIndices, SYSTEMATIC_RESAMPLING);
//      break;
    default:
      throw std::range_error("resample: no match for resampling function");
  }
}




