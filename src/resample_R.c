/*****************************************************************/
/* C-side interface to R                                         */
/*****************************************************************/

#include "resample.h"
#include "resample_R.h"


void rep2id_R(int *rep, int *sum, int *id) 
{
  rep2id(rep, *sum, id);
}


void is_increasing_R(int *n, const double *v, int *returned) 
{
  *returned = is_increasing(*n, v); 
}


void cumulative_sum_R(int *n, double *v) 
{
  cumulative_sum(*n, v);
}


void inverse_cdf_weights_R(int *nW, double *adWeights, int *nU, double *adUniforms, int *anIndices)
{
  inverse_cdf_weights(*nW, adWeights, *nU, adUniforms, anIndices);
}


void ess_R(int *n, double *weights, double *returned)
{
  *returned = ess(*n, weights);
}


void cov2_R(int *n, double *weights, double *returned)
{
  *returned = cov2(*n, weights);
}


void entropy_R(int *n, double *weights, double *returned) {
  *returned = entropy(*n, weights);
}


void resample_R(int *nW, double *adWeights, int *nI, int *anIndices,
                   int *nResamplingFunction)
{
  resample(*nW, adWeights, *nI, anIndices, *nResamplingFunction);
}


void stratified_resample_R(int *nW, double *adWeights, int *nI, int *anIndices)
{
  stratified_resample(*nW, adWeights, *nI, anIndices);
}


void multinomial_resample_R( int *nW, double *adWeights, int *nI, int *anIndices)
{
  multinomial_resample(*nW, adWeights, *nI, anIndices);
}


void systematic_resample_R(int *nW, double *adWeights, int *nI, int *anIndices)
{
  systematic_resample(*nW, adWeights, *nI, anIndices);
}

void residual_resample_R(int *nW, double *adWeights, int *nI, int *anIndices, 
                            int *nResidualResampleFunction)
{
  residual_resample(*nW, adWeights, *nI, anIndices, *nResidualResampleFunction);
}



