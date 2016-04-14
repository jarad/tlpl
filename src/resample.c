#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include "utility.h"
#include "resample.h"

#define NONE 1
#define ESS 2
#define COV 3
#define ENTROPY 4

#define NOT_SORTED 0
#define SORTED 1

#define NO_FXN_ERROR 0
#define FXN_ERROR 1


/***********************************************************************/
/* Utility functions                                                   */
/***********************************************************************/


// used in qsort and stolen from http://en.allexperts.com/q/C-1587/Qsort-function.htm
int compare_doubles (const void *X, const void *Y)
{
  double x = *((double *)X);
  double y = *((double *)Y);

  if (x > y)
  {
    return 1;
  }
  else
  {
    if (x < y)
    {
      return -1;
    }
    else
    {
      return 0;
    }
  }
}


int rep2id(int *rep, int sum, int *id)
{
  // This implementation seems poor. 
  // No error checking to assure we stay within the bounds of rep and id
  int i, j=0;

  i=0;
  while (i<sum) 
  {
    if (rep[j]>0) // If this particle is resampled (again)
    {
      rep[j]--; 
      id[i] = j;
      i++;
    } 
    else          // If all resamples of this particle are exhausted
    {
      j++;
    }
  }

  return NO_FXN_ERROR;
}



int inverse_cdf_weights(int nW, 
                        const double *adWeights, 
                        int nU, 
                        double *adUniforms,
                        int *anIndices,
                        int isSorted)
{
  if (!isSorted) 
    qsort(adUniforms, nU, sizeof(double), compare_doubles);

  double * adCumSum = malloc(nW * sizeof(double));
  if (adCumSum == NULL) { error("C: adCumSum not dynamically allocated"); } 

  memcpy(adCumSum, adWeights, nW * sizeof(double));

  cumulative_sum(nW, adCumSum);

  int i, j=0, found;
  for (i=0; i<nU; i++) 
  {
    found=0;
    while (!found) 
    {
      if (adUniforms[i] > adCumSum[j])
      {
        j++;
      }
      else 
      {
        found=1;
      }
    }
    anIndices[i] = j;
  } 

  free(adCumSum); adCumSum = NULL;

  return NO_FXN_ERROR;   
}


/***********************************************************************/
/* Particle weight nonuniformity                                       */
/***********************************************************************/

double cov2(int n, double *weights) 
{
  int i;
  double mean=0, var=0, tmp;

  // Calculate mean
  for (i=0; i<n; i++) mean += weights[i];
  mean /= n;

  // Calculate variance
  for (i=0; i<n; i++) 
  {
    tmp = weights[i]-mean;
    var += tmp*tmp;
  }
  var /= (n-1);
   
  // Return cov^2
  return var/(mean*mean);
}





int doResample(int n, double *weights, int nNonuniformity, double dThreshold)
{
  switch(nNonuniformity)
  {
    case NONE: // "none" means always resample
      return 1;
    case ESS: 
      return ess(    n, weights) < dThreshold ? 1 : 0;
    case COV: 
      return cov2(   n, weights) > dThreshold ? 1 : 0; // notice the greater than sign
    case ENTROPY: 
      return entropy(n, weights) < dThreshold ? 1 : 0;
    default:
      error("C: doResample: Nonuniformity measure not found.\n");
  }
  error("doResample exited switch without a proper response");
  return FXN_ERROR;
}


/***********************************************************************/
/* Resampling functions                                                */
/***********************************************************************/

int resample(int nW, double *adWeights, int nI, int *anIndices, 
             int nResamplingFunction)
{
  switch(nResamplingFunction)
  {
    case STRATIFIED_RESAMPLING:
      stratified_resample(nW, adWeights, nI, anIndices);
      break;
    case MULTINOMIAL_RESAMPLING:
      multinomial_resample(nW, adWeights, nI, anIndices);
      break;
    case SYSTEMATIC_RESAMPLING:
      systematic_resample(nW, adWeights, nI, anIndices);
      break;
    case RESIDUAL_RESAMPLING_THEN_STRATIFIED:
      residual_resample(nW, adWeights, nI, anIndices, STRATIFIED_RESAMPLING);
      break;
    case RESIDUAL_RESAMPLING_THEN_MULTINOMIAL: 
      residual_resample(nW, adWeights, nI, anIndices, MULTINOMIAL_RESAMPLING);
      break;
    case RESIDUAL_RESAMPLING_THEN_SYSTEMATIC : 
      residual_resample(nW, adWeights, nI, anIndices, SYSTEMATIC_RESAMPLING);
      break;
    default:
      REprintf("C: resample: no match for resampling function\n");
  }

  for (int i=0; i<nW; i++) adWeights[i] = 1.0 / nW;

  return NO_FXN_ERROR;
}



int stratified_resample(int nW, double *adWeights, int nI, int *anIndices)
{
  int i;
  double adUniforms[nI];

  GetRNGstate();
  for (i=0;i<nI;i++) adUniforms[i] = runif((double) i/nI, (double) (i+1)/nI);
  PutRNGstate();

  inverse_cdf_weights(nW, adWeights, nI, adUniforms, anIndices, SORTED);

  return NO_FXN_ERROR;
}




int multinomial_resample(int nW, double *adWeights, int nI, int *anIndices) 
{
  int i;
  double adUniforms[nI];

  GetRNGstate();
  for (i=0; i<nI; i++) adUniforms[i] = runif(0,1);
  PutRNGstate();

  inverse_cdf_weights(nW, adWeights, nI, adUniforms, anIndices, NOT_SORTED);

  return NO_FXN_ERROR;
}







int systematic_resample(int nW, double *adWeights, int nI, int *anIndices)
{
  int i;
  double adUniforms[nI];
  GetRNGstate();
  adUniforms[0] = runif(0, (float) 1.0 / nI);
  PutRNGstate();
  for (i=1; i<nI; i++) adUniforms[i] =  adUniforms[i-1] + (float) 1.0 / nI;

  inverse_cdf_weights(nW, adWeights, nI, adUniforms, anIndices, SORTED);

  return NO_FXN_ERROR;
}



int residual_resample(int nW, double *adWeights, int nI, int *anIndices,
                       int nResidualResampleFunction)
{
  // Particles are deterministically resampled floor(weights*nSamples) times
  int i, anDeterministicReps[nW], nDeterministicReps=0;
  double adExpectedSamples[nW];
  for (i=0; i<nW; i++) 
  {        
    // Expected samples 
    adExpectedSamples[i]    = adWeights[i]* nI;

    // Truncate to get deterministically resampled particles
    anDeterministicReps[i]  = adExpectedSamples[i];
     
    // Increment number of deterministic reps
    nDeterministicReps     += anDeterministicReps[i];

    // Remaining weight for use in random resampling
    adWeights[i]            = adExpectedSamples[i]-anDeterministicReps[i];
  }
  if (nDeterministicReps > nI) 
    REprintf("C: residual_resample: too many deterministic reps\n");
   
  rep2id(anDeterministicReps, nDeterministicReps, anIndices);


  // Particles are then randomly sampled with remaining weight
  nI -= nDeterministicReps;

  // Renormalize weights
  double sum=0;
  for (i=0; i<nW; i++) sum += adWeights[i];
  for (i=0; i<nW; i++) adWeights[i] /= sum;

  switch (nResidualResampleFunction) 
  {
    case STRATIFIED_RESAMPLING:
      stratified_resample( nW, adWeights, nI, &anIndices[nDeterministicReps]);
      break;
    case MULTINOMIAL_RESAMPLING:
      multinomial_resample(nW, adWeights, nI, &anIndices[nDeterministicReps]);
      break;
    case SYSTEMATIC_RESAMPLING:
      systematic_resample( nW, adWeights, nI, &anIndices[nDeterministicReps]);
      break;
    default:
      REprintf("C: residual_resample: no match for residual resampling function\n");
  }
       
  return NO_FXN_ERROR;
}


