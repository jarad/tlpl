#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rmath.h>

#include "Sckm.h"
#include "SckmParticle.h"
#include "SckmSwarm.h"



SckmSwarm **newSckmSwarms(Sckm *sckm, int _nParticles, int _nObs,
                        int *_state,
                        double *_probA, double *_probB, double *_rateA, double *_rateB)
{
    SckmSwarm **swarms;
    swarms = (SckmSwarm **) malloc((_nObs+1) * sizeof(SckmSwarm *));

    if (swarms == NULL) 
    {
      error("Memory for swarms could not be allocated.\n");
    }

    int i, nRxnOffset = _nParticles * sckm->r, nSpeciesOffset = _nParticles * sckm->s;

    for (i=0; i<=_nObs; i++) 
    {
        swarms[i] = newSckmSwarm(sckm, _nParticles, _state, _probA, _probB, _rateA, _rateB);

        _state += nSpeciesOffset;
        _probA += nRxnOffset;
        _probB += nRxnOffset;
        _rateA += nRxnOffset;
        _rateB += nRxnOffset;

    }

    return swarms;
}


void deleteSckmSwarms(SckmSwarm **swarms, int nObs)
{
    for (int i=0; i<=nObs; i++) deleteSckmSwarm(swarms[i]);
    free(swarms); swarms = NULL;
}


SckmSwarm *newSckmSwarm(Sckm *sckm, int _nParticles,
                        int *_state,
                        double *_probA, double *_probB, double *_rateA, double *_rateB)
{
    SckmSwarm *swarm;
    swarm = (SckmSwarm *) malloc(sizeof(SckmSwarm)); 

    if (swarm == NULL) 
    {
      error("Memory for swarm could not be allocated.\n");
    }

    swarm->nParticles = _nParticles;
    swarm->nStates    = sckm->s;
    swarm->nRxns      = sckm->r;

    // Allocate dWeights and make uniform
    swarm->adWeights = (double *) malloc(_nParticles * sizeof(double));

    if (swarm->adWeights == NULL) 
    {
      error("Memory for swarm->adWeights could not be allocated.\n");
    } 

    memset(swarm->adWeights, 0, _nParticles);

    swarm->logWeights        = 1; // log weights
    swarm->normalizedWeights = 0; // weights are unnormalized

    // Associate particle pointers
    swarm->pParticle = (SckmParticle **) malloc(_nParticles * sizeof(SckmParticle *));

    if (swarm->pParticle == NULL) 
    {
      error("Memory for swarm->pParticle could not be allocated.\n");
    }


    for (int i=0; i< _nParticles; i++) 
    {
        swarm->pParticle[i] = (SckmParticle *) malloc(sizeof(SckmParticle));

        if (swarm->pParticle[i] == NULL) 
        {
          error("Memory for swarm->pParticle[i] could not be allocated.\n");
        }

        swarm->pParticle[i]->state = _state; _state += sckm->s;
        swarm->pParticle[i]->probA = _probA; _probA += sckm->r;
        swarm->pParticle[i]->probB = _probB; _probB += sckm->r;
        swarm->pParticle[i]->rateA = _rateA; _rateA += sckm->r;
        swarm->pParticle[i]->rateB = _rateB; _rateB += sckm->r;
    }

    return swarm;
}

void deleteSckmSwarm(SckmSwarm *swarm)
{
    free(swarm->adWeights); swarm->adWeights = NULL;
    for (int i=0; i< swarm->nParticles; i++)
    {
      free(swarm->pParticle[i]); swarm->pParticle[i] = NULL;
    }

    free(swarm->pParticle); swarm->pParticle = NULL;
    free(swarm); swarm = NULL;
}

int renormalize(SckmSwarm *swarm)
{
    if (swarm->normalizedWeights) return 0;

    int i, n=swarm->nParticles;
    double *w;
    w = swarm->adWeights; 

    if (swarm->logWeights) 
    {
        double max=w[0];
        for (i=1; i<n; i++) max = fmax2(max, w[i]);
        for (i=0; i<n; i++) w[i] = exp(w[i]-max);
    }
    swarm -> logWeights = 0;

    double sum=0;
    for (i=0; i<n; i++) sum += w[i];
    for (i=0; i<n; i++) w[i] /= sum;

    return 0;
}

