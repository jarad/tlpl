#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector hazard_part(List sckm, IntegerVector const X)
{
  int nr = sckm["r"],
      ns = sckm["s"];  
  SEXP SEXP_lmult = sckm["lmult"],
       SEXP_Pre   = sckm["Pre"];
  IntegerMatrix Pre(SEXP_Pre);
  NumericVector hazard_part(clone(SEXP_lmult));
  Rcout << hazard_part << std::endl;
  
  for (int r=0; r<nr; r++) 
    for (int s=0; s<ns; s++)
      hazard_part[r] += Rf_lchoose(X[s], Pre(r,s)); 
      
  return exp(hazard_part);
}

//NumericVector hazard(NumericVector const theta, double tau,
//                     NumericVector hazard_part)
//{
//  return NumericVector = theta*hazard_part*tau;
//}


// [[Rcpp::export]]
void update_species(List sckm, IntegerVector const rxn_count)
{
  SEXP X_ = sckm["X"], S_ = sckm["stoich"];
  IntegerVector X(X_);
  IntegerMatrix S(S_);
  
  Rcout << X << "\n" << S << rxn_count << "\n";
  
  for (int r = 0; r < rxn_count.length(); r++)
    X = X + S(_,r) * rxn_count[r];
  
  Rcout << X << std::endl;
}




// --------------------------------------------------------------------------------------
// Discrete time
// --------------------------------------------------------------------------------------

//IntegerVector tau_leap_one_step(List sckm, NumericVector const hazard) {
//  int nr = sckm["r"];
//  
//  // Simulate number of reactions
//  IntegerVector number_of_reactions(nr);
//  for (int r=0; r<nr; r++) 
//    number_of_reactions[r] = rpois(1,hazard[r]);
//    
//  update_species
//}


//int tau_leap_one_step(Sckm *sckm, 
//                  const double *adHazard,                 
//                  int nWhileMax,                          
//                  int *anRxnCount, int *anX) // return: number of reactions and updated species
//{
//  int nSpecies=sckm->s;
//  int i, whileCount=0, anTempX[nSpecies];
//  while (1) 
//  {
//    memcpy(anTempX, anX, nSpecies*sizeof(int));
//
//    // Get number of reactions
//    GetRNGstate();
//    for (i=0; i<sckm->r; i++) anRxnCount[i] = rpois(adHazard[i]);
//    PutRNGstate();
//
//    update_species(sckm, anRxnCount, anTempX);
//
//    if (!anyNegative(nSpecies, anTempX)) 
//    {
//      memcpy(anX, anTempX, nSpecies*sizeof(int));
//      return 0;
//    }
//
//    // Limit how long the simulation tries to find a non-negative update
//    whileCount++;
//    if (whileCount>nWhileMax)  
//    {
//      error("C:tau_leap_one_step: Too many unsuccessful simulation iterations.");
//      return 1;
//    }
//  }
//} 

//int tau_leap(Sckm *sckm, const double *adTheta,
//         const double *adTau, int nSteps, 
//         int nWhileMax,
//         int *anRxnCount, int *anX)
//{
//  int nRxns = sckm->r, nSpecies = sckm->s;
//  int i, *ipLast, *ipCurrent;
//  ipLast     = anX;             // Points to last state
//  ipCurrent  = ipLast+nSpecies; // Points to current state
//
//  double adHazardPart[nRxns], adHazard[nRxns];
//  for (i=0; i<nSteps; i++)
//  {
//    memcpy(ipCurrent, ipLast, nSpecies*sizeof(int));
//    hazard(sckm, adTheta, ipCurrent, adTau[i], adHazardPart, adHazard);
//    tau_leap_one_step(sckm, adHazard, nWhileMax, &anRxnCount[i*nRxns], ipCurrent);
//    ipLast += nSpecies; ipCurrent += nSpecies;
//  }
//  return 0;
//}
//
//
//
//
//// --------------------------------------------------------------------------------------
//// Continuous time
//// --------------------------------------------------------------------------------------
//
//
//int next_to_fire(int nRxns, double *adCuSum) 
//{
//  int i, next=0;
//  double dUniform, dSum=adCuSum[nRxns-1]; 
//    
//  GetRNGstate();
//  dUniform = runif(0,1);
//  PutRNGstate();
//
//  for (i=0; i<nRxns; i++) 
//  {
//    adCuSum[i] /= dSum;
//    if (dUniform < adCuSum[i]) return next;
//    next++;
//  }   
//  // print error if this fails
//
//  error("next_to_fire failed to find a reaction\n");
//  return -1;
//}
//
//
//void gillespie_one_step_R(int *nSpecies, int *nRxns, int *anPre, int *anPost, double *adlMult,
//                             const double *adTheta, double *dT,  int *anRxnCount, int *anX)
//{
//  Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost, adlMult);
//  gillespie_one_step(sckm, adTheta, *dT, anRxnCount, anX);
//  deleteSckm(sckm);
//}
//
//
//int gillespie_one_step(Sckm *sckm, const double *adTheta, double dT, int *anRxnCount, int *anX)
//{
//  int nSpecies = sckm->s, nRxns = sckm->r;
//  int i, nRxnID, anX0[nSpecies];
//  double dCurrentTime=0, adHazardPart[nRxns], adHazard[nRxns];
//  while (1) 
//  {
//    memcpy(anX0, anX, nSpecies*sizeof(int));
//    hazard(sckm, adTheta, anX, 1, adHazardPart, adHazard);
//        
//    // Calculate cumulative hazard
//    for (i=1; i<nRxns; i++) adHazard[i] += adHazard[i-1];
//    if (adHazard[nRxns-1] < 0.0001) return 0;               // make this a function of dT?
//
//    dCurrentTime += rexp(1/adHazard[nRxns-1]);
//    if (dCurrentTime > dT) return 0;                        // stopping condition
//       
//    nRxnID = next_to_fire(nRxns, adHazard);
//    anRxnCount[nRxnID]++;
//
//    for (i=0; i<nSpecies; i++) anX[i] += sckm->Stoich[nSpecies * nRxnID + i]; 
//  }
//  return 0;               
//}
//
//
//
//void gillespie_R(int *nSpecies, int *nRxns, int *anPre, int *anPost, double *adlMult,
//               const double *adTheta,
//               double *adT, int *nSteps, int *anX)
//{
//  Sckm *sckm = newSckm(*nSpecies, *nRxns, anPre, anPost, adlMult);
//  gillespie(sckm, adTheta, adT, *nSteps, anX);
//  deleteSckm(sckm);
//}
//
//
//
//
//int gillespie(Sckm *sckm, const double *adTheta, double *adT, int nSteps, int *anX)
//{
//  int i, nSO=0, nr = sckm->r, ns = sckm->s;
//  int anRxnCount[nr];
//  for (i=0; i<nSteps;i++)
//  {
//    memcpy(&anX[nSO+ns], &anX[nSO], ns*sizeof(int));
//    nSO += ns;
//
//    gillespie_one_step(sckm, adTheta, adT[i], anRxnCount,  &anX[nSO]);
//  }
//  return 0;
//}


