#include <Rcpp.h>

using namespace Rcpp;

class Sckm
{
    IntegerVector X;
    IntegerMatrix Pre, Post, Stoichiometry;
    NumericVector theta, lmult, hazard_part, hazard;
    
  public:
    int n_reactions()    { return Stoichiometry.ncol(); }
    int n_species()      { return Stoichiometry.nrow(); }
    IntegerVector X_()   { return X;                    }
    IntegerMatrix Pre_() { return Pre;                  }
    IntegerMatrix S_()   { return Stoichiometry;        }
  
//    Sckm(IntegerMatrix, IntegerMatrix, IntegerVector, NumericVector);
    Sckm(List);
    
    void printX();
    
    void check();
    
    bool any_negative() { return is_true(any(X<0)); }
    
    // Simulation functions
    void update(int);
    void update(IntegerVector);
    double gillespie_step();
};

Sckm::Sckm(List sckm) {
  SEXP X_             = sckm["X"],
       Pre_           = sckm["Pre"], 
       Post_          = sckm["Post"], 
       theta_         = sckm["theta"],
       mult_          = sckm["mult"];
       
  X             = clone(X_);
  Pre           = clone(Pre_);
  Post          = clone(Post_);
  theta         = clone(theta_);
  lmult         = log(mult_);
  
  // Construct stoichiometry matrix
  int nr = Pre.ncol(), ns = Pre.nrow();
  IntegerMatrix S(nr, ns);
  for (int r=0; r<nr; r++) 
    for (int s=0; s<ns; s++) 
      S(s,r) = Post(r,s) - Pre(r,s);
  Stoichiometry = clone(S);
  
  // Set hazard_part and hazard
  NumericVector HP(nr), H(nr);
  for (int r=0; r<nr; r++) 
    for (int s=0; s<ns; s++)
      HP[r] += Rf_lchoose(X[s], Pre(r,s));
      
  hazard_part = exp(HP);
  hazard      = hazard_part*theta;
}

//Sckm::Sckm(IntegerMatrix Pre, IntegerMatrix Post, IntegerVector X, NumericVector mult) {
//  int n_reactions = Pre.nrow(), 
//      n_species   = Pre.ncol();
//  
//  lmult_ = log(mult);
//  
//  X_     = clone(X);
//  Pre_   = clone(Pre);
//  Post_  = clone(Post);
//}

void Sckm::update(int reaction_id) {
  X = X + Stoichiometry(_,reaction_id);
}

void Sckm::update(IntegerVector reaction_count) {
  for (int r = 0; r < n_reactions(); r++)
    X = X + Stoichiometry(_,r) * reaction_count[r];
}

void Sckm::printX() {
  Rcout << X << std::endl;
}

////////////////////////////////////////////////////////////////////
// Simulation functions
////////////////////////////////////////////////////////////////////

// Run one Gillespie step
double Sckm::gillespie_step() {
  int nr = n_reactions();
  NumericVector hazard_cumsum = cumsum(hazard);
  
  // Determine which reaction fired
  double total_hazard = hazard_cumsum[nr-1];
  double tmp = Rf_runif(0,1) * total_hazard;
  X = X + Stoichiometry(_, sum(tmp<hazard_cumsum)+1);
  
  // report time increment
  return exp(total_hazard);
}
