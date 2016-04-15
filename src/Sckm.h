#include <Rcpp.h>
// should use RcppArmadillo but I can get this to work

using namespace Rcpp;

class Sckm
{
    IntegerVector X_;
    IntegerMatrix Pre_, Post_, Stoichiometry;
    NumericVector lmult_;
    
  public:
    int n_reactions()   { return Pre_.nrow(); }
    int n_species()     { return Pre_.ncol(); }
    IntegerVector X()   { return X_;                   }
    IntegerMatrix Pre() { return Pre_;                 }
    IntegerMatrix S()   { return Stoichiometry;        }
  
    Sckm(IntegerMatrix, IntegerMatrix, IntegerVector, NumericVector);
    Sckm(List);
    
    void update(IntegerVector);
    void print();
    
    bool any_negative() { return any(X<0); }
    
    void check();
};

Sckm::Sckm(List sckm) {
  SEXP Pre_   = sckm["Pre"], 
       Post_  = sckm["Post"], 
       X_     = sckm["X"],
       mult_  = sckm["mult"];
       
  IntegerVector X(X_);
  IntegerMatrix Pre(Pre_), Post(Post_);
  NumericVector mult(mult_);
  
  Sckm(Pre, Post, X, mult);
}

Sckm::Sckm(IntegerMatrix Pre, IntegerMatrix Post, IntegerVector X, NumericVector mult) {
  int n_reactions = Pre.nrow(), 
      n_species   = Pre.ncol();
  
  lmult_ = log(mult);
  
  X_     = clone(X);
  Pre_   = clone(Pre);
  Post_  = clone(Post);
  
  // Fill in S
  IntegerMatrix S(n_species, n_reactions);
  for (int r=0; r<n_reactions; r++) 
    for (int s=0; s<n_species; s++) 
      S(s,r) = Post_(r,s) - Pre_(r,s);
    
  Stoichiometry = S;
}

void Sckm::update(IntegerVector reaction_count) {
  for (int r = 0; r < n_reactions(); r++)
    X_ = X_ + Stoichiometry(_,r) * reaction_count[r];
}

void Sckm::print() {
  Rcout << Pre_ << std::endl;
}
