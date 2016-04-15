#include <Sckm.h>

using namespace Rcpp;

// [[Rcpp::export]]
void test(List sckm_list) {
  
  Sckm sckm(sckm_list);
  
  int nr = sckm.n_reactions(); 
  IntegerVector reaction_count = as<IntegerVector>(rpois(nr,1));
  
  sckm.printX();
  sckm.update(reaction_count);
  sckm.printX();
}
