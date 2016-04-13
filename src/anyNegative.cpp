#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
int anyNegative(IntegerVector v)                   
{
    for (int i=0; i<v.size(); i++) 
    {
        if (v[i]<0) return 1;
    }
    return 0;
}

