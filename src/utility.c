
#include <R.h>
#include <Rmath.h>
#include "utility.h"

/* Utility function to determine if any element in a vector is negative */
int anyNegative(const int n, const int *v)                   
{
    int i;
    for (i=0; i<n; i++) 
    {
        if (v[i]<0) return 1;
    }
    return 0;
}





