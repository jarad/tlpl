
int compare_doubles(const void *, const void *);
int is_increasing(int , const double *);
int cumulative_sum(int , double *);
int rep2id(int *, int , int *);
int inverse_cdf_weights(int , double *, int , double *, int *);

double ess(int , double *);
double cov2(int , double *);
double entropy(int , double *);

int doResample(int n, double *weights, int nNonuniformity, double dThreshold);
int resample(int , double *, int , int *, int );

int multinomial_resample(int, double *, int, int *);
int stratified_resample(int , double *, int , int *);
int systematic_resample(int , double *, int , int *);
int residual_resample(int , double *, int , int *, int );


