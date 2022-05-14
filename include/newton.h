#ifndef NEWTON

#define NEWTON

#include <lapacke.h>
#include <math.h>
#include <string.h>

#endif

typedef struct {

    unsigned int iter;
    double error;
    int flag;

} newtonResults;

double normVec(double *F, const unsigned int N);
void system_func(double *F, const double *x, const double *A, const unsigned int n);

newtonResults newtonLU(double *x0, double *A, double *J,
                        double *F, const double tolr, const double tolA,
                        const unsigned int maxIt, const unsigned int n);

void finDiffJac(const double *F0, const double *x0, const double *A, 
                double *J,  const unsigned int N);
