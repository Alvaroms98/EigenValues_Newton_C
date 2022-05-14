#include "newton.h"


newtonResults newtonLU(double *x0, double *A, double *J,
                        double *F, const double tolr, const double tolA,
                        const unsigned int maxIt, const unsigned int n){
    
    const unsigned int N = n + 1;

    // Function evaluation 
    system_func(F, x0, A, n);

    // initial error
    double r0 = normVec(F, N);
    double err = 1;

    // Vector allocation for pivot in A = P * L * U decomposition
    int *P = (int*) malloc( N * sizeof(int) );

    // Vector allocation of x_{n+1}
    double *x = (double*) malloc( N * sizeof(double) );

    unsigned int iter = 0;
    while ( (iter < maxIt) && (err > tolr*r0 + tolA) ) {
        iter++;

        // Compute Jacobian matrix
        finDiffJac(F, x0, A, J, N);

        // Solve linear system A * x = b (solution is saved in b)
        // Thats because we need to copy F -> x
        memcpy(x, F, N * sizeof(double) );

        // Change the value of F -> F'(x) * y = -F(x)
        for (int i = 0; i < N; i++) {
            x[i] = -1.0 * x[i];
        }
        LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, 1, J, N, P, x, 1);

        // At this point x = x_{n+1} - x_{n} -> x_{n+1} = x + x_{n}
        for (int i = 0; i < N; i++) {
            x[i] += x0[i];
        }

        // function evaluation for the next iteration
        system_func(F, x, A, n);

        // Error calculation
        err = normVec(F, N);

        // Update x0
        memcpy(x0, x, N * sizeof(double));
    }

    free(P);
    free(x);

    newtonResults res;
    res.iter = iter;
    res.error = err;
    res.flag = (iter == maxIt) ? 1 : 0;

    return res;
}


double normVec(double *F, const unsigned int N) {
    double sol = 0.0;
    for (int i = 0; i < N; i++) {
        sol += F[i] * F[i];
    }
    sol = sqrt(sol);
    return sol;
}

void system_func(double *F, const double *x, const double *A, const unsigned int n) {
    const unsigned int N = n + 1;

    for (int i = 0; i < N; i++) {
        double tmp = 0.0;

        if (i < n) {
            for (int j = 0; j < n; j++) {
                tmp += A[i*n + j] * x[j];
            }
            F[i] = tmp - x[N-1] * x[i];

        } else {
            for (int j = 0; j < n; j++) {
                tmp += x[j] * x[j];
            }
            F[i] = sqrt(tmp) - 1;
        }
    }
}

void finDiffJac(const double *F0, const double *x0, const double *A, double *J,  const unsigned int N) {
    double h = 1e-7;

    // Allocate memory for system evaluation
    double *F = (double*) malloc( N * sizeof(double) ); 
    double *x = (double*) malloc( N * sizeof(double) );

    
    // Columns of the jacobian
    for (int j = 0; j < N; j++) {

        memcpy(x, x0, N*sizeof(double));
        x[j] += h;
        system_func(F, x, A, N-1);

        // rows
        for (int i = 0; i < N; i++){
            J[i*N + j] = (F[i] - F0[i]) / h;
        }
    }
    
    free(F);
    free(x);
}