#include "utils.h"
#include "newton.h"

int main(int argc, char *argv[]) {

    char *inputMatrixFile;
    char *inputX0;
    int flagX0 = 0;

    if (argc < 2) {
        printf("Usage: %s <Matrix_input_file.csv> [X0_input_file.csv]\n", argv[0]);
        return EXIT_SUCCESS;
    } else if (argc == 2) {
        inputMatrixFile = argv[1];
    } else if (argc == 3) {
        inputMatrixFile = argv[1];
        inputX0 = argv[2];
        flagX0 = 1;
    } else {
        printf("Usage: %s <Matrix_input_file.csv> [X0_input_file.csv]\n", argv[0]);
        return EXIT_FAILURE;
    }
    
    // Set problem size by reading de number of lines of the .csv
    int n = getFileRows(inputMatrixFile);
    printf("Problem size: %dx%d\n",n,n);

    // Allocate memory for matrix A and x0
    double *A = (double*) malloc( n * n * sizeof(double) );
    double *x0 = (double*) calloc( (n+1), sizeof(double) );

    // Read matrix from file
    int flag = readMatrix(inputMatrixFile, A, n);
    if (flag > 0) {
        printf("Error reading the input matrix\n");
        return EXIT_FAILURE;
    }

    // Read initial values (if provided)
    if (flagX0) {
        int n2 = getFileRows(inputX0);
        if ( (n+1) != n2) {
            printf("Error: Matrix size and the length of the initial vector are different\n");
            return EXIT_FAILURE;
        } else {
            readX0(inputX0, x0, n2);
        }
    } else { // If x0 isn't provided, deafult value is 1.0 (0.0 cause divergence)
        for (int i = 0; i < (n+1); i++) {
            x0[i] = 1.0;
        }
    }

    // Allocate memory for Jacobian and System Function
    double *J = (double*) malloc( (n+1) * (n+1) * sizeof(double) );
    double *F = (double*) malloc( (n+1) * sizeof(double) );

    
    // Check everything is correct
    printf("\nMatrix A:\n");
    printMatrix(A, n);
    printf("\nInitial solution:\n");
    printVector(x0,n+1);

    /* CALCULATION OF NONLINEAR SYSTEM */
    double tolr = 1e-8;
    double tolA = 1e-12;
    unsigned int maxIt = 1000;

    newtonResults result = newtonLU(x0,A,J,F,tolr,tolA,maxIt,n);

    /* Print Results */
    printf("\n\t**** NONLINEAR SYSTEM RESULTS ****\n");
    if (!result.flag) {
        printf("\nResult has reached convergence\n");
        printf("Iterations needed: %d\tError: %g\n", result.iter, result.error);
    } else {
        printf("Result has not reached convergence\n");
    }
    printf("\nEigen Vector:\n");
    printVector(x0, n);
    printf("Eigen Value: %lf\n", x0[n]);


    /* LAPACK eigenvalue calculation */
    double *WR = (double*) malloc( n * sizeof(double));
    double *WL = (double*) malloc( n * sizeof(double));
    double *VR = (double*) malloc( n * n * sizeof(double));

    LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, A, n, WR, WL, NULL, n, VR, n);
    printf("\n\n\t**** LAPACK RESULTS ****\n");

    for (int i = 0; i < n; i++) {
        printf("\nEigen Value %d/%d: %f\n", i+1, n, WR[i]);
        printf("Eigen Vector %d/%d: ", i+1, n);
        for (int j = 0; j < n; j++){
            if (j < (n-1))
                printf("%lf,", VR[j * n + i]);
            else
                printf("%lf\n", VR[j * n + i]);
        }
    }

    // Free memory and terminate program
    free(A);
    free(x0);
    free(J);
    free(F);
    free(WR);
    free(WL);
    free(VR);
    return EXIT_SUCCESS;
}