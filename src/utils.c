#include "utils.h"

int readMatrix(const char *matrixFilename, double *A, const unsigned int n) {
    FILE *file = fopen(matrixFilename, "r");

    if (file == NULL) {
        printf("Error opening matrix input file\n");
        return EXIT_FAILURE;
    }


    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if ( j == (n-1) )
                fscanf(file, "%lf", &A[i*n + j]);
            else
                fscanf(file, "%lf,", &A[i*n + j]);
        }
    }

    fclose(file);
    return EXIT_SUCCESS;
}

int readX0(const char *x0FileName, double *x0, const unsigned int n) {
    FILE *file = fopen(x0FileName, "r");

    if (file == NULL) {
        printf("Error opening x0 input file\n");
        return EXIT_FAILURE;
    }

    for (int i = 0; i < n; i++) {
        if ( i == (n-1) )
            fscanf(file, "%lf", &x0[i]);
        else
            fscanf(file, "%lf,", &x0[i]);
    }

    fclose(file);
    return EXIT_SUCCESS;
}

int getFileRows(const char *filename) {
    FILE *file = fopen(filename, "r");

    if (file == NULL) {
        printf("Error opening matrix input file\n");
        return EXIT_FAILURE;
    }

    int rows = 1;
    char c;

    while ((c = fgetc(file)) != EOF) {
        if (c == '\n')
            rows++;
    }
    fclose(file);

    return rows;
}

void printMatrix(const double *A, const unsigned int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (j < (n-1))
                printf("%lf,", A[i * n + j]);
            else
                printf("%lf\n", A[i * n + j]);
        }
    }
}

void printVector(const double *x, const unsigned int n) {
    for (int i = 0; i < n; i++) {
        if (i < (n-1))
            printf("%lf,", x[i]);
        else
            printf("%lf\n", x[i]);
    }
}