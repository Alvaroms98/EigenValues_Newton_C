#ifndef UTILS

#define UTILS
#include <stdio.h>
#include <stdlib.h>

#endif

int getFileRows(const char *filename);
int readMatrix(const char *matrixFilename, double *A, const unsigned int n);
void printMatrix(const double *A, const unsigned int n);
int readX0(const char *x0FileName, double *x0, const unsigned int n);
void printVector(const double *x, const unsigned int n);
