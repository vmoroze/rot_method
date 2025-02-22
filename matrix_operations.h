#ifndef _MATRIX_OPERATIONS
#define _MATRIX_OPERATIONS

#include <iostream>
#include <fstream>
#include <cmath>

double _max(double a, double b);
double _mod(double a, double b);
int createMatrixFromFile(char *filename, double *matrix, int len);
int createMatrixWithoutFile(double *matrix, int len, int option);
int createMatrixWithoutFile2(double *matrix, int len, int option);
int createColumnB(double *matrix, double *B, int len);
int printMatrix(double *matrix, int len, int endOfOutput, int mode);
double checkDiscrepancy(double *temp, double *A, double *vector, double *B, int len);
double checkInaccuracy(double *vector, int len);
double Discrepancy1(double tr, double *X, int n);
double Discrepancy2(double len, double *X, int n);
int multiplyMatrixVector(double *matrix, double *vector, double *result, int len);
double trace(double *a, int n);
double lenOfMatrix(double *a, int n);

#endif