#ifndef _MAIN_ALGORITHM
#define _MAIN_ALGORITHM

#include "matrix_operations.h"
#include <cmath>

int FindEigenvalues(int n, double *a, double left, double right, double *x, int *nX, double e, int *iter);
int Rot(int n, double *a);
int n_(int n, double *a, double l);

#endif