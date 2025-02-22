#include <iostream>
#include <fstream>
#include <time.h>
#include "matrix_operations.h"
#include "eigen_alg.h"
#define eps 1e-30

int main(int argc, char *argv[]){
    /*
        n - dim of matrix
	    m - num of output row and column
		e - accuracy of finding eigenvalues
	    k - option of input matrix (0 - from file)
	    filename - file name which consist matrix
    */

    int n, m, k, i, err;
	int iter;
	int nX;
	double e, left, right;
    char *filename;
	double *massA, *massB, *massX;
	clock_t t;
	double lenForD, trForD;

    try{
		n = atoi(argv[1]);
		m = atoi(argv[2]);
		e = strtod(argv[3], NULL);
		k = atoi(argv[4]);
	} catch(...){
		std::cout << "incorrect input" << std::endl;
		return -1;
	}

	std::cout << std::endl;
	printf("Please, enter left and right borders: ");
	scanf("%lf%lf", &left, &right);

	if(right - left < eps){
		std::cout << "incorrect input" << std::endl;
		return -1;
	}

    try{
		massA = new double[n * n];
		massX = new double[n];
	} catch(...){
		std::cout << "some trouble with memory" << std::endl;
		return -2;
	}

    if(k == 0){
		try{
			filename = argv[5];
		} catch(...){
			std::cout << "incorrect name of file" << std::endl;
			return -3;
		}
		if(createMatrixFromFile(filename, massA, n) != 1){
			delete [] massA;
			delete [] massX;
			return -4;
		}
	} else{
		if(createMatrixWithoutFile2(massA, n, k) != 1){
			delete [] massA;
            delete [] massX;
			return -5;
		}
	}

    std::cout << std::endl;
	std::cout << "Matrix A:" << std::endl;
	if(printMatrix(massA, n, m, 0) != 1){
		delete [] massA;
		delete [] massX;
		return -6;
	}
	std::cout << std::endl;

	lenForD = lenOfMatrix(massA, n);
	trForD = trace(massA, n);

    t = clock();
	err = FindEigenvalues(n, massA, left, right, massX, &nX, e, &iter);
	t = clock() - t;

    if(err != 0){
		std::cout << "incorrect matrix" << std::endl;
		delete [] massA;
		delete [] massX;
		return -7;
	}

	std::cout << "Number of values:\t" << nX << std::endl;
	
	if(nX > 0){
    	std::cout << "Values:" << std::endl;
		if(nX > m){
			printMatrix(massX, n, m, 1);
		}
		else{
			printMatrix(massX, n, nX, 1);
		}
		std::cout << std::endl;
	}

	std::cout << "Time:\t\t\t" << t << std::endl;

    std::cout << "Iterations:\t\t" << iter << std::endl;
	std::cout << std::endl;

	if(nX == n){
		std::cout << "All eigenvalues are located on the specified segment, so we can calculate the discrepancy:" << std::endl;
		std::cout << "Discrepancy I:\t\t" << Discrepancy1(trForD, massX, n) << std::endl;
		std::cout << "Discrepancy II:\t\t" << Discrepancy2(lenForD, massX, n) << std::endl;
		std::cout << std::endl;
	}

    delete [] massA;
	delete [] massX;
	return 1;
}