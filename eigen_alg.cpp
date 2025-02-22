#include "eigen_alg.h"
#define eps 1e-30

//см стр. 94
int FindEigenvalues(int n, double *a, double left, double right, double *x, int *nXout, double e, int *iterout){
	int i, j;
	int iter, count, nX;
	int bef;
	double curL, curR, curM;

	Rot(n, a);

	right += eps;
	left -= eps;

	nX = n_(n, a, right) - n_(n, a, left);

	if (nX == 0){
		*nXout = 0;
		*iterout = 0;
		return 1;
	}

	iter = i = 0;
	bef = n_(n, a, left);
	curL = left;
	curR = right;

	while(i < nX){
		while(curR - curL > e){
			curM = 0.5 * (curL + curR);

			if(n_(n, a, curM) < i + 1 + bef){
				curL = curM;
			}
			else{
				curR = curM;
			}

			iter++;
		}

		curM = 0.5 * (curL + curR);
		count = n_(n, a, curR) - n_(n, a, curL);

		for(j = 0; j < count; j++){
			x[i + j] = curM;
		}

		i += count;

		curL = curM;
		curR = right;
	}

	*nXout = nX;
	*iterout = iter;
    return 0;
}

int n_(int n, double *a, double lmb){			//см. стр. 23 и 96
	int res, i;
	double l;

	l = a[0*n + 0] - lmb;
	res = l < 0. ? 1 : 0;

	for(i = 1; i < n; i++){
		if(fabs(l) < eps){
			l = 1e-10;
		}
		l = a[i*n + i] - lmb - a[i*n + i - 1] * a[(i - 1)*n + i] / l;

		if (l < 0){
			res++;
		}
	}

	return res;
}

int Rot(int n, double *a){			//см. стр. 52, 62 и 67
	int i, j, k;
	double x, y, r, s;
	double a_ii, a_ij, a_ji, a_jj;
	double cosPhi, sinPhi;

	for (i = 1; i < n - 1; i++){
		for (j = i + 1; j < n; j++){
			x = a[i*n + i - 1];
			y = a[j*n + i - 1];

			if(fabs(y) < eps){
				continue;
			}
			r = sqrt(x*x + y*y);
			if(r < eps){
				continue;
			}
			cosPhi = x / r;
			sinPhi = -y / r;

			a[i*n + i - 1] = a[(i - 1)*n + i] = r;
			a[j*n + i - 1] = a[(i - 1)*n + j] = 0.;

			for(k = i + 1; k < n; k++){
				if (k == j){
					continue;
				}
				x = a[i*n + k];
				y = a[j*n + k];
				a[k*n + i] = a[i*n + k] = x*cosPhi - y*sinPhi;
				a[k*n + j] = a[j*n + k] = x*sinPhi + y*cosPhi;
			}

			x = a[i*n + i];
			y = a[j*n + j];
			r = a[i*n + j];
			s = a[j*n + i];

			a_ii = x*cosPhi - s*sinPhi;
			a_ji = x*sinPhi + s*cosPhi;
			a_ij = r*cosPhi - y*sinPhi;
			a_jj = r*sinPhi + y*cosPhi;

			a[i*n + i] = a_ii*cosPhi - a_ij*sinPhi;
			a[j*n + i] = a_ii*sinPhi + a_ij*cosPhi;
			a[i*n + j] = a[j*n + i];
			a[j*n + j] = a_ji*sinPhi + a_jj*cosPhi;
		}
	}

	return 0;
}