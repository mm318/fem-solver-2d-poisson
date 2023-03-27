#ifndef _DGB_H_
#define _DGB_H_


// helper functions
int i4_min(int i1, int i2);
int i4_max(int i1, int i2);


// matrix functions
int dgb_fa(int n, int ml, int mu, double *a, int *pivot);
double* dgb_sl(int n, int ml, int mu, double *a_lu, int *pivot, double *b, int job);
double* dgb_mxv(int m, int n, int ml, int mu, double *a, double *x);


#endif
