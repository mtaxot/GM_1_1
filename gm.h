#include<stdio.h>
#include<math.h>
#include<string.h>
#include<stdlib.h>

/**
 ** generate the preprocessed data.
 ** input dim n, output dim n.
 **/
double* pre_process(double *input, int n);

/**
 ** generate the origin data.
 ** input dim n, output dim n.
 **/
double* gen_origin(double *input, int n);


/**
 ** generate the matrix B.
 ** input dim 1* n, output dim (n-1) * 2.
 **/
double** gen_B_mat(double *pre_processed, int n);

/**
 ** generate the matrix Y.
 ** input dim 1* (n-1), output dim (n-1) * 1.
 ** for convenience we ouput 1* (n-1)
 **/
double** gen_Y_mat(double *origin, int n);

/**
 ** est param a and u.
 ** input B and Y, output [a, u]
 **/
double* estau(double **matB, double **matY, int matrow, int matcols, int ycols);

/**
 ** matrix transpose
 **/
double** transpose_mat(double **mat, int rows, int cols);


/**
 ** gauss matrix det
 **/
double** gauss(double **A, int n);

/**
 ** gauss matrix multiply.
 **/
double** mat_multiply(double **A, double **B, int arows, int acols, int brows, int bcols);


/**
 ** this is the final model given a and u, k
 ** k is the index of the sequence start from 0,
 ** the output is x1, so we need to calc the origin data after this.
 **/
double predict_x1(double a, double u, double x0_1, int k);

double check_precision(double *origin, double *pred, int n);


void cat_vec(double *vec, int len);

void cat_mat(double **mat, int rows, int cols);