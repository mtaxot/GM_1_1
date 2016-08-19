#include "gm.h"


double* pre_process(double *input, int n){
	double *output = (double*)malloc(n * sizeof(double));
	int i;
	double sum = 0;
	for(i = 0; i < n; i++){
		sum += input[i];
		output[i] = sum;
	}
	return output;
}

double* gen_origin(double *input, int n){
	double *output = (double*)malloc(n * sizeof(double));
	int i;
	double orig = 0;
	for(i = 0; i < n; i++){
		output[i] = input[i]-orig;
		orig = input[i];
	}
	return output;
}

double** gen_B_mat(double *pre_processed, int n){
	double **mat = (double**)malloc((n - 1)*sizeof(double*));
	int i;
	for(i = 0; i < n - 1; i++){
		mat[i] = malloc(2 * sizeof(double));
		mat[i][0] = -(pre_processed[i] + pre_processed[i + 1]) / 2.0;
		mat[i][1] = 1.0;
	}
	return mat;
}

double** gen_Y_mat(double *origin, int n){
	double **output = (double**)malloc((n - 1) * sizeof(double*));
	int i;
	for(i = 0; i < n - 1; i++){
		output[i] = (double*)malloc(1 * sizeof(double));
	}
	
	
	for(i = 0; i < n - 1; i++){
		output[i][0] = origin[i + 1];
	}
	return output;
}

double** transpose_mat(double **mat, int rows, int cols){
	int new_rows = cols;
	int new_cols = rows;
	double** new_mat = (double**)malloc(new_rows*sizeof(double*));
	int i, j;
	for(i = 0; i < new_rows; i ++){
		new_mat[i] = (double*)malloc(new_cols * sizeof(double));
	}
	
	for(i = 0; i < rows; i++){
		for(j = 0; j < cols; j++){
			new_mat[j][i] = mat[i][j];
		}
	}
	
	return new_mat;
	
}

double** gauss(double **A, int n)
{
    int i, j, k;
    double max, temp;
    double **t = (double**)malloc(n * sizeof(double*));
	for(i = 0; i < n; i++){
		t[i] = (double*)malloc(n * sizeof(double));
	}
	
	double **B = (double**)malloc(n * sizeof(double*));
	for(i = 0; i < n; i++){
		B[i] = (double*)malloc(n * sizeof(double));
	}
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            t[i][j] = A[i][j];
        }
    }
    for (i = 0; i < n; i++){
        for (j = 0; j < n; j++){
            B[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    for (i = 0; i < n; i++){
        max = t[i][i];
        k = i;
        for (j = i + 1; j < n; j++){
            if (fabs(t[j][i]) > fabs(max)){
                max = t[j][i];
                k = j;
            }
        }
        if (k != i){
            for (j = 0; j < n; j++){
                temp = t[i][j];  
                t[i][j] = t[k][j];  
                t[k][j] = temp;  
                temp = B[i][j];  
                B[i][j] = B[k][j];  
                B[k][j] = temp;  
            }
        }
        if (t[i][i] == 0){
            printf("There is no inverse matrix!\n");
			for(i = 0; i < n; i++){
				free(t[i]);
			}
			free(t);
            return NULL;
        }
        temp = t[i][i];
        for (j = 0; j < n; j++){
            t[i][j] = t[i][j] / temp;
            B[i][j] = B[i][j] / temp;
        }
        for (j = 0; j < n; j++){
            if (j != i){
                temp = t[j][i];
                for (k = 0; k < n; k++){
                    t[j][k] = t[j][k] - t[i][k] * temp;
                    B[j][k] = B[j][k] - B[i][k] * temp;
                }
            }
        }
    }
	for(i = 0; i < n; i++){
		free(t[i]);
	}
	free(t);
    return B;
}

double** mat_multiply(double **a, double **b, int arows, int acols, int brows, int bcols){
	if(acols != brows){
		return NULL;
	}

	int i, j, k;
	double **c = (double**)malloc(arows * sizeof(double*));
	for(i = 0; i < arows; i++){
		c[i] = (double*)malloc(bcols * sizeof(double));
	}

	for(i = 0; i < arows; i++){
		
		for(k = 0; k < bcols; k++){
			double temp = 0;
			for(j = 0; j < acols; j++){
				temp += a[i][j] * b[j][k];
			}
			c[i][k] = temp;
		}
	}
	return c;
}


double* estau(double **matB, double **matY, int matrow, int matcols, int ycols){
	int i, j;
	double **matB_transpose = transpose_mat(matB, matrow, matcols);
	/*
	printf("BT:\n");
	for(i = 0; i < matcols; i++){
		for(j = 0; j < matrow; j++){
			printf("%f ", matB_transpose[i][j]);
		}
		printf("\n");
	}
	getchar();
	*/
	double **c = mat_multiply(matB_transpose, matB, matcols, matrow, matrow, matcols);
	if(c == NULL){
		return NULL;
	}
	/*
	for(i = 0; i < matcols; i++){
		for(j = 0; j < matcols; j++){
			printf("%f ", c[i][j]);
		}
		printf("\n");
	}
	getchar();
	*/
	double **c1 = gauss(c, matcols);
	if(c1 == NULL){
		return NULL;
	}
	
	double **d = mat_multiply(c1, matB_transpose, matcols, matcols, matcols, matrow);
	
	
	double **au = mat_multiply(d, matY, matcols, matrow, ycols, 1);
	
	double *au_vec = (double*)malloc(matcols * sizeof(double));
	
	for(i = 0; i < matcols; i++){
		au_vec[i] = au[i][0];
	}
	for(i = 0; i < matcols; i++){
		free(matB_transpose[i]);
	}
	free(matB_transpose);
	for(i = 0; i < matcols; i++){
		free(c[i]);
	}
	free(c);
	for(i = 0; i < matcols; i++){
		free(c1[i]);
	}
	free(c1);
	for(i = 0; i < matcols; i++){
		free(d[i]);
	}
	free(d);

	for(i = 0; i < matcols; i++){
		free(au[i]);
	}
	free(au);
	
	return au_vec;
}

double predict_x1(double a, double u, double x0_0, int k){
	double res = (x0_0 - u / a) * exp(-a * k) + u / a;
	return res;
}
double check_precision(double *origin, double *pred, int n){
	double *delta = (double*)malloc(n * sizeof(double));
	int i;
	for(i = 0; i < n; i++){
		delta[i] = fabs(origin[i] - pred[i]);
	}
	double max = -999999999;
	double min = 999999999;
	for(i = 0; i < n; i++){
		if(delta[i] > max){
			max = delta[i];
		}
		if(delta[i] < min){
			min = delta[i];
		}
	}
	double *y = (double*)malloc(n * sizeof(double));
	for(i = 0; i < n; i++){
		y[i] = (min + 0.5 * max) / (delta[i] + 0.5 * max);
	}
	
	double sum = 0;
	for(i = 0; i < n; i++){
		sum += y[i];
	}
	
	free(delta);
	free(y);
	return sum / n;//sum / n >=0.6 better.
}



void cat_vec(double *vec, int len){
	int i;
	for(i = 0; i < len; i++){
		printf("%f ",vec[i]);
	}
	printf("\n");
}

void cat_mat(double **mat, int rows, int cols){
	int i,j;
	for(i = 0; i < rows; i++){
		for(i = 0; i < cols; i++){
			printf("%f ",mat[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}