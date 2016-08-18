#include "gm.h"

int main(int argc, char** argv){
	int i,j;
	double input[6] = {2.67, 3.13, 3.25, 3.36, 3.56, 3.72};
	double *pre_processed = pre_process(input, 6);
	printf("X0:\n");
	for(i = 0; i < 6; i++){
		printf("%f ",input[i]);
	}
	printf("\n\nX1:\n");
	for(i = 0; i < 6; i++){
		printf("%f ",pre_processed[i]);
	}
	printf("\n\nB:\n");
	double **B = gen_B_mat(pre_processed, 6);
	for(i = 0; i < 5; i++){
		for(j = 0; j < 2; j++){
			printf("%f ",B[i][j]);
		}
		printf("\n");
	}
	printf("\nY:\n");
	double **Y = gen_Y_mat(input, 6);
	for(i = 0; i < 5; i++){
		printf("%f ",Y[i][0]);
	}
	printf("\n\n");
	double *au = estau(B, Y, 5, 2, 5);
	if(au != NULL){
		printf("a=%f u=%f\n",au[0], au[1]);
	}

	printf("\nmodel gen value x1:\n");
	for(i = 0; i < 6; i++){
		printf("%f ", predict_x1(au[0], au[1], input[0], i));
	}
	printf("\n");
	return 0;
}