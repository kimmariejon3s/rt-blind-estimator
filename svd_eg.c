#include <stdio.h>
#include <stdlib.h>
#include <svdlib.h>

/* Values for 250 Hz, Bach wav, for Matlab and my code */
#ifdef MATLAB
#define MIN	78
#define MAX 	2574
#else
#define MIN	131
#define MAX	2381
#endif	
#define ROW 	MAX - MIN + 1 
#define COL 	2		/* Fixed don't change!! */

int main(void) {
	int ii, i, j;
	SVDVerbosity = 1;

	SVDRec svd_mat = svdNewSVDRec();
	DMat dmat_a_mat = svdNewDMat(ROW, COL);

	/* Define the matrix */
	for (i = 0; i < ROW; i++) {
		dmat_a_mat->value[i][0] = 1.0;
		dmat_a_mat->value[i][1] = ((double) (i + MIN)) / MAX; 
	}

	svd_mat = svdLAS2A(svdConvertDtoS(dmat_a_mat), 0);

	/* Have svd, now print svd vals and get pinv */
	printf("S:\n");
	for (i = 0; i < COL; i++) {
		printf("s %d: %le\n", i, svd_mat->S[i]);
		svd_mat->S[i] = 1.0 / svd_mat->S[i];
	}

	printf("\nV:\n");
	for (i = 0; i < COL; i++)
		for (j = 0; j < COL; j++) {
			printf("V %d,%d: %lf\n", i, j, svd_mat->Vt->value[i][j]);
			svd_mat->Vt->value[i][j] *= svd_mat->S[i];
	}	

	printf("\nU:\n");
	for (i = 0; i < COL; i++) {
		for (j = 0; j < ROW; j++)
			printf("U %d,%d: %lf\n", i, j, svd_mat->Ut->value[i][j]);
	}

	printf("pinv():\n");
	for (j = 0; j < ROW; j++) {
		dmat_a_mat->value[j][0] = svd_mat->Ut->value[0][j] * 
			svd_mat->Vt->value[0][0];

		dmat_a_mat->value[j][1] = svd_mat->Ut->value[1][j] * 
			svd_mat->Vt->value[1][1];
	
		dmat_a_mat->value[j][0] += svd_mat->Ut->value[1][j] * 
				svd_mat->Vt->value[1][0]; 

		dmat_a_mat->value[j][1] += svd_mat->Ut->value[0][j] *
				svd_mat->Vt->value[0][1];	

		if (j < 5) {
			printf("%d,0: %lf\n", j, dmat_a_mat->value[j][0]);
			printf("%d,1: %lf\n", j, dmat_a_mat->value[j][1]);
		}
	}

	return 0;
}
