// This is the main DLL file.

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "svdlib.h"
#include "svdutil.h"

__declspec(dllexport) void ProcessData(char* datafile, char* resultnameblock) {

	clock_t start = clock();
	
	double las2end[2] = {-1.0e-30, 1.0e-30};
	double kappa = 1e-6;
	int iterations = 0;
	int readFormat = SVD_F_SB;
	int writeFormat = SVD_F_DB;

	SMat A = NULL;
	SVDRec R = NULL;
	printf("Loading the matrix...\n");
	A = svdLoadSparseMatrix(datafile, readFormat);

	if (!A) printf("failed to read sparse matrix.  Did you specify the correct file type with the -r argument?\n");
	int dimensions = A->rows;
	if (A->cols < dimensions) dimensions = A->cols;
	printf("Computing the SVD...\n");
	R = svdLAS2(A, dimensions, iterations, las2end, kappa);
	if (R == NULL) {
		printf("Erro processing data\n");
		return;
	}
	printf("MULTIPLICATIONS BY A      = %6ld\n", (SVDCount[SVD_MXV] - R->d) / 2 + R->d);
	printf("MULTIPLICATIONS BY A^T    = %6ld\n", (SVDCount[SVD_MXV] - R->d) / 2);

	R->Ut->rows = R->d;
	R->Vt->rows = R->d;
	char fName[260];
	sprintf(fName, "%s-Ut", resultnameblock);
	svdWriteDenseMatrix(R->Ut, fName, writeFormat);
	sprintf(fName, "%s-S", resultnameblock);
	svdWriteDenseArray(R->S, R->d, fName, FALSE);
	sprintf(fName, "%s-Vt", resultnameblock);
	svdWriteDenseMatrix(R->Vt, fName, writeFormat);

	clock_t end = clock();
	printf("Time %2.3f sec.\n\n", (double)(end - start)/CLOCKS_PER_SEC);
}

