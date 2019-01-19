
// C code MPI cannon matrix - matrix multiplication
// authors: Aravind I and Kapil D.
// date: 11/13/2018
	
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define ROOT 0
#define N 4096
#define index(i,j,rowsize) (int)i*rowsize+j


typedef struct
{
	uint32_t size; //number of elements
    uint32_t allocSize; //size of allocated data
    float * data;
}vector;

vector * readfile(char * filename, int rows, int cols){
		
	FILE * mat_file;
	//char float_in[sizeof(float)];
	
	//initialize return
	vector * retVal = (vector *) malloc(sizeof(vector));
    memset(retVal,0,sizeof(vector));
	
	//initialize matrix
	int size = rows * cols;
	float * array = (float *) calloc(size,sizeof(float));
	
	//open file
	mat_file = fopen(filename,"rb");
	if(mat_file == NULL){
		printf("error: invalid filename \"%s\"", filename);
		exit(-1);
	}
	
	//read data into vector
	fread((void*)array,sizeof(float),size,mat_file);
	
	//close file
	fclose(mat_file);
	
	//set return values
	retVal->data = array;
    retVal->size = size;
    retVal->allocSize = size;
    return retVal;
	
}

float * allocateFloatMatrix(int rows, int cols) {
	//initialize matrix
	int size = rows * cols;
	float * array = (float *)calloc(size, sizeof(float));

	return array;
}

double * allocateDoubleMatrix(int rows, int cols) {
	//initialize matrix
	int size = rows * cols;
	double * array = (double *)calloc(size, sizeof(double));

	return array;
}


/* This function performs a serial matrix-matrix multiplication c = a*b */

matrixMultiply(int n, float *a, float *b, double *c)
{

	int i, j, k;

	for (i=0; i<n; i++)

		for (j=0; j<n; j++)

			for (k=0; k<n; k++)

				c[i*n+j] += a[i*n+k]*b[k*n+j];

}

/* This function performs a serial matrix-matrix multiplication of transposed b matrix c = a*b' */
matrixMultiplyT(int n, float *a, float *b, double *c)
{

	int i, j, k;

	for (i=0; i<n; i++)

		for (j=0; j<n; j++)
	
			for (k=0; k<n; k++)

				c[i*n+j] += a[i*n+k]*b[j*n+k];
}
/* find the transpose of matrix*/
matrixTranspose(int n, float *b, float *c)
{
		int i, j;
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
			{
				c[j*n + i] = b[i*n + j];
			}
}

int main(int argc, char *argv[])
{
	//initialize
	int i, j, k, l, mpi_rank, mpi_size, mpi_rowsize, mpi_colsize, subN, sqrtP;
	int row_rank, col_rank, row, col, destR, destC, src, srcR, srcC;
	// declare variables to store time of parallelism
	double execTime, execStart, execEnd;

	MPI_Comm rowComm, colComm;
	
	MPI_Init(&argc, &argv);

	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

	if (argc < 2) {
		printf("error: need one argument for filename");
		exit(-1);
	}

	// declare variables of type vector to read matrices
	vector * mat1;
	vector * mat2;
	float * vector1;
	float * vector2;
	//create data on root
	// read matrices 
	if (mpi_rank == ROOT)
	{
		mat1 = readfile(argv[1], N, N);
		mat2 = readfile(argv[2], N, N);
		vector1 = mat1->data;
		vector2 = mat2->data;
	}
	// find the sub matrix dimention
	sqrtP = (int) sqrt(mpi_size);
	subN = N/sqrtP;

	//allocate memory for N/4xN rows, and NxN/4 columns
	float * row_mat, *col_mat, *row_mat_rec, *col_mat_rec, *col_matT;
	double *can_res, *can_out;

	//allocate memory for the buffers
	row_mat = allocateFloatMatrix(subN, subN);
	row_mat_rec = allocateFloatMatrix(subN, subN);
	col_mat = allocateFloatMatrix(subN, subN);
	col_matT = allocateFloatMatrix(subN, subN);
	col_mat_rec = allocateFloatMatrix(subN, subN);
	can_res = allocateDoubleMatrix(subN, subN);
	can_out = allocateDoubleMatrix(N, N);

	//create and commit datatypes
	MPI_Datatype arrtype, resized_arrtype, arrtypeD, resized_arrtypeD;

	int sizes[2] = { N,N };
	int subsizes[2] = { subN,subN };
	int starts[2] = { 0,0 };

	MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT, &arrtype);
	MPI_Type_create_resized(arrtype, 0, subN * sizeof(float), &resized_arrtype);
	MPI_Type_commit(&resized_arrtype);
	
	MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_DOUBLE, &arrtypeD);
	MPI_Type_create_resized(arrtypeD, 0, subN * sizeof(double), &resized_arrtypeD);
	MPI_Type_commit(&resized_arrtypeD);

	//calculate send counts and displacements
	int * counts, * displs;
	counts = (int *) malloc(mpi_size * sizeof(int));
	displs = (int *) malloc(mpi_size * sizeof(int));

	for(i = 0; i < mpi_size; i++)
	{		
		counts[i] = 1;
		displs[i] = N*(i/sqrtP) + (i%sqrtP);
	}	
	
	//start timer, compute dot product and record execution time 
	execStart = MPI_Wtime();
	
	//scatterv subarrays
	MPI_Scatterv(vector1, counts, displs, resized_arrtype, row_mat, subN*subN, MPI_FLOAT, ROOT, MPI_COMM_WORLD);
	MPI_Scatterv(vector2, counts, displs, resized_arrtype, col_mat, subN*subN, MPI_FLOAT, ROOT, MPI_COMM_WORLD);
	
	//get row comm and rank
	row = mpi_rank/sqrtP;
	MPI_Comm_split(MPI_COMM_WORLD,row,mpi_rank,&rowComm);
	MPI_Comm_rank(rowComm,&row_rank);
	MPI_Comm_size(rowComm, &mpi_rowsize);
	//get col comm and rank
	col = mpi_rank%sqrtP;
	MPI_Comm_split(MPI_COMM_WORLD,col,mpi_rank,&colComm);
	MPI_Comm_rank(colComm,&col_rank);
	MPI_Comm_size(colComm, &mpi_colsize);
	
	//MPI_Barrier(MPI_COMM_WORLD);
	
	// find the source and destination in row communicator - to left shift rows by row number
	destR = row_rank-row;
    if (destR < 0) {
        destR = destR + mpi_rowsize;
    }
	
	srcR = row_rank+row;
	if (srcR > (mpi_rowsize-1)) {
        srcR = srcR - mpi_rowsize;
    }
	
	// find the source and destination in column communicator  - to north shift columns by column number
	destC = col_rank - col;
    if (destC < 0) {
        destC = destC + mpi_colsize;
    }
	srcC = col_rank+col;
	if (srcC > (mpi_colsize-1)) {
        srcC = srcC - mpi_colsize;
    }
	// left shift rows by row number
	MPI_Sendrecv(row_mat, subN*subN, MPI_FLOAT, destR, 0, row_mat_rec, subN*subN, MPI_FLOAT, srcR, MPI_ANY_TAG, rowComm, MPI_STATUS_IGNORE);

	//north shift columns by column number
	MPI_Sendrecv(col_mat, subN*subN, MPI_FLOAT, destC, 1, col_mat_rec, subN*subN, MPI_FLOAT, srcC, MPI_ANY_TAG, colComm, MPI_STATUS_IGNORE);

	 
	
	for (l=0; l<sqrtP; l++)
	{
		memcpy(row_mat, row_mat_rec, sizeof(float)*subN*subN);
		memcpy(col_mat, col_mat_rec, sizeof(float)*subN*subN);
		
		// Finding the transpose of matrix B
		matrixTranspose(subN, col_mat, col_matT);
		//perform a partial matrix-vector multiplication on each process
		matrixMultiplyT(subN, row_mat, col_matT, can_res);
		
		// find the source and destination in row communicator  - to left shift all rows once
	    if (row_rank != 0) {
			destR = row_rank - 1;
		} else {
			destR = mpi_rowsize - 1;
		}
		
	    srcR = row_rank + 1;
		if (srcR == mpi_rowsize) {
			srcR = 0;
		}
		// find the source and destination in column communicator  - to north shift all columns once
		if (col_rank != 0) {
		destC = col_rank - 1;
		} else {
			destC = mpi_colsize - 1;
		}
		
	    srcC = col_rank + 1;
		if (srcC == mpi_colsize) {
			srcC = 0;
		}
		
		//left shift all rows once
		MPI_Sendrecv(row_mat, subN*subN, MPI_FLOAT, destR, 2, row_mat_rec, subN*subN, MPI_FLOAT, srcR, MPI_ANY_TAG, rowComm, MPI_STATUS_IGNORE);
		
		//north shift all columns once
		MPI_Sendrecv(col_mat, subN*subN, MPI_FLOAT, destC, 3, col_mat_rec, subN*subN, MPI_FLOAT, srcC, MPI_ANY_TAG, colComm, MPI_STATUS_IGNORE);
	}

	// gather the matrix multiplication results from all procs
	MPI_Gatherv(can_res, subN*subN, MPI_DOUBLE, can_out, counts, displs, resized_arrtypeD, ROOT, MPI_COMM_WORLD);
	
		//stop timer
	execEnd = MPI_Wtime();
	execTime = execEnd - execStart;

    //free datatypes
	MPI_Type_free(&resized_arrtype);
	MPI_Type_free(&resized_arrtypeD);
	if (mpi_rank == ROOT)
	{
		printf("Execution time for dot product: %f seconds\n", execTime);
		printf("Result: %f, %f, %f \n ", can_out[0], can_out[2047*N + 2047], can_out[4095*N + 4095]);
		free(vector1);
		free(vector2);
	}
	free(row_mat);
	free(col_mat);
	free(row_mat_rec);
	free(col_mat_rec);
	free(col_matT);
	free(can_res);
	free(can_out);

	//shut down MPI
	MPI_Finalize();

	return 0;
}
