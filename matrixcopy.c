
#include <stdio.h>
#include"mpi.h"
#include<stdlib.h>

void printMatrix(int row, int col, float **matrix);
void generateMatrix(int matrixSize, float **matrix);
float **transposeMatrix(int matrixSize, float **matrix);
float **alloc_2d_float(int rows, int cols);
float **calloc_2d_float(int rows, int cols);
void calculateProduct(int r1, int c1, float** m1, int r2, int c2, float**m2, float** res, int res_start_position);

int main( int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  if (world_size < 2) {
    fprintf(stderr, "World size must be greater than 1");
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  
  MPI_Request request=MPI_REQUEST_NULL;
  int matrixSize = 100;
  int modulo = matrixSize % world_size;
  if(modulo != 0) {
    matrixSize = matrixSize + (matrixSize - ((matrixSize / world_size) * world_size));
  }
  int rowsPerCore = matrixSize / world_size; 
  float** a;
  float** b;
  int i;
  if (world_rank == 0) {    
    a  = alloc_2d_float(matrixSize, matrixSize);
    b  = alloc_2d_float(matrixSize,matrixSize);
    generateMatrix(matrixSize ,a);
    generateMatrix(matrixSize, b); 
    //printMatrix(matrixSize, matrixSize, a);
    //printMatrix(matrixSize, matrixSize, b);
    b = transposeMatrix(matrixSize, b);
    for(i = 0; i < world_size; i ++) {
      MPI_Isend(&a[rowsPerCore * i][0], rowsPerCore * matrixSize, MPI_FLOAT, i, 0, MPI_COMM_WORLD, &request);
      MPI_Isend(&b[rowsPerCore * i][0], rowsPerCore * matrixSize, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &request);
    }
    free(a);
    free(b);
  } 
  
  
  double t1, t2; 
  
  t1 = MPI_Wtime(); 
  
  a  = alloc_2d_float(rowsPerCore, matrixSize);
  b  = alloc_2d_float(rowsPerCore,matrixSize);
  
  MPI_Recv(&(a[0][0]), rowsPerCore * matrixSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(&(b[0][0]), rowsPerCore * matrixSize, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  float **res;
  res = calloc_2d_float(rowsPerCore, matrixSize);
  calculateProduct(rowsPerCore, matrixSize, a, rowsPerCore, matrixSize, b, res, world_rank * rowsPerCore);
  
  for (i = 1; i < world_size; i++) {
    MPI_Isend(&b[0][0], rowsPerCore * matrixSize, MPI_FLOAT, getMod((world_rank + i), world_size), 0, MPI_COMM_WORLD, &request);
    MPI_Recv(&(b[0][0]), rowsPerCore * matrixSize, MPI_FLOAT, getMod((world_rank - i), world_size), 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    calculateProduct(rowsPerCore, matrixSize, a, rowsPerCore, matrixSize, b, res, getMod(world_rank * rowsPerCore + rowsPerCore * i, matrixSize));
  }
  
  MPI_Isend(&(res[0][0]), rowsPerCore * matrixSize, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &request);
  free(res);
  
  free(a);
  free(b);

  t2 = MPI_Wtime();
  printf("Elapsed time is %f\n", t2 - t1);
  
  if(world_rank == 0) {
    float **total_res;
    total_res = alloc_2d_float(matrixSize, matrixSize);
    for(i = 0; i < world_size; i++) {
      MPI_Recv(&total_res[rowsPerCore * i][0], rowsPerCore * matrixSize, MPI_FLOAT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    // printMatrix(matrixSize, matrixSize, total_res);
  }
  
  MPI_Finalize();
  return 0;
}

int getMod(int a, int b){ return (a%b+b)%b; }

void printMatrix(int row, int col, float **matrix) {
  int i, j;
  for (i=0;i<row;i++)
  {
    for(j=0;j<col;j++)
    {
      printf("%g ",matrix[i][j]);
      fflush(stdout);
    }
    printf("\n");
    fflush(stdout);
  }
  printf("\n\n\n\n");
  fflush(stdout);
}

float **alloc_2d_float(int rows, int cols) {
  float *data = (float *)malloc(rows*cols*sizeof(float));
  float **array= (float **)malloc(rows*sizeof(float*));
  int i;
  for (i=0; i<rows; i++)
    array[i] = &(data[cols*i]);
  
  return array;
}

float **calloc_2d_float(int rows, int cols) {
  float *data = (float *)calloc(rows*cols,sizeof(float));
  float **array= (float **)calloc(rows,sizeof(float*));
  int i;
  for (i=0; i<rows; i++)
    array[i] = &(data[cols*i]);
  
  return array;
}



float randomFloat(float a, float b) {
  float random = ((float) rand() + 1) / ((float) RAND_MAX + 1);
  float diff = b - a;
  float r = random * diff;
  return a + r;
}

void generateMatrix(int matrixSize, float** matrix) {
  int i,j;
  for(i = 0; i < matrixSize; i++) {
    for(j = 0; j < matrixSize; j++) {
      matrix[i][j] = randomFloat(-1, 1);
    }
  }
}

float** transposeMatrix(int matrixSize, float** matrix) {
  float** result;
  result = alloc_2d_float(matrixSize,matrixSize);
  int i, j;
  for(i = 0; i < matrixSize; i++) {
    for(j = 0; j< matrixSize; j++) {
      result[j][i] = matrix[i][j];
    }
  }
  return result;
}

void calculateProduct(int r1, int c1, float** m1, int r2, int c2, float**m2, float** res, int res_start_position) {
  int i, j, k;
  for(i = 0;i < r1; i++) {
    for(j = 0; j < r2; j++) {
      for(k = 0; k < c1; k++) {
	// This is assuming matrix m2 has been transposed!
	res[i][res_start_position + j] += m1[i][k] * m2[j][k];  
      }
    }
  }  
}
