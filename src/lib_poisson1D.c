/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include <cblas.h>

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv){
  //AB = malloc ((*lab)*(*la+(*kv))*sizeof(double));
  AB[0] = 0;
  AB[1] = 0;
  AB[2] = 2;
  AB[3] = -1;
  for (int i = 4; i < (*lab)*(*la)-4; i+=4){
    AB[i] = 0;
    AB[i+1] = -1;
    AB[i+2] = 2;
    AB[i+3] = -1;
  }
  AB[(*lab)*(*la)-3] = -1;
  AB[(*lab)*(*la)-2] = 2;
  AB[(*lab)*(*la)-1] = 0;
  // TODO: Fill AB with the tridiagonal Poisson operator
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv){
  for (int i = 0; i < (*lab)*(*la); i = 1 + i + (*lab)){
    AB[i] = 1;
  }
  // TODO: Fill AB with the identity matrix
  // Only the main diagonal should have 1, all other entries are 0
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  RHS[0]= (*BC0);
  RHS[(*la)-1] = (*BC1);
  // TODO: Compute RHS vector
}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  for (int i = 0; i < (*la); i++){
    EX_SOL[i] = (*BC0) + X[i]*((*BC1)-(*BC0));
  }
  // TODO: Compute the exact analytical solution at each grid point
  // This depends on the source term f(x) used in set_dense_RHS_DBC_1D
}  

void set_grid_points_1D(double* x, int* la){
  double space = 1./((*la)+1);
  for (int i = 0; i < (*la); i++){
    x[i] = (i+1) * space;
  }
  // TODO: Generate uniformly spaced grid points in [0,1]
}

double relative_forward_error(double* x, double* y, int* la){
    double sum = 0;
    for (int i = 0; i< (*la); i++){
      sum += pow((x[i]-y[i]),2);
    }
    sum = sqrt(sum);
    double sum2 = 0;
    for (int i = 0; i < (*la); i++){
      sum2 += pow(x[i],2);
    }
    sum2 = sqrt(sum2);
    return sum/sum2;
  }
  // TODO: Compute the relative error using BLAS functions (dnrm2, daxpy or manual loop)

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
  // TODO: Return the correct index formula for column-major band storage
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info){
  for (int i = 0; i < sizeof(ipiv); i++){
    ipiv[i]=i+1;
  }
  for (int i = 0; i < sizeof(ipiv); i++){
    printf("%d ", ipiv[i]);
    printf("\n");
  }

  AB[2] = AB[2];
  AB[3] = AB[3]/AB[2];
  for (int i = 4; i < (*la)*(*lab); i+=4){
    AB[i+1] = AB[i+1];
    AB[i+2] = AB[i+2] - AB[i-1]*AB[i+1];
    AB[i+3] = AB[i+3]/AB[i+2];
  }
  // TODO: Implement specialized LU factorization for tridiagonal matrices
  return 0;
}
