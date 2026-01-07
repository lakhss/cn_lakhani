/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"
#include "math.h"

#define PI       3.14159265358979323846

void eig_poisson1D(double* eigval, int *la){
  for (int i = 0; i < (*la); i++){
    eigval[i] = 2 - 2*cos((i*PI)/((*la)+1));
  }
  // TODO: Compute all eigenvalues for the 1D Poisson operator
}

double eigmax_poisson1D(int *la){
  return 2 - 2*cos(((*la)*PI)/((*la)+1));;
  // TODO: Compute and return the maximum eigenvalue for the 1D Poisson operator
}

double eigmin_poisson1D(int *la){
  // TODO: Compute and return the minimum eigenvalue for the 1D Poisson operator
  return 2 - 2*cos((PI)/((*la)+1));
}

double richardson_alpha_opt(int *la){
  return 2/(eigmax_poisson1D(la)+eigmin_poisson1D(la));
}

/**
 * Solve linear system Ax=b using Richardson iteration with fixed relaxation parameter alpha.
 * The iteration is: x^(k+1) = x^(k) + alpha*(b - A*x^(k))
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iteration
  // 1. Compute residual r = b - A*x (use dgbmv for matrix-vector product)
  // 2. Update x = x + alpha*r (use daxpy)
  // 3. Check convergence: ||r||_2 < tol (use dnrm2)
  // 4. Store residual norm in resvec and repeat
  //1 
  int beta = 0;
  int i = 0;
  while (cblas_dnrm2(*la, RHS, 1) >= *tol){
  cblas_dgbmv (CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, *alpha_rich, AB, *(kl)+(*(ku))+1, X, 1, beta, RHS, 1);
  cblas_daxpy(*la,*alpha_rich, RHS, 1, X, 1);
  resvec[i] = cblas_dnrm2(*la, RHS, 1);
  }

}

/**
 * Extract MB for Jacobi method from tridiagonal matrix.
 * Such as the Jacobi iterative process is: x^(k+1) = x^(k) + D^(-1)*(b - A*x^(k))
 */
void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal elements from AB and store in MB
  // MB should contain only the diagonal of A
  for (int i = 0; i < (*la)*(*lab); i+=3){
    MB[i] = 0;
    MB[i+1] = AB[i+1];
    MB[i+2] = 0;
  }
}

/**
 * Extract MB for Gauss-Seidel method from tridiagonal matrix.
 * Such as the Gauss-Seidel iterative process is: x^(k+1) = x^(k) + (D-E)^(-1)*(b - A*x^(k))
 */
void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  // TODO: Extract diagonal and lower diagonal from AB
  // MB should contain the lower triangular part (including diagonal) of A
  for (int i = 0; i < (*la)*(*lab); i+=3){
    MB[i] = 0;
    MB[i+1] = AB[i+1];
    MB[i+2] = AB[i+2];
  }
}

/**
 * Solve linear system Ax=b using preconditioned Richardson iteration.
 * The iteration is: x^(k+1) = x^(k) + M^(-1)*(b - A*x^(k))
 * where M is either D for Jacobi or (D-E) for Gauss-Seidel.
 * Stops when ||b - A*x^(k)||_2  / ||b||_2 < tol or when reaching maxit iterations.
 */
void richardson_MB(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  // TODO: Implement Richardson iterative method
}

