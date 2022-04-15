#pragma once

#include <complex>

#ifdef __APPLE__
  #include <Accelerate/Accelerate.h>
  typedef __CLPK_doublecomplex CpxDbl;
#elif __linux__
  #include <lapacke.h>
  extern "C" {
    #include <cblas.h>
  }
  typedef lapack_complex_double CpxDbl;
#endif

// Compute the singular value decomposition of a bidiagonal matrix
int dbdsqrCpp(char uplo, int n, int ncvt, int nru, int ncc, 
              double *D, double *E, double *VT, int ldvt, 
              double *U, int ldu, double *C, int ldc);
// Compute LU factorization of A
int dgetrfCpp(int m, int n, double *A, int lda, int *ipiv);
// Compute Cholesky factorization of A
int dpotrfCpp(char uplo, int n, double *A, int lda);
// Solve linear system of equations using LU factorization from dgetrf
int dgetrsCpp(char trans, int n, int nrhs, double *A, int lda, int *ipiv, 
              double *B, int ldb);
// Solve linear system of equations using Cholesky factorization from dpotrf
int dpotrsCpp(char uplo, int n, int nrhs, double *A, int lda, double *b,
              int ldb);
// Invert a matrix given the LU factorization from dgetrf
int dgetriCpp(int n, double *A, int lda, int *ipiv);
// Compute eigenvalues and left and right eigenvectors
int dgeevCpp(char jobvl, char jobvr, int n, double *A, int lda, double *WR, 
             double *WI, double *VL, int ldvl, double *VR, int ldvr);
// Return the reciprocal condition number using norm of A and LU factorization
double dgeconCpp(char norm, int n, double *A, int lda, double anorm);
// Return the matrix norm of A
double dlangeCpp(char norm, int m, int n, double *A, int lda);
// Solve the Sylvester system of equations
int dtrsylCpp(char trana, char tranb, int isgn, int m, int n, double *A, 
              int lda, double *B, int ldb, double *C, int ldc, double scale);
// Schur decomposition of a Hessenberg matrix
int dhseqrCpp(char job, char compz, int n, int ilo, int ihi, double *H, int ldh,
              double *WR, double *WI, double *Z, int ldz);
// Reduce matrix to upper Hessenger form
int dgehrdCpp(int n, int ilo, int ihi, double *A, int lda, double *tau);
// Generate the orthogonal matrix Q implictly created in dgehrd
int dorghrCpp(int n, int ilo, int ihi, double *A, int lda, double *tau);
// Balance a matrix (for better eigenvalue conditioning)
int dgebalCpp(char job, int n, double *A, int lda, int &ilo, int &ihi,
              double *scale);

// Compute LU factorization of complex A
int zgetrfCpp(int m, int n, std::complex<double> *A, int lda, int *ipiv);
// Solve complex linear system of equations using LU
int zgetrsCpp(char trans, int n, int nrhs, std::complex<double> *A, int lda, 
              int *ipiv, std::complex<double> *B, int ldb);
