#include "LapackWrapper.h"
#include <memory>

using CpxPtr = CpxDbl*;

int dbdsqrCpp(char uplo, int n, int ncvt, int nru, int ncc, 
              double *D, double *E, double *VT, int ldvt, 
              double *U, int ldu, double *C, int ldc) {
  int info;
  std::unique_ptr<double[]> work(new double[4*n]);
  dbdsqr_(&uplo, &n, &ncvt, &nru, &ncc, D, E, VT, &ldvt, U, &ldu, C, &ldc,
          work.get(), &info);
  return info;
}

int dgetrfCpp(int m, int n, double *A, int lda, int *ipiv) {
  int info;
  dgetrf_(&m, &n, A, &lda, ipiv, &info);
  return info;
}

int dgetrsCpp(char trans, int n, int nrhs, double *A, int lda, int *ipiv, 
              double *B, int ldb) {
  int info;
  dgetrs_(&trans, &n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
  return info;
}

int dpotrfCpp(char uplo, int n, double *A, int lda) {
  int info;
  dpotrf_(&uplo, &n, A, &lda, &info);
  return info;
}

int dpotrsCpp(char uplo, int n, int nrhs, double *A, int lda, double *b, 
              int ldb) {
  int info;
  dpotrs_(&uplo, &n, &nrhs, A, &lda, b, &ldb, &info);
  return info;
}

int dgetriCpp(int n, double *A, int lda, int *ipiv) {
  int info;
  int lwork = -1;
  double optWork;
  dgetri_(&n, A, &lda, ipiv, &optWork, &lwork, &info);
  lwork = optWork;
  std::unique_ptr<double[]> work(new double[lwork]);
  dgetri_(&n, A, &lda, ipiv, work.get(), &lwork, &info);
  return info;
}

int dgeevCpp(char jobvl, char jobvr, int n, double *A, int lda, double *WR, 
             double *WI, double *VL, int ldvl, double *VR, int ldvr) {
  int info;
  int lwork = -1;
  double optWork;
  // find optimal lwork
  dgeev_(&jobvl, &jobvr, &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, &optWork,
         &lwork, &info);
  lwork = int(optWork);
  std::unique_ptr<double[]> work(new double[lwork]);
  dgeev_(&jobvl, &jobvr, &n, A, &lda, WR, WI, VL, &ldvl, VR, &ldvr, work.get(), 
         &lwork, &info);
  return info;
}

int zgetrfCpp(int m, int n, std::complex<double> *A, int lda, int *ipiv) {
  int info;
  zgetrf_(&m, &n, CpxPtr(A), &lda, ipiv, &info);
  return info;
}

int zgetrsCpp(char trans, int n, int nrhs, std::complex<double> *A, int lda, 
              int *ipiv, std::complex<double> *B, int ldb) {
  int info;
  zgetrs_(&trans, &n, &nrhs, CpxPtr(A), &lda, ipiv, CpxPtr(B), &ldb, &info);
  return info;
}

double dgeconCpp(char norm, int n, double *A, int lda, double anorm) {
  int info;
  double rcond;
  std::unique_ptr<double[]> work(new double[4*n]);
  std::unique_ptr<int[]> iwork(new int[n]);
  dgecon_(&norm, &n, A, &lda, &anorm, &rcond, work.get(), iwork.get(), &info);
  return rcond;
}

double dlangeCpp(char norm, int m, int n, double *A, int lda) {
  std::unique_ptr<double[]> work(new double[m]);
  return dlange_(&norm, &m, &n, A, &lda, work.get());
}

int dtrsylCpp(char trana, char tranb, int isgn, int m, int n, double *A, 
              int lda, double *B, int ldb, double *C, int ldc, double scale) {
  int info;
  dtrsyl_(&trana, &tranb, &isgn, &m, &n, A, &lda, B, &ldb, C, &ldc, &scale, 
          &info);
  return info;
}

int dhseqrCpp(char job, char compz, int n, int ilo, int ihi, double *H, int ldh,
              double *WR, double *WI, double *Z, int ldz) {
  int info;
  int lwork = -1;
  double optWork;
  // find optimal lwork
  dhseqr_(&job, &compz, &n, &ilo, &ihi, H, &ldh, WR, WI, Z, &ldz, &optWork, 
          &lwork, &info);
  lwork = int(optWork);
  std::unique_ptr<double[]> work(new double[lwork]);
  dhseqr_(&job, &compz, &n, &ilo, &ihi, H, &ldh, WR, WI, Z, &ldz, work.get(), 
          &lwork, &info);
  return info;
}

int dgehrdCpp(int n, int ilo, int ihi, double *A, int lda, double *tau) {
  int info;
  int lwork = -1;
  double optWork;
  // find optimal lwork
  dgehrd_(&n, &ilo, &ihi, A, &lda, tau, &optWork, &lwork, &info);
  lwork = int(optWork);
  std::unique_ptr<double[]> work(new double[lwork]);
  dgehrd_(&n, &ilo, &ihi, A, &lda, tau, work.get(), &lwork, &info);
  return info;
}

int dorghrCpp(int n, int ilo, int ihi, double *A, int lda, double *tau) {
  int info;
  int lwork = -1;
  double optWork;
  // find optimal lwork
  dorghr_(&n, &ilo, &ihi, A, &lda, tau, &optWork, &lwork, &info);
  lwork = int(optWork);
  std::unique_ptr<double[]> work(new double[lwork]);
  dorghr_(&n, &ilo, &ihi, A, &lda, tau, work.get(), &lwork, &info);
  return info;
}

int dgebalCpp(char job, int n, double *A, int lda, int &ilo, int &ihi,
              double *scale) {
  int info;
  dgebal_(&job, &n, A, &lda, &ilo, &ihi, scale, &info);
  return info;
}
