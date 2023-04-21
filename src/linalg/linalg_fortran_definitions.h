/*---------------------------------------------------------------------*/
/*! \file

\brief build options and definition of fortran functions

\level 1


*/
/*---------------------------------------------------------------------*/

#ifndef LINALG_FORTRAN_DEFINITIONS_H
#define LINALG_FORTRAN_DEFINITIONS_H

// append underscores, if necessary. Important for linking to fortran routines
#undef CCA_APPEND_U
#define CCA_APPEND_U (1)

#ifdef CCA_APPEND_U
// required to use lapack functions
#define dsytrf dsytrf_
#define dsytri dsytri_
#define dhgeqz dhgeqz_
#define dgghrd dgghrd_
#define dgeqp3 dgeqp3_
#define dggbal dggbal_

#endif

// fortran routines from the lapack package
#ifdef __cplusplus
extern "C"
{
#endif

  void dsytrf(
      char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
  void dsytri(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *info);
  void dhgeqz(char *job, char *compq, char *compz, int *n, int *ilo, int *ihi, double *h, int *ldh,
      double *t, int *ldt, double *alphar, double *alphai, double *beta, double *q, int *ldq,
      double *z, int *ldz, double *work, int *lwork, int *info);
  void dgghrd(char *compq, char *compz, int *n, int *ilo, int *ihi, double *a, int *lda, double *b,
      int *ldb, double *q, int *ldq, double *z, int *lzd, int *info);
  void dgeqp3(int *m, int *n, double *a, int *lda, int *jpvt, double *tau, double *work, int *lwork,
      int *info);
  void dggbal(const char *job, const int *n, double *A, const int *lda, double *B, const int *ldb,
      int *ilo, int *ihi, double *lscale, double *rscale, double *work, int *info);

#ifdef __cplusplus
}

#endif

#endif
