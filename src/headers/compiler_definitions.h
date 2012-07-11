/*!---------------------------------------------------------------------
\file compiler_definitions.h
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/

#ifndef COMPILER_DEFINITIONS_H
#define COMPILER_DEFINITIONS_H


/*----------------------------------------------------------------------*
 | special definitions for special compilers.....                       |
 *----------------------------------------------------------------------*/
/* append underslashs, if necessary. Important for linking to fortran routines! */

 #undef CCA_APPEND_U

/* append underslash for gnu's linux compiler gcc and g77 */
/* refer to src/fortran for the respective routines */

/* LINUX_MUENCH is still in use!  */
#ifdef LINUX_MUENCH
#define CCA_APPEND_U (1)
#endif

#ifdef CCA_APPEND_U
#define c1ab                c1ab_
#define c1inv3              c1inv3_
#define c1inv6              c1inv6_
#define c1invf              c1invf_
#define c1jacb              c1jacb_
#define colsol              colsol_
/* required for lapack access. do not remove! */
#define dgesv               dgesv_
#define dgetrf              dgetrf_
#define dgetri              dgetri_
#define dgetrs              dgetrs_
#define dsyev               dsyev_
#define dsyevd              dsyevd_
#define dsygv               dsygv_
#define dsytrf              dsytrf_
#define dsytri              dsytri_
#define dsytrs              dsytrs_
#define dveczero            dveczero_

#define fortranpow          fortranpow_
#define iluk                iluk_
#define iveczero            iveczero_
#define lusol               lusol_
#define mlpcgupdupdvec      mlpcgupdupdvec_
#define mlpcgupdvec         mlpcgupdvec_
#define mlpcgvecvec         mlpcgvecvec_
#define mlpcgveczero        mlpcgveczero_
#define mumps_interface     mumps_interface_
#define qat2v2              qat2v2_
#define s8jacb              s8jacb_
#define v2call              v2call_
#define v3call              v3call_
#define v3call_struct       v3call_struct_

#endif

/* methods located in src/fortran and (theoretically) still in usable  (gjb 7/12) */
/* however, the corresponding source files are currently not in CmakeMakeLists and hence not compiled */
void colsol(double *a, double *v, int *maxa, int *nn, int *nrr, int *nrc, int *nwa, int *nqm, int *nr1, int *nr2, int *kkk, double *det, int *isc, int *nsch, int *ipr, int *info);
void iluk(int *n, double *a, int *ja, int *ia, int *lfil, double *alu, int *jlu, int *ju, int *levs, int *iwk, double *w, int *jw, int *ierr);
void lusol(int *n, double *y, double *x, double *alu, int *jlu, int *ju);
void mlpcgveczero(double *x, int *n);
void mlpcgvecvec(double *x, double *y, double *sum, int *n);
void mlpcgupdupdvec(double *a, double *y, double *facy, double *x, double *facx, int *init, int *n);
void mlpcgupdvec(double *y, double *x, double *fac, int *init, int *n);
void dveczero(double *x, int *n);
void iveczero(int *x, int *n);
void fortranpow(double *V,double *R,double *RE);

/* fortran routines from the lapack package, used in src/linalg/linalg_utils.cpp */
#ifdef __cplusplus
extern "C"
{
#endif

void dsytrf(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);
void dsytri(char *uplo, int *n, double *a, int *lda, int *ipiv, double *work, int *info);
void dgetrf(int *m,int *n, double *a, int *lda, int *ipiv, int* info);
void dgetri(int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

#ifdef __cplusplus
}

#endif

#endif /* COMPILER_DEFINITIONS_H */
