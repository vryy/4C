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
/* required for lapack access */
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

/*----------------------------------------------------------------------*
 | basic data types, do not use INT, DOUBLE or CHAR in baci !           |
 | only still here to allow the compilation of shell8 and some C stuff  |
 *----------------------------------------------------------------------*/
#ifdef INT
#undef INT
#endif
typedef int       INT;
#ifdef DOUBLE
#undef DOUBLE
#endif
typedef double    DOUBLE;
#ifdef CHAR
#undef CHAR
#endif
typedef char      CHAR;

/* methods located in src/fortran and (theoretically) still in usable  (gjb 7/12) */
/* however, the corresponding source files are currently not in CmakeMakeLists and hence not compiled */
void colsol(DOUBLE *a, DOUBLE *v, INT *maxa, INT *nn, INT *nrr, INT *nrc, INT *nwa, INT *nqm, INT *nr1, INT *nr2, INT *kkk, DOUBLE *det, INT *isc, INT *nsch, INT *ipr, INT *info);
void iluk(INT *n, DOUBLE *a, INT *ja, INT *ia, INT *lfil, DOUBLE *alu, INT *jlu, INT *ju, INT *levs, INT *iwk, DOUBLE *w, INT *jw, INT *ierr);
void lusol(INT *n, DOUBLE *y, DOUBLE *x, DOUBLE *alu, INT *jlu, INT *ju);
void mlpcgveczero(DOUBLE *x, INT *n);
void mlpcgvecvec(DOUBLE *x, DOUBLE *y, DOUBLE *sum, INT *n);
void mlpcgupdupdvec(DOUBLE *a, DOUBLE *y, DOUBLE *facy, DOUBLE *x, DOUBLE *facx, INT *init, INT *n);
void mlpcgupdvec(DOUBLE *y, DOUBLE *x, DOUBLE *fac, INT *init, INT *n);
void dveczero(DOUBLE *x, INT *n);
void iveczero(INT *x, INT *n);
void fortranpow(DOUBLE *V,DOUBLE *R,DOUBLE *RE);


#endif /* COMPILER_DEFINITIONS_H */
