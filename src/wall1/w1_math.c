/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*----------------------------------------------------------------------*
 | inversion of arbitrary  N x N - matrix                    ah 9/02    |
 *----------------------------------------------------------------------*/

void w1_inverse_matrix
(
  INT N,          /* I: size of matrix: N x N                  */
  DOUBLE **A,     /* I: Original Matrix (will be destroyed!!!) */
  DOUBLE **Y      /* O: Inverse Matrix                         */
)
{
  DOUBLE d;
  INT i,j;

  static ARRAY    col_a;
  static DOUBLE  *col;
  static ARRAY    indx_a;
  static INT     *indx;
  static ARRAY    Acopy_a;
  static DOUBLE **Acopy;
  static ARRAY    Ycopy_a;
  static DOUBLE **Ycopy;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("w1_inverse_matrix");
#endif
  /*----------------------------------------------------------------------*/
  col   = amdef("col"  ,&col_a  ,N+1,1,"DV");
  indx  = amdef("indx" ,&indx_a ,N+1,1,"IV");
  Acopy = amdef("Acopy",&Acopy_a,N+1,N+1,"DA");
  Ycopy = amdef("Ycopy",&Ycopy_a,N+1,N+1,"DA");
  /*------------- matrix shifting: 0->1, because loops starting with 1 ---*/
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      Acopy[i+1][j+1] = A[i][j];
    }
  }
  /*--------------------------------- LU-decomposition of Acopy = L*U  ---*/
  /*-------------------------- have a look at "numerical recipes in C"  --*/
  w1_ludcmp(Acopy,N,indx,&d);
  for(j=1;j<=N;j++)
  {
    d *= Acopy[j][j];/*------ determinant -> check if marix is singular---*/
    if(fabs(d)< 1.0E-15)dserror("inversion of singular matrix");
    for(i=1;i<=N;i++) col[i]=0.0;
    col[j]=1.0;
    w1_lubksb(Acopy,N,indx,col);/* solve A*X=B for LU-decomposed matrix A ---*/
    for(i=1;i<=N;i++) Ycopy[i][j]=col[i];
  }
  /*---------------------- matrix back shifting to usual notaion: 1->0 ---*/
  for(i=0;i<N;i++)
  {
    for(j=0;j<N;j++)
    {
      Y[i][j] = Ycopy[i+1][j+1];
    }
  }
  /*----------------------------------------------------------------------*/
  amdel(&col_a);
  amdel(&indx_a);
  amdel(&Acopy_a);
  amdel(&Ycopy_a);
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of w1_inverse_matrix */


/*----------------------------------------------------------------------*
 | LU-Decomposition of  NxN - matrix ->"numerical recipes in C"  ah 9/02|
 *----------------------------------------------------------------------*/

#define NRANSI
#define TINY 1.0e-20;

void w1_ludcmp(DOUBLE **a, /* I: Original Matrix (will be rearanged)     */
               INT n,      /* I: size of matrix: N x N                   */
               INT *indx,  /* O:row premutation pf pivoting              */
               DOUBLE *d)  /* O:+1/-1,depending on even or odd row permut*/
{
  INT i,imax,j,k;
  DOUBLE big,dum,sum,temp;
  DOUBLE *vv;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("w1_ludcmp");
#endif
  /*----------------------------------------------------------------------*/

  vv=w1_vector(1,n);
  *d=1.0;
  for (i=1;i<=n;i++) {
        big=0.0;
        for (j=1;j<=n;j++)
              if ((temp=fabs(a[i][j])) > big) big=temp;
        if (big == 0.0) dserror("Singular matrix in routine w1_ludcmp\n");
        vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
        for (i=1;i<j;i++) {
              sum=a[i][j];
              for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
              a[i][j]=sum;
        }
        big=0.0;
        for (i=j;i<=n;i++) {
              sum=a[i][j];
              for (k=1;k<j;k++)
                    sum -= a[i][k]*a[k][j];
              a[i][j]=sum;
              if ( (dum=vv[i]*fabs(sum)) >= big) {
                    big=dum;
                    imax=i;
              }
        }
        if (j != imax) {
              for (k=1;k<=n;k++) {
                    dum=a[imax][k];
                    a[imax][k]=a[j][k];
                    a[j][k]=dum;
              }
              *d = -(*d);
              vv[imax]=vv[j];
        }
        indx[j]=imax;
        if (a[j][j] == 0.0) a[j][j]=TINY;
        if (j != n) {
              dum=1.0/(a[j][j]);
              for (i=j+1;i<=n;i++) a[i][j] *= dum;
        }
  }
  w1_free_vector(vv,1);
  /*------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
}/* end of w1_ludcmp */
#undef TINY
#undef NRANSI
/*----------------------------------------------------------------------*
 | solving AX=B, A is LU decomposed  ->"numerical recipes in C"  ah 9/02|
 *----------------------------------------------------------------------*/

void w1_lubksb(DOUBLE **a,      /* I: LU-decomposed matrix                 */
            INT n,           /* I: size of matrix: N x N                */
            INT *indx,       /* I:row premutation pf pivoting           */
            DOUBLE *b)       /* In:B will be changed-> Out:Solution X   */
{
  INT i,ii=0,ip,j;
  DOUBLE sum;
  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("w1_lubksb");
#endif
  /*----------------------------------------------------------------------*/
  for (i=1;i<=n;i++) {
        ip=indx[i];
        sum=b[ip];
        b[ip]=b[i];
        if (ii)
              for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
        else if (sum) ii=i;
        b[i]=sum;
  }
  for (i=n;i>=1;i--) {
        sum=b[i];
        for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
        b[i]=sum/a[i][i];
  }
  /*---------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
}/* end of w1_lubksb */

#define NR_END 1
#define FREE_ARG char*

DOUBLE *w1_vector(INT nl, INT nh)
/* allocate a DOUBLE vector with subscript range v[nl..nh] */
{
	DOUBLE *v;

	v=(DOUBLE *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(DOUBLE)));
	if (!v) dserror("allocation failure in w1_vector()");
	return v-nl+NR_END;
}/* end of w1_vector */

void w1_free_vector(DOUBLE *v, INT nl)
/* free a DOUBLE vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}/* end of w1_free_vector */

#undef NR_END
#undef FREE_ARG



#endif
#endif
