#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*----------------------------------------------------------------------*
 | inversion of arbitrary  N x N - matrix                    ah 9/02    |
 *----------------------------------------------------------------------*/

void w1_inverse_matrix
(
  int N,          /* I: size of matrix: N x N                  */
  double **A,     /* I: Original Matrix (will be destroyed!!!) */
  double **Y      /* O: Inverse Matrix                         */
)
{
  double d;
  int i,j;
  
  static ARRAY    col_a;       
  static double  *col;
  static ARRAY    indx_a;       
  static int     *indx;
  static ARRAY    Acopy_a;       
  static double **Acopy;
  static ARRAY    Ycopy_a;       
  static double **Ycopy;
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

void w1_ludcmp(double **a, /* I: Original Matrix (will be rearanged)     */
               int n,      /* I: size of matrix: N x N                   */
               int *indx,  /* O:row premutation pf pivoting              */
               double *d)  /* O:+1/-1,depending on even or odd row permut*/
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;
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
  w1_free_vector(vv,1,n);
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

void w1_lubksb(double **a,      /* I: LU-decomposed matrix                 */
            int n,           /* I: size of matrix: N x N                */
            int *indx,       /* I:row premutation pf pivoting           */
            double *b)       /* In:B will be changed-> Out:Solution X   */
{
  int i,ii=0,ip,j;
  double sum;
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

double *w1_vector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) dserror("allocation failure in w1_vector()");
	return v-nl+NR_END;
}/* end of w1_vector */

void w1_free_vector(double *v, int nl, int nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}/* end of w1_free_vector */

#undef NR_END
#undef FREE_ARG



#endif
