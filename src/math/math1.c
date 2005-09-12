/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

------------------------------------------------------------------------*/
#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 |  copies a matrix of identical size                     m.gee 6/01    |
 *----------------------------------------------------------------------*/
void math_array_copy(DOUBLE **from, INT n, INT m, DOUBLE **to)
{
INT i,size;
DOUBLE *ptr_from;
DOUBLE *ptr_to;
#ifdef DEBUG
dstrc_enter("math_array_copy");
#endif
/*----------------------------------------------------------------------*/
/* the arrays to be copied here have to be of the following style:

   a[0]    = [ hole array as a vector ]
   a[1..n] = ptr to start of this row in the vector above
*/
size=n*m;
ptr_from = from[0];
ptr_to   = to[0];
for (i=0; i<size; i++) *(ptr_to++) = *(ptr_from++);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_array_copy */
/*----------------------------------------------------------------------*
 |  calcs the inverse of an unsym 3x3 matrix and determinant m.gee 6/01 |
 *----------------------------------------------------------------------*/
void math_inv3(DOUBLE **a, DOUBLE *det)
{
INT i,j;
DOUBLE b00,b01,b02,b10,b11,b12,b20,b21,b22;
DOUBLE detinv;
#ifdef DEBUG
dstrc_enter("math_inv3");
#endif
/*----------------------------------------------------------------------*/
/* the array to be inved here has to be of the following style:

   a[0]    = [ hole array as a vector ]
   a[1..n] = ptr to start of this row in the vector above
*/
b00 = a[0][0];
b01 = a[0][1];
b02 = a[0][2];
b10 = a[1][0];
b11 = a[1][1];
b12 = a[1][2];
b20 = a[2][0];
b21 = a[2][1];
b22 = a[2][2];

a[0][0] =   b11*b22 - b21*b12;
a[1][0] = - b10*b22 + b20*b12;
a[2][0] =   b10*b21 - b20*b11;
a[0][1] = - b01*b22 + b21*b02;
a[1][1] =   b00*b22 - b20*b02;
a[2][1] = - b00*b21 + b20*b01;
a[0][2] =   b01*b12 - b11*b02;
a[1][2] = - b00*b12 + b10*b02;
a[2][2] =   b00*b11 - b10*b01;

*det = b00*a[0][0]+b01*a[1][0]+b02*a[2][0];
detinv = 1.0/(*det);
for (i=0; i<3; i++)
for (j=0; j<3; j++) a[i][j] *= detinv;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_inv3 */
/*----------------------------------------------------------------------*
 |  transpose of an an unsym nxn matrix                      m.gee 6/01 |
 *----------------------------------------------------------------------*/
void math_tran(DOUBLE **a, INT n)
{
INT i,j;
DOUBLE change;
#ifdef DEBUG
dstrc_enter("math_tran");
#endif
/*----------------------------------------------------------------------*/
/* the array to be transed here has to be of the following style:

   a[0]    = [ hole array as a vector ]
   a[1..n] = ptr to start of this row in the vector above
*/
for (i=0; i<n; i++)
{
   for (j=i+1; j<n; j++)
   {
      change = a[j][i];
      a[j][i] = a[i][j];
      a[i][j] = change;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_tran */



/*----------------------------------------------------------------------*
 |  make a vector unit lenght and return orig lenght         m.gee 6/01 |
 *----------------------------------------------------------------------*/
void math_unvc(
    DOUBLE       *enorm,
    DOUBLE       *vec,
    INT           n)
{

  INT    i;
  DOUBLE skalar;

#ifdef DEBUG
  dstrc_enter("math_unvc");
#endif


  skalar=0.0;
  for (i=0; i<n; i++) skalar += vec[i]*vec[i];
  *enorm = sqrt(skalar);
  if (*enorm < EPS13) dserror("Vector of lenght < EPS13 appeared");
  for (i=0; i<n; i++) vec[i] /= (*enorm);


#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of math_unvc */





/*----------------------------------------------------------------------*
 |  r(I) = A(I,K)*b(K) -----  r = A*b                        m.gee 6/01 |
 |  or                                                                  |
 |  r(I) += A(I,K)*b(K)*factor                                          |
 *----------------------------------------------------------------------*/
void math_matvecdense(DOUBLE  *r,
                         DOUBLE **A,
                         DOUBLE  *b,
                         INT      ni,
                         INT      nk,
                         INT      init,
                         DOUBLE   factor)
{
INT i,k;
DOUBLE sum;
#ifdef DEBUG
dstrc_enter("math_matvecdense");
#endif
/*----------------------------------------------------------------------*/
if (init==0)
{
   for (i=0; i<ni; i++) r[i]=0.0;
}
for (i=0; i<ni; i++)
{
   sum=0.0;
   for (k=0; k<nk; k++) sum += A[i][k]*b[k];
   r[i] += sum*factor;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_matvecdense */
/*----------------------------------------------------------------------*
 |  r(I) = A(K,I)*b(K) -----  r = A*b                        m.gee 6/01 |
 |  or                                                                  |
 |  r(I) += A(K,I)*b(K)*factor                                          |
 *----------------------------------------------------------------------*/
void math_mattrnvecdense(DOUBLE  *r,
                         DOUBLE **A,
                         DOUBLE  *b,
                         INT      ni,
                         INT      nk,
                         INT      init,
                         DOUBLE   factor)
{
INT i,k;
DOUBLE *dptr;
DOUBLE sum;
#ifdef DEBUG
dstrc_enter("math_mattrnvecdense");
#endif
/*----------------------------------------------------------------------*/
if (init==0)
{
   dptr = r;
   for (i=0; i<ni; i++) *(dptr++) = 0.0;
}
for (i=0; i<ni; i++)
{
   sum=0.0;
   for (k=0; k<nk; k++) sum += A[k][i]*b[k];
   r[i] += sum*factor;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_mattrnvecdense */
/*----------------------------------------------------------------------*
 |  R[i][j] = A[i][k]*B[k][j] -----  R=A*B                   m.gee 7/01 |
 | if init==0 result is inited to 0.0                                   |
 |        !=0 result is added  to R multiplied by factor                |
 |  R[i][j] += A[i][k]*B[k][j]*factor                                   |
 |                                                                      |
 | important: this routine works only with arrays that are dynamically  |
 |            allocated the way the AM routines do it                   |
 |            (you can use BLAS as well for this, but you have to care  |
               for the fact, that rows and columns are changed by BLAS )|
 |                                                                      |
 *----------------------------------------------------------------------*/
void math_matmatdense(DOUBLE **R,
                         DOUBLE **A,
                         DOUBLE **B,
                         INT      ni,
                         INT      nk,
                         INT      nj,
                         INT      init,
                         DOUBLE   factor)
{
INT i,j,k;
DOUBLE sum;
#ifdef DEBUG
dstrc_enter("math_matmatdense");
#endif
/*----------------------------------------------------------------------*/
if (init==0)
{
   for (i=0; i<ni; i++)
   {
      for (j=0; j<nj; j++)
      {
         sum=0.0;
         for (k=0; k<nk; k++) sum += A[i][k]*B[k][j];
         R[i][j] = sum;
      }
   }
}
else
{
   for (i=0; i<ni; i++)
   {
      for (j=0; j<nj; j++)
      {
         sum=0.0;
         for (k=0; k<nk; k++) sum += A[i][k]*B[k][j];
         R[i][j] += sum*factor;
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_matmatdense */
/*----------------------------------------------------------------------*
 |  R[i][j] = A[k][i]*B[k][j] -----  R=A^t * B               m.gee 7/01 |
 | if init==0 result is inited to 0.0                                   |
 |        !=0 result is added  to R multiplied by factor                |
 |  R[i][j] += A[k][i]*B[k][j]*factor                                   |
 |                                                                      |
 | important: this routine works only with arrays that are dynamically  |
 |            allocated the way the AM routines do it                   |
 |            (you can use BLAS as well for this, but you have to care  |
               for the fact, that rows and columns are changed by BLAS )|
 |                                                                      |
 *----------------------------------------------------------------------*/
void math_mattrnmatdense(DOUBLE **R,
                            DOUBLE **A,
                            DOUBLE **B,
                            INT      ni,
                            INT      nk,
                            INT      nj,
                            INT      init,
                            DOUBLE   factor)
{
INT i,j,k;
DOUBLE sum;
#ifdef DEBUG
dstrc_enter("math_mattrnmatdense");
#endif
/*----------------------------------------------------------------------*/
if (init==0)
{
   for (i=0; i<ni; i++)
   {
      for (j=0; j<nj; j++)
      {
         sum=0.0;
         for (k=0; k<nk; k++) sum += A[k][i]*B[k][j];
         R[i][j] = sum;
      }
   }
}
else
{
   for (i=0; i<ni; i++)
   {
      for (j=0; j<nj; j++)
      {
         sum=0.0;
         for (k=0; k<nk; k++) sum += A[k][i]*B[k][j];
         R[i][j] += sum*factor;
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_mattrnmatdense */
/*----------------------------------------------------------------------*
 |  R[i][j] = A[i][k]*B[j][k] -----  R=A * B^t               m.gee 7/01 |
 | if init==0 result is inited to 0.0                                   |
 |        !=0 result is added  to R multiplied by factor                |
 |  R[i][j] += A[i][k]*B[j][k]*factor                                   |
 |                                                                      |
 | important: this routine works only with arrays that are dynamically  |
 |            allocated the way the AM routines do it                   |
 |            (you can use BLAS as well for this, but you have to care  |
               for the fact, that rows and columns are changed by BLAS )|
 |                                                                      |
 *----------------------------------------------------------------------*/
void math_matmattrndense(DOUBLE **R,
                            DOUBLE **A,
                            DOUBLE **B,
                            INT      ni,
                            INT      nk,
                            INT      nj,
                            INT      init,
                            DOUBLE   factor)
{
INT i,j,k;
DOUBLE sum;
#ifdef DEBUG
dstrc_enter("math_matmattrndense");
#endif
/*----------------------------------------------------------------------*/
if (init==0)
{
   for (i=0; i<ni; i++)
   {
      for (j=0; j<nj; j++)
      {
         sum=0.0;
         for (k=0; k<nk; k++) sum += A[i][k]*B[j][k];
         R[i][j] = sum;
      }
   }
}
else
{
   for (i=0; i<ni; i++)
   {
      for (j=0; j<nj; j++)
      {
         sum=0.0;
         for (k=0; k<nk; k++) sum += A[i][k]*B[j][k];
         R[i][j] += sum*factor;
      }
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_matmattrndense */

/* used by math_sym_inv and math_unsym_inv */
static DOUBLE   Inverse[4000000];

/*----------------------------------------------------------------------*
 |  do inverse of sym matrix A[n][n]                         m.gee 7/01 |
 |                                                                      |
 | important: this routine works only with arrays that are dynamically  |
 |            allocated the way the AM routines do it                   |
 |            uses LAPACKS dsytrf and dsytri routines                   |
 |                                                                      |
 *----------------------------------------------------------------------*/
void math_sym_inv(DOUBLE **A, INT dim)
{
INT      i,j,k,n;
char     uplo[5];
INT      ipiv[500];
INT      lwork;
DOUBLE   work[5000];
/*DOUBLE   Inverse[250000];*/
INT      info;

#ifdef DEBUG
dstrc_enter("math_sym_inv");
#endif
lwork   = 5000;
info    = 0;
k       = 0;
strncpy(uplo,"L ",2);
/*----------------------------------------------------------------------*/
if (dim>=500) dserror("This routine only sym dense matrices upto size 500");

for (i=0; i<dim; i++)
for (j=0; j<dim; j++)
{
   Inverse[k] = A[i][j];
   k++;
}

n = dim;
dsytrf(uplo,&dim,&(Inverse[0]),&n,&(ipiv[0]),&(work[0]),&lwork,&info);

if (info)
   dserror("Inversion of dense matrix failed");

dsytri(uplo,&dim,&(Inverse[0]),&n,&(ipiv[0]),&(work[0]),&info);

if (info)
   dserror("Inversion of dense matrix failed");

/*--------------------------------------------------- symmetrize result */
for (i=0; i<dim; i++)
for (j=0; j<dim; j++)
{
   A[i][j] = Inverse[info];
   info++;
}
for (i=0; i<dim; i++) for (j=0; j<i; j++) A[i][j] = A[j][i];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_sym_inv */
/*----------------------------------------------------------------------*
 |  do inverse of unsymmetric matrix A[m][n]                 genk 05/02 |
 |                                                                      |
 | important: this routine works only with arrays that are dynamically  |
 |            allocated the way the AM routines do it                   |
 |            uses LAPACKS dgetrf and dgetri routines                   |
 |                                                                      |
 *----------------------------------------------------------------------*/
void math_unsym_inv(DOUBLE **A, INT dimr, INT dimc)
{
INT      i,j;
INT      ipiv[2000];
INT      lwork=2000;
DOUBLE   work[2000];
INT      info=0;

#ifdef DEBUG
dstrc_enter("math_unsym_inv");
#endif
/*----------------------------------------------------------------------*/
if (dimr>=2000 || dimc>=2000)
   dserror("This routine only unsym dense matrices upto size 2000");

for (i=0; i<dimr; i++)
for (j=0; j<dimc; j++)
{
   Inverse[info] = A[i][j];
   info++;
}

dgetrf(&dimr,
       &dimc,
       &(Inverse[0]),
       &dimr,
       &(ipiv[0]),
       &info);

if (info) dserror("Inversion of dense matrix failed");

dgetri(&dimc,
       &(Inverse[0]),
       &dimc,
       &(ipiv[0]),
       &(work[0]),
       &lwork,
       &info);

if (info) dserror("Inversion of dense matrix failed");

/*---------------------------------------------------- copy result to A */
for (i=0; i<dimr; i++)
for (j=0; j<dimc; j++)
{
   A[i][j] = Inverse[info];
   info++;
}
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_sym_inv */

/*----------------------------------------------------------------------*
 |  do inverse of unsymmetric matrix A[m][n]                 genk 05/02 |
 |                                                                      |
 | important: this routine works only with arrays that are dynamically  |
 |            allocated the way the AM routines do it                   |
 |            uses LAPACKS dgetrf and dgetri routines                   |
 |                                                                      |
 *----------------------------------------------------------------------*/
void math_unsym_inv6x6(DOUBLE **A)
{
INT      i,j;
INT      ipiv[6];
INT      lwork=6;
INT      dimr=6;
INT      dimc=6;
DOUBLE   work[6];
DOUBLE   Inverse[36];
INT      info=0;

#ifdef DEBUG
dstrc_enter("math_unsym_inv6x6");
#endif
/*----------------------------------------------------------------------*/

for (i=0; i<6; i++)
for (j=0; j<6; j++)
{
   Inverse[info] = A[i][j];
   info++;
}

dgetrf(&dimr,
       &dimc,
       &(Inverse[0]),
       &dimr,
       &(ipiv[0]),
       &info);

if (info) dserror("Inversion of dense matrix failed");

dgetri(&dimc,
       &(Inverse[0]),
       &dimc,
       &(ipiv[0]),
       &(work[0]),
       &lwork,
       &info);

if (info) dserror("Inversion of dense matrix failed");

/*---------------------------------------------------- copy result to A */
for (i=0; i<6; i++)
for (j=0; j<6; j++)
{
   A[i][j] = Inverse[info];
   info++;
}
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_sym_inv6x6 */

/*----------------------------------------------------------------------*
 |  spatproduct  spt = a*(b x c)                            m.gee 10/01 |
 *----------------------------------------------------------------------*/
void math_sppr(DOUBLE *spat, DOUBLE *a, DOUBLE *b, DOUBLE *c)
{
INT i;
DOUBLE d[3];
#ifdef DEBUG
dstrc_enter("math_sppr");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------------------------- cross product */
d[0] = b[1]*c[2] - b[2]*c[1];
d[1] = b[2]*c[0] - b[0]*c[2];
d[2] = b[0]*c[1] - b[1]*c[0];
/*------------------------------------------------------ scalar product */
*spat = 0.0;
for (i=0; i<3; i++) *spat += a[i]*d[i];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_sppr */

/*-----------------------------------------------------------------------*
 | add up two matriced A = A+B                                           |
 *-----------------------------------------------------------------------*/
void math_addab(DOUBLE **a, DOUBLE **b, INT dim1, INT dim2,DOUBLE fact)
{
int i,j;
#ifdef DEBUG
dstrc_enter("math_addab");
#endif

for (i=0;i<dim1;i++)
{
   for (j=0;j<dim2;j++)
   {
       a[i][j]+=fact*b[i][j];
   }
}

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of math_addab */

/*!---------------------------------------------------------------------
\brief extract digits from integer number

<pre>                                                         genk 04/02
</pre>
\param  num	 INT   (i)    integer number
\param *it	 INT   (o)    integer on position "thousand"
\param *ih       INT   (o)    integer on position "hundred"
\param *id       INT   (o)    integer on position "ten"
\param *id       INT   (o)    integer on position "one"
\return void

------------------------------------------------------------------------*/
void math_intextract(
                    INT num,
                    INT *it,
		    INT *ih,
		    INT *id,
		    INT *io
	            )
{
INT nit, nih, nid, nio;

#ifdef DEBUG
dstrc_enter("intextract");
#endif

nit = num/1000;
nih = (num-nit*1000)/100;
nid = (num-nit*1000-nih*100)/10;
nio = num -nit*1000-nih*100-nid*10;

*it=nit;
*ih=nih;
*id=nid;
*io=nio;

 /*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of intextract*/

