#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 |  copies a matrix of identical size                     m.gee 6/01    |
 *----------------------------------------------------------------------*/
void math_array_copy(double **from, int n, int m, double **to)
{
int i,size;
double *ptr_from;
double *ptr_to;
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
void math_inv3(double **a, double *det)
{
int i,j;
double b00,b01,b02,b10,b11,b12,b20,b21,b22;
double detinv;
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
void math_tran(double **a, int n)
{
int i,j;
double change;
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
void math_unvc(double *enorm,double *vec, int n)
{
int i,j;
double skalar;
double tol = 1.0E-12;
#ifdef DEBUG 
dstrc_enter("math_unvc");
#endif
/*----------------------------------------------------------------------*/
skalar=0.0;
for (i=0; i<n; i++) skalar += vec[i]*vec[i];
*enorm = sqrt(skalar);
if (*enorm <= tol) dserror("Vector of lenght zero appeared");
for (i=0; i<n; i++) vec[i] /= (*enorm);
/*----------------------------------------------------------------------*/
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
void math_matvecdense(double  *r,
                         double **A,
                         double  *b,
                         int      ni,
                         int      nk,
                         int      init,
                         double   factor)
{
int i,k;
double sum;
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
void math_mattrnvecdense(double  *r,
                         double **A,
                         double  *b,
                         int      ni,
                         int      nk,
                         int      init,
                         double   factor)
{
int i,k;
double *dptr;
double sum;
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
void math_matmatdense(double **R,
                         double **A,
                         double **B,
                         int      ni,
                         int      nk,
                         int      nj,
                         int      init,
                         double   factor)
{
int i,j,k;
double sum;
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
void math_mattrnmatdense(double **R,
                            double **A,
                            double **B,
                            int      ni,
                            int      nk,
                            int      nj,
                            int      init,
                            double   factor)
{
int i,j,k;
double sum;
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
void math_matmattrndense(double **R,
                            double **A,
                            double **B,
                            int      ni,
                            int      nk,
                            int      nj,
                            int      init,
                            double   factor)
{
int i,j,k;
double sum;
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
/*----------------------------------------------------------------------*
 |  do inverse of sym matrix A[n][n]                         m.gee 7/01 |
 |                                                                      |
 | important: this routine works only with arrays that are dynamically  |
 |            allocated the way the AM routines do it                   |
 |            uses LAPACKS dsytrf and dsytri routines                   |
 |                                                                      |
 *----------------------------------------------------------------------*/
void math_sym_inv(double **A, int dim)
{
int      i,j;
char     uplo = 'L';
int      ipiv[2000];
int      lwork=2000;
double   work[2000];
double   Inverse[4000000];
int      info=0;

#ifdef DEBUG 
dstrc_enter("math_sym_inv");
#endif
/*----------------------------------------------------------------------*/
if (dim>=2000) dserror("This routine only sym dense matrices upto size 2000");

for (i=0; i<dim; i++) 
for (j=0; j<dim; j++)
{
   Inverse[info] = A[i][j];
   info++;
}

dsytrf(&uplo,
       &dim,
       &(Inverse[0]),
       &dim,
       &(ipiv[0]),
       &(work[0]),
       &lwork,
       &info);

if (info) dserror("Inversion of dense matrix failed");

dsytri(&uplo,
       &dim,
       &(Inverse[0]),  
       &dim,
       &(ipiv[0]),
       &(work[0]),
       &info);

if (info) dserror("Inversion of dense matrix failed");

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
 |            uses LAPACKS dsytrf and dsytri routines                   |
 |                                                                      |
 *----------------------------------------------------------------------*/
void math_unsym_inv(double **A, int dimr, int dimc)
{
int      i,j;
char     uplo = 'L';
int      ipiv[2000];
int      lwork=2000;
double   work[2000];
double   Inverse[4000000];
int      info=0;

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
 |  spatproduct  spt = a*(b x c)                            m.gee 10/01 |
 *----------------------------------------------------------------------*/
void math_sppr(double *spat, double *a, double *b, double *c)
{
int i;
double d[3];
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
 
 
