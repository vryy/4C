#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 |  transform forces to local/global coordinates          m.gee 12/01   |
C!    BERECHNEN DER PHYSIKALISCHEN KRAEFTE UND MOMENTE IN              +
C!    UNTERSCHIEDLICHEN TANGENTIALEN KOORDINATENSYSTEMEN               +
 *----------------------------------------------------------------------*/
void s8_tforce(double **force,
             int      ngauss,
             double **akov,
             double **akon,
             ELEMENT *ele)
{
int    i,j;
int    option;
double sn[3][3],sm[3][3];
double b[3][3];/*-------------------- basis vectors of orthogonal basis */
double cnorm;
double anorm;
double sum;

const double caeins = 0.99999999999;
const double one    = 1.0;
const double zero   = 0.0;

#ifdef DEBUG 
dstrc_enter("s8_tforce");
#endif
/*----------------------------------------------------------------------*/
sn[0][0] = force[0][ngauss];
sn[0][1] = force[2][ngauss];
sn[0][2] = force[3][ngauss];
sn[1][0] = force[8][ngauss];
sn[1][1] = force[1][ngauss];
sn[1][2] = force[4][ngauss];
sn[2][0] = force[16][ngauss];
sn[2][1] = force[17][ngauss];
sn[2][2] = force[9][ngauss];

sm[0][0] = force[5][ngauss];
sm[0][1] = force[7][ngauss];
sm[0][2] = force[10][ngauss];
sm[1][0] = force[14][ngauss];
sm[1][1] = force[6][ngauss];
sm[1][2] = force[11][ngauss];
sm[2][0] = force[12][ngauss];
sm[2][1] = force[13][ngauss];
sm[2][2] = force[15][ngauss];
/*----------------------------------------------------------------------*/
if (ele->e.s8->forcetyp==s8_xyz)       option=0;
if (ele->e.s8->forcetyp==s8_rst)       option=1;
if (ele->e.s8->forcetyp==s8_rst_ortho) option=2;
/*----------------------------------------------------------------------*/
switch(option)
{
case 0:/*--------------------------------- global coordinate orientated */
      for (i=0; i<3; i++)
      for (j=0; j<3; j++) b[i][j] = 0.0;
      
      b[0][2] = akov[1][0]*akov[2][1] - akov[2][0]*akov[1][1];
      b[1][2] = akov[2][0]*akov[0][1] - akov[0][0]*akov[2][1];
      b[2][2] = akov[0][0]*akov[1][1] - akov[1][0]*akov[0][1];
      cnorm = sqrt(b[0][2]*b[0][2] + b[1][2]*b[1][2] + b[2][2]*b[2][2]);
      b[0][2] /= cnorm;
      b[1][2] /= cnorm;
      b[2][2] /= cnorm;
      if (FABS(b[1][2]) < caeins)
      {
         b[0][0] = b[2][2];
         b[1][0] = zero;
         b[2][0] = -b[0][2];
         anorm = sqrt(b[0][0]*b[0][0] + b[2][0]*b[2][0]);
         b[0][0] /= anorm;
         b[2][0] /= anorm;
         
         b[0][1] =  b[1][2]*b[2][0];
         b[1][1] =  b[2][2]*b[0][0] - b[0][2]*b[2][0];
         b[2][1] = -b[1][2]*b[0][0]; 
      }
      else
      {
         b[0][0] = one;
         b[1][0] = zero;
         b[2][0] = zero;
         b[0][1] = zero;
         b[1][1] = zero;
         if (b[1][2] > zero) b[2][1] = -one;
         else                b[2][1] =  one;
      }
      /*------------------------------------- make transformation of sn */
      s8_tettr(sn,akov,b);
      /*------------------------------------- make transformation of sm */
      s8_tettr(sm,akov,b);
      /*----------------------------------------------------------------*/
break;
case 1:/*---------------------------------------------local coordinates */
       for (i=0; i<3; i++) for (j=0; j<3; j++) b[i][j] = akov[i][j];
       /*-------------------------------------------------- normalize b */
       for (j=0; j<3; j++)
       {
          sum = 0.0;
          for (i=0; i<3; i++)  sum += b[i][j]*b[i][j];
          sum = sqrt(sum);
          sum = 1.0/sum;
          for (i=0; i<3; i++)  b[i][j] *= sum;
       }
      /*------------------------------------- make transformation of sn */
      s8_tettr(sn,akov,b);
      /*------------------------------------- make transformation of sm */
      s8_tettr(sm,akov,b);
      /*----------------------------------------------------------------*/
break;
case 2:
      /*---------------------------------- local orthogonal coordinates */
      for (i=0; i<3; i++)
      for (j=0; j<3; j++) b[i][j] = 0.0;
      
      b[0][2] = akov[1][0]*akov[2][1] - akov[2][0]*akov[1][1];
      b[1][2] = akov[2][0]*akov[0][1] - akov[0][0]*akov[2][1];
      b[2][2] = akov[0][0]*akov[1][1] - akov[1][0]*akov[0][1];
      cnorm = sqrt(b[0][2]*b[0][2] + b[1][2]*b[1][2] + b[2][2]*b[2][2]);
      b[0][2] /= cnorm;
      b[1][2] /= cnorm;
      b[2][2] /= cnorm;
      anorm = sqrt(akov[0][0]*akov[0][0] + akov[1][0]*akov[1][0] + akov[2][0]*akov[2][0]);
      b[0][0] = akov[0][0]/anorm;
      b[1][0] = akov[1][0]/anorm;
      b[2][0] = akov[2][0]/anorm;
      
      b[0][1] = akov[1][2]*akov[2][0] - akov[2][2]*akov[1][0];
      b[1][1] = akov[2][2]*akov[0][0] - akov[0][2]*akov[2][0];
      b[2][1] = akov[0][2]*akov[1][0] - akov[1][2]*akov[0][0];

      /*------------------------------------- make transformation of sn */
      s8_tettr(sn,akov,b);
      /*------------------------------------- make transformation of sm */
      s8_tettr(sm,akov,b);
      /*----------------------------------------------------------------*/
break;
default:
   dserror("Unknown type of option");
break;
}
/*----------------------------------------------------------------------*/
force[0][ngauss]  = sn[0][0];
force[2][ngauss]  = sn[0][1];
force[3][ngauss]  = sn[0][2];
force[8][ngauss]  = sn[1][0];
force[1][ngauss]  = sn[1][1];
force[4][ngauss]  = sn[1][2];
force[16][ngauss] = sn[2][0];
force[17][ngauss] = sn[2][1];
force[9][ngauss]  = sn[2][2];

force[5][ngauss]  = sm[0][0];
force[7][ngauss]  = sm[0][1];
force[10][ngauss] = sm[0][2];
force[14][ngauss] = sm[1][0];
force[6][ngauss]  = sm[1][1];
force[11][ngauss] = sm[1][2];
force[12][ngauss] = sm[2][0];
force[13][ngauss] = sm[2][1];
force[15][ngauss] = sm[2][2];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_tforce */




/*----------------------------------------------------------------------*
 |  transform forces to local/global coordinates          m.gee 12/01   |
C.......................................................................
C!    TENSORTRANSFORMATION                                             .
C!    XX(K,L)=T(I,K)*X(I,J)*T(J,L)                                     .
C!    T(I,J)=A(I)*B(J)                                                 .
C!    X(I,J)     - TENSOR 2.TER STUFE                                  .
C!    A UND B    - BASISVEKTOREN   A ALTE SYS./ B NEUE SYS.            .
C!    T          - TRANSFORMATIONSMATRIX                               .
C.......................................................................
 NOTE:
   This is NOT a multi-purpose routine, as it takes
   a static 2D array as first argument, a dynamic as second and a
   static as third argument. This routine will fail with any other 
   combination of static or dynamic 2D arrays
 *----------------------------------------------------------------------*/
void s8_tettr(double x[3][3], double **a, double b[3][3])
{
int    i,j,k;
double t[3][3];
double h[3][3];
double sum;
#ifdef DEBUG 
dstrc_enter("s8_tettr");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
for (j=0; j<3; j++)
{
   t[i][j] = 0.0;
   for (k=0; k<3; k++) t[i][j] += a[k][i]*b[k][j];
}
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
for (j=0; j<3; j++)
{
   sum = 0.0;
   for (k=0; k<3; k++) sum += x[i][k]*t[k][j];
   h[i][j] = sum;
}
/*----------------------------------------------------------------------*/
for (i=0; i<3; i++)
for (j=0; j<3; j++)
{
   sum = 0.0;
   for (k=0; k<3; k++) sum += t[k][i]*h[k][j];
   x[i][j] = sum;
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_tettr */

