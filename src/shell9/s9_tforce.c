/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - s9_tforce:  transforms the forces (stress resultants -> "Schnittgroessen")
               to local/global coordinates ("XYZ", "RST", "RST_ortho")
 - s9_tstress: transforms the stresses to local/global coordinates
               ("XYZ", "RST", "RST_ortho")
 - s9_tettr:   transformation routine


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*!
\addtogroup SHELL9
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief transform forces to local/global coordinates

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine transforms forces (stress resultants -> "Schnittgroessen")
to local/global coordinate systems ("XYZ", "RST", "RST_ortho")
</pre>
\param DOUBLE   **force   (i/o)  forces to be transformed
\param INT        ngauss   (i)   ID of actual GP
\param DOUBLE  ***akov     (i)   metrics at reference surface of each kinematic layer
\param INT        klay     (i)   actual kinematic layer
\param ELEMENT   *ele      (i)   my element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: is not called, as no forces are calculated yet

*----------------------------------------------------------------------*/
void s9_tforce(DOUBLE   **force,
               INT        ngauss,
               DOUBLE  ***akov,
               INT        klay,    /* actual kin layer */
               ELEMENT   *ele)
{
INT    i,j;
INT    option;
DOUBLE sn[3][3],sm[3][3];
DOUBLE b[3][3];/*-------------------- basis vectors of orthogonal basis */
DOUBLE cnorm;
DOUBLE anorm;
DOUBLE sum;

const DOUBLE caeins = 0.99999999999;
const DOUBLE one    = 1.0;
const DOUBLE zero   = 0.0;

static ARRAY    akovKL_a;     static DOUBLE **akovKL;     /* kovariant basis vectors at Int point of actual kin. Layer */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_tforce");
#endif
/*----------------------------------------------------------------------*/
/*- init the basis vectors of actual kin. Layer ------------------------*/
akovKL     = amdef("akovKL"  ,&akovKL_a,3,3,"DA");
for (i=0; i<3; i++) for (j=0; j<3; j++) akovKL[i][j] = akov[i][j][klay];

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
if (ele->e.s9->forcetyp==s9_xyz)       option=0;
if (ele->e.s9->forcetyp==s9_rst)       option=1;
if (ele->e.s9->forcetyp==s9_rst_ortho) option=2;
/*----------------------------------------------------------------------*/
switch(option)
{
case 0:/*--------------------------------- global coordinate orientated */
      for (i=0; i<3; i++)
      for (j=0; j<3; j++) b[i][j] = 0.0;

      b[0][2] = akovKL[1][0]*akovKL[2][1] - akovKL[2][0]*akovKL[1][1];
      b[1][2] = akovKL[2][0]*akovKL[0][1] - akovKL[0][0]*akovKL[2][1];
      b[2][2] = akovKL[0][0]*akovKL[1][1] - akovKL[1][0]*akovKL[0][1];
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
      s9_tettr(sn,akovKL,b);
      /*------------------------------------- make transformation of sm */
      s9_tettr(sm,akovKL,b);
      /*----------------------------------------------------------------*/
break;
case 1:/*---------------------------------------------local coordinates */
       for (i=0; i<3; i++) for (j=0; j<3; j++) b[i][j] = akovKL[i][j];
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
      s9_tettr(sn,akovKL,b);
      /*------------------------------------- make transformation of sm */
      s9_tettr(sm,akovKL,b);
      /*----------------------------------------------------------------*/
break;
case 2:
      /*---------------------------------- local orthogonal coordinates */
      for (i=0; i<3; i++)
      for (j=0; j<3; j++) b[i][j] = 0.0;

      b[0][2] = akovKL[1][0]*akovKL[2][1] - akovKL[2][0]*akovKL[1][1];
      b[1][2] = akovKL[2][0]*akovKL[0][1] - akovKL[0][0]*akovKL[2][1];
      b[2][2] = akovKL[0][0]*akovKL[1][1] - akovKL[1][0]*akovKL[0][1];
      cnorm = sqrt(b[0][2]*b[0][2] + b[1][2]*b[1][2] + b[2][2]*b[2][2]);
      b[0][2] /= cnorm;
      b[1][2] /= cnorm;
      b[2][2] /= cnorm;
      anorm = sqrt(akovKL[0][0]*akovKL[0][0] + akovKL[1][0]*akovKL[1][0] + akovKL[2][0]*akovKL[2][0]);
      b[0][0] = akovKL[0][0]/anorm;
      b[1][0] = akovKL[1][0]/anorm;
      b[2][0] = akovKL[2][0]/anorm;

      b[0][1] = akovKL[1][2]*akovKL[2][0] - akovKL[2][2]*akovKL[1][0];
      b[1][1] = akovKL[2][2]*akovKL[0][0] - akovKL[0][2]*akovKL[2][0];
      b[2][1] = akovKL[0][2]*akovKL[1][0] - akovKL[1][2]*akovKL[0][0];

      /*------------------------------------- make transformation of sn */
      s9_tettr(sn,akovKL,b);
      /*------------------------------------- make transformation of sm */
      s9_tettr(sm,akovKL,b);
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

/*- uninit the basis vectors of actual kin. Layer ----------------------*/
amdel(&akovKL_a);
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_tforce */


/*!----------------------------------------------------------------------
\brief calculate the physical stresses in different coordinate systems

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine calculates the physical stresses in different coordinate
systems ("XYZ", "RST", "RST_ortho"). It is equal to the routine 'S9SFOR'
in the old CARAT.
</pre>
\param DOUBLE   **force    (o)   transformed stresses
\param DOUBLE    *stress   (i)   stresses to be transformed
\param INT        ngauss   (i)   ID of actual GP
\param DOUBLE   **akov     (i)   basic vectors at actual GP
\param ELEMENT   *ele      (i)   my element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by:   s9_stress()    [s9_stress.c]

*----------------------------------------------------------------------*/
void s9_tstress(DOUBLE   **force,
                DOUBLE    *stress,
                INT        ngauss,
                DOUBLE   **akov,
                ELEMENT   *ele)
{
INT    i,j;
INT    option;
DOUBLE s[3][3];   /*--------------------------stresses to be transformed to*/
DOUBLE b[3][3];   /*-------------------- basis vectors of orthogonal basis */

const DOUBLE caeins = 0.99999999999;
const DOUBLE one    = 1.0;
const DOUBLE zero   = 0.0;

DOUBLE cnorm,anorm;
DOUBLE sum;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_tstress");
#endif
/*----------------------------------------------------------------------*/
s[0][0] = stress[0];
s[0][1] = stress[1];
s[0][2] = stress[3];
s[1][0] = stress[1];
s[1][1] = stress[2];
s[1][2] = stress[4];
s[2][0] = stress[3];
s[2][1] = stress[4];
s[2][2] = stress[5];
/*----------------------------------------------------------------------*/
if (ele->e.s9->forcetyp==s9_xyz)       option=0;
if (ele->e.s9->forcetyp==s9_rst)       option=1;
if (ele->e.s9->forcetyp==s9_rst_ortho) option=2;
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

break;
default:
   dserror("Unknown type of option");
break;
}
/*------------------------------------- make transformation of stresses */
s9_tettr(s,akov,b);
/*----------------------------------------------------------------------*/
force[0][ngauss]  = s[0][0];
force[1][ngauss]  = s[0][1];
force[2][ngauss]  = s[1][1];
force[3][ngauss]  = s[0][2];
force[4][ngauss]  = s[1][2];
force[5][ngauss]  = s[2][2];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_tstress */


/*!----------------------------------------------------------------------
\brief transformation to local/global coordinates

<pre>                                                         m.gee 12/01
.........................................................................
.        TENSORTRANSFORMATION                                           .
.        XX(K,L)=T(I,K)*X(I,J)*T(J,L)                                   .
.        T(I,J)=A(I)*B(J)                                               .
.        X(I,J)     - TENSOR 2.TER STUFE                                .
.        A UND B    - BASISVEKTOREN   A ALTE SYS./ B NEUE SYS.          .
.        T          - TRANSFORMATIONSMATRIX                             .
.........................................................................
 NOTE:
   This is NOT a multi-purpose routine, as it takes
   a static 2D array as first argument, a dynamic as second and a
   static as third argument. This routine will fail with any other
   combination of static or dynamic 2D arrays
</pre>
\param  DOUBLE x[3][3]   (i/o) matrix to be transformed(i)-> transformed matrix(o)
\param  DOUBLE **a        (i)  basis vectors old system
\param  DOUBLE b[3][3]    (i)  basis vectors new system

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by:  s9_tforce(); s9_tstress()  [s9_tforce.c]

*----------------------------------------------------------------------*/
void s9_tettr(DOUBLE x[3][3], DOUBLE **a, DOUBLE b[3][3])
{
INT    i,j,k;
DOUBLE t[3][3];
DOUBLE h[3][3];
DOUBLE sum;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_tettr");
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
} /* end of s9_tettr */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
