/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_keu' which calculates the elastic and
       initaial displacement stiffness (total Lagr) for a wall element

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"

/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | elastic and initial displacement stiffness (total lagrange)  ah 06/02|
 |----------------------------------------------------------------------|
 | keu     -->  ELEMENT STIFFNESS MATRIX keu                             |
 | Boplin  -->  DERIVATIVE OPERATOR                                     |
 | D       -->  Tangent matrix                                          |
 | F       -->  Deformation gradient                                    |
 | fac     -->  INTEGRATION FACTOR: THICK*DET J * WR * WS               |
 | ND      -->  TOTAL NUMBER DEGREES OF FREEDOM OF ELEMENT              |
 | NEPS    -->  ACTUAL NUMBER OF STRAIN COMPONENTS   =4                 |
 |                                                                      |
 *----------------------------------------------------------------------*/
void w1_keu(double  **keu, 
            double  **boplin, 
            double  **D,
            double   *F,
            double    fac, 
            int       nd,
            int       neps)
{
int       i, j, r, s, t, l, m;
double    bfcfb;
double    fb[4], cfb[4], fcfb[4];
double    Fmatrix[4][4];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_keu");
#endif
/*------------------------------ write Vector F as a Matrix Fmatrix-----*/
for(i=0;i<neps;i++)
{
  for(j=0;j<neps;j++)
  {
    Fmatrix[i][j]=0;
  }
}
Fmatrix[0][0] = F[0];
Fmatrix[0][2] = 0.5 * F[2];
Fmatrix[0][3] = 0.5 * F[2];
Fmatrix[1][1] = F[1];
Fmatrix[1][2] = 0.5 * F[3];
Fmatrix[1][3] = 0.5 * F[3];
Fmatrix[2][1] = F[2];
Fmatrix[2][2] = 0.5 * F[0];
Fmatrix[2][3] = 0.5 * F[0];
Fmatrix[3][0] = F[3];
Fmatrix[3][2] = 0.5 * F[1];
Fmatrix[3][3] = 0.5 * F[1];

/*----------------------------------------------------------------------*/
for (j=0; j<nd; j++)
{
   for (r=0; r<neps; r++)
   {
      fcfb[r]=0;
      for (s=0; s<neps; s++)
      {
         cfb[s]=0;
         for (t=0; t<neps; t++)
         {
           fb[t]=0;
           for (l=0; l<neps; l++)
           { 
             fb[t] +=  Fmatrix[l][t] * boplin[l][j];
           }
           cfb[s] += D[s][t] * fb[t];
         }
       fcfb[r] +=  Fmatrix[r][s] * cfb[s];
       }
   }
   for (i=0; i<nd; i++)
   {
     bfcfb=0;
     for (m=0; m<neps; m++)
     {
       bfcfb += boplin[m][i] * fcfb[m];
     }
     keu[i][j] += bfcfb * fac;
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_keu */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
