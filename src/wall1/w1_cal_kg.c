/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_kg' which calculates the geometric
       stiffness part (total Lagr) for a wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0711 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"


#ifdef GEMM
extern ALLDYNA      *alldyn;
#endif

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | geometric stiffness part (total lagrange)                    ah 06/02|
 |----------------------------------------------------------------------|
 | kg      -->  ELEMENT STIFFNESS MATRIX kg                             |
 | Boplin  -->  DERIVATIVE OPERATOR                                     |
 | Stress  -->  STRESS MATRIX                                           |
 | fac     -->  INTEGRATION FACTOR: THICK*DET J * WR * WS               |
 | ND      -->  TOTAL NUMBER DEGREES OF FREEDOM OF ELEMENT              |
 | NEPS    -->  ACTUAL NUMBER OF STRAIN COMPONENTS   =4                 |
 |                                                                      |
 *----------------------------------------------------------------------*/
void w1_kg(ELEMENT *ele,
           DOUBLE  **kg,
           DOUBLE  **boplin,
           DOUBLE  **stress,
           DOUBLE    fac,
           INT       nd,
           INT       neps,
	   INT       ip
	   )
{
INT i, j, r, m;
DOUBLE int_stress[4][4];
#ifdef GEMM
STRUCT_DYNAMIC *sdyn;
DOUBLE alpha_f, xsi;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_kg");
#endif
/*----------------------------------------------------------------------*/

#ifdef GEMM
sdyn = alldyn[0].sdyn;
alpha_f = sdyn->alpha_f;
xsi     = sdyn->xsi;
#endif

for(i=0; i<4; i++)
  for(j=0; j<4; j++)
  int_stress[i][j] = 0.0;

for(i=0; i<4; i++)
  for(j=0; j<4; j++)
  {
#ifdef GEMM
  int_stress[i][j] = (1.0-alpha_f+xsi) * stress[i][j] + (alpha_f-xsi) * ele->e.w1->PK_history.a.d3[ip][i][j];
#else
  int_stress[i][j] = stress[i][j];
#endif
  }
/*------------------------------------------------------------New format*/
/*---------------------------------------------- perform B^T * SIGMA * B*/
for(i=0; i<nd; i++)
   for(j=0; j<nd; j++)
      for(r=0; r<neps; r++)
         for(m=0; m<neps; m++)
            kg[i][j] += boplin[r][i]*int_stress[r][m]*boplin[m][j]*fac;


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_kg */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
