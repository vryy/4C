/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_fint' which evaluates the internal element
       forces for large def (total Lagr) for a wall element

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
#include "wall1_prototypes.h"

#ifdef GEMM
extern ALLDYNA      *alldyn;
#endif

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | evaluate internal element forces for large def (total Lagr) ah 06/02  |
 *----------------------------------------------------------------------*/
void w1_fint(ELEMENT *ele,        /*active element pointer*/
             DOUBLE **stress,     /* 2.PK stresses        */
             DOUBLE **int_b_bar,
             DOUBLE  *fint,       /* internal forces      */
             DOUBLE   fac,        /* detJ*wr*ws*thickness */
             INT      nd,         /* Element-DOF          */
             INT      ip          /*Integration point     */
	       )

{
/*----------------------------------------------------------------------*/
INT i, j;
DOUBLE st[4];
DOUBLE int_stress[4][4];
#ifdef GEMM
STRUCT_DYNAMIC *sdyn;
DOUBLE alpha_f, xsi;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_fint");
#endif
/*--------------------------------------------- F*S = Vektor fs[i] -----*/

#ifdef GEMM
sdyn = alldyn[0].sdyn;
alpha_f = sdyn->alpha_f;
xsi     = sdyn->xsi;
#endif

/*----------------------------------------------------------B_bar format*/
for(i=0;i<4;i++)
 for(j=0;j<4;j++)
    int_stress[i][j] = 0.0;

for(i=0; i<4; i++)
  for(j=0; j<4; j++){
#ifdef GEMM
  int_stress[i][j] = (1.0-alpha_f+xsi) * stress[i][j] + (alpha_f-xsi) * ele->e.w1->PK_history.a.d3[ip][i][j];
#else
  int_stress[i][j] = stress[i][j];
#endif
  }

st[0] = fac * int_stress[0][0];
st[1] = fac * int_stress[1][1];
st[2] = fac * int_stress[0][2];
st[3] = fac * int_stress[0][2];

for(i=0; i<nd; i++)
  for(j=0; j<4; j++)
    fint[i] += int_b_bar[j][i]*st[j];
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_fint */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
