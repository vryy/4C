/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_b_bar' which calculates the
       b_bar operator matrix for a wall element at point r,s

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


#ifdef GEMM
extern ALLDYNA      *alldyn;
#endif

/*!
\addtogroup WALL1
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 | calc b_bar operator at Gaussian point r,s                            |
 *----------------------------------------------------------------------*/
void w1_b_barop(ELEMENT *ele, DOUBLE **b_bar,DOUBLE **int_b_bar, DOUBLE **boplin,
                DOUBLE *F, INT numeps, INT nd, INT ip)
{
INT i, j, k;
DOUBLE Fmatrix[4][4];
#ifdef GEMM
STRUCT_DYNAMIC *sdyn;
DOUBLE alpha_f;
#endif
/*------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("w1_b_barop");
#endif
/*------------------------------------------------------------*/
#ifdef GEMM
sdyn = alldyn[0].sdyn;
alpha_f = sdyn->alpha_f;
#endif

for(i=0;i<numeps;i++)
 for(j=0;j<numeps;j++)
     Fmatrix[i][j]=0;


/*---------------------------write Vector F as a matrix Fmatrix*/

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
/*-------------------------------------------------B_bar operator*/
  for(i=0; i<numeps; i++)
    for(j=0; j<nd; j++)
      for(k=0; k<numeps; k++)
      b_bar[i][j] += Fmatrix[k][i]*boplin[k][j];
/*----------------------------------------------------------------*/

/*-------------------- Interpolated B_bar operator(for E-M scheme)*/
  for(i=0; i<numeps; i++)
    for(j=0; j<nd; j++) {
#ifdef GEMM
  int_b_bar[i][j] = (1.0-alpha_f) * b_bar[i][j] + alpha_f * ele->e.w1->b_bar_history.a.d3[ip][i][j];
#else
  int_b_bar[i][j] = b_bar[i][j];
#endif
    }
/*-----------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_boplin */
/*-----------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/

#endif
