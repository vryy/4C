/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_history' which is used to update 
       B_bar operator and 2_nd PK at each integration point ip

<pre>
Maintainer: Andrea Hund
            hund@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hund/
            0771 - 685-6122
</pre>

*----------------------------------------------------------------------*/
#ifdef GEMM
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*! 
\addtogroup WALL1 
*//*! @{ (documentation module open)*/
/*----------------------------------------------------------------------*
 | Update B_bar operator and 2_nd PK stresses at each INT. point ip     |
 *----------------------------------------------------------------------*/
void w1_history(ELEMENT *ele, DOUBLE **b_bar, DOUBLE **boplin,
                DOUBLE **stress, DOUBLE *F, INT numeps, INT nd, INT ip)
{

DOUBLE Fmatrix[4][4];
INT i, j, k, m;
#ifdef DEBUG 
dstrc_enter("w1_history");
#endif
for(i=0;i<numeps;i++)
 for(j=0;j<numeps;j++)
    Fmatrix[i][j]=0;

/*-----------------------------write Vector F as a matrix Fmatrix*/  

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
      
/*-------------------------------------------update b_bar_history*/ 
 for(i=0; i<numeps; i++)
   for(j=0; j<nd; j++)
   ele->e.w1->b_bar_history.a.d3[ip][i][j] = 0.0;


 for(i=0; i<numeps; i++)
   for(j=0; j<nd; j++)
   ele->e.w1->b_bar_history.a.d3[ip][i][j] = b_bar[i][j];


/*--------------------------------------------update 2PK stresses*/
for(i=0; i<4; i++)
  for(j=0; j<4; j++)
  ele->e.w1->PK_history.a.d3[ip][i][j] = 0.0;


for(i=0; i<4; i++)
  for(j=0; j<4; j++)
  ele->e.w1->PK_history.a.d3[ip][i][j] = stress[i][j];
/*---------------------------------------------------------------*/   
#ifdef DEBUG 
dstrc_exit();
#endif
return;      
} /* end of w1_history */
#endif /*GEMM*/
/*! @} (documentation module close)*/
