/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_fint' which evaluates the internal element
       forces for large def (total Lagr) for a wall element

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
void w1_fint(
	     ELEMENT *ele,        /*active element pointer*/
	     double **stress,     /* 2.PK stresses        */ 
             double  *F,          /* Deformation gradient */ 
             double **int_b_bar,      
             double  *fint,       /* internal forces      */ 
             double   fac,        /* detJ*wr*ws*thickness */ 
             int      nd,         /* Element-DOF          */	     
	     int      ip          /*Integration point     */
	     )
	              
{
/*----------------------------------------------------------------------*/
int i, j, k;
double st[4];
double int_stress[4][4];
#ifdef GEMM
STRUCT_DYNAMIC *sdyn;
double alpha_f, xsi;       
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
