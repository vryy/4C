/*!----------------------------------------------------------------------
\file
\brief contains the routine 'w1_keu' which calculates the elastic and
       initaial displacement stiffness (total Lagr) for a wall element

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
void w1_keu(    
	    DOUBLE  **keu, 
            DOUBLE  **b_bar,
            DOUBLE  **int_b_bar,
	    DOUBLE  **D,
            DOUBLE   *F,
            DOUBLE    fac, 
            INT       nd,
            INT       neps)
{

INT       i, j, k, m;          /*New format*/
#ifdef GEMM
STRUCT_DYNAMIC *sdyn;    
DOUBLE    alpha_f, xsi;
#endif
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("w1_keu");
#endif

#ifdef GEMM
sdyn = alldyn[0].sdyn;
alpha_f = sdyn->alpha_f;
xsi     = sdyn->xsi;
#endif
/*------------------------------------------------------------new format*/
  for(i=0; i<nd; i++)
    for(j=0; j<nd; j++)
      for(k=0; k<neps; k++)
        for(m=0; m<neps; m++)
#ifdef GEMM        
	keu[i][j] +=  ((1.0-alpha_f+xsi)/(1.0-alpha_f)) * (int_b_bar[k][i]*D[k][m]*b_bar[m][j]*fac);
#else
        keu[i][j] +=  int_b_bar[k][i]*D[k][m]*b_bar[m][j]*fac;
#endif	
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of w1_keu */
/*----------------------------------------------------------------------*/
#endif /*D_WALL1*/
/*! @} (documentation module close)*/
