/*!-----------------------------------------------------------------------------
\file
\brief contains the routine 'w1_strain_energy' which calculates the strain_energy
of an element

*-----------------------------------------------------------------------------*/
#ifdef GEMM
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*--------------------------------------------------------------------------*
 | Calculation of strain energy for an element                              |
 *-------------------------------------------------------------------------*/
void w1_strain_energy(ELEMENT *ele, DOUBLE **stress, DOUBLE *strain, DOUBLE fac)

{

DOUBLE strain_matrix[4][4];
INT i, j, k, m;

for(i=0; i<4; i++)
  for(j=0; j<4; j++)
   strain_matrix[i][j] = 0.0;


/*------------------------------------------write strain as a matrix strain_matrix*/  
strain_matrix[0][0]=strain[0]; 
strain_matrix[0][2]=strain[2]; 
strain_matrix[1][1]=strain[1]; 
strain_matrix[1][3]=strain[2]; 
strain_matrix[2][0]=strain[2]; 
strain_matrix[2][2]=strain[1]; 
strain_matrix[3][1]=strain[2]; 
strain_matrix[3][3]=strain[0];   


/*-----------------------------------calculate element contribution to strain energy*/
/*--------------------------St.Venant's Kirchhoff Material(Quadratic energy storage)*/
for(i=0;i<4;i++)
 for(j=0;j<4;j++)
    ele->e.w1->strain_energy = ele->e.w1->strain_energy 
                             + 0.25 * fac * stress[i][j] * strain_matrix[i][j];
			     
/*----------------------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
}/* end of w1_strain_energy */
#endif /*GEMM*/
/*! @} (documentation module close)*/
