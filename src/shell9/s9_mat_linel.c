/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - 's9_mat_linel' which calculates the constitutive matrix for a St. Venant-
                  Kirchhoff Material
 - 's9_mat_stress1' which calculates the PK_II-Stresses from the green
                  lagrange strains

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates the constitutive matrix for St. Venant-Kirchhoff-Material                                      

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine calculates the constitutive matrix for a St. Venant-Kirchhoff-
Material. The C-Matrix is defined according to the kovariant basis vectors
of the shell body. 
NOTE: Sorting of stresses and strains differ in shell8 and shell9
      s8: [11,12,13,22,23,33] <=> s9: [11,12,22,13,23,33]
</pre>
\param  STVENANT  *mat     (i) Material properties
\param  DOUBLE   **g       (i) kovariant basis vector
\param  DOUBLE   **CC      (o) konsitutive matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_call_mat()   [s9_call_mat.c]

*----------------------------------------------------------------------*/
void s9_mat_linel(STVENANT *mat, DOUBLE **g, DOUBLE **CC)
{
INT i,j,k,l;
DOUBLE xsi=1.0; /*----- shear correction coefficient not yet introduced */
DOUBLE C[3][3][3][3]; /*--------------------------- constitutive tensor */
DOUBLE l1,l2;/*----------------------------------------- lame constants */
DOUBLE emod,nue;/*--------------------------------------- mat constants */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_mat_linel");
#endif
/*----------------------------------------------------------------------*/
emod = mat->youngs;
nue  = mat->possionratio;
l1 = (emod*nue) / ((1.0+nue)*(1.0-2.0*nue));
l2 = emod/ (2.0*(1.0+nue));
/*---------this is not very fast, but corresponds nicely with theory... */
for (i=0; i<3; i++)
for (j=0; j<3; j++)
for (k=0; k<3; k++)
for (l=0; l<3; l++)
C[i][j][k][l] = l1*g[i][j]*g[k][l] + l2*( g[i][k]*g[j][l]+g[i][l]*g[k][j] );
/*----------------------------------------------------------------------*/
CC[0][0] = C[0][0][0][0];
CC[0][1] = C[0][0][1][0];
CC[0][2] = C[0][0][1][1];
CC[0][3] = C[0][0][2][0];
CC[0][4] = C[0][0][2][1];
CC[0][5] = C[0][0][2][2];

CC[1][0] = C[1][0][0][0];
CC[1][1] = C[1][0][1][0];
CC[1][2] = C[1][0][1][1];
CC[1][3] = C[1][0][2][0];
CC[1][4] = C[1][0][2][1];
CC[1][5] = C[1][0][2][2];

CC[2][0] = C[1][1][0][0];
CC[2][1] = C[1][1][1][0];
CC[2][2] = C[1][1][1][1];
CC[2][3] = C[1][1][2][0];
CC[2][4] = C[1][1][2][1];
CC[2][5] = C[1][1][2][2];
 
CC[3][0] = C[2][0][0][0];
CC[3][1] = C[2][0][1][0];
CC[3][2] = C[2][0][1][1];
CC[3][3] = C[2][0][2][0]/xsi;
CC[3][4] = C[2][0][2][1]/xsi;
CC[3][5] = C[2][0][2][2];

CC[4][0] = C[2][1][0][0];
CC[4][1] = C[2][1][1][0];
CC[4][2] = C[2][1][1][1];
CC[4][3] = C[2][1][2][0]/xsi;
CC[4][4] = C[2][1][2][1]/xsi;
CC[4][5] = C[2][1][2][2];

CC[5][0] = C[2][2][0][0];
CC[5][1] = C[2][2][1][0];
CC[5][2] = C[2][2][1][1];
CC[5][3] = C[2][2][2][0];
CC[5][4] = C[2][2][2][1];
CC[5][5] = C[2][2][2][2];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_mat_linel */

/*!----------------------------------------------------------------------
\brief calculates the PK II stresses                                      

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine calculates the PK II stresses from green-lagrange strains
and the constitutive matrix
NOTE: Sorting of stresses and strains differ in shell8 and shell9
      s8: [11,12,13,22,23,33] <=> s9: [11,12,22,13,23,33]
</pre>
\param  DOUBLE   *stress  (o) PK II stresses
\param  DOUBLE   *strain  (i) green-lagrange strains
\param  DOUBLE  **C       (i) konsitutive matrix

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9_call_mat()   [s9_call_mat.c]

*----------------------------------------------------------------------*/
void s9_mat_stress1(DOUBLE stress[6], DOUBLE strain[6], DOUBLE **C)
{
DOUBLE   E[6];
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_mat_linel");
#endif
/*----------------------------------------------------------------------*/
E[0] = strain[0];
E[2] = strain[2];
E[5] = strain[5];
E[1] = strain[1] * 2.0;
E[3] = strain[3] * 2.0;
E[4] = strain[4] * 2.0;

math_matvecdense(stress,C,E,6,6,0,1.0);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s9_mat_linel */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
