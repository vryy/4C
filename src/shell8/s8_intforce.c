#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | make internal forces                                  m.gee 11/01    |
 *----------------------------------------------------------------------*/
void s8_intforce(double *intforce, double *stress_r, double **bop,
                 int iel, int numdf, int nstress_r, double weight)
{
int nd;

#ifdef DEBUG 
dstrc_enter("s8_intforce");
#endif
/*----------------------------------------------------------------------*/
nd = iel*numdf;
/*
make intforce[nd] = transposed(bop[nstress_r][nd]) * stress_r[nstress_r]
*/
math_mattrnvecdense(intforce,bop,stress_r,nd,nstress_r,1,weight);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_intforce */
#endif
