#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | calculate Ke += Bt * D * B                             m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8_BtDB(double **estif, double **bop, double **D, int iel,
                int numdf, double weight, double **work)
{
int i,j,k;
int dim;
double sum;
#ifdef DEBUG 
dstrc_enter("s8_BtDB");
#endif
/*----------------------------------------------------------------------*/
dim = iel*numdf;
/*------------------------------------ make multiplication work = D * B */
math_matmatdense(work,D,bop,12,12,dim,0,1.0);
/*--------------------------- make multiplication estif += bop^t * work */
math_mattrnmatdense(estif,bop,work,dim,12,dim,1,weight);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_BtDB */
 
