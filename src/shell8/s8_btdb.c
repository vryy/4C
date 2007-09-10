/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/gee/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | calculate Ke += Bt * D * B                             m.gee 6/01    |
 *----------------------------------------------------------------------*/
void s8_BtDB(DOUBLE **estif, DOUBLE **bop, DOUBLE **D, INT iel,
                INT numdf, DOUBLE weight, DOUBLE **work)
{
INT dim;
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
#endif

#endif
