/*!----------------------------------------------------------------------
\file
\brief contains the routine 
 - s9_BtDB: which calculates the stiffness matrix Ke

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief calculates the stiffness matrix Ke                                       

<pre>                     m.gee 6/01              modified by    sh 10/02
This routine calculates the stiffness matrix Ke
                 Ke += Bt * D * B
</pre>
\param  DOUBLE **estif   (i/o) element stiffness matrix to be modified
\param  DOUBLE **bop      (i)  B-Operator matrix
\param  DOUBLE **D        (i)  over thickness integrated material matrix
\param  INT      iel      (i)  number of nodes to this element
\param  INT      numdf    (i)  number of degrees of freedom to one node
\param  DOUBLE   weight   (i)  weight at GP
\param  DOUBLE **work     (i)  matrix to store interim values

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9static_keug() [s9_static_keug.c]

*----------------------------------------------------------------------*/
void s9_BtDB(DOUBLE **estif, DOUBLE **bop,    DOUBLE **D,   INT iel,
             INT      numdf, DOUBLE   weight, DOUBLE **work)
{
INT i,j,k;
INT dim;
DOUBLE sum;
#ifdef DEBUG 
dstrc_enter("s9_BtDB");
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
} /* end of s9_BtDB */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
 
