#ifdef D_WALL1
#include "../headers/standardtypes.h"
#include "wall1.h"
#include "wall1_prototypes.h"
/*!---------------------------------------------------------------------
\brief degnenerateed shape functions and their natural derivatives

<pre>                                                         he  05/03

In this routine the degnerated shape function and their natural first
derivative with respect to r is evaluated for T R I A N G L E S
and  RECTANGLES. They are used for the integrals over the element
edges.

</pre>
\param  *funct     double   (o)    shape functions
\param **deriv     double   (o)    1st natural deriv. of shape funct.
\param **deriv2    double   (o)    2nd natural deriv. of shape funct.
\param   r 	   double   (i)    coordinate
\param   typ 	   DIS_TYP  (i)    element type
\param   icode	   int	    (i)    evaluation flag
\return void

------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | BESTIMMUNG DER Formfunktionen f. Linienintegr. der Trian. he 04/03   |
 *----------------------------------------------------------------------*/
void w1_degrectri(DOUBLE     *funct,
                  DOUBLE    **deriv,
                  DOUBLE      r,
                  DIS_TYP     typ,
                  INT         option
		      )
{
DOUBLE         rr,rp,rm,r2;
const double   q12 = ONE/TWO;

#ifdef DEBUG
dstrc_enter("w1_degrectri");
#endif

/*----------------------------------------------------------------------*/
/* if option ==0 only funtion evaluation, if option==1 also derivatives */
/*----------------------------------------------------------------------*/
rr = r*r;
rp = ONE+r;
rm = ONE-r;
r2 = ONE-rr;

switch(typ)
{
/*------------------------------------------------ rectangular elements */
case tri3:/*-------------- linear rectangular interpolation */
   funct[0] = q12*rm;
   funct[1] = q12*rp;
   if (option==1)              /*--- check for derivative evaluation ---*/
   {
      deriv[0][0]= -q12;
      deriv[0][1]=  q12;
   }
break;
default:
   dserror("unknown typ of interpolation");
break;
} /* end of switch typ */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of w1_degrectri */

/*!---------------------------------------------------------------------
\brief jacobian matrix for integration over the element edges

<pre>                                                        he  05/03

In this routine the jacobian matrix and its determinant is calculated

</pre>
\param  *ele       ELEMENT  (i)    actual element
\param **deriv     double   (i)    natural deriv. of shape funcs
\param **xjm       double   (o)    jacobian matrix
\param  *det       double   (o)    determinant of jacobian matrix
\param   iel       int      (i)    num. of nodes of act. ele
\param  *iedgnod   int      (i)    edgenodes
\return void

------------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | BESTIMMUNG DER JACOBI-MATRIX f. Linienintegr. der Trian.  he 04/03   |
 *----------------------------------------------------------------------*/
void w1_edgejaco(ELEMENT    *ele,
                 DOUBLE    **deriv,
                 DOUBLE    **xjm,
                 DOUBLE     *det,
                 INT         iel,
                 INT        *iedgnod
	           )
{
INT k;
INT node;

#ifdef DEBUG
dstrc_enter("w1_edgejaco");
#endif

/*---------------------------------- determine jacobian at point r,s ---*/
xjm[0][0] = ZERO ;
xjm[0][1] = ZERO ;

for (k=0; k<iel; k++) /* loop all nodes of the element */
{
     node=iedgnod[k];
     xjm[0][0] += deriv[0][k] * ele->node[node]->x[0];
     xjm[0][1] += deriv[0][k] * ele->node[node]->x[1];
} /* end loop over iel */

/*------------------------------------------ determinant of jacobian ---*/
*det = sqrt(xjm[0][0]* xjm[0][0] + xjm[0][1]* xjm[0][1]);

if(*det<ZERO)
dserror("NEGATIVE JACOBIAN DETERMINANT ALONG EDGE\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of w1_edgejaco */

#endif
