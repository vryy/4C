/*!----------------------------------------------------------------------
\file
\brief contains the routine 'ale2_jaco' which calculates the Jacobian 
matrix for a 2d ale element

*----------------------------------------------------------------------*/
#ifdef D_ALE
#include "../headers/standardtypes.h"
#include "ale2.h"

/*! 
\addtogroup Ale 
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief  calculate the Jacobian matrix  

<pre>                                                              mn 06/02 
This routine calculates the Jacobian matrix  at a point r,s for 
a 2D ale element.

</pre>
\param *deriv   DOUBLE     (i)   derivatives of the shape functions
\param **xjm    DOUBLE     (o)   the Jacobian matrix
\param *det     DOUBLE     (i)   determinant of the Jacobian matrix
\param *ele     ELEMENT    (i)   the element
\param iel      INT        (i)   number of nodes of the element

\warning There is nothing special to this routine
\return void                                               
\sa caling: ---; called by: ale2_static_ke()

*----------------------------------------------------------------------*/
void ale2_jaco(DOUBLE    **deriv,
             DOUBLE    **xjm,
             DOUBLE     *det,
             ELEMENT    *ele,
             INT         iel)
{
/*----------------------------------------------------------------------*/
INT k;
#ifdef DEBUG 
dstrc_enter("ale2_jaco");
#endif
/*---------------------------------- determine jacobian at point r,s ---*/       
   xjm[0][0] = 0.0 ;
   xjm[0][1] = 0.0 ;
   xjm[1][0] = 0.0 ;
   xjm[1][1] = 0.0 ;
   
   for (k=0; k<iel; k++)
   {
        xjm[0][0] += deriv[0][k] * ele->node[k]->x[0] ;
        xjm[0][1] += deriv[0][k] * ele->node[k]->x[1] ;
        xjm[1][0] += deriv[1][k] * ele->node[k]->x[0] ;
        xjm[1][1] += deriv[1][k] * ele->node[k]->x[1] ;
   }
/*------------------------------------------ determinant of jacobian ---*/        
     *det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];
   
      if (*det<0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");         
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of ale2_jaco */
#endif
/*! @} (documentation module close)*/
