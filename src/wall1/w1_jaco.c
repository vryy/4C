#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s                    al 9/01    |
 *----------------------------------------------------------------------*/
void w1_jaco(double     *funct,
                double    **deriv,
                double    **xjm,
                double     *det,
                ELEMENT    *ele,
                int         iel)
{
/*----------------------------------------------------------------------*/
int k;
#ifdef DEBUG 
dstrc_enter("w1_jaco");
#endif
/*---------------------------------------------- initialize jacobian ---*/
/*      CALL MXIR8 (XJM,NUMDFE,NUMDFE,ZERO) */                                      
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
} /* end of w1_jaco */
