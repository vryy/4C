#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s,t                  al 9/01    |
 *----------------------------------------------------------------------*/
void b1_jaco(double    **deriv,
                double    **xjm,
                double     *det,
                ELEMENT    *ele,
                int         iel)
{
/*----------------------------------------------------------------------*/
int i,j,k,l;
double dum;
#ifdef DEBUG 
dstrc_enter("b1_jaco");
#endif
/*-------------------------------- determine jacobian at point r,s,t ---*/       
  for (i=0; i<3; i++)
  {
    for (j=0; j<3; j++)
    {
      dum=0.0;
      for (l=0; l<iel; l++)
      {
      dum += deriv[i][l]*ele->node[l]->x[j];
      }
      xjm[i][j]=dum;
    }
  }
/*------------------------------------------ determinant of jacobian ---*/        
  *det = xjm[0][0]*xjm[1][1]*xjm[2][2]+
         xjm[0][1]*xjm[1][2]*xjm[2][0]+
         xjm[0][2]*xjm[1][0]*xjm[2][1]-
         xjm[0][2]*xjm[1][1]*xjm[2][0]-
         xjm[0][0]*xjm[1][2]*xjm[2][1]-
         xjm[0][1]*xjm[1][0]*xjm[2][2];
   
  if (*det<0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");         
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of b1_jaco */
/*----------------------------------------------------------------------*/
