/* NOT USED ANY MORE 
   FUNCTION IS IN f3_calfuncderiv.c
*/

#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
/*----------------------------------------------------------------------*
 | calculate operator matrix at point r,s, t                genk 05/02  |
 *----------------------------------------------------------------------*/
void f3_jaco(double     *funct,
             double    **deriv,
             double    **xjm,
             double     *det,
             ELEMENT    *ele,
             int         iel)

{
/*----------------------------------------------------------------------*/
int i,j,l;
double dum;

#ifdef DEBUG 
dstrc_enter("f3_jaco");
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
       
if(*det<ZERO)
{   
   printf("\n");
   printf("GLOBAL ELEMENT %i\n",ele->Id);
   dserror("NEGATIVE JACOBIAN DETERMINANT\n");
}
      
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of f3_jaco */
    


#endif
