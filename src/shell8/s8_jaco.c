#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | calculate jacobian matrix                             m.gee 10/01    |
 *----------------------------------------------------------------------*/
void s8jaco(double    *funct,
               double   **deriv,
               double   **x,
               double   **xjm,
               double    *hte,
               double   **a3ref,
               double     e3,
               int        iel,
               double    *det,
               double    *deta,
               int        init)
{
int          i,j,k;
double       x1r;
double       x2r;
double       x3r;
double       x1s;
double       x2s;
double       x3s;

static ARRAY gkov_a;  static double **gkov;
static ARRAY gkon_a;  static double **gkon;
static ARRAY gmkov_a; static double **gmkov;
static ARRAY gmkon_a; static double **gmkon;

#ifdef DEBUG 
dstrc_enter("s8jaco");
#endif
/*----------------------------------------------------------------------*/
/*---------------------------------------------------------- init phase */
if (init==1)
{
   gkov     = amdef("gkov"  ,&gkov_a,3,3,"DA");         
   gkon     = amdef("gkon"  ,&gkon_a,3,3,"DA");         
   gmkov    = amdef("gmkov" ,&gmkov_a,3,3,"DA");        
   gmkon    = amdef("gmkon" ,&gmkon_a,3,3,"DA");        
   goto end;
}
else if (init==-1)
{
   amdel(&gkov_a);   
   amdel(&gkon_a);   
   amdel(&gmkov_a);  
   amdel(&gmkon_a);  
   goto end;
}
/*--------------------------------------------------- calculate metrics */
s8mtr(x,a3ref,e3,gkov,gkon,gmkov,gmkon,det,funct,deriv,hte,iel,1.0,"g1g2g3");
xjm[0][0]=gkov[0][0];
xjm[0][1]=gkov[1][0];
xjm[0][2]=gkov[2][0];
xjm[1][0]=gkov[0][1];
xjm[1][1]=gkov[1][1];
xjm[1][2]=gkov[2][1];
xjm[2][0]=gkov[0][2];
xjm[2][1]=gkov[1][2];
xjm[2][2]=gkov[2][2];

x1r=xjm[0][0];
x2r=xjm[0][1];
x3r=xjm[0][2];
x1s=xjm[1][0];
x2s=xjm[1][1];
x3s=xjm[1][2];

*deta = DSQR(x1r*x2s - x2r*x1s) + DSQR(x3r*x1s - x3s*x1r) + DSQR(x2r*x3s - x3r*x2s);
*deta = sqrt(*deta);

if (*deta <= EPS14) dserror("Element Area equal 0.0 or negativ detected");
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8jaco */
