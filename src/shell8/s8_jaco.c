#ifdef D_SHELL8
#include "../headers/standardtypes.h"
#include "shell8.h"
/*----------------------------------------------------------------------*
 | calculate jacobian matrix                             m.gee 10/01    |
 *----------------------------------------------------------------------*/
void s8jaco(DOUBLE    *funct,
               DOUBLE   **deriv,
               DOUBLE   **x,
               DOUBLE   **xjm,
               DOUBLE    *hte,
               DOUBLE   **a3ref,
               DOUBLE     e3,
               INT        iel,
               DOUBLE    *det,
               DOUBLE    *deta,
               INT        init)
{
DOUBLE       x1r;
DOUBLE       x2r;
DOUBLE       x3r;
DOUBLE       x1s;
DOUBLE       x2s;
DOUBLE       x3s;

static ARRAY gkov_a;  static DOUBLE **gkov;
static ARRAY gkon_a;  static DOUBLE **gkon;
static ARRAY gmkov_a; static DOUBLE **gmkov;
static ARRAY gmkon_a; static DOUBLE **gmkon;

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
#endif
