/*!---------------------------------------------------------------------
\file
\brief contains init phase of the shell contact routines

---------------------------------------------------------------------*/
#ifdef S8CONTACT
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "s8contact.h"
#include "shell8.h"

/*! 
\addtogroup CONTACTS8 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief the contact main structure

<pre>                                                         m.gee 2/03    
defined in s8_contact_init.c
</pre>

*----------------------------------------------------------------------*/
extern struct _SHELLCONTACT shellcontact;
/*!---------------------------------------------------------------------
\brief make value of gap function                                              

<pre>                                                        m.gee 2/03 
</pre>
\param actcnode    SHELLNODE*    (i)   the active node
\param ssurf       int*          (i)   indicates top or bottom of slave node
\param msurf       int*          (i)   indicates top or bottom of closest master node
\param actele      ELEMENT*      (i)   element to be projected on 
\param xires       double*       (o)   local coodinates of projection point
\param g           double*       (o)   gap function
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_gapfunction(SHELLNODE  *actcnode,
                            int        *ssurf,
                            int        *msurf,
                            ELEMENT    *actele,
                            double      xi[],
                            double     *g)
{
int          i,j,k;
int          iel;
double       h2;
double       a3r[3][4];
double       a3c[3][4];
double       xr[3][4];
double       xc[3][4];
double       funct[4];
double       deriv[2][4];
double       deriv2[4];
double       thetam,thetas;
double       xs[3],xbar[3];
double       diff[3];
double       tau[3][2];
double       nue[3];
double       dummy;
#ifdef DEBUG 
dstrc_enter("s8_contact_gapfunction");
#endif
/*------------------------------------------ build geometry for element */
iel = actele->numnp;
for (k=0; k<iel; k++)
{
   h2 = actele->e.s8->thick_node.a.dv[k] / 2.0;
   a3r[0][k] = actele->e.s8->a3ref.a.da[0][k] * h2;
   a3r[1][k] = actele->e.s8->a3ref.a.da[1][k] * h2;
   a3r[2][k] = actele->e.s8->a3ref.a.da[2][k] * h2;

   a3c[0][k] = a3r[0][k] + actele->node[k]->sol.a.da[0][3];
   a3c[1][k] = a3r[1][k] + actele->node[k]->sol.a.da[0][4];
   a3c[2][k] = a3r[2][k] + actele->node[k]->sol.a.da[0][5];

   xr[0][k]  = actele->node[k]->x[0];   
   xr[1][k]  = actele->node[k]->x[1];   
   xr[2][k]  = actele->node[k]->x[2];   
   
   xc[0][k]  = xr[0][k] + actele->node[k]->sol.a.da[0][0];
   xc[1][k]  = xr[1][k] + actele->node[k]->sol.a.da[0][1];
   xc[2][k]  = xr[2][k] + actele->node[k]->sol.a.da[0][2];
}
/*---------------- check whether top or bottom projection shall be done */
if (*msurf == 1) thetam =  1.0;
else             thetam = -1.0;
if (*ssurf == 1) thetas =  1.0;
else             thetas = -1.0;
/*----------------------------------- current coodinates of slave point */
xs[0] = actcnode->xc[0] + thetas * actcnode->xc[3];
xs[1] = actcnode->xc[1] + thetas * actcnode->xc[4];
xs[2] = actcnode->xc[2] + thetas * actcnode->xc[5];
/*----------------------------------------- evaluate functions point xi */
s8_contact_functderiv(funct,deriv,deriv2,xi[0],xi[1]);
/*------------------ make base vectors at the given point thetam and xi */
/*-------------------------------------------------- make tau1 and tau2 */
for (j=0; j<2; j++)
   for (i=0; i<3; i++)
   {
      tau[i][j] = 0.0;
      for (k=0; k<iel; k++)
         tau[i][j] += deriv[j][k] * (xc[i][k] + thetam * a3c[i][k]);
   }
/*-------------------- construct vector nue orthogonal to tau1 and tau2 */
nue[0] = tau[1][0]*tau[2][1] - tau[2][0]*tau[1][1];
nue[1] = tau[2][0]*tau[0][1] - tau[0][0]*tau[2][1];
nue[2] = tau[0][0]*tau[1][1] - tau[1][0]*tau[0][1];
/*--------------------------------------------- make vector unit length */
math_unvc(&dummy,nue,3);
/* if we are on the top surface of the master element, this already is outward */
/* if we are on the bottom surface of the element this nue is inward, switch it */
if (thetam==-1)
for (i=0; i<3; i++) nue[i] = -nue[i];
/*----------------------------------- make global projection point xbar */
for (i=0; i<3; i++) xbar[i] = 0.0;
for (k=0; k<iel; k++) 
   for (i=0; i<3; i++)
      xbar[i] += funct[k] * (xc[i][k] + thetam*a3c[i][k]);
/*---------------------------------------------- make difference vector */
for (i=0; i<3; i++)
   diff[i] = xs[i] - xbar[i];
/*--------------------------------- make projection onto outward normal */
*g = 0.0;
for (i=0; i<3; i++)
   (*g) += diff[i] * nue[i];
/*------------------------------------ gap is positive when penetrating */
*g = -(*g);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8_contact_gapfunction */


 













/*! @} (documentation module close)*/

#endif
