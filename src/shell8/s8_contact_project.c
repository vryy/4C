/*!---------------------------------------------------------------------
\file
\brief contains init phase of the shell contact routines

---------------------------------------------------------------------*/
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "s8contact.h"
#include "shell8.h"
/*! 
\addtogroup CONTACT 
*//*! @{ (documentation module open)*/
#ifdef S8CONTACT
/*!----------------------------------------------------------------------
\brief the contact main structure

<pre>                                                         m.gee 2/03    
defined in s8_contact_init.c
</pre>

*----------------------------------------------------------------------*/
extern struct _SHELLCONTACT shellcontact;
/*!---------------------------------------------------------------------
\brief make orthogonal projection onto element                                              

<pre>                                                        m.gee 2/03 
</pre>
\param actcnode    SHELLNODE*    (i)   the active node
\param ssurf       int*          (i)   indicates top or bottom of slave node
\param msurf       int*          (i)   indicates top or bottom of closest master node
\param actele      ELEMENT*      (i)   element to be projected on 
\param xires       double*       (o)   local coodinates of projection point
\param success     int*          (o)   indicates wether projection is inside element
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_orthproject(SHELLNODE  *actcnode,
                            int        *ssurf,
                            int        *msurf,
                            ELEMENT    *actele,
                            double      xires[],
                            double     *distance,
                            int        *success)
{
int          i,j,k,l,m,n,iter;
int          iel;
double       h2;
double       l2;
double       a3r[3][4];
double       a3c[3][4];
double       xr[3][4];
double       xc[3][4];
double       funct[4];
double       deriv[2][4];
double       deriv2[4];
double       thetam,thetas;
double       xs[3],xm[3];
double       xi[2],dxi[2];
double       xbar[3];
double       a3bar[3];
double       xbarxi[3],xbareta[3],xbarxieta[3];
double       a3barxi[3],a3bareta[3],a3barxieta[3];
double       F1,F2,F1xi,F1eta,F2xi,F2eta,detF;
double       sum[3];
int          second;
#ifdef DEBUG 
dstrc_enter("s8_contact_orthproject");
#endif
/*----------------------------------------------------------------------*/
*success = 0;
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
/*----------------------------------------------- start at point xi = 0 */
xi[0]  = 0.0;
xi[1]  = 0.0;
second = 0;
/*------------------------------------------------ start iteration loop */
start:
for (iter=0; iter<80; iter++)
{
   /*-------------------------- evaluate function F1 and F2 at point xi */
   s8_contact_functderiv(funct,deriv,deriv2,xi[0],xi[1]);
   /*------ interpolate xbar, a3bar, xbarxi, xbareta, a3barxi, a3bareta */
   for (i=0; i<3; i++)
   {
      xbar[i]       = 0.0;
      a3bar[i]      = 0.0;
      xbarxi[i]     = 0.0;
      xbareta[i]    = 0.0;
      a3barxi[i]    = 0.0;
      a3bareta[i]   = 0.0;
      xbarxieta[i]  = 0.0;
      a3barxieta[i] = 0.0;
   }
   for (k=0; k<iel; k++)
      for (i=0; i<3; i++)
      {
         xbar[i]       += funct[k]    * xc[i][k];
         xbarxi[i]     += deriv[0][k] * xc[i][k];
         xbareta[i]    += deriv[1][k] * xc[i][k];
         xbarxieta[i]  += deriv2[k]   * xc[i][k];

         a3bar[i]      += funct[k]    * a3c[i][k];
         a3barxi[i]    += deriv[0][k] * a3c[i][k];
         a3bareta[i]   += deriv[1][k] * a3c[i][k];
         a3barxieta[i] += deriv2[k]   * a3c[i][k];
      }
   /*---------------------------------------------------------- make F1 */
   F1 = 0.0;
   for (i=0; i<3; i++)
   {
      F1 +=          xbarxi[i]   * xs[i];
      F1 += thetam * a3barxi[i]  * xs[i];
      F1 -=          xbar[i]     * xbarxi[i];
      F1 -= thetam * xbar[i]     * a3barxi[i];
      F1 -= thetam * a3bar[i]    * xbarxi[i];
      F1 -=          a3bar[i]    * a3barxi[i];
   }
   /*---------------------------------------------------------- make F2 */
   F2 = 0.0;
   for (i=0; i<3; i++)
   {
      F2 +=          xbareta[i]  * xs[i];
      F2 += thetam * a3bareta[i] * xs[i];
      F2 -=          xbar[i]     * xbareta[i];
      F2 -= thetam * xbar[i]     * a3bareta[i];
      F2 -= thetam * a3bar[i]    * xbareta[i];
      F2 -=          a3bar[i]    * a3bareta[i];
   }
   /*-------------------------------------------------------- make F1xi */
   F1xi = 0.0;
   for (i=0; i<3; i++)
   {
      F1xi -=                xbarxi[i]  * xbarxi[i];
      F1xi -= 2.0 * thetam * xbarxi[i]  * a3barxi[i];
      F1xi -=                a3barxi[i] * a3barxi[i];
   }
   /*------------------------------------------------------- make F1eta */
   F1eta = 0.0;
   for (i=0; i<3; i++)
   {
      F1eta +=          xbarxieta[i]  * xs[i];
      F1eta += thetam * a3barxieta[i] * xs[i];
      F1eta -=          xbareta[i]    * xbarxi[i];
      F1eta -=          xbar[i]       * xbarxieta[i];
      F1eta -= thetam * xbareta[i]    * a3barxi[i];
      F1eta -= thetam * xbar[i]       * a3barxieta[i];
      F1eta -= thetam * a3bareta[i]   * xbarxi[i];
      F1eta -= thetam * a3bar[i]      * xbarxieta[i];
      F1eta -=          a3bareta[i]   * a3barxi[i];
      F1eta -=          a3bar[i]      * a3barxieta[i];
   }
   /*-------------------------------------------------------- make F2xi */
   F2xi = 0.0;
   for (i=0; i<3; i++)
   {
      F2xi +=          xbarxieta[i]  * xs[i];
      F2xi += thetam * a3barxieta[i] * xs[i];
      F2xi -=          xbarxi[i]     * xbareta[i];
      F2xi -=          xbar[i]       * xbarxieta[i];
      F2xi -= thetam * xbarxi[i]     * a3bareta[i];
      F2xi -= thetam * xbar[i]       * a3barxieta[i];
      F2xi -= thetam * a3barxi[i]    * xbareta[i];
      F2xi -= thetam * a3bar[i]      * xbarxieta[i];
      F2xi -=          a3barxi[i]    * a3bareta[i];
      F2xi -=          a3bar[i]      * a3barxieta[i];
   }
   /*------------------------------------------------------- make F2eta */
   F2eta = 0.0;
   for (i=0; i<3; i++)
   {
      F2eta -=                xbareta[i]  * xbareta[i];
      F2eta -= 2.0 * thetam * xbareta[i]  * a3bareta[i];
      F2eta -=                a3bareta[i] * a3bareta[i];
      
   }
   /*---------------------------------- solution of the equation system */
   detF = F1xi * F2eta - F1eta * F2xi;
/*   if (detF <= 0.0) printf("Determinant smaller zero detF = %20.10E\n",detF);*/
   detF = 1.0/detF;
   
   dxi[0] = detF * (F2eta*(-F1) + F1eta *   F2 );
   dxi[1] = detF * (F2xi *  F1  + F1xi  * (-F2));

   xi[0] += dxi[0];
   xi[1] += dxi[1];

   l2 = dxi[0]*dxi[0] + dxi[1]*dxi[1];
   l2 = sqrt(l2);
   
   if (l2 < EPS12) break;
   
}/* end of  for (iter=0; iter<80; iter++)  */
/*------------------------ check whether projection point is in element */
if ( -1.01  < xi[0] && xi[0] < 1.01) *success = 1;
else                                 *success = 0;
if (*success)
if ( -1.01  < xi[1] && xi[1] < 1.01) *success = 1;
else                                 *success = 0;
/*--------------------------------------------------- check convergence */
/*------------------------------ we demand convergence near the element */
if (l2 > EPS12)
if (-2.0  < xi[0] && xi[0] < 2.0 &&
    -2.0  < xi[1] && xi[1] < 2.0)
{
/*   printf("Local projection not converged : L2 %20.10E\n",l2);
   printf("xi[0] %f xi[1] %f try %d\n",xi[0],xi[1],second);
   fflush(stdout);*/
   second++;
   if (second<5) goto start;
/*   printf("Local projection finally failed: L2 %20.10E\n",l2);
   printf("try another starting point\n");
   fflush(stdout);*/
   if (second<10) 
   {
      xi[0] = 0.5;
      xi[1] = 0.5;
      goto start;
   }
   printf("Final failure                : L2 %20.10E\n",l2);
   printf("xi[0] %f xi[1] %f try %d\n",xi[0],xi[1],second);
   fflush(stdout);
}
/*-------------------------- give the coordinates up to calling routine */
xires[0] = xi[0];
xires[1] = xi[1];
/*------------------------------------------------------- make distance */
if (*success)
{
   for (i=0; i<3; i++)
   {
      xm[i] = 0.0;
      for (k=0; k<iel; k++)
         xm[i] += funct[k] * (xc[i][k] + thetam * a3c[i][k]);
   }
   *distance = (xs[0]-xm[0])*(xs[0]-xm[0])+(xs[1]-xm[1])*(xs[1]-xm[1])+(xs[2]-xm[2])*(xs[2]-xm[2]);
   *distance = sqrt(*distance);
}
else
   *distance = VERYLARGEREAL;
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8_contact_orthproject */
 



/*!---------------------------------------------------------------------
\brief evaluate shape functions at a given point                                              

<pre>                                                        m.gee 2/03 
</pre>
\param funct    double*          (o)   values of shape functions at xi
\param deriv    double**         (o)   derivatives of shape functions at xi
\param deriv2   double**         (o)   2nd derivatives of shape functions at xi
\param r        double           (i)   xi
\param s        double           (i)   eta 
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_functderiv(double     funct[], 
                           double    deriv[][4], 
                           double    deriv2[],
                           double      r, 
                           double      s)
{
int            i;
const double   q12 = 1.0/2.0;
const double   q14 = 1.0/4.0;
const double   q16 = 1.0/6.0;
double         rp,rm,sp,sm;
#ifdef DEBUG 
dstrc_enter("s8_contact_functderiv");
#endif
/*----------------------------------------------------------------------*/
rp = 1.0+r;
rm = 1.0-r;
sp = 1.0+s;
sm = 1.0-s;
funct[0] = q14*rp*sp;
funct[1] = q14*rm*sp;
funct[2] = q14*rm*sm;
funct[3] = q14*rp*sm;
deriv[0][0]= q14*sp;
deriv[0][1]=-q14*sp;
deriv[0][2]=-q14*sm;
deriv[0][3]= q14*sm;
deriv[1][0]= q14*rp;
deriv[1][1]= q14*rm;
deriv[1][2]=-q14*rm;
deriv[1][3]=-q14*rp;
deriv2[0]= q14;
deriv2[1]=-q14;
deriv2[2]= q14;
deriv2[3]=-q14;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_contact_functderiv */

/*!---------------------------------------------------------------------
\brief evaluate allmetrics at a given point                                              

<pre>                                                        m.gee 2/03 
</pre>
\param funct    double*          (o)   values of shape functions at xi
\param deriv    double**         (o)   derivatives of shape functions at xi
\param deriv2   double**         (o)   2nd derivatives of shape functions at xi
\param r        double           (i)   xi
\param s        double           (i)   eta 
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_metrics(double x[][4],
                        double a3[][4],
                        double e3,
                        double gkov[][3],
                        double gkon[][3],
                        double gmkov[][3],
                        double gmkon[][3],
                        double funct[],
                        double deriv[][4],
                        int    iel)
{
int            i,j,k,idim,ialpha,inode;
double         det;
#ifdef DEBUG 
dstrc_enter("s8_contact_metrics");
#endif
/*----------------------------------------------------------------------*/
/*------------------------------------ interpolation of kovariant g1,g2 */
for (ialpha=0; ialpha<2; ialpha++)
{
   for (idim=0; idim<3; idim++)
   {
      gkov[idim][ialpha] = 0.0;
      for (inode=0; inode<iel; inode++)
      {
         gkov[idim][ialpha] += deriv[ialpha][inode] * (x[idim][inode]+e3*a3[idim][inode]);
      }
   }
}
/*------------------------------------------------- interpolation of g3 */
/* in the shell formultation g3 is the director */
#if 0
for (idim=0; idim<3; idim++)
{
   gkov[idim][2]=0.0;
   for (inode=0; inode<iel; inode++)
   {
      gkov[idim][2] +=
      
         funct[inode] * a3[idim][inode];
   }
}
#endif
/* in contact g3 should be orthogonal to g1,g2 */
gkov[0][2] = gkov[1][0]*gkov[2][1] - gkov[2][0]*gkov[1][1];
gkov[1][2] = gkov[2][0]*gkov[0][1] - gkov[0][0]*gkov[2][1];
gkov[2][2] = gkov[0][0]*gkov[1][1] - gkov[1][0]*gkov[0][1];
/*--------- kontravariant basis vectors g1,g2,g3 (inverse of kovariant) */
for (i=0; i<3; i++)
for (j=0; j<3; j++) gkon[i][j] = gkov[i][j];
s8_contact_inv3(gkon,&det);
s8_contact_trans(gkon,3);
/*--------------------------------------- kovariant metric tensor gmkov */
for (i=0; i<3; i++)
{
   for (j=i; j<3; j++)
   {
      gmkov[i][j]=0.0;
      for (k=0; k<3; k++) 
      gmkov[i][j] += gkov[k][i]*gkov[k][j];
   }
}   
      gmkov[1][0] = gmkov[0][1];
      gmkov[2][0] = gmkov[0][2];
      gmkov[2][1] = gmkov[1][2];
/*----------------------------------------- kontravariant metric tensor */
for (i=0; i<3; i++)
for (j=0; j<3; j++) gmkon[i][j] = gmkov[i][j];
s8_contact_inv3(gmkon,&det);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_contact_metrics */



/*!---------------------------------------------------------------------
\brief invert 3x3 matrix                                              

<pre>                                                        m.gee 2/03 
</pre>
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_inv3(double a[][3], double *det)
{
int i,j;
double b00,b01,b02,b10,b11,b12,b20,b21,b22;
double detinv;
#ifdef DEBUG 
dstrc_enter("s8_contact_inv3");
#endif
/*----------------------------------------------------------------------*/
/* the array to be inved here has to be of the following style:

   a[0]    = [ hole array as a vector ]
   a[1..n] = ptr to start of this row in the vector above
*/   
b00 = a[0][0];
b01 = a[0][1];
b02 = a[0][2];
b10 = a[1][0];
b11 = a[1][1];
b12 = a[1][2];
b20 = a[2][0];
b21 = a[2][1];
b22 = a[2][2];

a[0][0] =   b11*b22 - b21*b12;
a[1][0] = - b10*b22 + b20*b12;
a[2][0] =   b10*b21 - b20*b11;
a[0][1] = - b01*b22 + b21*b02;
a[1][1] =   b00*b22 - b20*b02;
a[2][1] = - b00*b21 + b20*b01;
a[0][2] =   b01*b12 - b11*b02;
a[1][2] = - b00*b12 + b10*b02;
a[2][2] =   b00*b11 - b10*b01;

*det = b00*a[0][0]+b01*a[1][0]+b02*a[2][0];
detinv = 1.0/(*det);
for (i=0; i<3; i++)
for (j=0; j<3; j++) a[i][j] *= detinv;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_contact_inv3 */
/*!---------------------------------------------------------------------
\brief transpose matrix                                              

<pre>                                                        m.gee 2/03 
</pre>
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_trans(double a[][3], int n)
{
int i,j;
double change;
#ifdef DEBUG 
dstrc_enter("s8_contact_trans");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<n; i++)
{
   for (j=i+1; j<n; j++)
   {
      change = a[j][i];
      a[j][i] = a[i][j];
      a[i][j] = change;
   }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_contact_trans */
/*!---------------------------------------------------------------------
\brief evaluate da                                              

<pre>                                                        m.gee 2/03 
</pre>
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_deta(double gkov[][3], double *deta)
{
int i,j;
double change;
#ifdef DEBUG 
dstrc_enter("s8_contact_deta");
#endif
/*----------------------------------------------------------------------*/
*deta = (gkov[0][0]*gkov[1][1] - gkov[1][0]*gkov[0][1])*(gkov[0][0]*gkov[1][1] - gkov[1][0]*gkov[0][1])
      + (gkov[2][0]*gkov[0][1] - gkov[2][1]*gkov[0][0])*(gkov[2][0]*gkov[0][1] - gkov[2][1]*gkov[0][0]) 
      + (gkov[1][0]*gkov[2][1] - gkov[2][0]*gkov[1][1])*(gkov[1][0]*gkov[2][1] - gkov[2][0]*gkov[1][1]);
*deta = sqrt(*deta);
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of s8_contact_deta */



 


#endif
/*! @} (documentation module close)*/
