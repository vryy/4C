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
\brief calculate contact stiffness and forces                                              

<pre>                                                        m.gee 2/03 
calculate contact stiffness and forces 
</pre>
\param  actcnode    SHELLNODE*     (i/o)  the contacting slave node
\param  actele      ELEMENT*       (i)    the element that is contacted
\param  xi          double*        (i)    local coorinates of projection point
\param  ssurf       int            (i)    indicates whether top or bottom of slave contated
\param  msurf       int            (i)    indicates, whether top or bottom of element is contacted
\return void                                               

------------------------------------------------------------------------*/
void s8_contact_make(SHELLNODE *actcnode,
                     ELEMENT   *actele,
                     double    *xi,
                     int        ssurf,
                     int        msurf)
{
int              i,j,k,l,m,n;
int              numdf;
double           sum;

double           tn;

double           thetas,thetam;

double          *force;
double         **stiff;
int             *lm;
ELEMENT         *ele;

double           funct[4];
double           deriv[2][4];
double           deriv2[4];

double           xinode[2];
double           h2;
double           xr[3][4];
double           a3r[3][4];
double           xc[3][4];
double           a3c[3][4];

double           gkov[3][3];
double           gkon[3][3];
double           gmkov[3][3];
double           gmkon[3][3];

double           gkovc[3][3];
double           gkonc[3][3];
double           gmkovc[3][3];
double           gmkonc[3][3];
double           gkovcab[3];


double           nue[3];

double           deta,detas,wgt;

double           xs[3];
double           xbar[3];
double           diff[3];
double           g;
double           N[30];
double           N1[30];
double           N2[30];
double           T1[30];
double           T2[30];
double           D1[30];
double           D2[30];
double           A[2][2],detA;
double           Nbar1[30];
double           Nbar2[30];

double           theta[30];

/* for friction */
int              friction = 0;
int              slip = 0;
double           tT_kov[2];                /* frictional tractions (kovariant components) */
double           tTtrial_kov[2];           /* trial frictional tractions (kovariant components) */
double           tTtrial_kon[2];           /* trial frictional tractions (kontravar. components) */
double           tTabs;                    /* norm of the trial frictional tractions (same in kov and kontra) */
double           lT[2];
double           t[3];
double           PHI;
double           dxi[2];
double           kappa1;
double           kappa2;
double           p_kov[2];
double           p_kon[2];
double           p[3];
double           Akon[2][2];

double           P1[30];
double           P2[30];
double           Pbar1[30];
double           Pbar2[30];
double           T11[30];
double           T12[30];
double           T21[30];
double           T22[30];
double           Tbar11[30];
double           Tbar12[30];
double           Tbar21[30];
double           Tbar22[30];
double           N12[30];

double           gkovr[3][3];
double           gkonr[3][3];
double           gmkovr[3][3];
double           gmkonr[3][3];
double           gkovrab[3];

double           his_xi[2];
ELEMENT         *his_ele;
SHELLCONTACTFLAG his_flag;

double           fac1,fac2,fac3;
double           p1p1;
double           p1p2;
double           p2p1;
double           p2p2;
double           M11;
double           M12;
double           M21;
double           M22;

double           slipstiff[30][30];

double           xctmp[6];
double           xiout[2];
double           distance;
int              success;

#ifdef DEBUG 
dstrc_enter("s8_contact_make");
#endif
/*----------------------------------------------------------------------*/
thetas = (double)ssurf;
thetam = (double)msurf;
/*------------------------------------------------ make location matrix */
numdf = actcnode->node->numdf;
for (i=0; i<actele->numnp; i++)
   numdf += actele->node[i]->numdf;
/*------------------------------------------ allocate vector and matrix */
if (ssurf==1)
{
   force = amdef("cforce",&(actcnode->forcetop),numdf,1,"DV");
   stiff = amdef("cstiff",&(actcnode->stifftop),numdf,numdf,"DA");
           amzero(&(actcnode->stifftop));
   lm    = amdef("clm"   ,&(actcnode->lmtop)   ,numdf,1,"IV");
}
if (ssurf==-1)
{
   force = amdef("cforce",&(actcnode->forcebot),numdf,1,"DV");
   stiff = amdef("cstiff",&(actcnode->stiffbot),numdf,numdf,"DA");
           amzero(&(actcnode->stiffbot));
   lm    = amdef("clm"   ,&(actcnode->lmbot)   ,numdf,1,"IV");
}
/*------------------------------------------------ fill location matrix */
j=0;
for (i=0; i<actcnode->node->numdf; i++)
   lm[j++] = actcnode->node->dof[i];
for (i=0; i<actele->numnp; i++)
for (k=0; k< actele->node[i]->numdf; k++)
   lm[j++] = actele->node[i]->dof[k];
dsassert(j==numdf,"Number of dofs wrong");
/*------------ get the jacobians from the elements adjacent to actcnode */
wgt = 0.0;
for (i=0; i<actcnode->node->numele; i++)
{
   ele = actcnode->node->element[i];
   for (k=0; k<ele->numnp; k++)
   {
      h2 = ele->e.s8->thick_node.a.dv[k] / 2.0;
      a3r[0][k] = ele->e.s8->a3ref.a.da[0][k] * h2;
      a3r[1][k] = ele->e.s8->a3ref.a.da[1][k] * h2;
      a3r[2][k] = ele->e.s8->a3ref.a.da[2][k] * h2;

      xr[0][k]  = ele->node[k]->x[0];   
      xr[1][k]  = ele->node[k]->x[1];   
      xr[2][k]  = ele->node[k]->x[2];   
   }
   /* find which node is actcnode->node */
   for (j=0; j<ele->numnp; j++)
   if (ele->node[j] == actcnode->node) break;
   dsassert(j < 4,"Cannot find matching node on element");
   if (j==0)
   {
      xinode[0] = 1.0;
      xinode[1] = 1.0;
   }   
   if (j==1)
   {
      xinode[0] = -1.0;
      xinode[1] =  1.0;
   }   
   if (j==2)
   {
      xinode[0] = -1.0;
      xinode[1] = -1.0;
   }   
   if (j==3)
   {
      xinode[0] =  1.0;
      xinode[1] = -1.0;
   }
   /*--------------------------- evaluate shape functions at this point */
   s8_contact_functderiv(funct,deriv,deriv2,xinode[0],xinode[1]);
   /*-------------- get metrics of this point on the contacting surface */
   s8_contact_metrics(xr,a3r,thetas,gkov,gkon,gmkov,gmkon,funct,deriv,ele->numnp);
   s8_contact_deta(gkov,&detas);
   /*--------------------------------------------- add to summed weight */
   wgt += detas;
}
/*-------------------- init the master element and the projection point */
for (k=0; k<actele->numnp; k++)
{
   h2        = actele->e.s8->thick_node.a.dv[k] / 2.0;
   
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
/*----------------- evaluate shape functions at the projection point xi */
s8_contact_functderiv(funct,deriv,deriv2,xi[0],xi[1]);
/*----------------------- build current metrics at the projection point */
s8_contact_metrics(xc,a3c,thetam,gkovc,gkonc,gmkovc,gmkonc,funct,deriv,actele->numnp);
/*---------------------------- build outward normal at projection point */
nue[0] = gkovc[0][2];
nue[1] = gkovc[1][2];
nue[2] = gkovc[2][2];
math_unvc(&detas,nue,3);
/*--------- if it is bottom surface, it has to be switched in direction */
if (msurf==-1)
for (i=0; i<3; i++) nue[i] = -nue[i];
/*------------------------------------- build value of the gap function */
/*--------------------------------------- make global slave coordinates */
xs[0] = actcnode->xc[0] + thetas * actcnode->xc[3];
xs[1] = actcnode->xc[1] + thetas * actcnode->xc[4];
xs[2] = actcnode->xc[2] + thetas * actcnode->xc[5];
/*---------------------------------------- make global projection point */
for (i=0; i<3; i++) xbar[i] = 0.0;
for (k=0; k<actele->numnp; k++) 
   for (i=0; i<3; i++)
      xbar[i] += funct[k] * (xc[i][k] + thetam*a3c[i][k]);
/*-------------------- make vector from projection point to slave point */     
for (i=0; i<3; i++)
   diff[i] = xs[i] - xbar[i];
/*----------------------------------------------- project diff onto nue */   
g = 0.0;
for (i=0; i<3; i++) g += diff[i] * nue[i];
g = -g;
/*--------------------- build vectors needed for forces and stiffnesses */
/*------------------------------------------------------ build vector N */
j=0;
for (l=0; l<2; l++)
for (i=0; i<3; i++) 
   N[j++] = nue[i]; 
for (k=0; k<actele->numnp; k++)
for (l=0; l<2; l++)
for (i=0; i<3; i++)
   N[j++] = -funct[k] * nue[i];
dsassert(j == numdf,"Number of dofs wrong");
/*---------------------------------------------- build vector N1 and N2 */
j=0;
for (l=0; l<2; l++)
for (i=0; i<3; i++)
{
   N1[j] = 0.0;
   N2[j] = 0.0;
   j++;
}
for (k=0; k<actele->numnp; k++)
for (l=0; l<2; l++)
for (i=0; i<3; i++)
{
   N1[j] = -deriv[0][k] * nue[i];
   N2[j] = -deriv[1][k] * nue[i];
   j++;
}
dsassert(j == numdf,"Number of dofs wrong");
/*--------------------------------------------- build vectors T1 and T2 */
j=0;
for (l=0; l<2; l++)
for (i=0; i<3; i++)
{
   T1[j] = gkovc[i][0];
   T2[j] = gkovc[i][1];
   j++;
}
for (k=0; k<actele->numnp; k++)
for (l=0; l<2; l++)
for (i=0; i<3; i++)
{
   T1[j] = -funct[k]*gkovc[i][0];
   T2[j] = -funct[k]*gkovc[i][1];
   j++;
}
dsassert(j == numdf,"Number of dofs wrong");
/*---------- make derivatives of basis vectors g_alphabeta = g_betalpha */
/*--------------------------------------- g_alphaalpha = g_betabeta = 0 */
for (i=0; i<3; i++)
{
   gkovcab[i] = 0.0;
   for (k=0; k<actele->numnp; k++)
      gkovcab[i] += deriv2[k] * (xc[i][k] + thetam * a3c[i][k]);
}
/*------------------------------------- make metrik A[2][2] (kovariant) */
A[0][0] = gmkovc[0][0];
A[0][1] = gmkovc[0][1];
A[1][0] = gmkovc[1][0];
A[1][1] = gmkovc[1][1];
sum = 0.0;
for (i=0; i<3; i++) 
   sum += gkovcab[i] * nue[i];
A[0][1] += g*sum;
A[1][0] += g*sum;
detA = A[0][0]*A[1][1] - A[0][1]*A[1][0];
detA = 1.0/detA;
/*------------------------------ make metrik Akon[2][2] (kontravariant) */
Akon[0][0] =  detA * A[1][1];
Akon[0][1] = -detA * A[0][1];
Akon[1][0] = -detA * A[1][0];
Akon[1][1] =  detA * A[0][0];
/*--------------------------------------------- build vectors D1 and D2 */
for (i=0; i<numdf; i++)
{
   D1[i] = detA*(  A[1][1]*(T1[i]+g*N1[i]) - A[0][1]*(T2[i]+g*N2[i]) );
   D2[i] = detA*( -A[1][0]*(T1[i]+g*N1[i]) + A[1][1]*(T2[i]+g*N2[i]) );
}
/*--------------------------------------- build vectors Nbar1 and Nbar2 */
for (i=0; i<numdf; i++)
{
  Nbar1[i] = N1[i] - sum * D2[i];
  Nbar2[i] = N2[i] - sum * D1[i];
}
/*------------------------------------------- make normal contact force */
if (ssurf==1) 
{
   tn = actcnode->top_ln + EPSN * g;
   actcnode->top_tn = tn;
}
else
{
   tn = actcnode->bot_ln + EPSN * g;
   actcnode->bot_tn = tn;
}
/*--------------------------------------- make vector of contact forces */
for (i=0; i<numdf; i++)
{
   force[i] =  - wgt * tn * N[i];
}
/*-------------------------------------------- make contact stiffnesses */
for (i=0; i<numdf; i++)
for (j=0; j<numdf; j++)
{
   stiff[i][j] += wgt * EPSN * N[i] * N[j];
   stiff[i][j] += wgt * tn * g * ( gmkonc[0][0]*Nbar1[i]*Nbar1[j] + gmkonc[0][1]*(Nbar1[i]*Nbar2[j] + Nbar2[i]*Nbar1[j]) + gmkonc[1][1]*Nbar2[i]*Nbar2[j]);
   stiff[i][j] -= wgt * tn * (D1[i]*N1[j] + D2[i]*N2[j] + N1[i]*D1[j] + N2[i]*D2[j]);
   stiff[i][j] += wgt * tn * sum * (D1[i]*D2[j]+D2[i]*D1[j]);
}
/*==============================================build frictional forces */
#if 1
/*------------------------------------------------------------ make dxi */
if (ssurf==1)
{
   his_xi[0]     = actcnode->hisxitop[0];
   his_xi[1]     = actcnode->hisxitop[1];
   his_ele       = actcnode->histopele;
   his_flag      = actcnode->histopflag;
   lT[0]         = actcnode->top_lt[0];
   lT[1]         = actcnode->top_lt[1];
}
else
{
   his_xi[0]     = actcnode->hisxibot[0];
   his_xi[1]     = actcnode->hisxibot[1];
   his_ele       = actcnode->hisbotele;
   his_flag      = actcnode->hisbotflag;
   lT[0]         = actcnode->bot_lt[0];
   lT[1]         = actcnode->bot_lt[1];
}
/* check whether we had contact here the last time step, if so we can slide */
if (his_flag == s8_c_on)
{
   friction = 1;
   if (actele == his_ele) /* we are slippping inside the same element */
   {
      dxi[0] = xi[0] - his_xi[0];
      dxi[1] = xi[1] - his_xi[1];
   }
   else /* we are slipping, but have been in another element last time */
   {
      printf("next element!\n");
      /* store the current coodinates of actnode temporary */
      for (j=0; j<6; j++) xctmp[j] = actcnode->xc[j];
      /* put the current coodinates of the last converged step in place */
      for (j=0; j<6; j++) actcnode->xc[j] = actcnode->xc_his[j];
      /* now the nodes has coodinates of last time step, project this onto current element */
      /* this projection will NOT be inside the element, but will give sufficient xi to do friction */
      s8_contact_orthproject(actcnode,&ssurf,&msurf,actele,xiout,&distance,&success);
      /* put the current coodinates back into place */
      for (j=0; j<6; j++) actcnode->xc[j] = xctmp[j];
      /* make delta xi */
      dxi[0] = xi[0] - xiout[0];
      dxi[1] = xi[1] - xiout[1];
   }
}
else
{
   friction = 0;
}
/*------------------------------------------------- in case of friction */
if (friction)
{
/*------------------------- build reference metrics at projection point */
s8_contact_metrics(xr,a3r,thetam,gkovr,gkonr,gmkovr,gmkonr,funct,deriv,actele->numnp);
/* frictional forces in directions T1 and T2 */
tTtrial_kov[0] = lT[0] + EPST * ( gmkovr[0][0] * dxi[0] + gmkovr[0][1] * dxi[1] );
tTtrial_kov[1] = lT[1] + EPST * ( gmkovr[1][0] * dxi[0] + gmkovr[1][1] * dxi[1] );
/* make absolute value of frictional force wit respect to current base vectors */
/* warning: this may be in terms of gmkonr metrics */
tTabs = tTtrial_kov[0]*gmkonc[0][0]*tTtrial_kov[0] + 
        tTtrial_kov[1]*gmkonc[1][1]*tTtrial_kov[1] +
        tTtrial_kov[0]*gmkonc[0][1]*tTtrial_kov[1] +
        tTtrial_kov[1]*gmkonc[1][0]*tTtrial_kov[0] ;
tTabs = sqrt(tTabs);  
/* tTtrial_kov are kovariant components in kontravariant basis */
/* make tTtrial_kon kontravariant components in covariant basis */
tTtrial_kon[0] = tTtrial_kov[0]*gmkonc[0][0] + tTtrial_kov[1]*gmkonc[1][0];
tTtrial_kon[1] = tTtrial_kov[0]*gmkonc[0][1] + tTtrial_kov[1]*gmkonc[1][1];
/* test the frictional traction vector t in kov and kon components should be equal */
/*
for (i=0; i<3; i++)
   t[i] = tTtrial_kov[0]*gkonc[i][0] + tTtrial_kov[1]*gkonc[i][1];
for (i=0; i<3; i++)
   t[i] = tTtrial_kon[0]*gkovc[i][0] + tTtrial_kon[1]*gkovc[i][1];
*/
/*-------------------------------------------------- make coulomb's law */
PHI = tTabs - NU * tn;
/*-------------------------------------------- decide for slip or stick */
if (PHI < 0.0) 
slip = 0;
else           
slip = 1;
/*--------------------------------- make return of frictional tractions */
if (slip)
{
   tT_kov[0] = NU * tn * tTtrial_kov[0]/tTabs;
   tT_kov[1] = NU * tn * tTtrial_kov[1]/tTabs;
}
else
{
   tT_kov[0] = tTtrial_kov[0];
   tT_kov[1] = tTtrial_kov[1];
}
/*------------------------------- put frictional tractions back to node */
/* from there they will be copied to the history in the converged status */  
if (ssurf==1)
{
   actcnode->top_tT[0] = tT_kov[0];
   actcnode->top_tT[1] = tT_kov[1];
}
else
{
   actcnode->bot_tT[0] = tT_kov[0];
   actcnode->bot_tT[1] = tT_kov[1];
}
/*---------------------------------------------------- make slip forces */
for (i=0; i<numdf; i++)
{
   force[i] += wgt * tT_kov[0]*D1[i] + tT_kov[1]*D2[i];
}
/*---------- make derivatives of basis vectors g_alphabeta = g_betalpha */
/*--------------------------------------- g_alphaalpha = g_betabeta = 0 */
/*------------------------------------------ in reference configuration */
for (i=0; i<3; i++)
{
   gkovrab[i] = 0.0;
   for (k=0; k<actele->numnp; k++)
      gkovrab[i] += deriv2[k] * (xr[i][k] + thetam * a3r[i][k]);
}
/*---------------------------------------------- make kappa1 and kappa2 */
/* kappa1 = g1,2 * g1              */
/* kappa2 = g2,1 * g2 (g1,2=g2,1)  */
kappa1 = 0.0;
kappa2 = 0.0;
for (i=0; i<3; i++)
{
   kappa1 += gkovrab[i] * gkovr[i][0];
   kappa2 += gkovrab[i] * gkovr[i][1];
}
/*-------------------------------------------------------------- make p */
p_kov[0] = tTtrial_kov[0]/tTabs;
p_kov[1] = tTtrial_kov[1]/tTabs;
p_kon[0] = tTtrial_kon[0]/tTabs;
p_kon[1] = tTtrial_kon[1]/tTabs;
for (i=0; i<3; i++)
   p[i] = p_kov[0]*gkonc[i][0] + p_kov[1]*gkonc[i][1];
/*----------------------------------- make vectors P1, P2, Pbar1, Pbar2 */   
sum = 0.0;
for (i=0; i<3; i++)
   sum += gkovcab[i] * p[i];
j=0;
for (l=0; l<2; l++)
for (i=0; i<3; i++)
{
   P1[j]    = 0.0;
   P2[j]    = 0.0;
   Pbar1[j] = P1[j] - sum * D2[j];
   Pbar2[j] = P2[j] - sum * D1[j];
   j++;
}
for (k=0; k<actele->numnp; k++)
for (l=0; l<2; l++)
for (i=0; i<3; i++)
{
   P1[j]    = -deriv[0][k] * p[i];
   P2[j]    = -deriv[1][k] * p[i];
   Pbar1[j] = P1[j] - sum * D2[j];
   Pbar2[j] = P2[j] - sum * D1[j];
   j++;
}
/*---------------------------------------- make vectors T11 T12 T21 T22 */
j=0;
for (l=0; l<2; l++)
for (i=0; i<3; i++)
{
   T11[j]   = 0.0;
   T12[j]   = 0.0;
   T21[j]   = 0.0;
   T22[j]   = 0.0;
   j++;
}
for (k=0; k<actele->numnp; k++)
for (l=0; l<2; l++)
for (i=0; i<3; i++)
{
   T11[j]   = -deriv[0][k] * gkovc[i][0]; 
   T12[j]   = -deriv[1][k] * gkovc[i][0];
   T21[j]   = -deriv[0][k] * gkovc[i][1];
   T22[j]   = -deriv[1][k] *  gkovc[i][1];
   j++;
}
/*---------------------------- make vectors Tbar11 Tbar21 Tbar12 Tbar22 */
/* Tbar11 = T11 - gkovcab*gkovc1 * D2 */
sum = 0;
for (i=0; i<3; i++) sum += gkovcab[i]*gkovc[i][0];
for (j=0; j<numdf; j++)
   Tbar11[j] = T11[j] - sum * D2[j];
/* Tbar12 = T12 - gkovcab*gkovc1 * D1 */
for (j=0; j<numdf; j++)
   Tbar12[j] = T12[j] - sum * D1[j];
/* Tbar21 = T21 - gkovcab*gkovc2 * D2 */
sum = 0.0;
for (i=0; i<3; i++) sum += gkovcab[i]*gkovc[i][1];
for (j=0; j<numdf; j++)
   Tbar21[j] = T21[j] - sum * D2[j];
/* Tbar22 = T22 - gkovcab*gkovc2 * D1 */
for (j=0; j<numdf; j++)
   Tbar22[j] = T22[j] - sum * D1[j];
/*-------------------------------- make vectors N12 = N21 (N11=N22=0) */
j=0;
for (l=0; l<2; l++)
for (i=0; i<3; i++)
{
   N12[j]   = 0.0;
   j++;
}
for (k=0; k<actele->numnp; k++)
for (l=0; l<2; l++)
for (i=0; i<3; i++)
{
   N12[j]   = -deriv2[k] * nue[i];
   j++;
}
/*----------------------------------------------------------------------*/
for (i=0; i<numdf; i++)
for (j=0; j<numdf; j++) slipstiff[i][j] = 0.0;
/*------------------------------------------------------ some constants */
p1p1 = p_kon[0]*p_kov[0];
p1p2 = p_kon[0]*p_kov[1];
p2p1 = p_kon[1]*p_kov[0];
p2p2 = p_kon[1]*p_kov[1];

M11  = gmkovr[0][0];
M12  = gmkovr[0][1];
M21  = gmkovr[1][0];
M22  = gmkovr[1][1];
/*----------------------------------- build frictional stiffness matrix */
if (slip)
{
   fac1 =   wgt * NU * EPSN;
   fac2 =   wgt * EPST * NU * tn * (1.0/tTabs);
   fac3 = - wgt * NU * tn;
   for (i=0; i<numdf; i++)
   for (j=0; j<numdf; j++)
   {
      /* expression 1 */
      slipstiff[i][j] -= 
      fac1 * 
      (
      p_kov[0]*D1[i] + p_kov[1]*D2[i]
      ) * N[j];
      /* expression 2 */
      slipstiff[i][j] += 
      fac2 *
      (
      D1[i]*D1[j]*( (1.0-p1p1)*(M11+kappa1*dxi[1]                  )-p2p1*(M21+kappa1*dxi[0]+2.0*kappa2*dxi[1]) )+
      D2[i]*D1[j]*( (1.0-p2p2)*(M21+kappa1*dxi[0]+2.0*kappa2*dxi[1])-p1p2*(M11+kappa1*dxi[1]                  ) )+
      D1[i]*D2[j]*( (1.0-p1p1)*(M12+2.0*kappa1*dxi[0]+kappa2*dxi[1])-p2p1*(M22+kappa2*dxi[0]                  ) )+
      D2[i]*D2[j]*( (1.0-p2p2)*(M22+kappa2*dxi[0]                  )-p1p2*(M12+2.0*kappa1*dxi[0]+kappa2*dxi[1]) )
      );
      /* expression 3 */                  
      slipstiff[i][j] += 
      fac3 *
      (
      p1p1*D1[i]*Pbar1[j] +
      p1p2*D2[i]*Pbar1[j] +
      p2p1*D1[i]*Pbar2[j] +
      p2p2*D2[i]*Pbar2[j]
      );
   }
}
else /* stick */
{
   fac1 = wgt * EPST;
   for (i=0; i<numdf; i++)
   for (j=0; j<numdf; j++)
   {
      slipstiff[i][j] += 
      fac1 *
      ( 
      gmkovr[0][0]*D1[i]*D1[j] + 
      gmkovr[0][1]*D1[i]*D2[j] + gmkovr[1][0]*D2[i]*D1[j] +
      gmkovr[1][1]*D2[i]*D2[j] +
      2.0*kappa1*dxi[0]*D1[i]*D2[j] +
      kappa1*dxi[1]*D1[i]*D1[j] + kappa2*dxi[1]*D1[i]*D2[j] +
      kappa1*dxi[0]*D2[i]*D1[j] + kappa2*dxi[0]*D2[i]*D2[j] +
      2.0*kappa2*dxi[1]*D2[i]*D1[j] 
      );
   }
}
/* make part of frictional stiffness which is identical for slip and stick */
fac1 = tT_kov[0]*Akon[0][0] + tT_kov[1]*Akon[0][1];
sum  = 0.0;
for (i=0; i<3; i++)
   sum += gkovcab[i]*gkovc[i][0];
for (i=0; i<numdf; i++)
for (j=0; j<numdf; j++)
{
   slipstiff[i][j] += wgt * fac1 *
                      (
                        T11[i]*D1[j] + T12[i]*D2[j] + D1[i]*T11[j] + D2[i]*T12[j] 
                      - sum * (D1[i]*D2[j] + D2[i]*D1[j])
                      + Tbar11[i]*D1[j] + Tbar21[i]*D2[j]
                      + D1[i]*Tbar11[j] + D2[i]*Tbar21[j]
                      + g * (N12[i]*D2[j] + D2[i]*N12[j])
                      - N[i]*Nbar1[j] 
                      - T1[i] * ( gmkonc[0][0]*Tbar11[j] + gmkonc[0][1]*Tbar21[j] )
                      - T2[i] * ( gmkonc[1][0]*Tbar11[j] + gmkonc[1][1]*Tbar21[j] )
                      - Nbar1[i]*N[j]
                      - ( gmkonc[0][0]*Tbar11[i] + gmkonc[0][1]*Tbar12[i] ) * T1[j]
                      - ( gmkonc[1][0]*Tbar11[i] + gmkonc[1][1]*Tbar21[i] ) * T2[j]
                      );
}

fac2 = tT_kov[0]*Akon[0][1] + tT_kov[1]*Akon[1][1]; 
sum  = 0.0;
for (i=0; i<3; i++)
   sum += gkovcab[i]*gkovc[i][0];
for (i=0; i<numdf; i++)
for (j=0; j<numdf; j++)
{
   slipstiff[i][j] += wgt * fac2 *
                      (
                        T21[i]*D1[j] + T22[i]*D2[j] + D1[i]*T21[j] + D2[i]*T22[j]
                      - sum * (D1[i]*D2[j] + D2[i]*D1[j])
                      + Tbar12[i]*D1[j] + Tbar22[i]*D2[j]
                      + D1[i]*Tbar12[j] + D2[i]*Tbar22[j]
                      + g * (N12[i]*D1[j] + D1[i]*N12[j])
                      - N[i]*Nbar2[j]
                      - T1[i] * ( gmkonc[0][0]*Tbar12[j] + gmkonc[0][1]*Tbar22[j] )
                      - T2[i] * ( gmkonc[1][0]*Tbar12[j] + gmkonc[1][1]*Tbar22[j] )
                      - Nbar2[i]*N[j]
                      - ( gmkonc[0][0]*Tbar12[i] + gmkonc[0][1]*Tbar22[i] ) * T1[j]
                      - ( gmkonc[1][0]*Tbar12[i] + gmkonc[1][1]*Tbar22[i] ) * T2[j]
                      );
}
/*=================== add frictional stiffness to the general stiffness */
for (i=0; i<numdf; i++)
for (j=0; j<numdf; j++)
   stiff[i][j] += slipstiff[i][j];
/*----------------------------------------------------- end of friction */
} /* end of if (friction) */
#endif 
if (ssurf==1)
{
   if (friction)
   printf("node %5d element %5d topbot %2d %2d gap %15.10f lambdan %15.8f tn %15.8f slip %d\n",
           actcnode->node->Id,actele->Id,ssurf,msurf,g,actcnode->top_ln,actcnode->top_tn,slip);
   else
   printf("node %5d element %5d topbot %2d %2d gap %15.10f lambdan %15.8f tn %15.8f\n",
           actcnode->node->Id,actele->Id,ssurf,msurf,g,actcnode->top_ln,actcnode->top_tn);
}
else
{
   if (friction)
   printf("node %5d element %5d topbot %2d %2d gap %15.10f lambdan %15.8f tn %15.8f slip %d\n",
           actcnode->node->Id,actele->Id,ssurf,msurf,g,actcnode->bot_ln,actcnode->bot_tn,slip);
   else
   printf("node %5d element %5d topbot %2d %2d gap %15.10f lambdan %15.8f tn %15.8f\n",
           actcnode->node->Id,actele->Id,ssurf,msurf,g,actcnode->bot_ln,actcnode->bot_tn);
}
fflush(stdout);        
/*========= build the toggle diagonal matrix theta (stored as a vector) */
j=0;
for (i=0; i<3; i++) theta[j++] = 1.0;
for (i=0; i<3; i++) theta[j++] = thetas;
for (k=0; k<actele->numnp; k++)
{
   for (i=0; i<3; i++) theta[j++] = 1.0;
   for (i=0; i<3; i++) theta[j++] = thetam;
}
dsassert(j == numdf,"Number of dofs wrong");
/*-------------------------------------- multiply the forces with theta */
for (i=0; i<numdf; i++) force[i] *= theta[i];
/*---------------------- mutliply the stiffness matrix with theta twice */
for (i=0; i<numdf; i++)
for (j=0; j<numdf; j++)
{
   stiff[i][j] *= theta[i];
}
for (i=0; i<numdf; i++)
for (j=0; j<numdf; j++)
{
   stiff[i][j] *= theta[j];
}
/*----------------------------------------------------------------------*/
end:
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s8_contact_make */


/*! @} (documentation module close)*/
#endif
