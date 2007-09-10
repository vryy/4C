/*!----------------------------------------------------------------------
\file
\brief service routines for fluid2 element

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0711 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_FLUID2TU
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
#include "fluid2.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

static FLUID_DYNAMIC *fdyn;
/*!---------------------------------------------------------------------
\brief set all arrays for element calculation

<pre>                                                         he  12/02


</pre>
\param   *ele      ELEMENT	     (i)    actual element
\param   *elev     ELEMENT	     (i)    actual element for velocity
\param   *kapepsn    DOUBLE	     (o)    kapeps at time n
\param   *kapepsg    DOUBLE	     (o)    kapeps at time n+g
\param   *kapepspro  DOUBLE	     (o)    kapeps at time n+g
\param   *eddyg      DOUBLE	     (o)    eddy-visc. at time n+g
\param   *eddypro    DOUBLE	     (o)    eddy-visc. for prod. term
\param   *kappa      DOUBLE	     (o)    kappa at time n+g
\param   *kappan     DOUBLE	     (o)    kappa at time n
\param   *epsilon    DOUBLE	     (o)    epsilon at time n+g
\param  **evel       DOUBLE	     (o)    element velocities
\param  **xzye       DOUBLE        (o)   nodal coordinates
\param  *ipos                        (i)    node array positions
\return void

------------------------------------------------------------------------*/
void f2_calset_tu(
	            ELEMENT         *ele,
                  ELEMENT         *elev,
                  DOUBLE          *kapepsn,
	            DOUBLE          *kapepsg,
	            DOUBLE          *kapepspro,
                  DOUBLE          *eddyg,
                  DOUBLE          *eddypro,
	            DOUBLE          *kappa,
	            DOUBLE          *kappan,
	            DOUBLE          *epsilon,
	            DOUBLE         **evel,
	            DOUBLE         **xyze,
                    ARRAY_POSITION *ipos
                 )
{
INT i,j;            /* simply a counter                                 */
INT kap_eps;
NODE  *actnode;     /* actual node                                      */
FLUID_DATA      *data;

#ifdef DEBUG
dstrc_enter("f2_calset_tu");
#endif

fdyn = alldyn[genprob.numff].fdyn;
data = fdyn->data;

/*-------------------------------------------- set element coordinates */
for(i=0;i<ele->numnp;i++)
{
   xyze[0][i]=ele->node[i]->x[0];
   xyze[1][i]=ele->node[i]->x[1];
}

/*---------------------------------------------------------------------*
 | position of the different solutions:                                |
 | node->sol_incement: solution history used for calculations          |
 |       sol_increment[0][i]: solution at (n-1)                        |
 |	   sol_increment[1][i]: solution at (n)                          |
 |	   sol_increment[2][i]: solution at (n+g)                        |
 |	   sol_increment[3][i]: solution at (n+1)                        |
 |	                  i=0 : kappa                                    |
 |	                  i=1 : eddy-viscosity                           |
 |	                  i=2 : epsilon                                  |
 |	                  i=3 : charact. lenght                          |
 *---------------------------------------------------------------------*/
 if (fdyn->kapeps_flag==0) kap_eps=0;
 if (fdyn->kapeps_flag==1) kap_eps=2;

   for(i=0;i<ele->numnp;i++) /* loop nodes of element */
   {
     actnode=ele->node[i];

/*----------------------------------- set element kapeps (n+1,i)       */
      kapepsg[i]  =actnode->sol_increment.a.da[ipos->velnp][kap_eps];
/*----------------------------------- set element kapeps (n)           */
      kapepsn[i]  =actnode->sol_increment.a.da[ipos->veln][kap_eps];
/*----------------------------------- set eddy viscosity for timestep  */
      eddyg[i]   =actnode->sol_increment.a.da[ipos->velnp][1];
      eddypro[i] =actnode->sol_increment.a.da[ipos->eddy][1];

/*--------------------- for kappa equation: epsilon is needed for R_t  */
    if (fdyn->kapeps_flag==0)
    {
     epsilon[i] = actnode->sol_increment.a.da[ipos->velnp][2];

     if (fdyn->kappan==2)
     {
      actnode->sol_increment.a.da[ipos->eddy][0] =
                           actnode->sol_increment.a.da[ipos->velnp][0];
     }
    }

/*-------- for epsilon equation: kappan is needed for production term  */
    if (fdyn->kapeps_flag==1)
    {
     kappa[i]  =actnode->sol_increment.a.da[ipos->velnp][0];

     if (fdyn->kappan==2)
     {
      actnode->sol_increment.a.da[ipos->eddy][2] =
                          actnode->sol_increment.a.da[ipos->velnp][2];
      kappan[i] =actnode->sol_increment.a.da[ipos->eddy][0];
     }
    }

/*----------------------------------- set element kapeps (pro)         */
      kapepspro[i]=actnode->sol_increment.a.da[ipos->eddy][kap_eps];

/*----------------------- get velocities calculated form Navier-Stokes */
    for(j=0;j<2;j++)
    {
      actnode=elev->node[i];
      evel[j][i]=actnode->sol_increment.a.da[ipos->velnp][j];
    }

 } /* end of loop over nodes of element */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calset */


/*!---------------------------------------------------------------------
\brief routine to calculate kapeps at integration point

<pre>                                                        he  12/02

</pre>
\param   *kapepsint   DOUBLE        (o)   kapeps at integration point
\param   *funct       DOUBLE        (i)   shape functions
\param   *ekapeps    DOUBLE        (i)   kapeps at element nodes
\param    iel	    INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_kapepsi(
             DOUBLE  *kapepsint,
             DOUBLE  *funct,
	       DOUBLE  *kapeps,
             INT      iel
	     )
{
INT     j;

#ifdef DEBUG
dstrc_enter("f2_kapepsi");
#endif
/*---------------------------------------------------------------------*/

   *kapepsint=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      *kapepsint += funct[j]*kapeps[j];
   } /* end loop over j */


/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_kapepsi */


/*!---------------------------------------------------------------------
\brief routine to calculate eddy-visc. at integration point

<pre>                                                        he  12/02

</pre>
\param   *eddyint     DOUBLE        (o)   eddy at integration point
\param   *funct       DOUBLE        (i)   shape functions
\param   **eddy       DOUBLE        (i)   kapeps at element nodes
\param    iel	    INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_eddyi(
             DOUBLE  *eddyint,
             DOUBLE  *funct,
	       DOUBLE  *eddy,
             INT      iel
	     )
{
INT     j;

#ifdef DEBUG
dstrc_enter("f2_eddyi");
#endif

   *eddyint=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      *eddyint   += funct[j]*eddy[j];
   } /* end loop over j */


/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_eddyi */

/*!---------------------------------------------------------------------
\brief routine to calculate kappa for epsilon at integration point

<pre>                                                        he  12/02

</pre>
\param   *kappaint     DOUBLE       (o)   kappa at integration point
\param   *kappanint    DOUBLE       (o)   kappan at integration point
\param   *eps_proint   DOUBLE       (o)   eps_pro at integration point
\param   *funct        DOUBLE       (i)   shape functions
\param   *kappa        DOUBLE       (i)   kappa from kappa equation
\param   *kappan       DOUBLE       (i)   kappan (start-value) kappa equation
\param   *eps_pro      DOUBLE       (i)   epsilon for production term
\param    iel	    INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_kappai_tu(
             DOUBLE     *kappaint,
             DOUBLE     *kappanint,
             DOUBLE     *eps_proint,
             DOUBLE     *funct,
             DOUBLE     *kappa,
             DOUBLE     *kappan,
             DOUBLE     *eps_pro,
             INT         iel
	     )
{
INT     j;

#ifdef DEBUG
dstrc_enter("f2_kappai_tu");
#endif
/*---------------------------------------------------------------------*/

   *kappaint  =ZERO;
   *kappanint =ZERO;
   *eps_proint=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      *kappaint     += funct[j]*kappa[j];
      *kappanint    += funct[j]*kappan[j];
      *eps_proint   += funct[j]*eps_pro[j];
    } /* end loop over j */


/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_kappai_tu */

/*!---------------------------------------------------------------------
\brief routine to calculate C_u with R_t for LOW-REYNOLD's MODEL

<pre>                                                        he  01/03

</pre>
\param    kapepsint     DOUBLE       (i)   kappa at integration point
\param  **epsilon       DOUBLE       (i)   epsilon at nodes
\param   *funct         DOUBLE       (i)   shape functions
\param    visc	      DOUBLE       (i)   molecular viscosity
\param   *C_u           DOUBLE       (o)   factor
\param    iel	    INT            (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_C_kappa(
             DOUBLE      kapepsint,
             DOUBLE     *epsilon,
             DOUBLE     *funct,
             DOUBLE      visc,
             DOUBLE     *C_u,
             INT         iel
	     )
{
INT     i;
DOUBLE  R_t;           /* turbulent Reynold's number                    */
DOUBLE  epsilonint;

#ifdef DEBUG
dstrc_enter("f2_C_kappa");
#endif
/*---------------------------------------------------------------------*/

epsilonint=ZERO;
for (i=0;i<iel;i++) /* loop over all nodes i of the element */
{
  epsilonint  += funct[i]*epsilon[i];
} /* end loop over j */

 R_t = pow(kapepsint,2) / (epsilonint*visc);

*C_u = 0.09 * exp(-3.4/pow(1+R_t/50,2));
/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_C_kappa */

/*!---------------------------------------------------------------------
\brief routine to calculate C_2 with R_t for LOW-REYNOLD's MODEL

<pre>                                                        he  01/03

</pre>
\param    kapepsint     DOUBLE       (i)   epsilon at integration point
\param    kappaint      DOUBLE       (i)   kappa at integration point
\param    visc	      DOUBLE       (i)   molecular viscosity
\param   *C_2           DOUBLE       (o)   factor
\param    iel	    INT            (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_C_eps(
             DOUBLE      kapepsint,
             DOUBLE      kappaint,
             DOUBLE      visc,
             DOUBLE     *C_2,
             INT         iel
	     )
{
DOUBLE  R_t;           /* turbulent Reynold's number                    */

#ifdef DEBUG
dstrc_enter("f2_C_eps");
#endif
/*---------------------------------------------------------------------*/

 R_t = pow(kappaint,2) / (kapepsint*visc);

*C_2 = 1.92 * (1-0.3*exp(-R_t*R_t));
/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_C_eps */

/*!---------------------------------------------------------------------
\brief routine to calculate (vderxy)^2 for LOW-REYNOLD's MODEL

<pre>                                                        he  01/03

</pre>
\param   **vderxy2     DOUBLE       (i)   2nd deriv. for velocity at INT.
\param    vderxy_12    DOUBLE       (o)   (vderxy)^2 at integration point
\return void

------------------------------------------------------------------------*/
void f2_v(
             DOUBLE    **vderxy2,
             DOUBLE     *vderxy_12
	   )
{
DOUBLE  factor;
INT     i;

#ifdef DEBUG
dstrc_enter("f2_v");
#endif
/*----------------------------------------------------------------------

                [ grad (grad (u)) ] ^2

----------------------------------------------------------------------
  (u1,11)^2 + (u1,21)^2 + (u2,11)^2 + (u2,21)^2  + (u1,12)^2 + (u1,22)^2
+ (u2,12)^2 + (u2,22)^2   =
  (u1,11)^2 + 2*(u1,12)^2  +  (u1,22)^2
+ (u2,22)^2 + 2*(u2,21)^2  +  (u2,11)^2                                */

*vderxy_12 = ZERO;

 for (i=0; i<3; i++)
 {
  if(i==2) factor = 2.0;
  if(i!=2) factor = 1.0;
  *vderxy_12 += factor*(vderxy2[0][i]*vderxy2[0][i]+vderxy2[1][i]*vderxy2[1][i]);
 }

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_v */

/*!---------------------------------------------------------------------
\brief routine to calculate factors for kappa-equation

<pre>                                                        he  12/02

</pre>
\param    C_u         DOUBLE        (i)   factor
\param    eddyint     DOUBLE        (i)   eddy at integration point
\param   *factor      DOUBLE        (o)   factor
\param   *factor1     DOUBLE        (o)   factor
\param   *factor2     DOUBLE        (o)   factor
\param   *sig         DOUBLE        (o)   factor
\return void

------------------------------------------------------------------------*/
void f2_fac_kappa(
                  DOUBLE   C_u,
                  DOUBLE   eddyint,
	            DOUBLE  *factor,
	            DOUBLE  *factor1,
	            DOUBLE  *factor2,
	            DOUBLE  *sig
	           )
{

#ifdef DEBUG
dstrc_enter("f2_fac_kappa");
#endif
/*---------------------------------------------------------------------*/

  *factor =2.0*C_u/eddyint;
  *factor2=    C_u/eddyint;
  *factor1=1.0;
  *sig    =1.0;

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_fac_kappa */

/*!---------------------------------------------------------------------
\brief routine to calculate factors for epsilon-equation

<pre>                                                        he  12/02

</pre>
\param    C_2        DOUBLE        (i)   factor
\param    eps_proint DOUBLE        (i)   epsilon at integration point
\param    kappaint   DOUBLE        (i)   kappa  at integration point
\param    kappanint  DOUBLE        (i)   kappan at integration point
\param   *factor     DOUBLE        (o)   factor
\param   *factor1    DOUBLE        (o)   factor
\param   *factor2    DOUBLE        (o)   factor
\param   *sig        DOUBLE        (o)   factor
\return void

------------------------------------------------------------------------*/
void f2_fac_eps(
                  DOUBLE   C_2,
                  DOUBLE   eps_proint,
                  DOUBLE   kappaint,
                  DOUBLE   kappanint,
	            DOUBLE  *factor,
	            DOUBLE  *factor1,
	            DOUBLE  *factor2,
	            DOUBLE  *sig
	           )
{

#ifdef DEBUG
dstrc_enter("f2_fac_eps");
#endif
/*---------------------------------------------------------------------*/

  *factor =2.0*C_2/kappaint;
  *factor2=    C_2/kappaint;
  *factor1=1.44*eps_proint/kappanint;
  *sig    =1.3;

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_fac_eps */

/*!---------------------------------------------------------------------
\brief routine to calculate production term (don't be confused from
factor 0.5 used in the comments in file f2_caltimerhs for the calculation
of the production rhs-part)

<pre>                                                        he  12/02

</pre>
\param   **vderxy     DOUBLE        (i)   vderi. at element nodes
\param   *production  DOUBLE        (o)   production term
\return void

------------------------------------------------------------------------*/
void f2_production(
	            DOUBLE  **vderxy,
	            DOUBLE  *production
	           )
{

#ifdef DEBUG
dstrc_enter("f2_production");
#endif
/*---------------------------------------------------------------------
             production = grad(u):(grad(u) + (grad(u))^T)
-----------------------------------------------------------------------*/

*production  = 2*(pow(vderxy[0][0],2)+pow(vderxy[1][1],2)+vderxy[1][0]*vderxy[0][1])+
                  pow(vderxy[0][1],2)+pow(vderxy[1][0],2);

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_production */

/*!---------------------------------------------------------------------
\brief routine to calculate eddyviscosity for RANS at integration point

<pre>                                                        he  01/03

</pre>
\param   *eleke        ELEMENT      (i)   actual element for kapeps
\param   *eddyint      DOUBLE       (o)   eddy-visc. at integration point
\param   *funct        DOUBLE       (i)   shape functions
\param   *eddy         DOUBLE       (i)   eddy visc. on nodes
\param   *ipos                      (i)   node array positions
\param    iel	    INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_eddyirans(
             ELEMENT    *eleke,
             DOUBLE     *eddyint,
             DOUBLE     *funct,
             DOUBLE     *eddy,
             ARRAY_POSITION *ipos,
             INT         iel
	     )
{
INT     i,j;
NODE    *actnode;      /* the actual node                               */

#ifdef DEBUG
dstrc_enter("f2_eddyirans");
#endif

   for(i=0;i<iel;i++) /* loop nodes of element */
   {
      actnode=eleke->node[i];
/*----------------------------------- get kappa from nodes             */
      eddy[i] =actnode->sol_increment.a.da[ipos->velnp][1];
   }

   *eddyint =ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      *eddyint    += funct[j]*eddy[j];
   } /* end loop over j */


/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_eddyirans */

/*!---------------------------------------------------------------------
\brief routine to calculate kapeps derivatives at integration point

<pre>                                                         he   12/02

In this routine the derivatives of the kapeps w.r.t x/y are calculated

</pre>
\param   *kapepsderxy   DOUBLE        (o)   kapeps derivativs
\param  **derxy         DOUBLE        (i)   global derivatives
\param   *ekapeps       DOUBLE        (i)   kapeps at element nodes
\param    iel	      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_kapepsder(
             DOUBLE  *kapepsderxy,
             DOUBLE **derxy,
	       DOUBLE  *kapeps,
             INT      iel
	       )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f2_kapepsder");
#endif

for (i=0;i<2;i++) /* loop directions i */
{
   kapepsderxy[i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      kapepsderxy[i] += derxy[i][j]*kapeps[j];
   } /* end of loop over j */
} /* end of loop over i */


/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_kapepsder */

/*!---------------------------------------------------------------------
\brief routine to calculate 2nd kapeps derivatives at integration point

<pre>                                                         he  12/02

In this routine the 2nd derivatives of the kapeps
w.r.t x/y are calculated

</pre>
\param   *kapepsderxy2  DOUBLE        (o)   2nd kapeps derivativs
\param  **derxy2        DOUBLE        (i)   2nd global derivatives
\param   *kapepsn       DOUBLE        (i)   velocites at element nodes
\param    iel	      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_kapepsder2(
                   DOUBLE  *kapepsderxy2,
                   DOUBLE **derxy2,
	             DOUBLE  *kapepsn,
	             INT      iel
	           )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f2_vder2");
#endif

for (i=0;i<3;i++)
{
   kapepsderxy2[i]=ZERO;
   for (j=0;j<iel;j++)
   {
      kapepsderxy2[i] += derxy2[i][j]*kapepsn[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_kapepsder2 */

/*!---------------------------------------------------------------------
\brief routine to calculate 2nd kapeps derivatives at integration point

<pre>                                                         he  03/03

In this routine velint_dc for DISC. CAPT. is calc.

</pre>
\param   *velint        DOUBLE        (i)   2nd kapeps derivativs
\param   *velint_dc     DOUBLE        (o)   2nd global derivatives
\param   *kapepsderxy   DOUBLE        (i)   kapepsderiv.
\return void

------------------------------------------------------------------------*/
void f2_vel_dc(
                   DOUBLE  *velint,
                   DOUBLE  *velint_dc,
	             DOUBLE  *kapepsderxy
	         )
{
DOUBLE skalar,square;
INT    idum;

#ifdef DEBUG
dstrc_enter("f2_vel_dc");
#endif
/*---------------------------------------------------------------------*/

if(FABS(kapepsderxy[0])<0.001) kapepsderxy[0] = 0.0;
if(FABS(kapepsderxy[1])<0.001) kapepsderxy[1] = 0.0;

skalar = velint[0]*kapepsderxy[0] + velint[1]*kapepsderxy[1];
square = pow(kapepsderxy[0],2) + pow(kapepsderxy[1],2);

if(square != 0.0 && fdyn->dis_capt == 1)
{
 velint_dc[0] = skalar*kapepsderxy[0]/square;
 velint_dc[1] = skalar*kapepsderxy[1]/square;
}
else
{
 velint_dc[0] = 0.0;
 velint_dc[1] = 0.0;
}


/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_vel_dc */


/*!---------------------------------------------------------------------
\brief permutation of element stiffness matrix

<pre>                                                         he  12/02

routine to add galerkin and stabilisation parts of the elment
stiffness matrix !
this is necessary since we would like to use the existing assembly
routines for the stiffness matrix

</pre>
\param  **estif   DOUBLE	 (i/o) ele stiffnes matrix
\param  **emass   DOUBLE	 (i)   ele mass matrix
\param  **tmp     DOUBLE	 (-)   working array
\param    iel	  INT		 (i)   number of nodes in ele
\return void
------------------------------------------------------------------------*/
void f2_estifadd_tu(
		   DOUBLE         **estif,
		   DOUBLE         **emass,
		   DOUBLE         **tmp,
		   INT              iel
	          )
{
INT    i,j,icol,irow;     /* simply some counters                       */
INT    nkapepsdof;        /* number of kapeps dofs                      */
DOUBLE thsl;              /* factor for LHS (THETA*DT)                  */

#ifdef DEBUG
dstrc_enter("f2_estifadd_tu");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*----------------------------------------------------- set some values */
nkapepsdof  = iel;
thsl   = fdyn->thsl;

/*--------------------------------------------- copy estif to tmp-array *
                           and mutlitply stiffniss matrix with THETA*DT */
for (i=0;i<nkapepsdof;i++)
{
   for (j=0;j<nkapepsdof;j++)
   {
      tmp[i][j] = estif[i][j] * thsl;
   } /* end of loop over j */
} /* end of loop over i */
/*------------------------------- add mass matrix for instationary case */
if (fdyn->nis==0)
{
   for (i=0;i<nkapepsdof;i++)
   {
      for (j=0;j<nkapepsdof;j++)
      {
         tmp[i][j] += emass[i][j];
      } /* end of loop over j */
   } /* end of loop over i */
} /* endif (fdyn->nis==0) */

/*---------------------------------------- estif=emass+(theta*dt)*estif */
for (i=0;i<nkapepsdof;i++)
{
   for (j=0;j<nkapepsdof;j++)
   {
    estif[i][j] = tmp[i][j];
    } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_permestif */


#endif
#endif
