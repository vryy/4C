/*!----------------------------------------------------------------------
\file
\brief service routines for fluid2 element

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2TU
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
#include "fluid2.h"
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
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

static FLUID_DYNAMIC *fdyn;
/*!---------------------------------------------------------------------
\brief set all arrays for element calculation

<pre>                                                         he  02/03


</pre>
\param   *ele      ELEMENT	     (i)    actual element
\param   *elev     ELEMENT	     (i)    actual element for velocity
\param   *kapomen    DOUBLE	     (o)    kapome at time n
\param   *kapomeg    DOUBLE	     (o)    kapome at time n+g
\param   *kapomepro  DOUBLE	     (o)    kapome at time n+g
\param   *eddyg      DOUBLE	     (o)    eddy-visc. at time n+g
\param   *eddypro    DOUBLE	     (o)    eddy-visc. for prod. term
\param   *kappan     DOUBLE	     (o)    kappa at time n
\param   *omega      DOUBLE	     (o)    omega
\param  **evel       DOUBLE	     (o)    ele velocity
\param  **xzye       DOUBLE        (o)   nodal coordinates
\return void

------------------------------------------------------------------------*/
void f2_calset_tu_1(
	              ELEMENT         *ele,
                    ELEMENT         *elev,
                    DOUBLE          *kapomen,
	              DOUBLE          *kapomeg,
	              DOUBLE          *kapomepro,
                    DOUBLE          *eddyg,
                    DOUBLE          *eddypro,
	              DOUBLE          *kappan,
	              DOUBLE          *omega,
	              DOUBLE         **evel,
	              DOUBLE         **xyze
                   )
{
INT i,j;            /* simply a counter                                 */
INT kap_ome;
NODE  *actnode;     /* actual node                                      */
FLUID_DATA      *data;

#ifdef DEBUG
dstrc_enter("f2_calset_tu_1");
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
 |	                  i=2 : omega                                    |
 |	                  i=3 : charact. lenght                          |
 *---------------------------------------------------------------------*/
 if (fdyn->kapomega_flag==0) kap_ome=0;
 if (fdyn->kapomega_flag==1) kap_ome=2;

   for(i=0;i<ele->numnp;i++) /* loop nodes of element */
   {
     actnode=ele->node[i];

/*----------------------------------- set element kapome (n+1,i)       */
      kapomeg[i]  =actnode->sol_increment.a.da[3][kap_ome];
/*----------------------------------- set element kapome (n)           */
      kapomen[i]  =actnode->sol_increment.a.da[1][kap_ome];
/*----------------------------------- set eddy viscosity for timestep  */
      eddyg[i]    =actnode->sol_increment.a.da[3][1];
      eddypro[i]  =actnode->sol_increment.a.da[2][1];

/*------------------- for kappa equation: omega is needed for factors  */
    if (fdyn->kapomega_flag==0)
    {
     omega[i] = actnode->sol_increment.a.da[3][2];

     if (fdyn->kappan==2)
     {
      actnode->sol_increment.a.da[2][0] = actnode->sol_increment.a.da[3][0];
     }
    }
/*---------- for omega equation: kappan is needed for production term  */
    if (fdyn->kapomega_flag==1 && fdyn->kappan==2)
    {
      actnode->sol_increment.a.da[2][2] = actnode->sol_increment.a.da[3][2];
      kappan[i] = actnode->sol_increment.a.da[2][0];
    }
/*----------------------------------- set element kapome (pro)         */
      kapomepro[i]=actnode->sol_increment.a.da[2][kap_ome];

/*----------------------- get velocities calculated form Navier-Stokes */
    for(j=0;j<2;j++)
    {
      actnode=elev->node[i];
      evel[j][i]=actnode->sol_increment.a.da[3][j];
    }

 } /* end of loop over nodes of element */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calset */

/*!---------------------------------------------------------------------
\brief routine to calculate kapome at integration point

<pre>                                                        he  02/03

</pre>
\param   *kapomeint   DOUBLE        (o)   kapome at integration point
\param   *funct       DOUBLE        (i)   shape functions
\param   *ekapome    DOUBLE        (i)   kapome at element nodes
\param    iel	    INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_kapomei(
             DOUBLE  *kapomeint,
             DOUBLE  *funct,
	       DOUBLE  *kapome,
             INT      iel
	     )
{
INT     j;

#ifdef DEBUG
dstrc_enter("f2_kapomei");
#endif
/*---------------------------------------------------------------------*/

   *kapomeint=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      *kapomeint += funct[j]*kapome[j];
   } /* end loop over j */


/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_kapepsi */

/*!---------------------------------------------------------------------
\brief routine to calculate kapome derivatives at integration point

<pre>                                                         he   02/03

In this routine the derivatives of the kapeps w.r.t x/y are calculated

</pre>
\param   *kapomederxy   DOUBLE        (o)   kapome derivativs
\param  **derxy         DOUBLE        (i)   global derivatives
\param   *ekapome       DOUBLE        (i)   kapome at element nodes
\param    iel	      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_kapomeder(
             DOUBLE  *kapomederxy,
             DOUBLE **derxy,
	       DOUBLE  *kapome,
             INT      iel
	       )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f2_kapomeder");
#endif

for (i=0;i<2;i++) /* loop directions i */
{
   kapomederxy[i]=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      kapomederxy[i] += derxy[i][j]*kapome[j];
   } /* end of loop over j */
} /* end of loop over i */


/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_kapepsder */

/*!---------------------------------------------------------------------
\brief routine to calculate 2nd kapome derivatives at integration point

<pre>                                                         he  02/03

In this routine the 2nd derivatives of the kapome
w.r.t x/y are calculated

</pre>
\param   *kapomederxy2  DOUBLE        (o)   2nd kapome derivativs
\param  **derxy2        DOUBLE        (i)   2nd global derivatives
\param   *kapomen       DOUBLE        (i)   kapome at element nodes
\param    iel	      INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_kapomeder2(
                   DOUBLE  *kapomederxy2,
                   DOUBLE **derxy2,
	             DOUBLE  *kapomen,
	             INT      iel
	           )
{
INT     i,j;

#ifdef DEBUG
dstrc_enter("f2_kapomeder2");
#endif

for (i=0;i<3;i++)
{
   kapomederxy2[i]=ZERO;
   for (j=0;j<iel;j++)
   {
      kapomederxy2[i] += derxy2[i][j]*kapomen[j];
   } /* end of loop over j */
} /* end of loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_kapomeder2 */

/*!---------------------------------------------------------------------
\brief routine to calculate factors for kappa-equation

<pre>                                                        he  02/03

</pre>
\param    *kapomederxy  DOUBLE      (i)   kapome deriv. at integr. point
\param    *omegaderxy   DOUBLE      (i)   omega deriv. at integr. point
\param     omegaint      DOUBLE      (i)   omega at integr. point
\param    xi            DOUBLE      (o)   factor
\return void

------------------------------------------------------------------------*/
void f2_xi_kappa(
                  DOUBLE  *kapomederxy,
                  DOUBLE  *omegaderxy,
	            DOUBLE   omegaint,
	            DOUBLE  *xi
	           )
{

#ifdef DEBUG
dstrc_enter("f2_xi_kappa");
#endif
/*---------------------------------------------------------------------*/

 *xi = (kapomederxy[0]*omegaderxy[0]+
        kapomederxy[1]*omegaderxy[1])/pow(omegaint,3);

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_xi_kappa */

/*!---------------------------------------------------------------------
\brief routine to calculate factors for kappa-equation

<pre>                                                        he  02/03

</pre>
\param     xi           DOUBLE      (i)   xi factor
\param    eddyint       DOUBLE      (i)   eddyvisc. at integration point
\param    kapomeint     DOUBLE      (i)   kapome at integration point
\param    omegaint      DOUBLE      (i)   omega at integration point
\param    visc          DOUBLE      (i)   viscosity
\param    factor        DOUBLE      (o)   factor
\param    factor1       DOUBLE      (o)   factor
\param    factor2       DOUBLE      (o)   factor
\param    sig           DOUBLE      (o)   factor
\return void

------------------------------------------------------------------------*/
void f2_fac_kappa_1(
                  DOUBLE   xi,
                  DOUBLE   eddyint,
                  DOUBLE   kapomeint,
                  DOUBLE   omegaint,
                  DOUBLE   visc,
                  DOUBLE  *factor,
                  DOUBLE  *factor1,
                  DOUBLE  *factor2,
	            DOUBLE  *sig
	           )
{
DOUBLE f;

#ifdef DEBUG
dstrc_enter("f2_fac_kappa_1");
#endif


/*---------------------------------------------------------------------*/
if(xi<=0.0)  f=1.0;
if(xi> 0.0)  f=(1+680*xi*xi)/(1+400*xi*xi);
/*---------------------------------------------------------------------*/

*factor  = 2.0*0.09*f/eddyint;
*factor2 = 0.09*f/eddyint;
*factor1 = 1.0;
*sig     = 0.5;

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_fac_kappa_1 */

/*!---------------------------------------------------------------------
\brief routine to calculate factors for omega-equation

<pre>                                                        he  02/03

</pre>
\param   **vderxy       DOUBLE      (i)   vel. deriv. at integr. point
\param    kapomeint     DOUBLE      (i)   omega at integr. point
\param    xi            DOUBLE      (o)   factor
\return void

------------------------------------------------------------------------*/
void f2_xi_ome(
                  DOUBLE  **vderxy,
	            DOUBLE   kapomeint,
	            DOUBLE  *xi
	           )
{
DOUBLE wert,zaehler;

#ifdef DEBUG
dstrc_enter("f2_xi_ome");
#endif

/*-----------------------------------------------------------------------
                      zaehler =
1/2 (grad u - grad^T u) 1/2 (grad u - grad^T u) : 1/2 (grad u + grad^T u)
-------------------------------------------------------------------------*/

 wert    = vderxy[0][1]-vderxy[1][0];
 zaehler = -2.0/8.0*(wert*wert*(vderxy[0][0]+vderxy[1][1]));

 *xi = FABS(zaehler/pow(0.09*kapomeint,3));

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_xi_ome */


/*!---------------------------------------------------------------------
\brief routine to calculate factors for omega-equation

<pre>                                                        he  02/03

</pre>
\param     xi           DOUBLE      (i)   xi factor
\param    ome_proint    DOUBLE      (i)   omega at integration point
\param    kappanint     DOUBLE      (i)   kappan at integration point
\param    visc          DOUBLE      (i)   viscosity
\param    factor        DOUBLE      (o)   factor
\param    factor1       DOUBLE      (o)   factor
\param    factor2       DOUBLE      (o)   factor
\param    sig           DOUBLE      (o)   factor
\return void

------------------------------------------------------------------------*/
void f2_fac_ome(
                  DOUBLE   xi,
                  DOUBLE   ome_proint,
                  DOUBLE   kappanint,
                  DOUBLE   visc,
                  DOUBLE  *factor,
                  DOUBLE  *factor1,
                  DOUBLE  *factor2,
	            DOUBLE  *sig
	           )
{
DOUBLE f;

#ifdef DEBUG
dstrc_enter("f2_fac_ome");
#endif

/*---------------------------------------------------------------------*/
f=(1+70*xi)/(1+80*xi);

*factor  = 2.0*0.072*f;
*factor2 =     0.072*f;
*factor1 = 0.52*ome_proint/kappanint;
*sig     = 0.5;

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_fac_ome */

/*!---------------------------------------------------------------------
\brief routine to calculate 2nd kapeps derivatives at integration point

<pre>                                                         he  03/03

In this routine velint_dc for DISC. CAPT. is calc.

</pre>
\param   *velint        DOUBLE        (i)   2nd kapeps derivativs
\param   *velint_dc     DOUBLE        (o)   2nd global derivatives
\param   *kapomederxy   DOUBLE        (i)   kapomederiv.
\return void

------------------------------------------------------------------------*/
void f2_vel_dc_1(
                   DOUBLE  *velint,
                   DOUBLE  *velint_dc,
	             DOUBLE  *kapomederxy
	         )
{
DOUBLE skalar,square;
INT    idum;

#ifdef DEBUG
dstrc_enter("f2_vel_dc_1");
#endif
/*---------------------------------------------------------------------*/

if(FABS(kapomederxy[0])<0.001) kapomederxy[0] = 0.0;
if(FABS(kapomederxy[1])<0.001) kapomederxy[1] = 0.0;

skalar = velint[0]*kapomederxy[0] + velint[1]*kapomederxy[1];
square = pow(kapomederxy[0],2) + pow(kapomederxy[1],2);

if(square != 0.0 && fdyn->dis_capt == 1)
{
 velint_dc[0] = skalar*kapomederxy[0]/square;
 velint_dc[1] = skalar*kapomederxy[1]/square;
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
} /* end of f2_vel_dc_1 */

/*!---------------------------------------------------------------------
\brief routine to calculate kappa for omega at integration point

<pre>                                                        he  12/02

</pre>
\param   *kappanint    DOUBLE       (o)   kappan at integration point
\param   *ome_proint   DOUBLE       (o)   ome_pro at integration point
\param   *funct        DOUBLE       (i)   shape functions
\param   *kappan       DOUBLE       (i)   kappan (start-value) kappa equation
\param   *kapomepro    DOUBLE       (i)   omega for production term
\param    iel	    INT           (i)   number of nodes in this element
\return void

------------------------------------------------------------------------*/
void f2_kappain(
             DOUBLE     *kappanint,
             DOUBLE     *ome_proint,
             DOUBLE     *funct,
             DOUBLE     *kappan,
             DOUBLE     *kapomepro,
             INT         iel
	     )
{
INT     j;

#ifdef DEBUG
dstrc_enter("f2_kappain");
#endif
/*---------------------------------------------------------------------*/

   *kappanint =ZERO;
   *ome_proint=ZERO;
   for (j=0;j<iel;j++) /* loop over all nodes j of the element */
   {
      *kappanint    += funct[j]*kappan[j];
      *ome_proint   += funct[j]*kapomepro[j];
    } /* end loop over j */


/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_kappain */


#endif
