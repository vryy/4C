/*!----------------------------------------------------------------------
\file
\brief service routines for fluid2 element 

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
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
/*!--------------------------------------------------------------------- 
\brief set all arrays for element calculation

<pre>                                                         he  02/03

				      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *data     FLUID_DATA	     (i)
\param   *ele      ELEMENT	     (i)    actual element
\param   *elev     ELEMENT	     (i)    actual element for velocity
\param   *kapomen    double	     (o)    kapome at time n
\param   *kapomeg    double	     (o)    kapome at time n+g
\param   *kapomepro  double	     (o)    kapome at time n+g
\param   *eddyg      double	     (o)    eddy-visc. at time n+g
\param   *eddypro    double	     (o)    eddy-visc. for prod. term
\param   *kappan     double	     (o)    kappa at time n
\param   *omega      double	     (o)    omega
\param  **evel       double	     (o)    ele velocity
\param  **xzye       double        (o)   nodal coordinates
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calset_tu_1( 
                    FLUID_DYN_CALC  *dynvar, 
                    FLUID_DATA      *data,     
	              ELEMENT         *ele,     
                    ELEMENT         *elev,
                    double          *kapomen,    
	              double          *kapomeg,
	              double          *kapomepro,
                    double          *eddyg,
                    double          *eddypro,
	              double          *kappan,    
	              double          *omega,    
	              double         **evel,
	              double         **xyze
                   )
{
int i,j;            /* simply a counter                                 */
int kap_ome;
NODE  *actnode;     /* actual node                                      */
#ifdef DEBUG 
dstrc_enter("f2_calset_tu_1");
#endif

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
 if (dynvar->kapomega_flag==0) kap_ome=0;
 if (dynvar->kapomega_flag==1) kap_ome=2;
   
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
    if (dynvar->kapomega_flag==0)
    {
     omega[i] = actnode->sol_increment.a.da[3][2];
    
     if (dynvar->kappan==2)
     {
      actnode->sol_increment.a.da[2][0] = actnode->sol_increment.a.da[3][0];
     }
    }
/*---------- for omega equation: kappan is needed for production term  */
    if (dynvar->kapomega_flag==1 && dynvar->kappan==2) 
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
\param   *kapomeint   double        (o)   kapome at integration point
\param   *funct       double        (i)   shape functions
\param   *ekapome    double        (i)   kapome at element nodes
\param    iel	    int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_kapomei(
             double  *kapomeint,     
             double  *funct,    
	       double  *kapome,     
             int      iel       
	     ) 
{
int     j;

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
\param   *kapomederxy   double        (o)   kapome derivativs
\param  **derxy         double        (i)   global derivatives
\param   *ekapome       double        (i)   kapome at element nodes
\param    iel	      int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_kapomeder(
             double  *kapomederxy,     
             double **derxy,    
	       double  *kapome,    
             int      iel       
	       ) 
{
int     i,j;

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
\param   *kapomederxy2  double        (o)   2nd kapome derivativs
\param  **derxy2        double        (i)   2nd global derivatives
\param   *kapomen       double        (i)   kapome at element nodes
\param    iel	      int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_kapomeder2(
                   double  *kapomederxy2,    
                   double **derxy2,    
	             double  *kapomen,      
	             int      iel        
	           ) 
{
int     i,j;

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
\param    *kapomederxy  double      (i)   kapome deriv. at integr. point
\param    *omegaderxy   double      (i)   omega deriv. at integr. point
\param     omegaint      double      (i)   omega at integr. point 
\param    xi            double      (o)   factor
\return void                                                                       

------------------------------------------------------------------------*/
void f2_xi_kappa(
                  double  *kapomederxy,     
                  double  *omegaderxy,    
	            double   omegaint,     
	            double  *xi     
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
\param     xi           double      (i)   xi factor
\param    eddyint       double      (i)   eddyvisc. at integration point
\param    kapomeint     double      (i)   kapome at integration point
\param    omegaint      double      (i)   omega at integration point
\param    visc          double      (i)   viscosity
\param    factor        double      (o)   factor
\param    factor1       double      (o)   factor
\param    factor2       double      (o)   factor
\param    sig           double      (o)   factor
\return void                                                                       

------------------------------------------------------------------------*/
void f2_fac_kappa_1(
                  double   xi,     
                  double   eddyint,    
                  double   kapomeint,    
                  double   omegaint,    
                  double   visc,    
                  double  *factor,    
                  double  *factor1,    
                  double  *factor2,    
	            double  *sig    
	           ) 
{
double f;

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
\param   **vderxy       double      (i)   vel. deriv. at integr. point
\param    kapomeint     double      (i)   omega at integr. point 
\param    xi            double      (o)   factor
\return void                                                                       

------------------------------------------------------------------------*/
void f2_xi_ome(
                  double  **vderxy,     
	            double   kapomeint,     
	            double  *xi     
	           ) 
{
double wert,zaehler;

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
\param     xi           double      (i)   xi factor
\param    ome_proint    double      (i)   omega at integration point
\param    kappanint     double      (i)   kappan at integration point
\param    visc          double      (i)   viscosity
\param    factor        double      (o)   factor
\param    factor1       double      (o)   factor
\param    factor2       double      (o)   factor
\param    sig           double      (o)   factor
\return void                                                                       

------------------------------------------------------------------------*/
void f2_fac_ome(
                  double   xi,     
                  double   ome_proint,    
                  double   kappanint,    
                  double   visc,    
                  double  *factor,    
                  double  *factor1,    
                  double  *factor2,    
	            double  *sig    
	           ) 
{
double f;

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
\param   *dynvar    FLUID_DYN_CALC    (i)
\param   *velint        double        (i)   2nd kapeps derivativs
\param   *velint_dc     double        (o)   2nd global derivatives
\param   *kapomederxy   double        (i)   kapomederiv.
\return void                                                                       

------------------------------------------------------------------------*/
void f2_vel_dc_1(
		       FLUID_DYN_CALC  *dynvar,
                   double  *velint,    
                   double  *velint_dc,    
	             double  *kapomederxy      
	         ) 
{
double skalar,square; 
int    idum;

#ifdef DEBUG 
dstrc_enter("f2_vel_dc_1");
#endif
/*---------------------------------------------------------------------*/

if(FABS(kapomederxy[0])<0.001) kapomederxy[0] = 0.0;
if(FABS(kapomederxy[1])<0.001) kapomederxy[1] = 0.0;

skalar = velint[0]*kapomederxy[0] + velint[1]*kapomederxy[1];
square = pow(kapomederxy[0],2) + pow(kapomederxy[1],2);

if(square != 0.0 && dynvar->dis_capt == 1)
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
\param   *kappanint    double       (o)   kappan at integration point
\param   *ome_proint   double       (o)   ome_pro at integration point
\param   *funct        double       (i)   shape functions
\param   *kappan       double       (i)   kappan (start-value) kappa equation
\param   *kapomepro    double       (i)   omega for production term
\param    iel	    int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_kappain(	          
             double     *kappanint,     
             double     *ome_proint,     
             double     *funct,    
             double     *kappan,    
             double     *kapomepro,    
             int         iel       
	     ) 
{
int     j;

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
