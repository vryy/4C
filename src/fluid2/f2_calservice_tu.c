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

<pre>                                                         he  12/02

				      
</pre>
\param   *dynvar   FLUID_DYN_CALC  (i)
\param   *data     FLUID_DATA	     (i)
\param   *ele      ELEMENT	     (i)    actual element
\param   *elev     ELEMENT	     (i)    actual element for velocity
\param   *kapepsn    double	     (o)    kapeps at time n
\param   *kapepsg    double	     (o)    kapeps at time n+g
\param   *kapepspro  double	     (o)    kapeps at time n+g
\param   *eddyg      double	     (o)    eddy-visc. at time n+g
\param   *eddypro    double	     (o)    eddy-visc. for prod. term
\param   *kappa      double	     (o)    kappa at time n+g
\param   *kappan     double	     (o)    kappa at time n
\param   *epsilon    double	     (o)    epsilon at time n+g
\param  **evel       double	     (o)    element velocities
\param  **xzye       double        (o)   nodal coordinates
\return void                                                                       

------------------------------------------------------------------------*/
void f2_calset_tu( 
                  FLUID_DYN_CALC  *dynvar, 
                  FLUID_DATA      *data,     
	            ELEMENT         *ele,     
                  ELEMENT         *elev,
                  double          *kapepsn,    
	            double          *kapepsg,
	            double          *kapepspro,
                  double          *eddyg,
                  double          *eddypro,
	            double          *kappa,    
	            double          *kappan,    
	            double          *epsilon,    
	            double         **evel,
	            double         **xyze
                 )
{
int i,j;            /* simply a counter                                 */
int kap_eps;
NODE  *actnode;     /* actual node                                      */
#ifdef DEBUG 
dstrc_enter("f2_calset_tu");
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
 |	                  i=2 : epsilon                                  |
 |	                  i=3 : charact. lenght                          |
 *---------------------------------------------------------------------*/
 if (dynvar->kapeps_flag==0) kap_eps=0;
 if (dynvar->kapeps_flag==1) kap_eps=2;
   
   for(i=0;i<ele->numnp;i++) /* loop nodes of element */
   {
     actnode=ele->node[i];

/*----------------------------------- set element kapeps (n+1,i)       */      
      kapepsg[i]  =actnode->sol_increment.a.da[3][kap_eps];
/*----------------------------------- set element kapeps (n)           */      
      kapepsn[i]  =actnode->sol_increment.a.da[1][kap_eps];
/*----------------------------------- set eddy viscosity for timestep  */
      eddyg[i]   =actnode->sol_increment.a.da[3][1];
      eddypro[i] =actnode->sol_increment.a.da[2][1];

/*--------------------- for kappa equation: epsilon is needed for R_t  */
    if (dynvar->kapeps_flag==0)
    {
     epsilon[i] = actnode->sol_increment.a.da[3][2];

     if (dynvar->kappan==2)
     {
      actnode->sol_increment.a.da[2][0] = actnode->sol_increment.a.da[3][0];
     }
    }

/*-------- for epsilon equation: kappan is needed for production term  */
    if (dynvar->kapeps_flag==1) 
    { 
     kappa[i]  =actnode->sol_increment.a.da[3][0];

     if (dynvar->kappan==2)
     {
      actnode->sol_increment.a.da[2][2] = actnode->sol_increment.a.da[3][2];
      kappan[i] =actnode->sol_increment.a.da[2][0];
     }
    }

/*----------------------------------- set element kapeps (pro)         */      
      kapepspro[i]=actnode->sol_increment.a.da[2][kap_eps];

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
\brief routine to calculate kapeps at integration point

<pre>                                                        he  12/02
				      
</pre>
\param   *kapepsint   double        (o)   kapeps at integration point
\param   *funct       double        (i)   shape functions
\param   *ekapeps    double        (i)   kapeps at element nodes
\param    iel	    int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_kapepsi(
             double  *kapepsint,     
             double  *funct,    
	       double  *kapeps,     
             int      iel       
	     ) 
{
int     j;

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
\param   *eddyint     double        (o)   eddy at integration point
\param   *funct       double        (i)   shape functions
\param   **eddy       double        (i)   kapeps at element nodes
\param    iel	    int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_eddyi(
             double  *eddyint,     
             double  *funct,    
	       double  *eddy,     
             int      iel       
	     ) 
{
int     j;

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
\param   *kappaint     double       (o)   kappa at integration point
\param   *kappanint    double       (o)   kappan at integration point
\param   *eps_proint   double       (o)   eps_pro at integration point
\param   *funct        double       (i)   shape functions
\param   *kappa        double       (i)   kappa from kappa equation
\param   *kappan       double       (i)   kappan (start-value) kappa equation
\param   *eps_pro      double       (i)   epsilon for production term
\param    iel	    int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_kappai_tu(	          
             double     *kappaint,     
             double     *kappanint,     
             double     *eps_proint,     
             double     *funct,    
             double     *kappa,    
             double     *kappan,    
             double     *eps_pro,    
             int         iel       
	     ) 
{
int     j;

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
\param    kapepsint     double       (i)   kappa at integration point
\param  **epsilon       double       (i)   epsilon at nodes
\param   *funct         double       (i)   shape functions
\param    visc	      double       (i)   molecular viscosity
\param   *C_u           double       (o)   factor
\param    iel	    int            (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_C_kappa(	          
             double      kapepsint,     
             double     *epsilon,
             double     *funct, 
             double      visc,    
             double     *C_u,
             int         iel       
	     ) 
{
int     i;
double  R_t;           /* turbulent Reynold's number                    */
double  epsilonint;       

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
\param    kapepsint     double       (i)   epsilon at integration point
\param    kappaint      double       (i)   kappa at integration point
\param    visc	      double       (i)   molecular viscosity
\param   *C_2           double       (o)   factor
\param    iel	    int            (i)   number of nodes in this element
\return void                                                                       
				      
------------------------------------------------------------------------*/
void f2_C_eps(	          
             double      kapepsint,     
             double      kappaint,
             double      visc,    
             double     *C_2,
             int         iel       
	     ) 
{
double  R_t;           /* turbulent Reynold's number                    */

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
\param   **vderxy2     double       (i)   2nd deriv. for velocity at int.
\param    vderxy_12    double       (o)   (vderxy)^2 at integration point
\return void                                                                       

------------------------------------------------------------------------*/
void f2_v(	          
             double    **vderxy2,
             double     *vderxy_12 
	   ) 
{
double  factor;
int     i;

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
\param    C_u         double        (i)   factor
\param    eddyint     double        (i)   eddy at integration point
\param   *factor      double        (o)   factor
\param   *factor1     double        (o)   factor
\param   *factor2     double        (o)   factor
\param   *sig         double        (o)   factor
\return void                                                                       

------------------------------------------------------------------------*/
void f2_fac_kappa(
                  double   C_u,     
                  double   eddyint,    
	            double  *factor,     
	            double  *factor1,     
	            double  *factor2,     
	            double  *sig     
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
\param    C_2        double        (i)   factor
\param    eps_proint double        (i)   epsilon at integration point
\param    kappaint   double        (i)   kappa  at integration point
\param    kappanint  double        (i)   kappan at integration point
\param   *factor     double        (o)   factor
\param   *factor1    double        (o)   factor
\param   *factor2    double        (o)   factor
\param   *sig        double        (o)   factor
\return void                                                                       

------------------------------------------------------------------------*/
void f2_fac_eps(
                  double   C_2,     
                  double   eps_proint,    
                  double   kappaint,    
                  double   kappanint,    
	            double  *factor,     
	            double  *factor1,     
	            double  *factor2,     
	            double  *sig     
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
\param   **vderxy     double        (i)   vderi. at element nodes
\param   *production  double        (o)   production term
\return void                                                                       

------------------------------------------------------------------------*/
void f2_production(
	            double  **vderxy,     
	            double  *production    
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
\param   *eddyint      double       (o)   eddy-visc. at integration point
\param   *funct        double       (i)   shape functions
\param   *eddy         double       (i)   eddy visc. on nodes
\param    iel	    int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_eddyirans(	          
             ELEMENT    *eleke,     
             double     *eddyint,     
             double     *funct,    
             double     *eddy,    
             int         iel       
	     ) 
{
int     i,j;
NODE    *actnode;      /* the actual node                               */

#ifdef DEBUG 
dstrc_enter("f2_eddyirans");
#endif

   for(i=0;i<iel;i++) /* loop nodes of element */
   {
      actnode=eleke->node[i];
/*----------------------------------- get kappa from nodes             */      
      eddy[i] =actnode->sol_increment.a.da[3][1];
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
\param   *kapepsderxy   double        (o)   kapeps derivativs
\param  **derxy         double        (i)   global derivatives
\param   *ekapeps       double        (i)   kapeps at element nodes
\param    iel	      int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_kapepsder(
             double  *kapepsderxy,     
             double **derxy,    
	       double  *kapeps,    
             int      iel       
	       ) 
{
int     i,j;

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
\param   *kapepsderxy2  double        (o)   2nd kapeps derivativs
\param  **derxy2        double        (i)   2nd global derivatives
\param   *kapepsn       double        (i)   velocites at element nodes
\param    iel	      int           (i)   number of nodes in this element
\return void                                                                       

------------------------------------------------------------------------*/
void f2_kapepsder2(
                   double  *kapepsderxy2,    
                   double **derxy2,    
	             double  *kapepsn,      
	             int      iel        
	           ) 
{
int     i,j;

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
\param   *dynvar    FLUID_DYN_CALC    (i)
\param   *velint        double        (i)   2nd kapeps derivativs
\param   *velint_dc     double        (o)   2nd global derivatives
\param   *kapepsderxy   double        (i)   kapepsderiv.
\return void                                                                       

------------------------------------------------------------------------*/
void f2_vel_dc(
		       FLUID_DYN_CALC  *dynvar,
                   double  *velint,    
                   double  *velint_dc,    
	             double  *kapepsderxy      
	         ) 
{
double skalar,square; 
int    idum;

#ifdef DEBUG 
dstrc_enter("f2_vel_dc");
#endif
/*---------------------------------------------------------------------*/

if(FABS(kapepsderxy[0])<0.001) kapepsderxy[0] = 0.0;
if(FABS(kapepsderxy[1])<0.001) kapepsderxy[1] = 0.0;

skalar = velint[0]*kapepsderxy[0] + velint[1]*kapepsderxy[1];
square = pow(kapepsderxy[0],2) + pow(kapepsderxy[1],2);

if(square != 0.0 && dynvar->dis_capt == 1)
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
\param  **estif   double	 (i/o) ele stiffnes matrix
\param  **emass   double	 (i)   ele mass matrix
\param  **tmp     double	 (-)   working array		
\param    iel	  int		 (i)   number of nodes in ele
\param   *dynvar  FLUID_DYN_CALC
\return void                                                                       
------------------------------------------------------------------------*/
void f2_estifadd_tu(                  
		   double         **estif,   
		   double         **emass, 
		   double         **tmp,   
		   int              iel,   
		   FLUID_DYN_CALC  *dynvar		   		    
	          ) 
{
int    i,j,icol,irow;     /* simply some counters                       */
int    nkapepsdof;        /* number of kapeps dofs                      */
double thsl;              /* factor for LHS (THETA*DT)                  */

#ifdef DEBUG 
dstrc_enter("f2_estifadd_tu");
#endif

/*----------------------------------------------------- set some values */
nkapepsdof  = iel;
thsl   = dynvar->thsl;

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
if (dynvar->nis==0)
{
   for (i=0;i<nkapepsdof;i++)
   {
      for (j=0;j<nkapepsdof;j++)
      {
         tmp[i][j] += emass[i][j];
      } /* end of loop over j */
   } /* end of loop over i */
} /* endif (dynvar->nis==0) */

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
