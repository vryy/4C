/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2_tu.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
/*----------------------------------------------------------------------*
 | integration loop for one fluid element                               |
 |                                                                      |
 |                                                          he    12/02 |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------
\brief integration loop for one fluid2_tu element

<pre>                                                       he   12/02

In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid2 element is calculated
      
</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *ele	 ELEMENT	   (i)    actual element
\param  *elev	 ELEMENT	   (i)    actual rans-element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     double	   (o)    element stiffness matrix
\param **emass     double	   (o)    element mass matrix
\param  *etforce   double	   (o)    element time force vector
\param  *eiforce   double	   (o)    element iter force vector
\param  *eproforce double	   (o)    element production force vector
\param  *funct     double	   (-)    natural shape functions
\param **deriv     double	   (-)	  deriv. of nat. shape funcs
\param **deriv2    double	   (-)    2nd deriv. of nat. shape f.
\param **xjm	 double	   (-)    jacobian matrix
\param **xyze      double          (-)    nodal coordinates
\param **derxy     double	   (-)	  global derivatives
\param **derxy2    double	   (-)    2nd global derivatives
\param  *kapomen   double	   (i)    kapome at time n
\param  *kapomeg   double	   (i)    kapome at time n+g
\param  *eddyg     double	   (i)    eddy-viscosity at time n+1
\param  *eddypro   double	   (i)    eddy-viscosity for prod. term
\param  *kappan    double	   (i)    kappa at time n
\param  *omega     double	   (-)    omega at nodes
\param  *kapomepro double	   (-)    kapome at nodes
\param  *omegaderxy   double	   (-)    omega derivates
\param  *kapomederxy  double	   (-)    kapome derivates
\param  *kapomederxy2 double	   (-)    second kapome derivatives
\param  *velint    double	   (-)    vel at integration point
\param  *velint_dc    double	   (-)    vel at integration point for DISC. CAPT.
\param **evel	 double	   (i)    velocity at nodes
\param **vderxy	 double	   (-)    velocity derivates
\param **vderxy2	 double	   (-)    velocity 2nd derivates
\param **wa1	 double	   (-)    working array
\param **wa2	 double	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f2_calint_tu_1(
               FLUID_DATA      *data,     
	         ELEMENT         *ele,     
	         ELEMENT         *elev, 
               FLUID_DYN_CALC  *dynvar, 
               double         **estif,   
	         double         **emass,   
	         double          *etforce, 
	         double          *eiforce, 
	         double          *eproforce, 
	         double          *funct,   
	         double         **deriv,   
	         double         **deriv2,  
	         double         **xjm,     
	         double         **xyze,
	         double         **derxy,   
	         double         **derxy2,  
	         double          *kapomen,   
	         double          *kapomeg,   
               double          *eddyg, 
               double          *eddypro, 
               double          *kappan, 
               double          *omega,    
               double          *kapomepro,    
               double          *omegaderxy,  
               double          *kapomederxy,  
               double          *kapomederxy2, 
	         double          *velint,  
	         double          *velint_dc,  
               double         **evel,
               double         **vderxy,
               double         **vderxy2,
               double         **wa1,     
	         double         **wa2      
	        )
{ 
int       i,j;        /* simply a counter                               */
int       iel;        /* number of nodes                                */
int       ntyp;       /* element type: 1 - quad; 2 - tri                */
int       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
int       nir,nis;    /* number of integration nodesin r,s direction    */
int       actmat;     /* material number of the element                 */
int       icode=2;    /* flag for eveluation of shape functions         */     
int       ihoel=0;    /* flag for higher order elements                 */
int       lr, ls;     /* counter for integration                        */
double    dens;       /* density                                        */
double    visc;       /* viscosity                                      */
double    xi;
double    fac;
double    factor,factor1;
double    factor2,sig,production;
double    facr, facs; /* integration weights                            */
double    det;        /* determinant of jacobian matrix                 */
double    e1,e2;      /* natural coordinates of integr. point           */
double    eddyint;
double    eddynint;
double    kappanint;
double    omegaint;
double    kapomeint;
double    kapomenint;
double    ome_proint;
DIS_TYP   typ;	    /* element type                                   */

#ifdef DEBUG 
dstrc_enter("f2_calint_tu_1");
#endif

/*----------------------------------------------------- initialisation */
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
ntyp = ele->e.f2_tu->ntyp; 
typ  = ele->distyp;

/*------- get integraton data and check if elements are "higher order" */
switch (ntyp)
{
case 1:  /* --> quad - element */
   icode   = 3;
   ihoel   = 1;
   /* initialise integration */
   nir = ele->e.f2_tu->nGP[0];
   nis = ele->e.f2_tu->nGP[1];
break;
case 2: /* --> tri - element */  
   if (iel>3)
   {
      icode   = 3;
      ihoel   = 1;
   }
   /* initialise integration */
   nir  = ele->e.f2_tu->nGP[0];
   nis  = 1;
   intc = ele->e.f2_tu->nGP[1];  
break;
default:
   dserror("ntyp unknown!");
} /* end switch(ntyp) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{    
   for (ls=0;ls<nis;ls++)
   {
/*--------------- get values of  shape functions and their derivatives */
      switch(ntyp)  
      {
      case 1:   /* --> quad - element */
	 e1   = data->qxg[lr][nir-1];
       facr = data->qwgt[lr][nir-1];
	 e2   = data->qxg[ls][nis-1];
	 facs = data->qwgt[ls][nis-1];
       f2_rec(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      case 2:   /* --> tri - element */              
	 e1   = data->txgr[lr][intc];
	 facr = data->twgt[lr][intc];
	 e2   = data->txgs[lr][intc];
	 facs = ONE;
	 f2_tri(funct,deriv,deriv2,e1,e2,typ,icode);
      break; 
      default:
         dserror("ntyp unknown!");
      } /* end switch(ntyp) */
/*-------------------------------------------- compute Jacobian matrix */
      f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;

/*------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);
      f2_gder2(ele,xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);

/*-------------------------------------- get vel. at integration point */               
      f2_veli(velint,funct,evel,iel);    

/*-------------------------------- get eddy-visc. at integration point */               
      f2_eddyi(&eddyint,funct,eddyg,iel);

/*---------------------------------- compute global derivates for vel. */
      f2_vder(vderxy,derxy,evel,iel);	

/*---------------------------- get kapome (n+1,i)  at integraton point */
      f2_kapomei(&kapomeint,funct,kapomeg,iel);

/*------------------ calculate stabilisation parameter for DISC. CAPT. */               
      f2_kapomeder(kapomederxy,derxy,kapomen,iel);
      f2_vel_dc_1(dynvar,velint,velint_dc,kapomederxy);

/*------------------------------- calculate factors for kappa equation */       
      if(dynvar->kapomega_flag==0)
      {
/*----------------- get kappa (n+1,i) derivatives at integration point */
       f2_kapomeder(kapomederxy,derxy,kapomeg,iel);
/*----------------- get omega (n+1,i) derivatives at integration point */
       f2_kapomeder(omegaderxy,derxy,omega,iel);
/*---------------------------- get kapome (n+1,i)  at integraton point */
       f2_kapomei(&omegaint,funct,omega,iel);
       f2_xi_kappa(kapomederxy,omegaderxy,omegaint,&xi);
       f2_fac_kappa_1(xi,eddyint,kapomeint,omegaint,visc,&factor,
                      &factor1,&factor2,&sig);
      }
/*------------------------------- calculate factors for omega equation */
      if(dynvar->kapomega_flag==1)
      {
       f2_kappain(&kappanint,&ome_proint,funct,kappan,kapomepro,iel); 
       f2_xi_ome(vderxy,kapomeint,&xi);
       f2_fac_ome(xi,ome_proint,kappanint,visc,&factor,&factor1,
                  &factor2,&sig);
      }

/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices                         |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/
 
 /*--------------------------------------------- compute matrix Kkapome */      
       f2_calkkapome(dynvar,estif,kapomeint,velint,eddyint,
                     funct,derxy,fac,visc,factor,sig,iel);
       
/*---------------------------------------------- compute matrix Mkapome */
	 if (dynvar->nis==0)	  	 	    
          f2_calmkapome(emass,funct,fac,iel);
  
/*----------------------------------------------------------------------*
 |         compute Stabilisation matrices                               |
 | NOTE:                                                                |
 |  Stabilisation matrices are all stored in one matrix "estif"         |
 |  Stabilisation mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
/*------------------------------------ stabilisation for matrix Kkapome */
       f2_calstabkkapome(ele,elev,dynvar,estif,kapomeint,velint,velint_dc,
                         eddyint,funct,derxy,derxy2,fac,visc,factor,sig,iel); 

/*------------------------------------ stabilisation for matrix Mkapome */
       if (dynvar->nis==0) 
          f2_calstabmkapome(ele,dynvar,emass,velint,velint_dc,funct,derxy,
                            fac,iel); 
/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 *----------------------------------------------------------------------*/ 

/*------------- calculate galerkin part of "Iter-RHS" (kapome dof) */
         f2_calgalifkapome(dynvar,eiforce,kapomeint,funct,fac,factor2,iel);

/*---------------- calculate stabilisation for "Iter-RHS" (kapome dof) */
         f2_calstabifkapome(dynvar,ele,eiforce,kapomeint,velint,velint_dc,
                            funct,derxy,fac,factor2,iel); 

/*----------------------------------------------------------------------*
 |         compute Production "Time" Force Vector                        |
 *----------------------------------------------------------------------*/ 
      if (dynvar->niturbu_pro!=0) 
      {
/*--------------------------------- get eddy-visc. at integration point */               
        f2_eddyi(&eddynint,funct,eddypro,iel);

/*------------ calculate production-term for "PROTime-RHS" (kapome-dofs)*/
        f2_production(vderxy,&production);

/*--------------- calculate galerkin part of "PROTime-RHS" (kapome-dofs)*/
        f2_calgalprofkapome(dynvar,eproforce,eddynint,funct,fac,factor1,
                            production,iel);

/*------------- calculate stabilisation for "PROTime-RHS" (kapome-dofs) */
        f2_calstabprofkapome(dynvar,ele,eproforce,eddynint,funct,fac,
                             factor1,production,velint,velint_dc,derxy,iel);
      }
/*----------------------------------------------------------------------*
 |         compute "Time" Force Vector                                  |
 *----------------------------------------------------------------------*/ 
      if (dynvar->niturbu_n!=0) 
      {
/*-------------------------------- get kapome (n) at integration point */
	    f2_kapomei(&kapomeint,funct,kapomen,iel);

/*--------------------- get kapomederivatives (n) at integration point */
          f2_kapomeder(kapomederxy,derxy,kapomen,iel);

/*-------------------- get kapomederivatives2 (n) at integration point */
          f2_kapomeder2(kapomederxy2,derxy2,kapomen,iel);

/*----------------- calculate galerkin part of "Time-RHS" (kapome-dofs)*/
          f2_calgaltfkapome(dynvar,etforce,kapomeint,velint,eddynint,funct,
                            derxy,vderxy,kapomederxy,visc,fac,factor,factor1,
                            factor2,sig,production,iel);

/*--------------- calculate stabilisation for "Time-RHS" (kapome-dofs) */
          f2_calstabtfkapome(dynvar,ele,etforce,kapomeint,velint,velint_dc,
                             eddynint,derxy,kapomederxy2,vderxy,kapomederxy,
                             visc,fac,factor,factor1,factor2,sig,production,iel);
     } /* endif dynvar->niturbu_n */
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calint_tu_1 */



#endif
