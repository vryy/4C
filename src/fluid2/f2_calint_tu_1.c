/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

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

static FLUID_DYNAMIC   *fdyn;
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
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param **emass     DOUBLE	   (o)    element mass matrix
\param  *etforce   DOUBLE	   (o)    element time force vector
\param  *eiforce   DOUBLE	   (o)    element iter force vector
\param  *eproforce DOUBLE	   (o)    element production force vector
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	 DOUBLE	   (-)    jacobian matrix
\param **xyze      DOUBLE          (-)    nodal coordinates
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param  *kapomen   DOUBLE	   (i)    kapome at time n
\param  *kapomeg   DOUBLE	   (i)    kapome at time n+g
\param  *eddyg     DOUBLE	   (i)    eddy-viscosity at time n+1
\param  *eddypro   DOUBLE	   (i)    eddy-viscosity for prod. term
\param  *kappan    DOUBLE	   (i)    kappa at time n
\param  *omega     DOUBLE	   (-)    omega at nodes
\param  *kapomepro DOUBLE	   (-)    kapome at nodes
\param  *omegaderxy   DOUBLE	   (-)    omega derivates
\param  *kapomederxy  DOUBLE	   (-)    kapome derivates
\param  *kapomederxy2 DOUBLE	   (-)    second kapome derivatives
\param  *velint    DOUBLE	   (-)    vel at integration point
\param  *velint_dc    DOUBLE	   (-)    vel at integration point for DISC. CAPT.
\param **evel	 DOUBLE	   (i)    velocity at nodes
\param **vderxy	 DOUBLE	   (-)    velocity derivates
\param **vderxy2	 DOUBLE	   (-)    velocity 2nd derivates
\param **wa1	 DOUBLE	   (-)    working array
\param **wa2	 DOUBLE	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f2_calint_tu_1(
               FLUID_DATA      *data,     
	         ELEMENT         *ele,     
	         ELEMENT         *elev, 
               DOUBLE         **estif,   
	         DOUBLE         **emass,   
	         DOUBLE          *etforce, 
	         DOUBLE          *eiforce, 
	         DOUBLE          *eproforce, 
	         DOUBLE          *funct,   
	         DOUBLE         **deriv,   
	         DOUBLE         **deriv2,  
	         DOUBLE         **xjm,     
	         DOUBLE         **xyze,
	         DOUBLE         **derxy,   
	         DOUBLE         **derxy2,  
	         DOUBLE          *kapomen,   
	         DOUBLE          *kapomeg,   
               DOUBLE          *eddyg, 
               DOUBLE          *eddypro, 
               DOUBLE          *kappan, 
               DOUBLE          *omega,    
               DOUBLE          *kapomepro,    
               DOUBLE          *omegaderxy,  
               DOUBLE          *kapomederxy,  
               DOUBLE          *kapomederxy2, 
	         DOUBLE          *velint,  
	         DOUBLE          *velint_dc,  
               DOUBLE         **evel,
               DOUBLE         **vderxy,
               DOUBLE         **vderxy2,
               DOUBLE         **wa1,     
	         DOUBLE         **wa2      
	        )
{ 
INT       i,j;        /* simply a counter                               */
INT       iel;        /* number of nodes                                */
INT       ntyp;       /* element type: 1 - quad; 2 - tri                */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration nodesin r,s direction    */
INT       actmat;     /* material number of the element                 */
INT       icode=2;    /* flag for eveluation of shape functions         */     
INT       ihoel=0;    /* flag for higher order elements                 */
INT       lr, ls;     /* counter for integration                        */
DOUBLE    dens;       /* density                                        */
DOUBLE    visc;       /* viscosity                                      */
DOUBLE    xi;
DOUBLE    fac;
DOUBLE    factor,factor1;
DOUBLE    factor2,sig,production;
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    eddyint;
DOUBLE    eddynint;
DOUBLE    kappanint;
DOUBLE    omegaint;
DOUBLE    kapomeint;
DOUBLE    kapomenint;
DOUBLE    ome_proint;
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

fdyn   = alldyn[genprob.numff].fdyn;

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
      f2_vel_dc_1(velint,velint_dc,kapomederxy);

/*------------------------------- calculate factors for kappa equation */       
      if(fdyn->kapomega_flag==0)
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
      if(fdyn->kapomega_flag==1)
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
       f2_calkkapome(estif,kapomeint,velint,eddyint,
                     funct,derxy,fac,visc,factor,sig,iel);
       
/*---------------------------------------------- compute matrix Mkapome */
	 if (fdyn->nis==0)	  	 	    
          f2_calmkapome(emass,funct,fac,iel);
  
/*----------------------------------------------------------------------*
 |         compute Stabilisation matrices                               |
 | NOTE:                                                                |
 |  Stabilisation matrices are all stored in one matrix "estif"         |
 |  Stabilisation mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
/*------------------------------------ stabilisation for matrix Kkapome */
       f2_calstabkkapome(ele,elev,estif,kapomeint,velint,velint_dc,
                         eddyint,funct,derxy,derxy2,fac,visc,factor,sig,iel); 

/*------------------------------------ stabilisation for matrix Mkapome */
       if (fdyn->nis==0) 
          f2_calstabmkapome(ele,emass,velint,velint_dc,funct,derxy,
                            fac,iel); 
/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 *----------------------------------------------------------------------*/ 

/*------------- calculate galerkin part of "Iter-RHS" (kapome dof) */
         f2_calgalifkapome(eiforce,kapomeint,funct,fac,factor2,iel);

/*---------------- calculate stabilisation for "Iter-RHS" (kapome dof) */
         f2_calstabifkapome(ele,eiforce,kapomeint,velint,velint_dc,
                            funct,derxy,fac,factor2,iel); 

/*----------------------------------------------------------------------*
 |         compute Production "Time" Force Vector                        |
 *----------------------------------------------------------------------*/ 
      if (fdyn->niturbu_pro!=0) 
      {
/*--------------------------------- get eddy-visc. at integration point */               
        f2_eddyi(&eddynint,funct,eddypro,iel);

/*------------ calculate production-term for "PROTime-RHS" (kapome-dofs)*/
        f2_production(vderxy,&production);

/*--------------- calculate galerkin part of "PROTime-RHS" (kapome-dofs)*/
        f2_calgalprofkapome(eproforce,eddynint,funct,fac,factor1,
                            production,iel);

/*------------- calculate stabilisation for "PROTime-RHS" (kapome-dofs) */
        f2_calstabprofkapome(ele,eproforce,eddynint,funct,fac,
                             factor1,production,velint,velint_dc,derxy,iel);
      }
/*----------------------------------------------------------------------*
 |         compute "Time" Force Vector                                  |
 *----------------------------------------------------------------------*/ 
      if (fdyn->niturbu_n!=0) 
      {
/*-------------------------------- get kapome (n) at integration point */
	    f2_kapomei(&kapomeint,funct,kapomen,iel);

/*--------------------- get kapomederivatives (n) at integration point */
          f2_kapomeder(kapomederxy,derxy,kapomen,iel);

/*-------------------- get kapomederivatives2 (n) at integration point */
          f2_kapomeder2(kapomederxy2,derxy2,kapomen,iel);

/*----------------- calculate galerkin part of "Time-RHS" (kapome-dofs)*/
          f2_calgaltfkapome(etforce,kapomeint,velint,eddynint,funct,
                            derxy,vderxy,kapomederxy,visc,fac,factor,factor1,
                            factor2,sig,production,iel);

/*--------------- calculate stabilisation for "Time-RHS" (kapome-dofs) */
          f2_calstabtfkapome(ele,etforce,kapomeint,velint,velint_dc,
                             eddynint,derxy,kapomederxy2,vderxy,kapomederxy,
                             visc,fac,factor,factor1,factor2,sig,production,iel);
     } /* endif fdyn->niturbu_n */
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calint_tu_1 */



#endif
