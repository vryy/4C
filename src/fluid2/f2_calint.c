/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*
 | integration loop for one fluid element                               |
 |                                                                      |
 |                                                           genk 03/02 |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------
\brief integration loop for one fluid2 element

<pre>                                                         genk 04/02

In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid2 element is calculated
      
</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *ele	   ELEMENT	   (i)    actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param **estif     double	   (o)    element stiffness matrix
\param **emass     double	   (o)    element mass matrix
\param  *etforce   double	   (o)    element time force vector
\param  *eiforce   double	   (o)    element iter force vector
\param  *funct     double	   (-)    natural shape functions
\param **deriv     double	   (-)	  deriv. of nat. shape funcs
\param **deriv2    double	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   double	   (-)    jacobian matrix
\param **derxy     double	   (-)	  global derivatives
\param **derxy2    double	   (-)    2nd global derivatives
\param **eveln     double	   (i)    ele vel. at time n
\param **evelng    double	   (i)    ele vel. at time n+g
\param  *epren     double	   (-)    ele pres. at time n
\param  *velint    double	   (-)    vel at integration point
\param  *vel2int   double	   (-)    vel at integration point
\param  *covint    double	   (-)    conv. vel. at integr. point
\param **vderxy    double	   (-)    global vel. derivatives
\param  *pderxy    double	   (-)    global pres. derivatives
\param **vderxy2   double	   (-)    2nd global vel. deriv.
\param **wa1	   double	   (-)    working array
\param **wa2	   double	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f2_calint(
               FLUID_DATA      *data,     
	       ELEMENT         *ele,     
	       FLUID_DYN_CALC  *dynvar, 
               double         **estif,   
	       double         **emass,   
	       double          *etforce, 
	       double          *eiforce, 
	       double          *funct,   
	       double         **deriv,   
	       double         **deriv2,  
	       double         **xjm,     
	       double         **derxy,   
	       double         **derxy2,  
	       double         **eveln,   
	       double         **evelng,  
	       double          *epren,   
	       double          *velint,  
	       double          *vel2int, 
	       double          *covint,  
	       double         **vderxy,  
	       double          *pderxy,  
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
int       ihoel=0;    /* flag for higher order elements                 */
int       icode=2;    /* flag for eveluation of shape functions         */     
int       lr, ls;     /* counter for integration                        */
double    dens;       /* density                                        */
double    visc;       /* viscosity                                      */
double    fac;
double    facr, facs; /* integration weights                            */
double    det;        /* determinant of jacobian matrix                 */
double    e1,e2;      /* natural coordinates of integr. point           */
double    preint;     /* pressure at integration point                  */
DIS_TYP   typ;	      /* element type                                   */

#ifdef DEBUG 
dstrc_enter("f2_calint");
#endif

/*----------------------------------------------------- initialisation */
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
ntyp = ele->e.f2->ntyp; 
typ  = ele->distyp;

/*------- get integraton data and check if elements are "higher order" */
switch (ntyp)
{
case 1:  /* --> quad - element */
   icode   = 3;
   ihoel   = 1;
   /* initialise integration */
   nir = ele->e.f2->nGP[0];
   nis = ele->e.f2->nGP[1];
break;
case 2: /* --> tri - element */  
   if (iel>3)
   {
      icode   = 3;
      ihoel   = 1;
   }
   /* initialise integration */
   nir  = ele->e.f2->nGP[0];
   nis  = 1;
   intc = ele->e.f2->nGP[1];  
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
      f2_jaco(funct,deriv,xjm,&det,ele,iel);
      fac = facr*facs*det;
/*------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);
/*------------------------- get velocities (n+g,i) at integraton point */
      f2_veli(velint,funct,evelng,iel);
/*-------------- get velocity (n+g,i) derivatives at integration point */
      f2_vder(vderxy,derxy,evelng,iel);
      
/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices                         |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/
      if(dynvar->nik>0)
      {
/*-------------------------------------------------- compute matrix Kvv */      
         f2_calkvv(dynvar,estif,velint,vderxy,funct,derxy,fac,visc,iel);
/*------------------------------------------ compute matrix Kvp and Kpv */
	 f2_calkvp(estif,funct,derxy,fac,iel);
/*-------------------------------------------------- compute matrix Mvv */
	 if (dynvar->nis==0)	  	 	    
            f2_calmvv(emass,funct,fac,iel);
      } /* endif (dynvar->nik>0) */
      
/*----------------------------------------------------------------------*
 |         compute Stabilisation matrices                               |
 | NOTE:                                                                |
 |  Stabilisation matrices are all stored in one matrix "estif"         |
 |  Stabilisation mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
      if (ele->e.f2->istabi>0)
      { 
 /*-------------- compute stabilisation parameter during ntegration loop*/
         if (ele->e.f2->iduring!=0)
            f2_calelesize2(ele,dynvar,funct,velint,wa1,visc,iel,ntyp);
/*------------------------------------ compute second global derivative */ 
         if (ihoel!=0)
            f2_gder2(ele,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
   
         if (dynvar->nie==0)
         {
/*---------------------------------------- stabilisation for matrix Kvv */
            f2_calstabkvv(ele,dynvar,estif,velint,vderxy,
                          funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kvp */
            f2_calstabkvp(ele,dynvar,estif,velint,
                          funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Mvv */
            if (dynvar->nis==0) 
               f2_calstabmvv(ele,dynvar,emass,velint, 
	                     funct,derxy,derxy2,fac,visc,iel,ihoel); 
            if (ele->e.f2->ipres!=0)	        
            {
/*---------------------------------------- stabilisation for matrix Kpv */
               f2_calstabkpv(dynvar,estif,velint,vderxy,
	                     funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Mpv */
	       if (dynvar->nis==0)
		  f2_calstabmpv(dynvar,emass,funct,derxy,fac,iel);
            } /* endif (ele->e.f2->ipres!=0) */
         } /* endif (dynvar->nie==0) */
/*---------------------------------------- stabilisation for matrix Kpp */
         if (ele->e.f2->ipres!=0)
	    f2_calstabkpp(dynvar,estif,derxy,fac,iel);  
      } /* endif (ele->e.f2->istabi>0) */
 
/*----------------------------------------------------------------------*
 |       compute "external" Force Vector                                |
 |   (at the moment there are no external forces implemented)           |
 |  but there can be due to self-weight /magnetism / etc. (b)           |
 *----------------------------------------------------------------------*/
/*-------------------- compute galerkin part of external RHS (vel dofs) *
         f2_calgalexfv();
/*--------------- compute stabilisation part of external RHS (vel dofs) *
         f2_calstabexfv();
/*--------------- compute stabilistaion part of external RHS (pre dofs) *
         f2_calgalexfp(); 
*/

/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 |      (for Newton iteration and for fixed-point iteration)            |
 *----------------------------------------------------------------------*/ 
      if (dynvar->nii!=0)
      {
/*-------------- get convective velocities (n+1,i) at integration point */
         f2_covi(vderxy,velint,covint);
/*-------------------- calculate galerkin part of "Iter-RHS" (vel dofs) */
         f2_calgalifv(dynvar,eiforce,covint,funct,fac,iel);
         if (ele->e.f2->istabi>0)
	 {
/*------------------- calculate stabilisation for "Iter-RHS" (vel dofs) */
            f2_calstabifv(dynvar,ele,eiforce,covint,velint,funct,
	                  derxy,derxy2,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (pre dofs) */
            if (ele->e.f2->ipres!=0)
               f2_calstabifp(dynvar,eiforce,covint,derxy,fac,iel);
	 } /* endif (ele->e.f2->istabi>0) */
      } /* endif (dynvar->nii!=0) */
    
/*----------------------------------------------------------------------*
 |         compute "Time" Force Vectors                                 |
 *----------------------------------------------------------------------*/
      if (dynvar->nif!=0)
      {
         if (dynvar->iprerhs>0)
	 {
/*------------------------------- get pressure (n) at integration point */
            f2_prei(&preint,funct,epren,iel);
/*------------------- get pressure derivatives (n) at integration point */
            f2_pder(pderxy,derxy,epren,iel);
	 } /* endif (dynvar->iprerhs>0) */
	 if (dynvar->isemim==0)
	 {
	    
	    /* in all but semi-implicit cases (n+gamma_bar) = (n)
	       --> hence we need the values according to u(n)
	       NOTE: since "time forces" are only calculated in the first
	       iteration step and in general U(n+1,0) =  U(n) - with only
	       exception being the dirichlet values - the stability 
	       parameter are not calculated for the field at (n) -
	       instead the ones from the field (n+1,0) are taken
	       (shouldn't be too much of a difference)!!!               */
	    
/*----------------------------- get velocities (n) at integration point */
	    f2_veli(velint,funct,eveln,iel);
/*------------------- get velocitiederivatives (n) at integration point */
            f2_vder(vderxy,derxy,eveln,iel);
/*------------- get 2nd velocities derivatives (n) at integration point */
	    if (ihoel!=0)
	       f2_vder2(vderxy2,derxy2,eveln,iel);	       
	 } /* endif (dynvar->isemim==0) */
	 if (dynvar->itwost!=0)
	 {
	    
	    /* for two-step methods we have values at two different times
	       involved in the computation of the time forces
	       --> velint  = U(n+g) from above;
	           vel2int = U(n) get now;                              */
	    
/*----------------------------- get velocities (n) at integration point */
	    f2_veli(vel2int,funct,eveln,iel);
/*------------------- get velocitiederivatives (n) at integration point */
            f2_vder(vderxy,derxy,eveln,iel);
/*------------- get 2nd velocities derivatives (n) at integration point */
	    if (ihoel!=0)
	       f2_vder2(vderxy2,derxy2,eveln,iel);            
	 } /* endif (dynvar->itwost!=0) */
	 if (dynvar->itwost==0)
	    vel2int=velint;
/*------------------ get convective velocities (n) at integration point */
         f2_covi(vderxy,velint,covint);        	    
/*--------------------- calculate galerkin part of "Time-RHS" (vel-dofs)*/
         f2_calgaltfv(dynvar,etforce,velint,vel2int,covint,
	              funct,derxy,vderxy,preint,visc,fac,iel);
/*-------------------- calculate galerkin part of "Time-RHS" (pre-dofs) */
         f2_calgaltfp(dynvar,&(etforce[2*iel]),funct,vderxy,fac,iel);
	 if (ele->e.f2->istabi>0)
	 {
/*------------------- calculate stabilisation for "Time-RHS" (vel-dofs) */
            f2_calstabtfv(dynvar,ele,etforce,velint,vel2int,
	                  covint,derxy,derxy2,vderxy,vderxy2,
		          pderxy,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Time-RHS" (pre-dofs) */
            if (ele->e.f2->ipres!=0)
	       f2_calstabtfp(dynvar,&(etforce[2*iel]),derxy,vderxy2,
	                     velint,covint,pderxy,visc,fac,ihoel,iel);
         } /* endif (ele->e.f2->istabi>0) */
      } /* endif  (dynvar->nif!=0) */
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calint */



#endif
