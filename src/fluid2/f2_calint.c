/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
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
\param  *hasext    int             (i)    element flag
\param **estif     double	   (o)    element stiffness matrix
\param **emass     double	   (o)    element mass matrix
\param  *etforce   double	   (o)    element time force vector
\param  *eiforce   double	   (o)    element iter force vector
\param **xyze      double          (-)    nodal coordinates
\param  *funct     double	   (-)    natural shape functions
\param **deriv     double	   (-)	  deriv. of nat. shape funcs
\param **deriv2    double	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   double	   (-)    jacobian matrix
\param **derxy     double	   (-)	  global derivatives
\param **derxy2    double	   (-)    2nd global derivatives
\param **eveln     double	   (i)    ele vel. at time n
\param **evelng    double	   (i)    ele vel. at time n+g
\param  *epren     double	   (-)    ele pres. at time n
\param  *edeadn    double	   (-)    ele dead load (selfweight) at n 
\param  *edeadng   double	   (-)    ele dead load (selfweight) at n+1
\param  *velint    double	   (-)    vel at integration point
\param  *vel2int   double	   (-)    vel at integration point
\param  *covint    double	   (-)    conv. vel. at integr. point
\param **vderxy    double	   (-)    global vel. derivatives
\param  *pderxy    double	   (-)    global pres. derivatives
\param **vderxy2   double	   (-)    2nd global vel. deriv.
\param **wa1	   double	   (-)    working array
\param **wa2	   double	   (-)    working array
\param   visc        double	   (-)    viscosity
\return void                                                   

------------------------------------------------------------------------*/
void f2_calint(
               FLUID_DATA      *data,     
	       ELEMENT         *ele,     
	       FLUID_DYN_CALC  *dynvar, 
               int             *hasext,
               double         **estif,   
	       double         **emass,   
	       double          *etforce, 
	       double          *eiforce, 
	       double         **xyze,
	       double          *funct,   
	       double         **deriv,   
	       double         **deriv2,  
	       double         **xjm,     
	       double         **derxy,   
	       double         **derxy2,  
	       double         **eveln,   
	       double         **evelng,  
	       double          *epren,   
	       double          *edeadn,
	       double          *edeadng,	        	       
	       double          *velint,  
	       double          *vel2int, 
	       double          *covint,  
	       double         **vderxy,  
	       double          *pderxy,  
	       double         **vderxy2, 
	       double         **wa1,     
	       double         **wa2,      
             double           visc      
	      )
{ 
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
      f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
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
         f2_calkvv(ele,dynvar,estif,velint,NULL,
	           vderxy,funct,derxy,fac,visc,iel);
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
            f2_calelesize2(ele,dynvar,xyze,funct,velint,wa1,visc,iel,ntyp);
/*------------------------------------ compute second global derivative */ 
         if (ihoel!=0)
            f2_gder2(ele,xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
   
         if (dynvar->nie==0)
         {
/*---------------------------------------- stabilisation for matrix Kvv */
            f2_calstabkvv(ele,dynvar,estif,velint,velint,NULL,vderxy,
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
               f2_calstabkpv(ele,dynvar,estif,velint,NULL,vderxy,
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
 |         compute "Iteration" Force Vectors                            |
 |      (for Newton iteration and for fixed-point iteration)            |
 |  at the moment only fixed-point-like iteration is used so this part  |
 |        was not checked yet!!! Thanx for your help!!!!                |
 *----------------------------------------------------------------------*/ 
      if (dynvar->nii!=0)
      {
         dserror("Iteration RHS not checked yet!\n");
/*-------------- get convective velocities (n+1,i) at integration point */
/*               covint = u*grad(u)                                     */
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
               f2_calstabifp(dynvar,&(eiforce[2*iel]),covint,derxy,fac,iel);
	 } /* endif (ele->e.f2->istabi>0) */
      } /* endif (dynvar->nii!=0) */

/*----------------------------------------------------------------------*
 |       compute "external" Force Vector (b)                            |
 |  dead load may vary over time, but stays constant over               |
 |  the whole domain --> no interpolation with shape funcs              |
 |  parts changing during the nonlinear iteration are added to          |
 |  Iteration Force Vector                                              |
 *----------------------------------------------------------------------*/
      if (*hasext!=0 && ele->e.f2->istabi>0)
      {
/*------- compute stabilisation part of external RHS (vel dofs) at (n+1)*/
         f2_calstabexfv(dynvar,ele,eiforce,derxy,derxy2,edeadng,
	                 velint,fac,visc,iel,ihoel,1); 
/*------ compute stabilisation part of external RHS (pre dofs) at (n+1) */
         if (ele->e.f2->ipres!=0)
            f2_calstabexfp(dynvar,&(eiforce[2*iel]),derxy,edeadng,fac,iel,1);  
      } /* endif (*hasext!=0 && ele->e.f2->istabi>0)  */

/*----------------------------------------------------------------------*
 |         compute "Time" Force Vector                                  |
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
         dsassert(dynvar->isemim==0,"only implicit methods implemented!\n");
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
	    dserror("TWOSTEP-method not checked yet!!!!");
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
         f2_calgaltfv(dynvar,etforce,vel2int,covint,
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
	                     vel2int,covint,pderxy,visc,fac,ihoel,iel);
         } /* endif (ele->e.f2->istabi>0) */
/*----------------------------------------------------------------------*
            | compute "external" Force Vector (b)                       |
            |  dead load may vary over time, but stays constant over    |
	    |  the whole domain --> no interpolation with shape funcs   |
            |  parts staying constant during nonlinear iteration are    |
            |  add to Time Force Vector                                 |
/*----------------------------------------------------------------------*/
         if (*hasext!=0)
	 {
/*--- compute galerkin part of external RHS (vel dofs) at (n) and (n+1) */
            f2_calgalexfv(dynvar,etforce,funct,edeadn,edeadng,fac,iel);
	    if (ele->e.f2->istabi>0)
	    {
/*-------- compute stabilisation part of external RHS (vel dofs) at (n) */
               f2_calstabexfv(dynvar,ele,etforce,derxy,derxy2,edeadn,
	                      velint,fac,visc,iel,ihoel,0);
/*--------------- compute stabilistaion part of external RHS (pre dofs) */
               if (ele->e.f2->ipres!=0)
                  f2_calstabexfp(dynvar,&(etforce[2*iel]),derxy,edeadn,fac,iel,0); 
            } /* endif (ele->e.f2->istabi>0) */
         } /* endif (*hasext!=0) */
      } /* endif  (dynvar->nif!=0) */
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
  /*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calint */

/*!---------------------------------------------------------------------
\brief integration loop for one fluid2 element for ALE

<pre>                                                         genk 10/02

In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid2 element is calculated
      
</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *ele	   ELEMENT	   (i)    actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *hasext    int             (i)    element flag
\param   imyrank   int             (i)    proc number
\param **estif     double	   (o)    element stiffness matrix
\param **emass     double	   (o)    element mass matrix
\param  *etforce   double	   (o)    element time force vector
\param  *eiforce   double	   (o)    element iter force vector
\param **xyze      double          (-)    nodal coordinates at n+theta
\param  *funct     double	   (-)    natural shape functions
\param **deriv     double	   (-)	  deriv. of nat. shape funcs
\param **deriv2    double	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   double	   (-)    jacobian matrix
\param **derxy     double	   (-)	  global derivatives
\param **derxy2    double	   (-)    2nd global derivatives
\param **eveln     double	   (i)    ele vel. at time n
\param **evelng    double	   (i)    ele vel. at time n+g
\param **ealecovn  double          (i)    ele ale-conv. vel. at time n
\param **ealecovng double          (i)    ele ale-conv. vel. at time n+1
\param **egridv    double          (i)    ele grid velocity
\param  *epren     double	   (-)    ele pres. at time n
\param  *edeadn    double	   (-)    ele dead load (selfweight) at n 
\param  *edeadng   double	   (-)    ele dead load (selfweight) at n+1
\param  *velint    double	   (-)    vel at integration point
\param  *vel2int   double	   (-)    vel at integration point
\param  *alecovint double          (-)    ale-conv. vel. at integr. point
\param  *gridvint  double          (-)    grid-velocity at integr. point
\param  *covint    double	   (-)    conv. vel. at integr. point
\param **vderxy    double	   (-)    global vel. derivatives
\param  *pderxy    double	   (-)    global pres. derivatives
\param **vderxy2   double	   (-)    2nd global vel. deriv.
\param  *ekappan   double          (i)    nodal curvature at n
\param  *ekappang  double          (i)    nodal curvature at n+g
\param **wa1	   double	   (-)    working array
\param **wa2	   double	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f2_calinta(
                FLUID_DATA      *data,     
	        ELEMENT         *ele,     
	        FLUID_DYN_CALC  *dynvar, 
                int             *hasext,
                int              imyrank,
                double         **estif,   
	        double         **emass,   
	        double          *etforce, 
	        double          *eiforce, 
	        double         **xyze,
	        double          *funct,   
	        double         **deriv,   
	        double         **deriv2,  
	        double         **xjm,     
	        double         **derxy,   
	        double         **derxy2,  
	        double         **eveln,   
	        double         **evelng, 
		double         **ealecovn,
		double         **ealecovng,
		double         **egridv, 
	        double          *epren,   
	        double          *edeadn,
	        double          *edeadng,	        	       
	        double          *velint,  
	        double          *vel2int,
		double          *alecovint,
		double          *gridvint, 
	        double          *covint,  
	        double         **vderxy,  
	        double          *pderxy,  
	        double         **vderxy2,
		double          *ekappan,
		double          *ekappang, 
	        double         **wa1,     
	        double         **wa2      
	      )
{ 
#ifdef D_FSI
int       i,k;        /* simply a counter                               */
int       iel;        /* number of nodes                                */
int       ntyp;       /* element type: 1 - quad; 2 - tri                */
int       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
int       nir,nis,nil;/* number of integration nodesin r,s direction    */
int       actmat;     /* material number of the element                 */
int       ihoel=0;    /* flag for higher order elements                 */
int       icode=2;    /* flag for eveluation of shape functions         */     
int       lr, ls;     /* counter for integration                        */
int       foundline;
int       ngline;
int       node;
double    dens;       /* density                                        */
double    visc;       /* viscosity                                      */
double    fac;
double    facr, facs; /* integration weights                            */
double    det;        /* determinant of jacobian matrix                 */
double    e1,e2;      /* natural coordinates of integr. point           */
double    preint;     /* pressure at integration point                  */
double    vn[2];      /* component of normal vector                     */
double    sigmaint;   /* surface tension                                */
double    gamma;      /* surface tension coeficient                     */
DIS_TYP   typ;	      /* element type                                   */
int       iedgnod[MAXNOD_F2];
int       line,ngnode;
GLINE    *gline[4];
FLUID_FREESURF_CONDITION *linefs[4];

#ifdef DEBUG 
dstrc_enter("f2_calinta");
#endif

/*----------------------------------------------------- initialisation */
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
gamma= mat[actmat].m.fluid->gamma;
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
      f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;
/*------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);
/*------------------------- get velocities (n+g,i) at integraton point */
      f2_veli(velint,funct,evelng,iel);
/*-------- get ale-convectcive velocities (n+g,i) at integration point */
      f2_veli(alecovint,funct,ealecovng,iel);
/*----------------------------- get grid velocity at integration point */
      f2_veli(gridvint,funct,egridv,iel);
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

         f2_calkvv(ele,dynvar,estif,velint,gridvint,
	           vderxy,funct,derxy,fac,visc,iel);
/*------------------------------------------ compute matrix Kvp and Kpv */
	 f2_calkvp(estif,funct,derxy,fac,iel);
/*-------------------------------------------------- compute matrix Kvg */
         if (ele->e.f2->fs_on==2)
	    f2_calkvg(dynvar,estif,vderxy,funct,fac,iel);
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
            f2_calelesize2(ele,dynvar,xyze,funct,velint,wa1,visc,iel,ntyp);
/*------------------------------------ compute second global derivative */ 
         if (ihoel!=0)
            f2_gder2(ele,xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
   
         if (dynvar->nie==0)
         {
/*---------------------------------------- stabilisation for matrix Kvv */
            f2_calstabkvv(ele,dynvar,estif,velint,alecovint,gridvint,vderxy,
                          funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kvp */
            f2_calstabkvp(ele,dynvar,estif,alecovint,
                          funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kvg */
            if (ele->e.f2->fs_on==2)
	       f2_calstabkvg(ele,dynvar,estif,vderxy,funct,derxy,derxy2,
	                     alecovint,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Mvv */
            if (dynvar->nis==0) 
               f2_calstabmvv(ele,dynvar,emass,alecovint, 
	                     funct,derxy,derxy2,fac,visc,iel,ihoel); 
            if (ele->e.f2->ipres!=0)	        
            {
/*---------------------------------------- stabilisation for matrix Kpv */
               f2_calstabkpv(ele,dynvar,estif,velint,gridvint,vderxy,
	                     funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kpg */
	       if (ele->e.f2->fs_on==2)
	          f2_calstabkpg(dynvar,estif,funct,vderxy,derxy,fac,iel);
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
 |         compute "Iteration" Force Vectors                            |
 |      (for Newton iteration and for fixed-point iteration)            |
 |  at the moment only fixed-point-like iteration is used so this part  |
 |        was not checked yet!!! Thanx for your help!!!!                |
 *----------------------------------------------------------------------*/ 
      if (dynvar->nii!=0)
      {
         dserror("Iteration RHS not checked for ALE!\n");
/*-------------- get convective velocities (n+1,i) at integration point */
/*               covint = c*grad(u)                                     */                  
         f2_covi(vderxy,alecovint,covint);
/*-------------------- calculate galerkin part of "Iter-RHS" (vel dofs) */
         f2_calgalifv(dynvar,eiforce,covint,funct,fac,iel);
         if (ele->e.f2->istabi>0)
	 {
/*------------------- calculate stabilisation for "Iter-RHS" (vel dofs) */
            f2_calstabifv(dynvar,ele,eiforce,covint,alecovint,funct,
	                  derxy,derxy2,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (pre dofs) */
            if (ele->e.f2->ipres!=0)
               f2_calstabifp(dynvar,&(eiforce[2*iel]),covint,derxy,fac,iel);
	 } /* endif (ele->e.f2->istabi>0) */
      } /* endif (dynvar->nii!=0) */

/*----------------------------------------------------------------------*
 |       compute "external" Force Vector (b)                            |
 |  dead load may vary over time, but stays constant over               |
 |  the whole domain --> no interpolation with shape funcs              |
 |  parts changing during the nonlinear iteration are added to          |
 |  Iteration Force Vector                                              |
 *----------------------------------------------------------------------*/
      if (*hasext!=0 && ele->e.f2->istabi>0)
      {
/*------- compute stabilisation part of external RHS (vel dofs) at (n+1)*/
         f2_calstabexfv(dynvar,ele,eiforce,derxy,derxy2,edeadng,
	                 alecovint,fac,visc,iel,ihoel,1); 
/*------ compute stabilisation part of external RHS (pre dofs) at (n+1) */
         if (ele->e.f2->ipres!=0)
            f2_calstabexfp(dynvar,&(eiforce[2*iel]),derxy,edeadng,fac,iel,1);  
      } /* endif (*hasext!=0 && ele->e.f2->istabi>0)  */
/*----------------------------------------------------------------------*
 |         compute "Time" Force Vector                                  |
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
	 dsassert(dynvar->isemim==0,"only implicit methods implemented!\n");
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
/*-------------- get ale-convective velocities (n) at integration point */
            f2_veli(alecovint,funct,ealecovn,iel);
/*------------------- get velocitiederivatives (n) at integration point */
            f2_vder(vderxy,derxy,eveln,iel);
/*------------- get 2nd velocities derivatives (n) at integration point */
	    if (ihoel!=0)
	       f2_vder2(vderxy2,derxy2,eveln,iel);	       
	 } /* endif (dynvar->isemim==0) */
/*------------------------------------------- check for TWOSTEP-methods */	 
         dsassert(dynvar->itwost==0,"TWOSTEP-methods not implemented for ALE!\n");
/*-------------------------- velocities at different times are necassary
         for TWOSTEP-methods. At the moment they are the same!          */	 
	 vel2int=velint;
/*------------------ get convective velocities (n) at integration point 
                     covint = c * grad(u)                               */
         f2_covi(vderxy,alecovint,covint);        	    
/*--------------------- calculate galerkin part of "Time-RHS" (vel-dofs)*/
         f2_calgaltfv(dynvar,etforce,vel2int,covint,
	              funct,derxy,vderxy,preint,visc,fac,iel);
/*-------------------- calculate galerkin part of "Time-RHS" (pre-dofs) */
         f2_calgaltfp(dynvar,&(etforce[2*iel]),funct,vderxy,fac,iel); 
	 if (ele->e.f2->istabi>0)
	 {
/*------------------- calculate stabilisation for "Time-RHS" (vel-dofs) */
            f2_calstabtfv(dynvar,ele,etforce,alecovint,vel2int,
	                  covint,derxy,derxy2,vderxy,vderxy2,
		          pderxy,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Time-RHS" (pre-dofs) */
            if (ele->e.f2->ipres!=0)
	       f2_calstabtfp(dynvar,&(etforce[2*iel]),derxy,vderxy2,
	                     vel2int,covint,pderxy,visc,fac,ihoel,iel);
         } /* endif (ele->e.f2->istabi>0) */
/*----------------------------------------------------------------------*
            | compute "external" Force Vector (b)                       |
            |  dead load may vary over time, but stays constant over    |
	    |  the whole domain --> no interpolation with shape funcs   |
            |  parts staying constant during nonlinear iteration are    |
            |  add to Time Force Vector                                 |
/*----------------------------------------------------------------------*/
         if (*hasext!=0)
	 {
/*--- compute galerkin part of external RHS (vel dofs) at (n) and (n+1) */
            f2_calgalexfv(dynvar,etforce,funct,edeadn,edeadng,fac,iel);
	    if (ele->e.f2->istabi>0)
	    {
/*-------- compute stabilisation part of external RHS (vel dofs) at (n) */
               f2_calstabexfv(dynvar,ele,etforce,derxy,derxy2,edeadn,
	                      alecovint,fac,visc,iel,ihoel,0);
/*--------------- compute stabilistaion part of external RHS (pre dofs) */
               if (ele->e.f2->ipres!=0)
                  f2_calstabexfp(dynvar,&(etforce[2*iel]),derxy,edeadn,fac,iel,0); 
            } /* endif (ele->e.f2->istabi>0) */
         } /* endif (*hasext!=0) */
      } /* endif  (dynvar->nif!=0) */
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */

/*----------------------------------------------------------------------*
 |             evaluate the integrals over the element edges            |
 | NOTE:                                                                |
 |   As long as the integrals over the edges contribute to the element  |
 |   stiffness matrix, there's no problem with interprocessor elements, |
 |   since parts which are not needed, are thrown away during assembly. |
 |   At the moment there's also no problem with RHS-parts due to        | 
 |   surface tension, since the RHS is not alreduced. This may change,  |
 |   if Volker writes back his version of building the RHS              |
 *----------------------------------------------------------------------*/
if (ele->e.f2->fs_on==2) 
{                                /* element with free surface condition */         
/*-------------------------------- check for presence of freesurface */
   foundline=0;
   /*---------------------------------- number of lines to this element */
   ngline=ele->g.gsurf->ngline;
   /*------- loop over lines, check for freesurface conditions on lines */
   for (i=0; i<ngline; i++)
   {
      gline[i] = ele->g.gsurf->gline[i];
      linefs[i] = gline[i]->freesurf;
#if 0 /* see NOTE */    
      if (gline[i]->proc!=imyrank) linefs[i]=NULL;
#endif
      if(linefs[i]==NULL) continue;
      foundline++;
   }
   if (foundline==0) goto end;  
   /*--------------------------------------- set number of gauss points */
   nil = IMAX(nir,2);
   /*---------------------------------- loop over lines at free surface */
   for (line=0; line<ngline; line++)
   {
      if (linefs[line]==NULL) continue;
      /*--------------------------------- check number of nodes on line */
      ngnode = gline[line]->ngnode;
      /*------------------------------------------------ get edge nodes */
      f2_iedg(iedgnod,ele,line,0);
      /*------------------------------ integration loop on actual gline */
      for (lr=0;lr<nil;lr++)
      {
      /*---------- get values of  shape functions and their derivatives */
     	 e1   = data->qxg[lr][nil-1];
   	 facr = data->qwgt[lr][nil-1];
   	 f2_degrectri(funct,deriv,e1,typ,1);
   	 /*------------------------------- compute jacobian determinant */
     	 f2_edgejaco(xyze,funct,deriv,xjm,&det,ngnode,iedgnod);
     	 fac = det*facr;
     	 /*--------------------------------- compute matrix Kgv and Kgg */
     	 f2_calkgedge(estif,funct,fac,iedgnod,iel,ngnode);
         if (dynvar->surftens>0)
	 {
            /*----------------------------------- compute normal vector */
            vn[0]=ZERO;
	    vn[1]=ZERO;
	    for(k=0;k<ngnode;k++) 
            {
               node=iedgnod[k];
               vn[0]+=deriv[0][k]*xyze[1][node];      
               vn[1]-=deriv[0][k]*xyze[0][node];	    
            }
	    /* REMARK:
	       jacobian determinant cancels out
	        --> identical with vnorm =length of normal vector       */
            /*----------------------------------------------------------*
	     |  surface tension effects at (n)                          |
	     *----------------------------------------------------------*/
     	    if (dynvar->fsstnif!=0)
	    {
	       /*------ calculate surface tension at gauss point at (n) */
	        sigmaint = gamma*f2_kappai(funct,ekappan,iedgnod,ngnode);
	       /*------------ calculate Time-RHS due to surface tension */ 
	        f2_calsurftenfv(etforce,funct,vn,sigmaint,
	                     dynvar->thsl,facr,ngnode,iedgnod);
            }
	    /*----------------------------------------------------------*
	     |  surface tension effects at (n+1)                        |
	     *----------------------------------------------------------*/
     	    if (dynvar->fsstnii!=0)
	    {
	       /*---- calculate surface tension at gauss point at (n+1) */	 	 
	        sigmaint = gamma*f2_kappai(funct,ekappang,iedgnod,ngnode);
	       /*------------ calculate Iter-RHS due to surface tension */
	        f2_calsurftenfv(etforce,funct,vn,sigmaint,
	                        dynvar->thsr,facr,ngnode,iedgnod);	    
	    }
	 }
      }
   }	       
}
else if (ele->e.f2->fs_on==1 && dynvar->surftens>0)
{
/*----------------------------------- check for presence of freesurface */
   foundline=0;
   /*---------------------------------- number of lines to this element */
   ngline=ele->g.gsurf->ngline;
   /*------- loop over lines, check for freesurface conditions on lines */
   for (i=0; i<ngline; i++)
   {
      gline[i] = ele->g.gsurf->gline[i];
      linefs[i] = gline[i]->freesurf;
#if 0 /* see NOTE */
      if (gline[i]->proc!=imyrank) linefs[i]=NULL;
#endif
      if(linefs[i]==NULL) continue;
      foundline++;
   }
   if (foundline==0) goto end;  
   /*--------------------------------------- set number of gauss points */
   nil = IMAX(nir,2);
   /*---------------------------------- loop over lines on free surface */
   for (line=0; line<ngline; line++)
   {
      if (linefs[line]==NULL) continue;
      /*--------------------------------- check number of nodes on line */
      ngnode = gline[line]->ngnode;	    
      /*------------------------------------------------ get edge nodes */
      f2_iedg(iedgnod,ele,line,0);
      /*------------------------------ integration loop on actual gline */
      for (lr=0;lr<nil;lr++)
      {
      /*---------- get values of  shape functions and their derivatives */
     	 e1   = data->qxg[lr][nil-1];
	 facr = data->qwgt[lr][nil-1];
	 f2_degrectri(funct,deriv,e1,typ,1);
         /*-------------------------------------- compute normal vector */
         vn[0]=ZERO;
	 vn[1]=ZERO;
	 for(k=0;k<ngnode;k++) 
         {
            node=iedgnod[k];
            vn[0]+=deriv[0][k]*xyze[1][node];      
            vn[1]-=deriv[0][k]*xyze[0][node];	    
         }
	 /* REMARK:
	    jacobian determinant cancels out --> identical with vnorm
	                                       =length of normal vector */
         /*-------------------------------------------------------------*
	  |  surface tension effects at (n)                             |
	  *-------------------------------------------------------------*/
     	 /*------------ calculate surface tension at gauss point at (n) */
	  sigmaint = gamma*f2_kappai(funct,ekappan,iedgnod,ngnode);
	 /*----------------- calculate Time-RHS due to surface tension */ 
	  f2_calsurftenfv(etforce,funct,vn,sigmaint,
	                  dynvar->thsl,facr,ngnode,iedgnod);
         /*-------------------------------------------------------------*
	  |  surface tension effects at (n+1)                           |
	  *-------------------------------------------------------------*/
     	 /*---------- calculate surface tension at gauss point at (n+1) */	 	 
	  sigmaint = gamma*f2_kappai(funct,ekappang,iedgnod,ngnode);
	 /*------------------ calculate Time-RHS due to surface tension */ 
	  f2_calsurftenfv(etforce,funct,vn,sigmaint,
	                  dynvar->thsr,facr,ngnode,iedgnod);
      }
   }
}
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

#else
dserror("FSI-functions not compiled in!\n");
#endif

return; 
} /* end of f2_calinta */

#endif
/*! @} (documentation module close)*/
