/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"

#ifdef DEBUG
void genkout(DOUBLE **matrix, DOUBLE *vector, char *title, INT ntitle,
             INT istart, INT iend, INT jstart, INT jend, INT flag);
#endif

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*!---------------------------------------------------------------------
\brief integration loop for one fluid2 element

<pre>                                                         genk 05/02

In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid2 element is calculated
      
</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *ele	   ELEMENT	   (i)    actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  *hasext    INT             (i)    element flag
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param **emass     DOUBLE	   (o)    element mass matrix
\param  *etforce   DOUBLE	   (o)    element time force vector
\param  *eiforce   DOUBLE	   (o)    element iter force vector
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param **eveln     DOUBLE	   (i)    ele vel. at time n
\param **evelng    DOUBLE	   (i)    ele vel. at time n+g
\param  *epren     DOUBLE	   (-)    ele pres. at time n
\param  *edeadn    DOUBLE	   (-)    ele dead load (selfweight) at n 
\param  *edeadng   DOUBLE	   (-)    ele dead load (selfweight) at n+1
\param  *velint    DOUBLE	   (-)    vel at integration point
\param  *vel2int   DOUBLE	   (-)    vel at integration point
\param  *covint    DOUBLE	   (-)    conv. vel. at integr. point
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param  *pderxy    DOUBLE	   (-)    global pres. derivatives
\param **vderxy2   DOUBLE	   (-)    2nd global vel. deriv.
\param **wa1	   DOUBLE	   (-)    working array
\param **wa2	   DOUBLE	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_calint(
               FLUID_DATA      *data, 
	       ELEMENT         *ele,
	       FLUID_DYN_CALC  *dynvar,
               INT             *hasext,
               DOUBLE         **estif,
	       DOUBLE         **emass,
	       DOUBLE          *etforce,
	       DOUBLE          *eiforce,
	       DOUBLE          *funct,
	       DOUBLE         **deriv,
	       DOUBLE         **deriv2,
	       DOUBLE         **xjm,
	       DOUBLE         **derxy,
	       DOUBLE         **derxy2,
	       DOUBLE         **eveln,
	       DOUBLE         **evelng,
	       DOUBLE          *epren,
	       DOUBLE          *edeadn,
	       DOUBLE          *edeadng,
	       DOUBLE          *velint,
	       DOUBLE          *vel2int,
	       DOUBLE          *covint,
	       DOUBLE         **vderxy,
	       DOUBLE          *pderxy,
	       DOUBLE         **vderxy2,
	       DOUBLE         **wa1,
	       DOUBLE         **wa2
	      )
{ 
INT      iel;	      /* number of nodes 			        */
INT      ntyp;        /* element type: 1 - hex; 2 - tet  	        */
INT      intc;        /* "integration case" for tet for further infos
         		 see f3_inpele.c and f3_intg.c                  */
INT      nir,nis,nit; /* number of integration nodes in r,s,t direction */
INT      actmat;      /* material number of the element                 */
INT      ihoel=0;     /* flag for higher order elements                 */
INT      icode=2;     /* flag for eveluation of shape functions         */   
INT      lr, ls, lt;  /* counter for integration                        */
DOUBLE   dens;        /* density                                        */
DOUBLE   visc;        /* viscosity                                      */
DOUBLE   fac;
DOUBLE   facr,facs,fact; /* integration weights                         */
DOUBLE   det;	      /* determinant of jacobian matrix                 */
DOUBLE   e1,e2,e3;    /* natural coordinates of integr. point           */
DOUBLE   preint;      /*pressure at integration point */
DIS_TYP  typ;         /* element type                                   */

STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_calint");
#endif		

/*----------------------------------------------------- initialisation -*/
gls    = ele->e.f3->stabi.gls;
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
ntyp = ele->e.f3->ntyp; 
typ  = ele->distyp;

/*------------------------------ check for proper stabilisation mode ---*/
if (ele->e.f3->stab_type != stab_gls) 
   dserror("routine with no or wrong stabilisation called");
   
/*------- get integraton data and check if elements are "higher order" -*/
switch (ntyp)
{
case 1:  /* --> hex - element */
   icode   = 3;
   ihoel   = 1;
   /* initialise integration */
   nir = ele->e.f3->nGP[0];
   nis = ele->e.f3->nGP[1];
   nit = ele->e.f3->nGP[2];
   break;
case 2: /* --> tet - element */  
   if (iel>4)
   {
      icode   = 3;
      ihoel   = 1;
   }
   /* initialise integration */
   nir  = ele->e.f3->nGP[0];
   nis  = 1;
   nit  = 1; 
   intc = ele->e.f3->nGP[1];  
   break;
default:
   dserror("ntyp unknown!");
} /* end switch (ntyp) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{
for (ls=0;ls<nis;ls++)
{
for (lt=0;lt<nit;lt++)
{
/*--------------- get values of  shape functions and their derivatives */
   switch(ntyp)  
   {
   case 1:   /* --> hex - element */
      e1   = data->qxg[lr][nir-1];
      facr = data->qwgt[lr][nir-1];
      e2   = data->qxg[ls][nis-1];
      facs = data->qwgt[ls][nis-1];
      e3   = data->qxg[lt][nit-1];
      fact = data->qwgt[lt][nit-1];
      f3_hex(funct,deriv,deriv2,e1,e2,e3,typ,icode);
      break;
   case 2:   /* --> tet - element */		  	
      e1   = data->txgr[lr][intc];
      facr = data->twgt[lr][intc];
      e2   = data->txgs[lr][intc];
      facs = ONE;
      e3   = data->txgt[lr][intc]; 
      fact = ONE;
      f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,icode); 
      break;
   default:
      dserror("ntyp unknown!");
   } /* end switch (ntyp) */
/*-------------------------------------------- compute Jacobian matrix */  
   f3_jaco(funct,deriv,xjm,&det,ele,iel);
   fac = facr*facs*fact*det;
/*------------------------------------------- compute global derivates */
   f3_gder(derxy,deriv,xjm,wa1,det,iel);
/*------------------------- get velocities (n+g,i) at integraton point */
   f3_veli(velint,funct,evelng,iel);
/*-------------- get velocity (n+g,i) derivatives at integration point */
   f3_vder(vderxy,derxy,evelng,iel);

/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices                         |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/
   if(dynvar->nik>0)
   {
/*-------------------------------------------------- compute matrix Kvv */
      f3_calkvv(dynvar,estif,velint,vderxy,funct,derxy,fac,visc,iel);
/*------------------------------------------ compute matrix Kvp and Kpv */
      f3_calkvp(estif,funct,derxy,fac,iel);
/*-------------------------------------------------- compute matrix Mvv */
      if (dynvar->nis==0)	  	 	    
         f3_calmvv(emass,funct,fac,iel);
   } /* endif (dynvar->nik>0) */
/*----------------------------------------------------------------------*
 |         compute Stabilisation matrices                               |
 | NOTE:                                                                |
 |  Stabilisation matrices are all stored in one matrix "estif"         |
 |  Stabilisation mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
   if (gls->istabi>0)
   { 
 /*-------------- compute stabilisation parameter during ntegration loop*/
      if (gls->iduring!=0)
         f3_calelesize2(ele,dynvar,velint,wa1,visc,iel,ntyp);
/*------------------------------------ compute second global derivative */ 
      if (ihoel!=0)
         f3_gder2(ele,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
   
      if (dynvar->nie==0)
      {
/*---------------------------------------- stabilisation for matrix Kvv */
         f3_calstabkvv(ele,dynvar,estif,velint,vderxy,
                          funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kvp */
         f3_calstabkvp(ele,dynvar,estif,velint,
                       funct,derxy,derxy2,fac,visc,iel,ihoel); 
/*---------------------------------------- stabilisation for matrix Mvv */
         if (dynvar->nis==0) 
            f3_calstabmvv(ele,dynvar,emass,velint,
	                     funct,derxy,derxy2,fac,visc,iel,ihoel);
         if (gls->ipres!=0)	        
         {
/*---------------------------------------- stabilisation for matrix Kpv */ 
            f3_calstabkpv(dynvar,estif,velint,vderxy,
	                  funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Mpv */
	    if (dynvar->nis==0)
	       f3_calstabmpv(dynvar,emass,funct,derxy,fac,iel);
         } /* endif (ele->e.f3->ipres!=0) */
      } /* endif (dynvar->nie==0) */
/*---------------------------------------- stabilisation for matrix Kpp */ 
      if (gls->ipres!=0)
	 f3_calstabkpp(dynvar,estif,derxy,fac,iel);  
   } /* endif (ele->e.f3->istabi>0) */ 

/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 |      (for Newton iteration and for fixed-point iteration)            |
 *----------------------------------------------------------------------*/
   if (dynvar->nii!=0)
   {
/*-------------- get convective velocities (n+1,i) at integration point */
      f3_covi(vderxy,velint,covint);
/*-------------------- calculate galerkin part of "Iter-RHS" (vel dofs) */
      f3_calgalifv(dynvar,eiforce,covint,funct,fac,iel);
      if (gls->istabi>0)
      {
/*------------------- calculate stabilisation for "Iter-RHS" (vel dofs) */
         f3_calstabifv(dynvar,ele,eiforce,covint,velint,funct,
	               derxy,derxy2,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (pre dofs) */
         if (gls->ipres!=0)
            f3_calstabifp(dynvar,&(eiforce[3*iel]),covint,derxy,fac,iel);
      } /* endif (ele->e.f3->istabi>0) */
   } /* endif (dynvar->nii!=0) */

/*----------------------------------------------------------------------*
 |       compute "external" Force Vector                                |
 |   (at the moment there are no external forces implemented)           |
 |  but there can be due to self-weight /magnetism / etc. (b)           |
 |  dead load may vary over time, but stays constant over               |
 |  the whole domain --> no interpolation with shape funcs              |
 |  parts changing during the nonlinear iteration are added to          |
 |  Iteration Force Vector                                              |
 *----------------------------------------------------------------------*/
   if (*hasext!=0 && gls->istabi>0)
   {
/*------- compute stabilisation part of external RHS (vel dofs) at (n+1)*/
      f3_calstabexfv(dynvar,ele,eiforce,derxy,derxy2,edeadng,
	              velint,fac,visc,iel,ihoel,1); 
/*------ compute stabilisation part of external RHS (pre dofs) at (n+1) */
      if (gls->ipres!=0)
          f3_calstabexfp(dynvar,&(eiforce[3*iel]),derxy,edeadng,fac,iel,1); 
   } /* endif (*hasext!=0 && ele->e.f3->istabi>0) */

/*----------------------------------------------------------------------*
 |         compute "Time" Force Vectors                                 |
 *----------------------------------------------------------------------*/
   if (dynvar->nif!=0)
   {
      if (dynvar->iprerhs>0)
      {
/*------------------------------- get pressure (n) at integration point */
         f3_prei(&preint,funct,epren,iel);
/*------------------- get pressure derivatives (n) at integration point */
         f3_pder(pderxy,derxy,epren,iel);
      } /* endif (dynvar->iprerhs>0) */

      /* in all but semi-implicit cases (n+gamma_bar) = (n)
	 --> hence we need the values according to u(n)
	 NOTE: since "time forces" are only calculated in the first
	 iteration step and in general U(n+1,0) =  U(n) - with only
	 exception being the dirichlet values - the stability 
	 parameter are not calculated for the field at (n) -
	 instead the ones from the field (n+1,0) are taken
	 (shouldn't be too much of a difference)!!!		  */
	    
/*----------------------------- get velocities (n) at integration point */
      f3_veli(velint,funct,eveln,iel);
/*------------------- get velocitiederivatives (n) at integration point */
      f3_vder(vderxy,derxy,eveln,iel);
/*------------- get 2nd velocities derivatives (n) at integration point */
      if (ihoel!=0)
	 f3_vder2(vderxy2,derxy2,eveln,iel);	       
/*---------------- due to historical reasons there exist two velocities */
      vel2int=velint;
/*------------------ get convective velocities (n) at integration point */
      f3_covi(vderxy,velint,covint);        	    
/*--------------------- calculate galerkin part of "Time-RHS" (vel-dofs)*/
      f3_calgaltfv(dynvar,etforce,velint,vel2int,covint,
	           funct,derxy,vderxy,preint,visc,fac,iel);
/*-------------------- calculate galerkin part of "Time-RHS" (pre-dofs) */
      f3_calgaltfp(dynvar,&(etforce[3*iel]),funct,vderxy,fac,iel);
      if (gls->istabi>0)
      {
/*------------------- calculate stabilisation for "Time-RHS" (vel-dofs) */
         f3_calstabtfv(dynvar,ele,etforce,velint,vel2int,
	               covint,derxy,derxy2,vderxy,vderxy2,
		       pderxy,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Time-RHS" (pre-dofs) */
         if (gls->ipres!=0)
	    f3_calstabtfp(dynvar,&(etforce[3*iel]),derxy,vderxy2,
	                  velint,covint,pderxy,visc,fac,ihoel,iel);
      } /* endif (ele->e.f3->istabi>0) */
/*----------------------------------------------------------------------*
            | compute "external" Force Vector                           |
            | (at the moment there are no external forces implemented)  |
            |  but there can be due to self-weight /magnetism / etc. (b)|
            |  dead load may vary over time, but stays constant over    |
	    |  the whole domain --> no interpolation with shape funcs   |
            |  parts staying constant during nonlinear iteration are    |
            |  add to Time Force Vector                                 |
 *----------------------------------------------------------------------*/
      if (*hasext!=0)
      {
/*--- compute galerkin part of external RHS (vel dofs) at (n) and (n+1) */
         f3_calgalexfv(dynvar,etforce,funct,edeadn,edeadng,fac,iel);
	 if (gls->istabi>0)
	 {
/*-------- compute stabilisation part of external RHS (vel dofs) at (n) */
            f3_calstabexfv(dynvar,ele,etforce,derxy,derxy2,edeadn,
	                   velint,fac,visc,iel,ihoel,0);
/*--------------- compute stabilistaion part of external RHS (pre dofs) */
            if (gls->ipres!=0)
               f3_calstabexfp(dynvar,&(etforce[3*iel]),derxy,edeadn,fac,iel,0); 
         } /* endif (ele->e.f3->istabi>0) */
      } /* endif (*hasext!=0) */
   } /* endif (dynvar->nif!=0)   */
}
}
} /* end of loop over integration points */

/* TEST - TEST - TEST - TEST - TEST - TEST - TEST - TEST - TEST - TEST */
/* output to the screen for testing */
#ifdef DEBUG

#if 0
if (ele->Id==0)
{
printf("\n");
printf("ELEMENT NUMBER ===== %d \n",ele->Id_loc);
printf("\n");
genkout(estif,NULL,"EKVV",4,0    ,3*iel,0    ,3*iel,0);
genkout(estif,NULL,"EKVP",4,0    ,3*iel,3*iel,4*iel,0);
genkout(estif,NULL,"EKPV",4,3*iel,4*iel,0    ,3*iel,0);
genkout(estif,NULL,"EKPP",4,3*iel,4*iel,3*iel,4*iel,0);
genkout(emass,NULL,"EMVV",4,0    ,3*iel,0    ,3*iel,0);
genkout(emass,NULL,"EMPV",4,3*iel,4*iel,0    ,3*iel,0);
genkout(NULL,etforce,"EFTV",4,0    ,3*iel,0,0,1);
genkout(NULL,etforce,"EFTP",4,3*iel,4*iel,0,0,1);
genkout(NULL,eiforce,"EFIV",4,0    ,3*iel,0,0,1);
genkout(NULL,eiforce,"EFIP",4,3*iel,4*iel,0,0,1);
exit(EXIT_FAILURE);
}
#endif

#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_calint */	

#endif      
