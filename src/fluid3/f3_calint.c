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
\brief integration loop for one fluid3 element

<pre>                                                         genk 05/02

In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid3 element is calculated
      
</pre>
\param  *data      FLUID_DATA	   (i)	  integration data
\param  *ele	   ELEMENT	   (i)    actual element
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
               INT             *hasext,
               DOUBLE         **estif,
               DOUBLE         **emass,
               DOUBLE          *etforce,
               DOUBLE          *eiforce,
               DOUBLE         **xyze,
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
               DOUBLE         **vderxy,
               DOUBLE          *pderxy,
               DOUBLE         **vderxy2,
               DOUBLE         **wa1,
               DOUBLE         **wa2
               )
{ 
INT      iel;	      /* number of nodes 			        */
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
DOUBLE   velint[3];
DOUBLE   covint[3];
DIS_TYP  typ;         /* element type                                   */

FLUID_DYNAMIC   *fdyn;

STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/

#ifdef DEBUG 
dstrc_enter("f3_calint");
#endif		

/*----------------------------------------------------- initialisation -*/
fdyn   = alldyn[genprob.numff].fdyn;
gls    = ele->e.f3->stabi.gls;
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
typ  = ele->distyp;

/*------------------------------ check for proper stabilisation mode ---*/
dsassert(ele->e.f3->stab_type == stab_gls,
         "routine with no or wrong stabilisation called");
   
/*------- get integraton data and check if elements are "higher order" -*/
intc = 0;
nir  = 0;
nis  = 0;
nit  = 0; 
switch (typ)
{
case hex8: case hex20: case hex27:  /* --> hex - element */
   icode   = 3;
   ihoel   = 1;
   /* initialise integration */
   nir = ele->e.f3->nGP[0];
   nis = ele->e.f3->nGP[1];
   nit = ele->e.f3->nGP[2];
   intc= 0;
   break;
case tet10: /* --> tet - element */  
   icode   = 3;
   ihoel   = 1;
/* do NOT break at this point!!! */
case tet4:    /* initialise integration */
   nir  = ele->e.f3->nGP[0];
   nis  = 1;
   nit  = 1; 
   intc = ele->e.f3->nGP[1];  
   break;
default:
   dserror("typ unknown!");
} /* end switch (typ) */

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
   switch(typ)  
   {
   case hex8: case hex20: case hex27:   /* --> hex - element */
      e1   = data->qxg[lr][nir-1];
      facr = data->qwgt[lr][nir-1];
      e2   = data->qxg[ls][nis-1];
      facs = data->qwgt[ls][nis-1];
      e3   = data->qxg[lt][nit-1];
      fact = data->qwgt[lt][nit-1];
      f3_hex(funct,deriv,deriv2,e1,e2,e3,typ,icode);
      break;
   case tet4: case tet10:   /* --> tet - element */		  	
      e1   = data->txgr[lr][intc];
      facr = data->twgt[lr][intc];
      e2   = data->txgs[lr][intc];
      facs = ONE;
      e3   = data->txgt[lr][intc]; 
      fact = ONE;
      f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,icode); 
      break;
   default:
      facr = facs = fact = 0.0;
      e1 = e2 = e3 = 0.0;
      dserror("typ unknown!");
   } /* end switch (typ) */
/*-------------------------------------------- compute Jacobian matrix */  
   f3_jaco(xyze,funct,deriv,xjm,&det,ele,iel);
   fac = facr*facs*fact*det;
/*------------------------------------------- compute global derivates */
   f3_gder(derxy,deriv,xjm,wa1,det,iel);
/*------------------------------------ compute second global derivative */ 
   if (ihoel!=0)
      f3_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
/*------------------------- get velocities (n+g,i) at integraton point */
   f3_veci(velint,funct,evelng,iel);
/*-------------- get velocity (n+g,i) derivatives at integration point */
   f3_vder(vderxy,derxy,evelng,iel);
 /*-------------- compute stabilisation parameter during ntegration loop*/
   if (gls->iduring!=0)
      f3_calelesize2(ele,velint,wa1,visc,iel,typ);

/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices                         |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/
/*-------------------------------------------------- compute matrix Mvv */
   if (fdyn->nis==0)	  	 	    
      f3_calmvv(emass,funct,fac,iel);
/*-------------------------------------------------- compute matrix Kvv */
   f3_calkvv(ele,estif,velint,NULL,vderxy,funct,
             derxy,fac,visc,iel);
/*------------------------------------------ compute matrix Kvp and Kpv */
   f3_calkvp(estif,funct,derxy,fac,iel);

/*----------------------------------------------------------------------*
 |         compute Stabilisation matrices                               |
 | NOTE:                                                                |
 |  Stabilisation matrices are all stored in one matrix "estif"         |
 |  Stabilisation mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
   
/*---------------------------------------- stabilisation for matrix Kvv */
   f3_calstabkvv(ele,gls,estif,velint,velint,NULL,vderxy,
                          funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kvp */
   f3_calstabkvp(ele,gls,estif,velint,
                 funct,derxy,derxy2,fac,visc,iel,ihoel); 
   if (fdyn->nis==0)
   {
/*---------------------------------------- stabilisation for matrix Mvv */
      f3_calstabmvv(ele,gls,emass,velint,
                   funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Mpv */
      f3_calstabmpv(gls,emass,funct,derxy,fac,iel);   
   }
/*---------------------------------------- stabilisation for matrix Kpv */ 
   f3_calstabkpv(ele,gls,estif,velint,NULL,vderxy,
                  funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kpp */ 
   f3_calstabkpp(gls,estif,derxy,fac,iel);  

/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 |      (for Newton iteration and for fixed-point iteration)            |
 *----------------------------------------------------------------------*/
   if (fdyn->nii!=0)
   {
/*-------------- get convective velocities (n+1,i) at integration point */
      f3_covi(vderxy,velint,covint);
/*-------------------- calculate galerkin part of "Iter-RHS" (vel dofs) */
      f3_calgalifv(eiforce,covint,funct,fac,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (vel dofs) */
      f3_calstabifv(gls,ele,eiforce,covint,velint,funct,
                     derxy,derxy2,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (pre dofs) */
      f3_calstabifp(gls,&(eiforce[3*iel]),covint,derxy,fac,iel);
   } /* endif (fdyn->nii!=0) */

/*----------------------------------------------------------------------*
 |       compute "external" Force Vector                                |
 |   (at the moment there are no external forces implemented)           |
 |  but there can be due to self-weight /magnetism / etc. (b)           |
 |  dead load may vary over time, but stays constant over               |
 |  the whole domain --> no interpolation with shape funcs              |
 |  parts changing during the nonlinear iteration are added to          |
 |  Iteration Force Vector                                              |
 *----------------------------------------------------------------------*/
   if (*hasext!=0)
   {
/*------- compute stabilisation part of external RHS (vel dofs) at (n+1)*/
      f3_calstabexfv(gls,ele,eiforce,derxy,derxy2,edeadng,
	             velint,fac,visc,iel,ihoel,1); 
/*------ compute stabilisation part of external RHS (pre dofs) at (n+1) */
      f3_calstabexfp(gls,&(eiforce[3*iel]),derxy,edeadng,fac,iel,1); 
   } /* endif (*hasext!=0) */

/*----------------------------------------------------------------------*
 |         compute "Time" Force Vectors                                 |
 *----------------------------------------------------------------------*/
   if (fdyn->nif!=0)
   {
/*------------------------------- get pressure (n) at integration point */
      preint=f3_scali(funct,epren,iel);
/*------------------- get pressure derivatives (n) at integration point */
      f3_pder(pderxy,derxy,epren,iel);	    
/*----------------------------- get velocities (n) at integration point */
      f3_veci(velint,funct,eveln,iel);
/*------------------- get velocitiederivatives (n) at integration point */
      f3_vder(vderxy,derxy,eveln,iel);
/*------------- get 2nd velocities derivatives (n) at integration point */
      if (ihoel!=0)
         f3_vder2(vderxy2,derxy2,eveln,iel);	       
/*------------------ get convective velocities (n) at integration point */
      f3_covi(vderxy,velint,covint);        	    
/*--------------------- calculate galerkin part of "Time-RHS" (vel-dofs)*/
      f3_calgaltfv(etforce,velint,covint,
                  funct,derxy,vderxy,preint,visc,fac,iel);
/*-------------------- calculate galerkin part of "Time-RHS" (pre-dofs) */
      f3_calgaltfp(&(etforce[3*iel]),funct,vderxy,fac,iel);
/*------------------- calculate stabilisation for "Time-RHS" (vel-dofs) */
      f3_calstabtfv(gls,ele,etforce,velint,velint,
                     covint,derxy,derxy2,vderxy,vderxy2,
                     pderxy,fac,visc,ihoel,iel);

/*------------------- calculate stabilisation for "Time-RHS" (pre-dofs) */
      f3_calstabtfp(ele,gls,&(etforce[3*iel]),derxy,vderxy2,
                   velint,covint,pderxy,visc,fac,ihoel,iel);

/*----------------------------------------------------------------------*
            | compute "external" Force Vector                           |
            |  but there can be due to self-weight /magnetism / etc. (b)|
            |  dead load may vary over time, but stays constant over    |
	    |  the whole domain --> no interpolation with shape funcs   |
            |  parts staying constant during nonlinear iteration are    |
            |  add to Time Force Vector                                 |
 *----------------------------------------------------------------------*/
      if (*hasext!=0)
      {
/*--- compute galerkin part of external RHS (vel dofs) at (n) and (n+1) */
         f3_calgalexfv(etforce,funct,edeadn,edeadng,fac,iel);

/*-------- compute stabilisation part of external RHS (vel dofs) at (n) */
         f3_calstabexfv(gls,ele,etforce,derxy,derxy2,edeadn,
                        velint,fac,visc,iel,ihoel,0);
/*--------------- compute stabilistaion part of external RHS (pre dofs) */
         f3_calstabexfp(gls,&(etforce[3*iel]),derxy,edeadn,fac,iel,0); 
      } /* endif (*hasext!=0) */
   } /* endif (fdyn->nif!=0)   */
}
}
} /* end of loop over integration points */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_calint */	

/*!---------------------------------------------------------------------
\brief integration loop for one fluid3 element

<pre>                                                         genk 05/02

In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid3 element is calculated
      
</pre>
\param  *data      FLUID_DATA       (i) integration data
\param  *ele       ELEMENT          (i) actual element
\param  *hasext    INT              (i) element flag
\param **estif     DOUBLE           (o) element stiffness matrix
\param **emass     DOUBLE           (o) element mass matrix
\param  *etforce   DOUBLE           (o) element time force vector
\param  *eiforce   DOUBLE           (o) element iter force vector
\param **xyze      DOUBLE           (-) nodal coordinates
\param  *funct     DOUBLE           (-) natural shape functions
\param **deriv     DOUBLE           (-)	deriv. of nat. shape funcs
\param **deriv2    DOUBLE           (-) 2nd deriv. of nat. shape f.
\param **xjm       DOUBLE           (-) jacobian matrix
\param **derxy     DOUBLE           (-)	global derivatives
\param **derxy2    DOUBLE           (-) 2nd global derivatives
\param **eveln     DOUBLE           (i) ele vel. at time n
\param **evelng    DOUBLE           (i) ele vel. at time n+g
\param  *epren     DOUBLE           (-) ele pres. at time n
\param  *edeadn    DOUBLE           (-) ele dead load (selfweight) at n 
\param  *edeadng   DOUBLE           (-) ele dead load (selfweight) at n+1
\param **vderxy    DOUBLE           (-) global vel. derivatives
\param  *pderxy    DOUBLE           (-) global pres. derivatives
\param **vderxy2   DOUBLE           (-) 2nd global vel. deriv.
\param **wa1       DOUBLE           (-) working array
\param **wa2       DOUBLE           (-) working array
                                  
                                  
\return void                                                   

------------------------------------------------------------------------*/
void f3_calinta(     	
                  FLUID_DATA      *data, 
                  ELEMENT         *ele,
                  INT             *hasext,
                  DOUBLE         **estif,
                  DOUBLE         **emass,
                  DOUBLE          *etforce,
                  DOUBLE          *eiforce,
                  DOUBLE         **xyze,
                  DOUBLE          *funct,
                  DOUBLE         **deriv,
                  DOUBLE         **deriv2,
                  DOUBLE         **xjm,
                  DOUBLE         **derxy,
                  DOUBLE         **derxy2,
                  DOUBLE         **eveln,
                  DOUBLE         **evelng,
                  DOUBLE         **ealecovn,
                  DOUBLE         **ealecovng,
                  DOUBLE         **egridv, 
                  DOUBLE          *epren,
                  DOUBLE          *edeadn,
                  DOUBLE          *edeadng,
                  DOUBLE         **vderxy,
                  DOUBLE          *pderxy,
                  DOUBLE         **vderxy2,
                  DOUBLE         **wa1,
                  DOUBLE         **wa2
	      )
{ 
#ifdef D_FSI
INT      i;
INT      iel;	        /* number of nodes 			        */
INT      intc;        /* "integration case" for tet for further infos
                          see f3_inpele.c and f3_intg.c                 */
INT      nir,nis,nit; /* number of integration nodes in r,s,t direction */
INT      nil;
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
DOUBLE   velint[3];
DOUBLE   alecovint[3];
DOUBLE   gridvint[3];
DOUBLE   covint[3];
DIS_TYP  typ;         /* element type                                   */
INT       iedgnod[MAXNOD_F3];
INT       surf,ngsurf,ngnode,foundsurf;
GSURF    *gsurf[6];
FLUID_FREESURF_CONDITION *surffs[4];
STAB_PAR_GLS *gls;	/* pointer to GLS stabilisation parameters	*/
DOUBLE  akov[3][3];
DOUBLE  h[3],da;
FLUID_DYNAMIC *fdyn;
#ifdef DEBUG 
dstrc_enter("f3_calinta");
#endif

/*----------------------------------------------------- initialisation */
fdyn   = alldyn[genprob.numff].fdyn;
gls    = ele->e.f3->stabi.gls;
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
typ  = ele->distyp;

/*------- get integraton data and check if elements are "higher order" */
switch (typ)
{
case hex8: case hex20: case hex27:  /* --> hex - element */
   icode   = 3;
   ihoel   = 1;
   /* initialise integration */
   nir = ele->e.f3->nGP[0];
   nis = ele->e.f3->nGP[1];
   nit = ele->e.f3->nGP[2];
   break;
case tet10: /* --> tet - element */  
   icode   = 3;
   ihoel   = 1;
/* do NOT break at this point!!! */
case tet4:    /* initialise integration */
   nir  = ele->e.f3->nGP[0];
   nis  = 1;
   nit  = 1; 
   intc = ele->e.f3->nGP[1];  
   break;
default:
   dserror("typ unknown!");
} /* end switch (typ) */

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
   switch(typ)  
   {
   case hex8: case hex20: case hex27:   /* --> hex - element */
      e1   = data->qxg[lr][nir-1];
      facr = data->qwgt[lr][nir-1];
      e2   = data->qxg[ls][nis-1];
      facs = data->qwgt[ls][nis-1];
      e3   = data->qxg[lt][nit-1];
      fact = data->qwgt[lt][nit-1];
      f3_hex(funct,deriv,deriv2,e1,e2,e3,typ,icode);
      break;
   case tet4: case tet10:   /* --> tet - element */		  	
      e1   = data->txgr[lr][intc];
      facr = data->twgt[lr][intc];
      e2   = data->txgs[lr][intc];
      facs = ONE;
      e3   = data->txgt[lr][intc]; 
      fact = ONE;
      f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,icode); 
      break;
   default:
      dserror("typ unknown!");
   } /* end switch (typ) */
/*-------------------------------------------- compute Jacobian matrix */  
   f3_jaco(xyze,funct,deriv,xjm,&det,ele,iel);
   fac = facr*facs*fact*det;
/*------------------------------------------- compute global derivates */
   f3_gder(derxy,deriv,xjm,wa1,det,iel);
/*------------------------------------ compute second global derivative */ 
   if (ihoel!=0)
      f3_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
/*------------------------- get velocities (n+g,i) at integraton point */
   f3_veci(velint,funct,evelng,iel);
/*-------- get ale-convectcive velocities (n+g,i) at integration point */
   f3_veci(alecovint,funct,ealecovng,iel);
/*----------------------------- get grid velocity at integration point */
   f3_veci(gridvint,funct,egridv,iel);
/*-------------- get velocity (n+g,i) derivatives at integration point */
   f3_vder(vderxy,derxy,evelng,iel);
 /*-------------- compute stabilisation parameter during integration loop*/
   if (gls->iduring!=0)
      f3_calelesize2(ele,velint,wa1,visc,iel,typ);

/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices                         |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/
/*-------------------------------------------------- compute matrix Kvv */
      f3_calkvv(ele,estif,velint,gridvint,vderxy,
                funct,derxy,fac,visc,iel);
/*------------------------------------------ compute matrix Kvp and Kpv */
      f3_calkvp(estif,funct,derxy,fac,iel);
/*-------------------------------------------------- compute matrix Mvv */
      f3_calmvv(emass,funct,fac,iel);
/*----------------------------------------------------------------------*
 |         compute Stabilisation matrices                               |
 | NOTE:                                                                |
 |  Stabilisation matrices are all stored in one matrix "estif"         |
 |  Stabilisation mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
   
/*---------------------------------------- stabilisation for matrix Kvv */
      f3_calstabkvv(ele,gls,estif,velint,alecovint,gridvint,vderxy,
                     funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kvp */
      f3_calstabkvp(ele,gls,estif,alecovint,
                     funct,derxy,derxy2,fac,visc,iel,ihoel); 
/*---------------------------------------- stabilisation for matrix Kvg*/ 
      if (ele->e.f3->fs_on==2)
         f3_calstabkvg(ele,gls,estif,vderxy,funct,derxy,derxy2,
                        alecovint,fac,visc,iel,ihoel);

/*---------------------------------------- stabilisation for matrix Mvv */
      f3_calstabmvv(ele,gls,emass,alecovint,
                   funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kpv */ 
      f3_calstabkpv(ele,gls,estif,velint,gridvint,vderxy,
                   funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kpg */
      if (ele->e.f3->fs_on==2)
         f3_calstabkpg(gls,estif,funct,vderxy,derxy,fac,iel);

/*---------------------------------------- stabilisation for matrix Mpv */
      f3_calstabmpv(gls,emass,funct,derxy,fac,iel);
/*---------------------------------------- stabilisation for matrix Kpp */ 
      f3_calstabkpp(gls,estif,derxy,fac,iel);  

/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 |      (for Newton iteration and for fixed-point iteration)            |
 *----------------------------------------------------------------------*/
   if (fdyn->nii!=0)
   {
/*-------------- get convective velocities (n+1,i) at integration point */
/*               covint = c*grad(u)                                     */                  
      f3_covi(vderxy,alecovint,covint);
/*-------------------- calculate galerkin part of "Iter-RHS" (vel dofs) */
      f3_calgalifv(eiforce,covint,funct,fac,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (vel dofs) */
      f3_calstabifv(gls,ele,eiforce,covint,alecovint,funct,
                   derxy,derxy2,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (pre dofs) */
      f3_calstabifp(gls,&(eiforce[3*iel]),covint,derxy,fac,iel);
   } /* endif (fdyn->nii!=0) */

/*----------------------------------------------------------------------*
 |       compute "external" Force Vector                                |
 |   (at the moment there are no external forces implemented)           |
 |  but there can be due to self-weight /magnetism / etc. (b)           |
 |  dead load may vary over time, but stays constant over               |
 |  the whole domain --> no interpolation with shape funcs              |
 |  parts changing during the nonlinear iteration are added to          |
 |  Iteration Force Vector                                              |
 *----------------------------------------------------------------------*/
   if (*hasext!=0)
   {
/*------- compute stabilisation part of external RHS (vel dofs) at (n+1)*/
      f3_calstabexfv(gls,ele,eiforce,derxy,derxy2,edeadng,
                     alecovint,fac,visc,iel,ihoel,1); 
/*------ compute stabilisation part of external RHS (pre dofs) at (n+1) */
      f3_calstabexfp(gls,&(eiforce[3*iel]),derxy,edeadng,fac,iel,1); 
   } /* endif (*hasext!=0) */

/*----------------------------------------------------------------------*
 |         compute "Time" Force Vectors                                 |
 *----------------------------------------------------------------------*/
   if (fdyn->nif!=0)
   {
/*------------------------------- get pressure (n) at integration point */
      preint=f3_scali(funct,epren,iel);
/*------------------- get pressure derivatives (n) at integration point */
      f3_pder(pderxy,derxy,epren,iel);	    
/*----------------------------- get velocities (n) at integration point */
      f3_veci(velint,funct,eveln,iel);
/*-------------- get ale-convective velocities (n) at integration point */
      f3_veci(alecovint,funct,ealecovn,iel);
/*------------------- get velocitiederivatives (n) at integration point */
      f3_vder(vderxy,derxy,eveln,iel);
/*------------- get 2nd velocities derivatives (n) at integration point */
      if (ihoel!=0)
         f3_vder2(vderxy2,derxy2,eveln,iel);	       
/*------------------ get convective velocities (n) at integration point 
                     covint = c * grad(u)                               */
      f3_covi(vderxy,alecovint,covint);        	    
/*--------------------- calculate galerkin part of "Time-RHS" (vel-dofs)*/
      f3_calgaltfv(etforce,velint,covint,
                  funct,derxy,vderxy,preint,visc,fac,iel);
/*-------------------- calculate galerkin part of "Time-RHS" (pre-dofs) */
      f3_calgaltfp(&(etforce[3*iel]),funct,vderxy,fac,iel);
/*------------------- calculate stabilisation for "Time-RHS" (vel-dofs) */
      f3_calstabtfv(gls,ele,etforce,alecovint,velint,
                     covint,derxy,derxy2,vderxy,vderxy2,
                     pderxy,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Time-RHS" (pre-dofs) */
      f3_calstabtfp(ele,gls,&(etforce[3*iel]),derxy,vderxy2,
                     velint,covint,pderxy,visc,fac,ihoel,iel);
/*----------------------------------------------------------------------*
            | compute "external" Force Vector                           |
            | (at the moment there are no external forces implemented)  |
            |  but there can be due to self-weight /magnetism / etc. (b)|
            |  dead load may vary over time, but stays constant over    |
            |  the whole domain --> no interpolation with shape funcs   |
            |  parts staying constant during nonlinear iteration are    |
            |  add to Time Force Vector                                 |
/-----------------------------------------------------------------------*/
      if (*hasext!=0)
      {
/*--- compute galerkin part of external RHS (vel dofs) at (n) and (n+1) */
         f3_calgalexfv(etforce,funct,edeadn,edeadng,fac,iel);
/*-------- compute stabilisation part of external RHS (vel dofs) at (n) */
         f3_calstabexfv(gls,ele,etforce,derxy,derxy2,edeadn,
                        alecovint,fac,visc,iel,ihoel,0);
/*--------------- compute stabilistaion part of external RHS (pre dofs) */
         f3_calstabexfp(gls,&(etforce[3*iel]),derxy,edeadn,fac,iel,0); 
      } /* endif (*hasext!=0) */
   } /* endif (fdyn->nif!=0)   */
}
}
} /* end of loop over integration points */

/*----------------------------------------------------------------------*
 |             evaluate the integrals over the element edges            |
 *----------------------------------------------------------------------*/

if (ele->e.f3->fs_on==2)  /* element at free surface (local lagrange)   */
{                                        
/*-------------------------------- check for presence of freesurface    */
   foundsurf=0;
   /*------------------------------- number of surfaces to this element */
   ngsurf=ele->g.gvol->ngsurf;
   /*------- loop over lines, check for freesurface conditions on lines */
   for (i=0; i<ngsurf; i++)
   {
      gsurf[i] = ele->g.gvol->gsurf[i];
      surffs[i] = gsurf[i]->freesurf;
      if(surffs[i]==NULL) continue;
      foundsurf++;
   }
   if (foundsurf==0) goto end;  
   /*--------------------------------------- set number of gauss points */
   nil = IMAX(nir,2);
   /*--------------------------------------------------- distyp of edge */
   switch (typ)
   {
   case hex8: typ=quad4; break;
   default: dserror("distyp not allowed for implicit free surface!\n");
   }
   /*------------------------------- loop over surfaces at free surface */
   for (surf=0; surf<ngsurf; surf++)
   {
      if (surffs[surf]==NULL) continue;
      /*--------------------------------- check number of nodes on surf */
      ngnode = gsurf[surf]->ngnode;
      /*------------------------------------------------ get edge nodes */
      f3_iedg(iedgnod,ele,surf);
      /*------------------------------ integration loop on actual gsurf */
      for (lr=0;lr<nil;lr++)
      {
      for (ls=0;ls<nil;ls++)
      {
      /*--------------- get values of  shape functions and their derivatives */
      switch(typ)  
      {
      case hex8: case hex20: case hex27:   /* --> hex - element */
         e1   = data->qxg[lr][nil-1];
         facr = data->qwgt[lr][nil-1];
         e2   = data->qxg[ls][nil-1];
         facs = data->qwgt[ls][nil-1];
         f3_rec(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      case tet4: case tet10:   /* --> tet - element */		  	
         dserror("implicit tri free surface not implemented yet!\n");              
         e1   = data->txgr[lr][intc];
         facr = data->twgt[lr][intc];
         e2   = data->txgs[lr][intc];
         facs = ONE;
         f3_tri(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      default: 
         dserror("typ unknown!");
      } /* end switch(typ) */    
      /*------------------------------------- metrics at gaussian point */
      f3_tvmr(xyze,akov,funct,deriv,iedgnod,ngnode);  
      /*------------------------- make h as cross product in ref config.
                                       to get area da on surface edge   */
      h[0] = akov[1][0]*akov[2][1] - akov[2][0]*akov[1][1];
      h[1] = akov[2][0]*akov[0][1] - akov[0][0]*akov[2][1];
      h[2] = akov[0][0]*akov[1][1] - akov[1][0]*akov[0][1];
      /*------------------------------------- make director unit lenght 
                                        and get midsurf area da from it */
      math_unvc(&da,h,3);

      fac = da*facr*facs;
      /*--------------------------------- compute matrix Kgv and Kgg */
      f3_calkgedge(estif,funct,fac,iedgnod,iel,ngnode);      
      } /* end of loop over integration points */
      } /* end of loop over integration points */
   } /* end of loop over glines */	       
} /* endif (ele->e.f3->fs_on==2) */

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
#else
dserror("FSI-functions not compiled in!\n");
#endif

return; 
} /* end of f3_calint */	


#endif      
