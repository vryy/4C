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
 |                                                           genk 03/02 |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------
\brief integration loop for one fluid2 element

<pre>                                                         genk 04/02

In this routine the element stiffness matrix, iteration-RHS and
time-RHS for one fluid2 element is calculated

NOTE:
Iteration-RHS has been renamed to eforce and etforce has been removed
since the part of the right hand side connected to old time step values
is now served within the mass-rhs procedure after element integration.
                                                             chfoe 02/05

</pre>
\param  *ele	   ELEMENT	   (i)    actual element
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param **emass     DOUBLE	   (o)    element mass matrix
\param  *eforce    DOUBLE	   (o)    element force vector
\param **xyze      DOUBLE          (-)    nodal coordinates
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param **evelng    DOUBLE	   (i)    ele vel. at time n+g
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param **wa1	   DOUBLE	   (-)    working array
\param **wa2	   DOUBLE	   (-)    working array
\param   visc      DOUBLE	   (-)    viscosity
\return void

------------------------------------------------------------------------*/
void f2_calint(
	       ELEMENT         *ele,
               DOUBLE         **estif,
               DOUBLE         **emass,
               DOUBLE          *eforce,
               DOUBLE         **xyze,
               DOUBLE          *funct,
               DOUBLE         **deriv,
               DOUBLE         **deriv2,
               DOUBLE         **xjm,
               DOUBLE         **derxy,
               DOUBLE         **derxy2,
               DOUBLE         **evelng,
               DOUBLE         **vderxy,
               DOUBLE         **wa1,
               DOUBLE         **wa2,
               DOUBLE           visc
	            )
{
INT       iel;        /* number of nodes                                */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration nodesin r,s direction    */
INT       actmat;     /* material number of the element                 */
INT       ihoel=0;    /* flag for higher order elements                 */
INT       icode=2;    /* flag for eveluation of shape functions         */
INT       lr, ls;     /* counter for integration                        */
DOUBLE    dens;       /* density                                        */
DOUBLE    fac;
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    velint[2];
DOUBLE    covint[2];
DIS_TYP   typ;	      /* element type                                   */
STAB_PAR_GLS *gls;    /* pointer to GLS stabilisation parameters        */
FLUID_DATA   *data;

#ifdef DEBUG
dstrc_enter("f2_calint");
#endif

/*----------------------------------------------------- initialisation */
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
typ  = ele->distyp;
gls  = ele->e.f2->stabi.gls;
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;
if (ele->e.f2->stab_type != stab_gls)
   dserror("routine with no or wrong stabilisation called");

/*------- get integraton data and check if elements are "higher order" */
switch (typ)
{
case quad4: case quad8: case quad9:  /* --> quad - element */
   icode   = 3;
   ihoel   = 1;
   /* initialise integration */
   nir = ele->e.f2->nGP[0];
   nis = ele->e.f2->nGP[1];
break;
case tri6: /* --> tri - element */
   icode   = 3;
   ihoel   = 1;
/* do NOT break at this point!!! */
case tri3:
   /* initialise integration */
   nir  = ele->e.f2->nGP[0];
   nis  = 1;
   intc = ele->e.f2->nGP[1];
break;
default:
   dserror("typ unknown!");
} /* end switch(typ) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{
   for (ls=0;ls<nis;ls++)
   {
/*--------------- get values of  shape functions and their derivatives */
      switch(typ)
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
         e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
         e2   = data->qxg[ls][nis-1];
         facs = data->qwgt[ls][nis-1];
         f2_rec(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      case tri3: case tri6:   /* --> tri - element */
         e1   = data->txgr[lr][intc];
         facr = data->twgt[lr][intc];
         e2   = data->txgs[lr][intc];
         facs = ONE;
         f2_tri(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      default:
         dserror("typ unknown!");
      } /* end switch(typ) */
/*--------------------------------------------- compute Jacobian matrix */
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;
/*-------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);
/*------------------------------------ compute second global derivative */
      if (ihoel!=0)
         f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
/*-------------------------- get velocities (n+g,i) at integraton point */
      f2_veci(velint,funct,evelng,iel);
/*--------------- get velocity (n+g,i) derivatives at integration point */
      f2_vder(vderxy,derxy,evelng,iel);
/*--------------- compute stabilisation parameter during ntegration loop*/
      if (gls->iduring!=0)
         f2_calelesize2(ele,xyze,funct,velint,visc,iel,typ);

/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices                         |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/
/*-------------------------------------------------- compute matrix Mvv */
      if (fdyn->nis==0)
         f2_calmvv(emass,funct,fac,iel);
/*-------------------------------------------------- compute matrix Kvv */
      f2_calkvv(ele,estif,velint,NULL,
	        vderxy,funct,derxy,fac,visc,iel);
/*------------------------------------------ compute matrix Kvp and Kpv */
      f2_calkvp(estif,funct,derxy,fac,iel);

/*----------------------------------------------------------------------*
 |         compute Stabilisation matrices                               |
 | NOTE:                                                                |
 |  Stabilisation matrices are all stored in one matrix "estif"         |
 |  Stabilisation mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
/*---------------------------------------- stabilisation for matrix Kvv */
      f2_calstabkvv(ele,gls,estif,velint,velint,NULL,vderxy,
                    funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kvp */
      f2_calstabkvp(ele,gls,estif,velint,
                    funct,derxy,derxy2,fac,visc,iel,ihoel);
      if (fdyn->nis==0)
      {
/*---------------------------------------- stabilisation for matrix Mvv */
         f2_calstabmvv(ele,gls,emass,velint,
	               funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Mpv */
         f2_calstabmpv(gls,emass,funct,derxy,fac,iel);
      }
/*---------------------------------------- stabilisation for matrix Kpv */
      f2_calstabkpv(ele,gls,estif,velint,NULL,vderxy,
                    funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kpp */
      f2_calstabkpp(gls,estif,derxy,fac,iel);

/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 |      (for Newton iteration)                                          |
 *----------------------------------------------------------------------*/
      if (fdyn->nii!=0)
      {
/*-------------- get convective velocities (n+1,i) at integration point */
/*               covint = u*grad(u)                                     */
         f2_covi(vderxy,velint,covint);
/*-------------------- calculate galerkin part of "Iter-RHS" (vel dofs) */
         f2_calgalifv(eforce,covint,funct,fac,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (vel dofs) */
         f2_calstabifv(gls,ele,eforce,covint,velint,funct,
	               derxy,derxy2,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (pre dofs) */
         f2_calstabifp(gls,&(eforce[2*iel]),covint,derxy,fac,iel);
      } /* endif (fdyn->nii!=0) */
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
\param  *ele       ELEMENT          (i)    actual element
\param   imyrank   INT              (i)    proc number
\param **estif     DOUBLE           (o)    element stiffness matrix
\param **emass     DOUBLE           (o)    element mass matrix
\param  *eforce    DOUBLE           (o)    element force vector
\param **xyze      DOUBLE           (-)    nodal coordinates at n+theta
\param  *funct     DOUBLE           (-)    natural shape functions
\param **deriv     DOUBLE           (-)	   deriv. of nat. shape funcs
\param **deriv2    DOUBLE           (-)    2nd deriv. of nat. shape f.
\param **xjm       DOUBLE           (-)    jacobian matrix
\param **derxy     DOUBLE           (-)	   global derivatives
\param **derxy2    DOUBLE           (-)    2nd global derivatives
\param **eveln     DOUBLE           (i)    ele vel. at time n
\param **evelng    DOUBLE           (i)    ele vel. at time n+g
\param **ealecovng DOUBLE           (i)    ele ale-conv. vel. at time n+1
\param **egridv    DOUBLE           (i)    ele grid velocity
\param **vderxy    DOUBLE           (-)    global vel. derivatives
\param  *ekappan   DOUBLE           (i)    nodal curvature at n
\param  *ekappang  DOUBLE           (i)    nodal curvature at n+g
\param **wa1       DOUBLE           (-)    working array
\param **wa2       DOUBLE           (-)    working array
\return void

------------------------------------------------------------------------*/
void f2_calinta(
                  ELEMENT         *ele,
                  INT              imyrank,
                  DOUBLE         **estif,
                  DOUBLE         **emass,
                  DOUBLE          *eforce,
                  DOUBLE         **xyze,
                  DOUBLE          *funct,
                  DOUBLE         **deriv,
                  DOUBLE         **deriv2,
                  DOUBLE         **xjm,
                  DOUBLE         **derxy,
                  DOUBLE         **derxy2,
                  DOUBLE         **eveln,
                  DOUBLE         **evelng,
                  DOUBLE         **ealecovng,
                  DOUBLE         **egridv,
                  DOUBLE         **vderxy,
                  DOUBLE          *ekappan,
                  DOUBLE          *ekappang,
                  DOUBLE          *ephin,
                  DOUBLE          *ephing,
                  DOUBLE         **evnng,
                  DOUBLE         **wa1,
                  DOUBLE         **wa2
               )
{
#ifdef D_FSI
INT       i,k;        /* simply a counter                               */
INT       iel;        /* number of nodes                                */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis,nil;/* number of integration nodesin r,s direction    */
INT       actmat;     /* material number of the element                 */
INT       ihoel=0;    /* flag for higher order elements                 */
INT       icode=2;    /* flag for eveluation of shape functions         */
INT       lr, ls;     /* counter for integration                        */
INT       foundline;
INT       ngline;
INT       node;
DOUBLE    dens;       /* density                                        */
DOUBLE    visc;       /* viscosity                                      */
DOUBLE    fac;
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    velint[2];
DOUBLE    vel2int[2];
DOUBLE    alecovint[2];
DOUBLE    gridvint[2];
DOUBLE    covint[2];
DOUBLE    vnint[2];
DOUBLE    vn[2];      /* component of normal vector                     */
DOUBLE    sigmaint;   /* surface tension                                */
DOUBLE    gamma;      /* surface tension coeficient                     */
DOUBLE    phiintn,phiintng,phiderxn,phiderxng;
DIS_TYP   typ;	      /* element type                                   */
INT       iedgnod[MAXNOD_F2];
INT       line,ngnode;
GLINE    *gline[4];
FLUID_FREESURF_CONDITION *linefs[4];
STAB_PAR_GLS *gls;      /* pointer to GLS stabilisation parameters      */
FLUID_DATA   *data;

#ifdef DEBUG
dstrc_enter("f2_calinta");
#endif

/*----------------------------------------------------- initialisation */
iel=ele->numnp;
actmat=ele->mat-1;
dens = mat[actmat].m.fluid->density;
visc = mat[actmat].m.fluid->viscosity;
gamma= mat[actmat].m.fluid->gamma;
typ  = ele->distyp;

gls    = ele->e.f2->stabi.gls;
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

if (ele->e.f2->stab_type != stab_gls)
   dserror("routine with no or wrong stabilisation called");

/*------- get integraton data and check if elements are "higher order" */
switch (typ)
{
case quad4: case quad8: case quad9:  /* --> quad - element */
   icode   = 3;
   ihoel   = 1;
   /* initialise integration */
   nir = ele->e.f2->nGP[0];
   nis = ele->e.f2->nGP[1];
break;
case tri6: /* --> tri - element */
      icode   = 3;
      ihoel   = 1;
/* do NOT break at this point!!! */
case tri3:    /* initialise integration */
   nir  = ele->e.f2->nGP[0];
   nis  = 1;
   intc = ele->e.f2->nGP[1];
break;
default:
   dserror("typ unknown!");
} /* end switch(typ) */

/*----------------------------------------------------------------------*
 |               start loop over integration points                     |
 *----------------------------------------------------------------------*/
for (lr=0;lr<nir;lr++)
{
   for (ls=0;ls<nis;ls++)
   {
/*--------------- get values of  shape functions and their derivatives */
      switch(typ)
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
         e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
         e2   = data->qxg[ls][nis-1];
         facs = data->qwgt[ls][nis-1];
         f2_rec(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      case tri3: case tri6:   /* --> tri - element */
         e1   = data->txgr[lr][intc];
         facr = data->twgt[lr][intc];
         e2   = data->txgs[lr][intc];
         facs = ONE;
         f2_tri(funct,deriv,deriv2,e1,e2,typ,icode);
      break;
      default:
         dserror("typ unknown!");
      } /* end switch(typ) */
/*--------------------------------------------- compute Jacobian matrix */
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;
/*-------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);
/*------------------------------------ compute second global derivative */
      if (ihoel!=0) f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
/*-------------------------- get velocities (n+g,i) at integraton point */
      f2_veci(velint,funct,evelng,iel);
/*--------- get ale-convectcive velocities (n+g,i) at integration point */
      f2_veci(alecovint,funct,ealecovng,iel);
/*------------------------------ get grid velocity at integration point */
      f2_veci(gridvint,funct,egridv,iel);
/*--------------- get velocity (n+g,i) derivatives at integration point */
      f2_vder(vderxy,derxy,evelng,iel);
/*--------------- compute stabilisation parameter during ntegration loop*/
      if (gls->iduring!=0)
         f2_calelesize2(ele,xyze,funct,velint,visc,iel,typ);

/*----------------------------------------------------------------------*
 |         compute "Standard Galerkin" matrices                         |
 | NOTE:                                                                |
 |  Standard Galerkin matrices are all stored in one matrix "estif"     |
 |  Standard Galerkin mass matrix is stored in "emass"                  |
 *----------------------------------------------------------------------*/
/*-------------------------------------------------- compute matrix Mvv */
      f2_calmvv(emass,funct,fac,iel);
/*-------------------------------------------------- compute matrix Kvv */
      f2_calkvv(ele,estif,velint,gridvint,
	        vderxy,funct,derxy,fac,visc,iel);
/*------------------------------------------ compute matrix Kvp and Kpv */
      f2_calkvp(estif,funct,derxy,fac,iel);
/*-------------------------------------------------- compute matrix Kvg */
      if (ele->e.f2->fs_on==2)
         f2_calkvg(estif,vderxy,funct,fac,iel);
/*-------------------------------------------------- compute matrix Kgg */
      if (ele->e.f2->fs_on==6)
         f2_calgfskgg(estif,funct,deriv,xjm,det,fac,iel);

/*----------------------------------------------------------------------*
 |         compute Stabilisation matrices                               |
 | NOTE:                                                                |
 |  Stabilisation matrices are all stored in one matrix "estif"         |
 |  Stabilisation mass matrices are all stored in one matrix "emass"    |
 *----------------------------------------------------------------------*/
/*---------------------------------------- stabilisation for matrix Kvv */
      f2_calstabkvv(ele,gls,estif,velint,alecovint,gridvint,vderxy,
                    funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kvp */
      f2_calstabkvp(ele,gls,estif,alecovint,
                    funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kvg */
      if (ele->e.f2->fs_on==2)
	 f2_calstabkvg(ele,gls,estif,vderxy,funct,derxy,derxy2,
                       alecovint,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Mvv */
      f2_calstabmvv(ele,gls,emass,alecovint,
                    funct,derxy,derxy2,fac,visc,iel,ihoel);
/*---------------------------------------- stabilisation for matrix Kpv */
      f2_calstabkpv(ele,gls,estif,velint,gridvint,vderxy,
                    funct,derxy,derxy2,fac,visc,iel,ihoel);

/*---------------------------------------- stabilisation for matrix Kpg */
      if (ele->e.f2->fs_on==2)
         f2_calstabkpg(gls,estif,funct,vderxy,derxy,fac,iel);
/*---------------------------------------- stabilisation for matrix Mpv */
      f2_calstabmpv(gls,emass,funct,derxy,fac,iel);
/*---------------------------------------- stabilisation for matrix Kpp */
      f2_calstabkpp(gls,estif,derxy,fac,iel);

/*----------------------------------------------------------------------*
 |         compute "Iteration" Force Vectors                            |
 |      (for Newton iteration)                                          |
 *----------------------------------------------------------------------*/
      if (fdyn->nii!=0)
      {
/*-------------- get convective velocities (n+1,i) at integration point */
/*               covint = c*grad(u)                                     */
         f2_covi(vderxy,alecovint,covint);
/*-------------------- calculate galerkin part of "Iter-RHS" (vel dofs) */
         f2_calgalifv(eforce,covint,funct,fac,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (vel dofs) */
         f2_calstabifv(gls,ele,eforce,covint,alecovint,funct,
                        derxy,derxy2,fac,visc,ihoel,iel);
/*------------------- calculate stabilisation for "Iter-RHS" (pre dofs) */
         f2_calstabifp(gls,&(eforce[2*iel]),covint,derxy,fac,iel);
      } /* endif (fdyn->nii!=0) */
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
if (ele->e.f2->fs_on==2)  /* element at free surface (local lagrange)   */
{
/*-------------------------------- check for presence of freesurface    */
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
     	 f2_edgejaco(xyze,deriv,xjm,&det,ngnode,iedgnod);
     	 fac = det*facr;
     	 /*--------------------------------- compute matrix Kgv and Kgg */
     	 f2_calkgedge(estif,funct,fac,iedgnod,iel,ngnode);
         if (fdyn->surftens>0)
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
     	    if (fdyn->fsstnif!=0)
	    {
            dserror("Surface tension rhs has to be checked by S. Genkinger");
            /* Im Rahmen der Umstellung auf die Verwendung der Massen-
               RHS wird etforce abgeschafft. Die korrekte Integration der
               Oberflaechenspannungseinfluesse muss ueberprueft werden! 
               chfoe 02/05 */
	       /*------ calculate surface tension at gauss point at (n) */
               sigmaint = gamma*f2_edgescali(funct,ekappan,iedgnod,ngnode);
	       /*------------ calculate Time-RHS due to surface tension */
               f2_calsurftenfv(eforce,funct,vn,sigmaint,
                              fdyn->thsl,facr,ngnode,iedgnod);
            } /* endif (fdyn->fsstnif!=0) */
	    /*----------------------------------------------------------*
	     |  surface tension effects at (n+1)                        |
	     *----------------------------------------------------------*/
            if (fdyn->fsstnii!=0)
            {
            dserror("Surface tension rhs has to be checked by S. Genkinger");
            /* Im Rahmen der Umstellung auf die Verwendung der Massen-
               RHS wird etforce abgeschafft. Die korrekte Integration der
               Oberflaechenspannungseinfluesse muss ueberprueft werden! 
               chfoe 02/05 */
               /*---- calculate surface tension at gauss point at (n+1) */
               sigmaint = gamma*f2_edgescali(funct,ekappang,iedgnod,ngnode);
               /*------------ calculate Iter-RHS due to surface tension */
               f2_calsurftenfv(eforce,funct,vn,sigmaint,
                               fdyn->thsr,facr,ngnode,iedgnod);
            } /* endif (fdyn->fsstnii!=0) */
         } /* endif (fdyn->surftens>0) */
      } /* end of loop over integration points */
   } /* end of loop over glines */
} /* endif (ele->e.f2->fs_on==2) */
else if (ele->e.f2->fs_on==5) /* element at free surface (height funct) */
{
/*-------------------------------- check for presence of freesurface    */
   foundline=0;
   /*---------------------------------- number of lines to this element */
   ngline=ele->g.gsurf->ngline;
   /*------- loop over lines, check for freesurface conditions on lines */
   for (i=0; i<ngline; i++)
   {
      gline[i] = ele->g.gsurf->gline[i];
      linefs[i] = gline[i]->freesurf;
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
         /*------- get values of  shape functions and their derivatives */
         e1   = data->qxg[lr][nil-1];
         facr = data->qwgt[lr][nil-1];
         f2_degrectri(funct,deriv,e1,typ,1);
         /*------------------------------- compute jacobian determinant */
         /*------------------ integration is performed along the x-axis */
         xjm[0][0] = ZERO ;
         for (k=0; k<ngnode; k++) /* loop all nodes of the element */
         {
            node=iedgnod[k];
            xjm[0][0] += deriv[0][k] * xyze[0][node] ;
         } /* end loop over iel */
         det = FABS(xjm[0][0]);
         fac = det*facr;

         /*--------------------------------- compute global derivatives */
         f2_edgegder(deriv,derxy,xjm,ngnode);
         /*------------------- get velocity at integration point at n+g */
         f2_edgeveci(velint,funct,evelng,ngnode,iedgnod);
         /*--------------------- get velocity at integration point at n */
         f2_edgeveci(vel2int,funct,eveln,ngnode,iedgnod);
         /*------------ get global derivative of height function at n+g */
         phiderxng = f2_phider(derxy,ephing,ngnode,iedgnod);
         /*-------------- get global derivative of height function at n */
         phiderxn = f2_phider(derxy,ephing,ngnode,iedgnod);
         /*-------------- get height function at integration point at n */
         phiintn   = f2_edgescali(funct,ephin,iedgnod,ngnode);
         /*------------ get height function at integration point at n+g */
         phiintng  = f2_edgescali(funct,ephing,iedgnod,ngnode);

         /*------- compute stabilisation parameter at integration point */
         f2_stabpar_hfsep(ele,xyze,NULL,NULL,velint,phiintn,phiintng,
                          phiderxng,e1,iedgnod,ngnode,typ);
         /*---------- compute Galerkin & stabilisation part of matrices */
         f2_calmat_vhf(emass,estif,funct,derxy,velint,phiderxng,
                           fac,iel,ngnode,iedgnod);
         /*------------------------------------------------ compute RHS */
         f2_calrhs_vhf(eforce,velint,vel2int,phiintn,funct,derxy,
                       phiderxng,phiderxn,fac,iel,ngnode,iedgnod);
      } /* end if loop over integration points */
   } /* end of loop over glines */
} /* endif (ele->e.f2->fs_on==5) */
else if (ele->e.f2->fs_on==1 && fdyn->surftens>0)
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
        /*------------ calculate surface tension at gauss point at (n)*/
         sigmaint = gamma*f2_edgescali(funct,ekappan,iedgnod,ngnode);
         /*----------------- calculate Time-RHS due to surface tension */
         dserror("Surface tension rhs has to be checked by S. Genkinger");
         /* Im Rahmen der Umstellung auf die Verwendung der Massen-
            RHS wird etforce abgeschafft. Die korrekte Integration der
            Oberflaechenspannungseinfluesse muss ueberprueft werden! 
            chfoe 02/05 */
         f2_calsurftenfv(eforce,funct,vn,sigmaint,
                        fdyn->thsl,facr,ngnode,iedgnod);
         /*-------------------------------------------------------------*
          |  surface tension effects at (n+1)                           |
          *-------------------------------------------------------------*/
         /*---------- calculate surface tension at gauss point at (n+1) */
         sigmaint = gamma*f2_edgescali(funct,ekappang,iedgnod,ngnode);
         /*------------------ calculate Time-RHS due to surface tension */
         dserror("Surface tension rhs has to be checked by S. Genkinger");
         /* Im Rahmen der Umstellung auf die Verwendung der Massen-
            RHS wird etforce abgeschafft. Die korrekte Integration der
            Oberflaechenspannungseinfluesse muss ueberprueft werden! 
            chfoe 02/05 */
         f2_calsurftenfv(eforce,funct,vn,sigmaint,
                         fdyn->thsr,facr,ngnode,iedgnod);
      } /* end of loop over integration points */
   } /* end of loop over glines */
}
else if(ele->e.f2->fs_on==6)
{
/*-------------------------------- check for presence of freesurface    */
   foundline=0;
   /*---------------------------------- number of lines to this element */
   ngline=ele->g.gsurf->ngline;
   /*------- loop over lines, check for freesurface conditions on lines */
   for (i=0; i<ngline; i++)
   {
      gline[i] = ele->g.gsurf->gline[i];
      linefs[i] = gline[i]->freesurf;
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
         f2_edgejaco(xyze,deriv,xjm,&det,ngnode,iedgnod);
         fac = det*facr;
         /*--------------------------------- compute normal at intpoint */
         f2_edgeveci(vnint,funct,evnng,ngnode,iedgnod);
         /*------------------------- compute matrix Kgv and Kgg at edge */
         f2_calgfskgedge(estif,funct,vnint,fac,iedgnod,iel,ngnode);
      }
   }
} /* endif (ele->e.f2->fs_on==1 && fdyn->surftens>0) */

/*----------------------------------------------------------------------*/
end:
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
