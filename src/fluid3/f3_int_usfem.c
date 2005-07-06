/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element using USFEM

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
/*! 
\addtogroup FLUID2 
*//*! @{ (documentation module open)*/
#ifdef D_FLUID3 
#include "../headers/standardtypes.h"
#include "../ale3/ale3.h"
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

/*!---------------------------------------------------------------------
\brief integration loop for one fluid2 element using USFEM

<pre>                                                         chfoe 04/04

In this routine the element 'stiffness' matrix and RHS for one 
fluid2 element is calculated
      
</pre>
\param  *ele	   ELEMENT	   (i)    actual element
\param  *hasext    INT             (i)    element flag
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param **emass     DOUBLE	   (o)    element mass matrix
\param  *etforce   DOUBLE	   (o)    element time force vector
\param  *eiforce   DOUBLE	   (o)    element iter force vector
\param **xyze      DOUBLE          (-)    nodal coordinates
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param **eveln     DOUBLE	   (i)    ele vel. at time n
\param **evelng    DOUBLE	   (i)    ele vel. at time n+g
\param **evhist    DOUBLE	   (i)    lin. combination of recent vel and acc
\param **egridv    DOUBLE	   (i)    grid velocity of element
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
void f3_int_usfem(
	              ELEMENT         *ele,
                      INT             *hasext,
                      DOUBLE         **estif,
	              DOUBLE          *force,
	              DOUBLE         **xyze,
	              DOUBLE          *funct,
	              DOUBLE         **deriv,
	              DOUBLE         **deriv2,
	              DOUBLE         **xjm,
	              DOUBLE         **derxy,
	              DOUBLE         **derxy2,
	              DOUBLE         **evelng,
	              DOUBLE         **evhist,
	              DOUBLE         **egridv,
	              DOUBLE          *epren,
	              DOUBLE          *edeadng,
	              DOUBLE         **vderxy,
                      DOUBLE         **vderxy2,
                      DOUBLE           visc,
	              DOUBLE         **wa1,
	              DOUBLE         **wa2
	             )
{ 
INT       i;          /* a couter                                       */
INT       iel;        /* number of nodes                                */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis,nit;/* number of integration nodesin r,s,t direction  */
INT       ihoel=0;    /* flag for higher order elements                 */
INT       icode=2;    /* flag for eveluation of shape functions         */     
INT       lr, ls, lt; /* counter for integration                        */
INT       is_ale;
DOUBLE    fac;        /* total integration vactor                       */
DOUBLE    facr, facs, fact; /* integration weights                      */
DOUBLE    det;        /* determinant of jacobian matrix at time (n+1)   */
DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
DOUBLE    gradp[3];   /* pressure gradient at integration point         */
DOUBLE    velint[3];  /* velocity vector at integration point           */
DOUBLE    histvec[3]; /* history data at integration point              */
DOUBLE    gridvelint[3]; /* grid velocity                               */
DIS_TYP   typ;	      /* element type                                   */

FLUID_DYNAMIC   *fdyn;
FLUID_DATA      *data;

#ifdef DEBUG 
dstrc_enter("f3_int_usfem");
#endif

/*--------------------------------------------------- initialisation ---*/
iel    = ele->numnp;
typ    = ele->distyp;
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

is_ale = ele->e.f3->is_ale;
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
/*------------- get values of  shape functions and their derivatives ---*/
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
   
   /*----------------------------------------- compute Jacobian matrix */
   f3_jaco(xyze,deriv,xjm,&det,ele,iel);
   fac = facr*facs*fact*det;
   
   /*---------------------------------------- compute global derivates */
   f3_gder(derxy,deriv,xjm,wa1,det,iel);

   /*--------------------------------- compute second global derivative */
   if (ihoel!=0)
   {
      f3_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
      f3_vder2(vderxy2,derxy2,evelng,iel);
   }

   /*---------------------- get velocities (n+g,i) at integraton point */
   f3_veci(velint,funct,evelng,iel);

   /*---------------- get history data (n,i) at integration point ---*/
   f3_veci(histvec,funct,evhist,iel);

   /*----------- get velocity (n+g,i) derivatives at integration point */
   f3_vder(vderxy,derxy,evelng,iel);

   /*--------------------- get grid velocity at integration point ---*/
   if(is_ale) f3_veci(gridvelint,funct,egridv,iel);
   
   /*------------------------------------- get pressure gradients ---*/
   gradp[0] = gradp[1] = gradp[2] = 0.0;
   
   for (i=0; i<iel; i++)
   {
      gradp[0] += derxy[0][i] * epren[i];
      gradp[1] += derxy[1][i] * epren[i];
      gradp[2] += derxy[2][i] * epren[i];
   }

   /*-------------- perform integration for entire matrix and rhs ---*/
   f3_calmat(estif,force,velint,histvec,gridvelint,vderxy,
             vderxy2,gradp,funct,derxy,derxy2,edeadng,fac,
                visc,iel,hasext,is_ale);
} /* end of loop over integration points lt*/
} /* end of loop over integration points ls */
} /* end of loop over integration points lr */
 
/*----------------------------------------------- to ensure assembly ---*/
*hasext = 1;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_int_usfem_ale */



/*!---------------------------------------------------------------------
\brief integration loop for one fluid2 element residual vector 
       basing on USFEM

<pre>                                                         chfoe 11/04

This routine calculates the elemental residual vector for one converged
fluid element. The field velint contains the converged Gauss point value
of the velocity. The elemental residual vector is used to obtain fluid
lift and drag forces and FSI coupling forces.
      
</pre>
\param  *ele	   ELEMENT	   (i)    actual element
\param  *hasext    INT             (i)    element flag
\param  *force     DOUBLE	   (o)    elemental force vector to be filled
\param **xyze      DOUBLE          (-)    nodal coordinates
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param **eveln     DOUBLE	   (i)    ele vel. at time n
\param **evelng    DOUBLE	   (i)    ele vel. at time n+g
\param **evhist    DOUBLE	   (i)    lin. combination of recent vel and acc
\param **ealecovng DOUBLE	   (i)    ALE convective velocity
\param  *epren     DOUBLE	   (-)    ele pres. at time n
\param  *edeadng   DOUBLE	   (-)    ele dead load (selfweight) at n+1
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param  *pderxy    DOUBLE	   (-)    global pres. derivatives
\param **vderxy2   DOUBLE	   (-)    2nd global vel. deriv.
\param   visc      DOUBLE          (i)    viscosity
\param **wa1	   DOUBLE	   (-)    working array
\param **wa2	   DOUBLE	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f3_int_res(
	        ELEMENT         *ele,
                INT             *hasext,
	        DOUBLE          *force,
	        DOUBLE         **xyze,
	        DOUBLE          *funct,
	        DOUBLE         **deriv,
	        DOUBLE         **deriv2,
	        DOUBLE         **xjm,
	        DOUBLE         **derxy,
	        DOUBLE         **derxy2,
	        DOUBLE         **evelng,
	        DOUBLE         **evhist,
	        DOUBLE         **ealecovng,
	        DOUBLE          *epren,
	        DOUBLE          *edeadng,
	        DOUBLE         **vderxy,
                DOUBLE         **vderxy2,
                DOUBLE           visc,
	        DOUBLE         **wa1,
	        DOUBLE         **wa2
	       )
{ 
INT       i;          /* a couter                                       */
INT       iel;        /* number of nodes                                */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       is_ale;     /* ALE or Euler element flag                      */
INT       nir,nis,nit;/* number of integration nodesin r,s,t direction  */
INT       ihoel=0;    /* flag for higher order elements                 */
INT       icode=2;    /* flag for eveluation of shape functions         */     
INT       lr, ls, lt; /* counter for integration                        */

DOUBLE    fac;        /* total integration vactor                       */
DOUBLE    facr, facs, fact; /* integration weights                      */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
DOUBLE    velint[3];  /* velocity vector at integration point           */
DOUBLE    histvec[3]; /* history data at integration point              */
DOUBLE    aleconv[3]; /* ALE convective velocity at Gauss point         */
DIS_TYP   typ;	      /* element type                                   */
DOUBLE    presint;    /* pressure at integration point                  */
DOUBLE    gradp[3];   /* pressure gradient                              */
FLUID_DYNAMIC   *fdyn;
FLUID_DATA      *data;


#ifdef DEBUG 
dstrc_enter("f3_int_res");
#endif

/*--------------------------------------------------- initialisation ---*/
iel    = ele->numnp;
typ    = ele->distyp;
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

is_ale = ele->e.f3->is_ale;

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
/*------------- get values of  shape functions and their derivatives ---*/
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
      
   /*----------------------------------------- compute Jacobian matrix */
   f3_jaco(xyze,deriv,xjm,&det,ele,iel);
   fac = facr*facs*fact*det;
   
   /*---------------------------------------- compute global derivates */
   f3_gder(derxy,deriv,xjm,wa1,det,iel);

   /*--------------------------------- compute second global derivative */
   if (ihoel!=0)
   {
      f3_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
      f3_vder2(vderxy2,derxy2,evelng,iel);
   }

   /*-------------------- get velocities (n+g,i) at integraton point ---*/
   f3_veci(velint,funct,evelng,iel);

   /*------------------- get history data (n,i) at integration point ---*/
   f3_veci(histvec,funct,evhist,iel);

   /*--------- get velocity (n+g,i) derivatives at integration point ---*/
   f3_vder(vderxy,derxy,evelng,iel);

   /*-------------- get ALE convective velocity at integration point ---*/
   if(is_ale) f3_veci(aleconv,funct,ealecovng,iel);
   
   /*---------------------------------------- get pressure gradients ---*/
   gradp[0] = gradp[1] = gradp[2] = 0.0;
   
   for (i=0; i<iel; i++)
   {
      gradp[0] += derxy[0][i] * epren[i];
      gradp[1] += derxy[1][i] * epren[i];
      gradp[2] += derxy[2][i] * epren[i];
   }

   /*----------------------------- get pressure at integration point ---*/
   presint = f3_scali(funct,epren,iel);
      
   /*----------------- perform integration for entire matrix and rhs ---*/
   f3_calresvec(force,velint,histvec,vderxy,vderxy2,funct,derxy,derxy2,
                edeadng,aleconv,&presint,gradp,fac,visc,iel,hasext,
                 is_ale);
} /* end of loop over integration points lt*/
} /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_int_res*/


#endif
/*! @} (documentation module close)*/
