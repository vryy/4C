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
#ifdef D_FLUID3_F

#include "../headers/standardtypes.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../fluid3_fast/f3f_prototypes.h"
#include "../fluid3/fluid3.h"

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
/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01  |
  | vector of material laws                                            |
  | defined in global_control.c                                        |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

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
void f3fint_usfem(
	              ELEMENT         *ele[LOOPL],
                      DOUBLE          *elecord,
                      DOUBLE          *tau,
                      INT             *hasext,
                      DOUBLE          *estif,
	              DOUBLE          *force,
	              DOUBLE          *funct,
	              DOUBLE          *deriv,
	              DOUBLE          *deriv2,
	              DOUBLE          *xjm,
	              DOUBLE          *derxy,
	              DOUBLE          *derxy2,
	              DOUBLE          *evelng,
	              DOUBLE          *evhist,
	              DOUBLE          *egridv,
	              DOUBLE          *epren,
	              DOUBLE          *edeadng,
                      DOUBLE          *velint,
                      DOUBLE          *histint,
                      DOUBLE          *gridvelint,
	              DOUBLE          *vderxy,
                      DOUBLE          *vderxy2,
                      DOUBLE          *pderxy,
	              DOUBLE          *wa1,
	              DOUBLE          *wa2,
                      INT              sizevec[6]
	             )
{ 
INT       l;          /* a couter                                       */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis,nit;/* number of integration nodesin r,s,t direction  */
INT       icode=2;    /* flag for eveluation of shape functions         */
INT       lr, ls, lt; /* counter for integration                        */
INT       actmat;     /* actual material number                         */
DOUBLE    fac[LOOPL]; /* total integration vactor                       */
DOUBLE    facr=0.0, facs=0.0, fact=0.0;
                      /* integration weights                            */
DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
DIS_TYP   typ;	      /* element type                                   */

FLUID_DYNAMIC   *fdyn;
FLUID_DATA      *data;

/*Fortran variables - passed to fortran subroutines as parameters.*/
INT      flagvec[2];
DOUBLE   paravec[4];

INT      inttyp;
DOUBLE   det[LOOPL];

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("f3fint_usfem");
#endif

/*--------------------------------------------------- initialisation ---*/
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;
actmat = ele[0]->mat-1;
typ    = ele[0]->distyp;

flagvec[0] = 0;                                  /* ihoel */
flagvec[1] = ele[0]->e.f3->is_ale;               /* isale */

paravec[0] = 1.0;
paravec[1] = mat[actmat].m.fluid->viscosity;
paravec[2] = fdyn->thsl;
paravec[3] = fdyn->dta;

intc = 0;
nir  = 0;
nis  = 0;
nit  = 0;


switch (typ)
{
case hex8: case hex20: case hex27:  /* hex - element */
   icode   = 3;
   flagvec[0]   = 1;
   /* initialise integration */
   nir = ele[0]->e.f3->nGP[0];
   nis = ele[0]->e.f3->nGP[1];
   nit = ele[0]->e.f3->nGP[2];
   intc= 0;
break;

case tet10: /* tet - element */
   icode   = 3;
   flagvec[0]   = 1;
   /* do NOT break at this point!!! */

case tet4:    /* initialise integration */
   nir  = ele[0]->e.f3->nGP[0];
   nis  = 1;
   nit  = 1;
   intc = ele[0]->e.f3->nGP[1];
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
   case hex8: case hex20: case hex27:   /* hex - element */
      e1   = data->qxg[lr][nir-1];
      e2   = data->qxg[ls][nis-1];
      e3   = data->qxg[lt][nit-1];
      facr = data->qwgt[lr][nir-1];
      facs = data->qwgt[ls][nis-1];
      fact = data->qwgt[lt][nit-1];

      if (typ==hex8)
         inttyp=8;
      else if(typ==hex20)
         inttyp=20;
      else if(typ==hex27)
         inttyp=27;

      f3fhex(funct, deriv, deriv2, &e1, &e2, &e3, &inttyp, &icode, sizevec);

   break;
   case tet4: case tet10:  /* tet - element */
      e1   = data->txgr[lr][intc];
      facr = data->twgt[lr][intc];
      e2   = data->txgs[lr][intc];
      facs = ONE;
      e3   = data->txgt[lr][intc];
      fact = ONE;

      if (typ==tet4)
         inttyp=4;
      else if(typ==tet10)
         inttyp=10;

      f3ftet(funct, deriv, deriv2, &e1, &e2, &e3, &inttyp, 
             &icode, sizevec);

   break;
   default:
      dserror("typ unknown!");
   } /* end switch (typ) */

   /*----------------------------------------- compute Jacobian matrix -*/
   f3fjaco(funct, deriv, xjm, det, elecord, sizevec);
   for(l=0;l<sizevec[4];l++)
      fac[l] = facr * facs * fact * det[l];

   /*---------------------------------------- compute global derivates -*/
   f3fgder(derxy, deriv, xjm, wa1, det, sizevec);

   /*--------------------------------- compute second global derivative */
   if (flagvec[0]!=0)
   {
      f3fgder2loop(elecord, xjm, wa1, wa2, derxy, derxy2, deriv2, sizevec);
      f3fvder2(vderxy2, derxy2, evelng, sizevec);
   }

   /*----------------------- get velocities (n+g,i) at integraton point */
   f3fveli(velint, funct, evelng, sizevec);
   
   /*------------------- get history data (n,i) at integration point ---*/
   f3fveli(histint, funct, evhist, sizevec);

   /*----------- get velocity (n+g,i) derivatives at integration point -*/
   f3fvder(vderxy, derxy, evelng, sizevec);

   /*--------------------- get grid velocity at integration point ---*/
   if(flagvec[1]!=0) f3fveli(gridvelint, funct, egridv, sizevec);
   
   /*------------------------------------- get pressure gradients ---*/
   f3fpder(pderxy, derxy, epren, sizevec);

   /*-------------- perform integration for entire matrix and rhs ---*/
   f3fcalmat(estif,force,velint,histint,gridvelint,vderxy,
             vderxy2,pderxy,funct,derxy,derxy2,edeadng,fac,tau,
             hasext,paravec,flagvec,sizevec);
} /* end of loop over integration points lt*/
} /* end of loop over integration points ls */
} /* end of loop over integration points lr */
 
/*----------------------------------------------- to ensure assembly ---*/
for(l=0;l<sizevec[4];l++)
   hasext[l] = 1;

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3fint_usfem_ale */



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
void f3fint_res(
	        ELEMENT         *ele[LOOPL],
                DOUBLE          *elecord,
                DOUBLE          *tau,
                INT             *hasext,
	        DOUBLE          *force,
	        DOUBLE          *funct,
	        DOUBLE          *deriv,
	        DOUBLE          *deriv2,
	        DOUBLE          *xjm,
	        DOUBLE          *derxy,
	        DOUBLE          *derxy2,
	        DOUBLE          *evelng,
	        DOUBLE          *evhist,
	        DOUBLE          *ealecovng,
	        DOUBLE          *epren,
	        DOUBLE          *edeadng,
                DOUBLE          *velint,
                DOUBLE          *histint,
                DOUBLE          *aleconvint,
	        DOUBLE          *vderxy,
                DOUBLE          *vderxy2,
                DOUBLE          *pderxy,
	        DOUBLE          *wa1,
	        DOUBLE          *wa2,
                INT              sizevec[6]
	       )
{ 
INT       l;          /* a couter                                       */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis,nit;/* number of integration nodesin r,s,t direction  */
INT       icode=2;    /* flag for eveluation of shape functions         */
INT       actmat;     /* actual material number                         */
INT       lr, ls, lt; /* counter for integration                        */

DOUBLE    fac[LOOPL];       /* total integration vector                 */
DOUBLE    facr=0.0, facs=0.0, fact=0.0;
                      /* integration weights                            */
DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
DOUBLE    presint[LOOPL];
DIS_TYP   typ;	      /* element type                                   */
FLUID_DYNAMIC   *fdyn;
FLUID_DATA      *data;

/*Fortran variables - passed to fortran subroutines as parameters.*/
INT      flagvec[2];
DOUBLE   paravec[4];

INT      inttyp;
DOUBLE   det[LOOPL];

#ifdef DEBUG 
dstrc_enter("f3fint_res");
#endif

/*--------------------------------------------------- initialisation ---*/
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;
actmat = ele[0]->mat-1;
typ    = ele[0]->distyp;

flagvec[0] = 0;                                  /* ihoel */
flagvec[1] = ele[0]->e.f3->is_ale;               /* isale */

paravec[0] = 1.0;
paravec[1] = mat[actmat].m.fluid->viscosity;
paravec[2] = fdyn->thsl;
paravec[3] = fdyn->dta;

intc = 0;
nir  = 0;
nis  = 0;
nit  = 0;


/*------- get integraton data and check if elements are "higher order" */
switch (typ)
{
case hex8: case hex20: case hex27:  /* --> hex - element */
   icode      = 3;
   flagvec[0] = 1;
   /* initialise integration */
   nir = ele[0]->e.f3->nGP[0];
   nis = ele[0]->e.f3->nGP[1];
   nit = ele[0]->e.f3->nGP[2];
   break;
case tet10: /* --> tet - element */
   icode      = 3;
   flagvec[0] = 1;
/* do NOT break at this point!!! */
case tet4:    /* initialise integration */
   nir  = ele[0]->e.f3->nGP[0];
   nis  = 1;
   nit  = 1;
   intc = ele[0]->e.f3->nGP[1];
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

      if (typ==hex8)
         inttyp=8;
      else if(typ==hex20)
         inttyp=20;
      else if(typ==hex27)
         inttyp=27;

      f3fhex(funct, deriv, deriv2, &e1, &e2, &e3, &inttyp, &icode, sizevec);

   break;
   case tet4: case tet10:   /* --> tet - element */
      e1   = data->txgr[lr][intc];
      facr = data->twgt[lr][intc];
      e2   = data->txgs[lr][intc];
      facs = ONE;
      e3   = data->txgt[lr][intc];
      fact = ONE;

      if (typ==tet4)
         inttyp=4;
      else if(typ==tet10)
         inttyp=10;

      f3ftet(funct, deriv, deriv2, &e1, &e2, &e3, &inttyp, 
             &icode, sizevec);
   break;
   default:
      dserror("typ unknown!");
   } /* end switch (typ) */

   /*----------------------------------------- compute Jacobian matrix -*/
   f3fjaco(funct, deriv, xjm, det, elecord, sizevec);
   for(l=0;l<sizevec[4];l++)
      fac[l] = facr * facs * fact * det[l];

   /*---------------------------------------- compute global derivates -*/
   f3fgder(derxy, deriv, xjm, wa1, det, sizevec);

   /*--------------------------------- compute second global derivative */
   if (flagvec[0]!=0)
   {
      f3fgder2loop(elecord, xjm, wa1, wa2, derxy, derxy2, deriv2, sizevec);
      f3fvder2(vderxy2, derxy2, evelng, sizevec);
   }

   /*----------------------- get velocities (n+g,i) at integraton point */
   f3fveli(velint, funct, evelng, sizevec);
   
   /*------------------- get history data (n,i) at integration point ---*/
   f3fveli(histint, funct, evhist, sizevec);

   /*----------- get velocity (n+g,i) derivatives at integration point -*/
   f3fvder(vderxy, derxy, evelng, sizevec);

   /*--------------------- get grid velocity at integration point ---*/
   if(flagvec[1]!=0) f3fveli(aleconvint, funct, ealecovng, sizevec);
   
   /*------------------------------------- get pressure gradients ---*/
   f3fpder(pderxy, derxy, epren, sizevec);

   /*----------------------------- get pressure at integration point ---*/
   f3fprei(presint,funct,epren,sizevec);
      
   /*----------------- perform integration for entire matrix and rhs ---*/
   f3fcalresvec(force,velint,histint,vderxy,vderxy2,funct,derxy,derxy2,
                edeadng,aleconvint,presint,pderxy,fac,hasext,tau,
                paravec,flagvec,sizevec);
} /* end of loop over integration points lt*/
} /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3fint_res*/


#endif
/*! @} (documentation module close)*/
