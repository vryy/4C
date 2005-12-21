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
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
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
The fully linearised 2D fluid element is called here!
      
</pre>
\param  *ele	   ELEMENT	   (i)    actual element
\param  *hasext    INT             (i)    element flag
\param **estif     DOUBLE	   (o)    element stiffness matrix
\param  *eforce    DOUBLE	   (o)    element iter force vector
\param **xyze      DOUBLE          (-)    nodal coordinates
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **deriv2    DOUBLE	   (-)    2nd deriv. of nat. shape f.
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **derxy2    DOUBLE	   (-)    2nd global derivatives
\param **evelng    DOUBLE	   (i)    ele vel. at time n+g
\param **evhist    DOUBLE	   (i)    lin. combination of recent vel and acc
\param **egridv    DOUBLE	   (i)    grid velocity of element
\param  *epren     DOUBLE	   (-)    ele pres. at time n
\param  *edeadng   DOUBLE	   (-)    ele dead load (selfweight) at n+1
\param  *velint    DOUBLE	   (-)    vel at integration point
\param  *vel2int   DOUBLE	   (-)    vel at integration point
\param  *covint    DOUBLE	   (-)    conv. vel. at integr. point
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param **vderxy2   DOUBLE	   (-)    2nd global vel. deriv.
\param **wa1	   DOUBLE	   (-)    working array
\param **wa2	   DOUBLE	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f2_int_usfem(
	              ELEMENT         *ele,
                      INT             *hasext,
                      DOUBLE         **estif,
	              DOUBLE          *eforce,
	              DOUBLE         **xyze,
	              DOUBLE          *funct,
	              DOUBLE         **deriv,
	              DOUBLE         **deriv2,
	              DOUBLE         **xjm,
	              DOUBLE         **derxy,
	              DOUBLE         **derxy2,
	              DOUBLE         **evelng,
	              DOUBLE         **eveln,
	              DOUBLE         **evhist,
	              DOUBLE         **egridv,
	              DOUBLE          *epren,
	              DOUBLE          *edeadng,
	              DOUBLE         **vderxy,
                      DOUBLE         **vderxy2,
                      DOUBLE           visc,
	              DOUBLE         **wa1,
	              DOUBLE         **wa2,
                      DOUBLE           estress[3][MAXNOD_F2]
	             )
{ 
INT       i;          /* a couter                                       */
INT       iel;        /* number of nodes                                */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration nodesin r,s direction    */
INT       ihoel=0;    /* flag for higher order elements                 */
INT       icode=2;    /* flag for eveluation of shape functions         */     
INT       lr, ls;     /* counter for integration                        */
INT       is_ale;
DOUBLE    fac;        /* total integration vactor                       */
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix at time (n+1)   */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    gradp[2];   /* pressure gradient at integration point         */
DOUBLE    velint[2];  /* velocity vector at integration point           */
DOUBLE    histvec[2]; /* history data at integration point              */
DOUBLE    gridvelint[2]; /* grid velocity                               */
DOUBLE    divuold;
DIS_TYP   typ;	      /* element type                                   */

FLUID_DYNAMIC   *fdyn;
FLUID_DATA      *data;

#ifdef DEBUG 
dstrc_enter("f2_int_usfem");
#endif

/*--------------------------------------------------- initialisation ---*/
iel    = ele->numnp;
typ    = ele->distyp;
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

is_ale = ele->e.f2->is_ale;

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
/*------------- get values of  shape functions and their derivatives ---*/
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
      
      /*------------------------ compute Jacobian matrix at time n+1 ---*/
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr * facs * det;
      
      /*----------------------------------- compute global derivates ---*/
      f2_gder(derxy,deriv,xjm,det,iel);

      /*---------------- get velocities (n+1,i) at integration point ---*/
      f2_veci(velint,funct,evelng,iel);

      /*---------------- get history data (n,i) at integration point ---*/
      f2_veci(histvec,funct,evhist,iel);
      
      /*--------------------- get grid velocity at integration point ---*/
      if(is_ale) f2_veci(gridvelint,funct,egridv,iel);
           
      /*-------- get velocity (n,i) derivatives at integration point ---*/
      f2_vder(vderxy,derxy,eveln,iel);
      divuold = vderxy[0][0] + vderxy[1][1];

      /*------ get velocity (n+1,i) derivatives at integration point ---*/
      f2_vder(vderxy,derxy,evelng,iel);

      /*--------------------------- compute second global derivative ---*/ 
      if(fdyn->stresspro == 0)/* no stress projection, do second derivs */
      {
         if (ihoel!=0)
         {
            f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
            f2_vder2(vderxy2,derxy2,evelng,iel);
         }
      }
      else      /* get second derivatives from global stress projection */
      {
         /* Eine Hilfskonstruktion: */
         f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);

         vderxy2[0][0] = vderxy2[0][1] = vderxy2[0][2]= 0.0;
         vderxy2[1][0] = vderxy2[1][1] = vderxy2[1][2]= 0.0;

         for (i=0; i<iel; i++)
         {
            vderxy2[0][0] += estress[0][i] * derxy[0][i];
            vderxy2[0][1] += estress[2][i] * derxy[1][i];
            vderxy2[0][2] += estress[0][i] * derxy[1][i];
            vderxy2[1][0] += estress[2][i] * derxy[0][i];
            vderxy2[1][1] += estress[1][i] * derxy[1][i];
            vderxy2[1][2] += estress[2][i] * derxy[1][i];
         }
      }

      /*------------------------------------- get pressure gradients ---*/
      gradp[0] = gradp[1] = 0.0;
      
      for (i=0; i<iel; i++)
      {
         gradp[0] += derxy[0][i] * epren[i];
         gradp[1] += derxy[1][i] * epren[i];
      }
      /*-------------- perform integration for entire matrix and rhs ---*/
      f2_calmat(estif,eforce,velint,histvec,gridvelint,vderxy,
                vderxy2,gradp,funct,derxy,derxy2,edeadng,fac,
                visc,iel,hasext,is_ale);
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

/*------------------------------------------- assure assembly of rhs ---*/
*hasext = 1;
return; 
} /* end of f2_int_usfem_ale */



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
\param **evelng    DOUBLE	   (i)    ele vel. at time n+g
\param **evhist    DOUBLE	   (i)    lin. combination of recent vel and acc
\param **ealecovng DOUBLE	   (i)    ALE convective velocity
\param  *epren     DOUBLE	   (-)    ele pres. at time n
\param  *edeadng   DOUBLE	   (-)    ele dead load (selfweight) at n+1
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param **vderxy2   DOUBLE	   (-)    2nd global vel. deriv.
\param   visc      DOUBLE          (i)    viscosity
\param **wa1	   DOUBLE	   (-)    working array
\param **wa2	   DOUBLE	   (-)    working array
\return void                                                   

------------------------------------------------------------------------*/
void f2_int_res(
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
	        DOUBLE         **wa2,
                DOUBLE           estress[3][MAXNOD_F2]
	       )
{ 
INT       i;          /* a couter                                       */
INT       iel;        /* number of nodes                                */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       is_ale;     /* ALE or Euler element flag                      */
INT       nir,nis;    /* number of integration nodesin r,s direction    */
INT       ihoel=0;    /* flag for higher order elements                 */
INT       icode=2;    /* flag for eveluation of shape functions         */     
INT       lr, ls;     /* counter for integration                        */

DOUBLE    fac;        /* total integration vactor                       */
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DOUBLE    velint[2];  /* velocity vector at integration point           */
DOUBLE    histvec[2]; /* history data at integration point              */
DOUBLE    aleconv[2]; /* ALE convective velocity at Gauss point         */
DIS_TYP   typ;	      /* element type                                   */
DOUBLE    presint;    /* pressure at integration point                  */
DOUBLE    gradp[2];   /* pressure gradient                              */
FLUID_DYNAMIC   *fdyn;
FLUID_DATA      *data;


#ifdef DEBUG 
dstrc_enter("f2_int_usfem");
#endif

/*--------------------------------------------------- initialisation ---*/
iel    = ele->numnp;
typ    = ele->distyp;
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

is_ale = ele->e.f2->is_ale;

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
/*------------- get values of  shape functions and their derivatives ---*/
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
      
      /*------------------------------------ compute Jacobian matrix ---*/
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;

      /*----------------------------------- compute global derivates ---*/
      f2_gder(derxy,deriv,xjm,det,iel);

      /*------------------------ get velocities at integration point ---*/
      f2_veci(velint,funct,evelng,iel);

      /*---------------------- get history data at integration point ---*/
      f2_veci(histvec,funct,evhist,iel);

      /*----------- get ALE convective velocity at integration point ---*/
      if(is_ale) f2_veci(aleconv,funct,ealecovng,iel);

      /*-------------- get velocity derivatives at integration point ---*/
      f2_vder(vderxy,derxy,evelng,iel);

      /*--------------------------- compute second global derivative ---*/ 
      if(fdyn->stresspro == 0)/* no stress projection, do second derivs */
      {
         if (ihoel!=0)
         {
            f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
            f2_vder2(vderxy2,derxy2,evelng,iel);
         }
      }
      else      /* get second derivatives from global stress projection */
      {
         /* Eine Hilfskonstruktion: */
         f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel); 
         
         vderxy2[0][0] = vderxy2[0][1] = vderxy2[0][2]= 0.0;
         vderxy2[1][0] = vderxy2[1][1] = vderxy2[1][2]= 0.0;

         for (i=0; i<iel; i++)
         {
            vderxy2[0][0] += estress[0][i] * derxy[0][i];
            vderxy2[0][1] += estress[2][i] * derxy[1][i];
            vderxy2[0][2] += estress[0][i] * derxy[1][i];
            vderxy2[1][0] += estress[2][i] * derxy[0][i];
            vderxy2[1][1] += estress[1][i] * derxy[1][i];
            vderxy2[1][2] += estress[2][i] * derxy[1][i];
         }
      }

      /*------------------------------------- get pressure gradients ---*/
      gradp[0] = gradp[1] = 0.0;
      
      for (i=0; i<iel; i++)
      {
         gradp[0] += derxy[0][i] * epren[i];
         gradp[1] += derxy[1][i] * epren[i];
      }
      /*-------------------------- get pressure at integration point ---*/
      presint = f2_scali(funct,epren,iel);
      
      /*-------------- perform integration for entire matrix and rhs ---*/
      f2_calresvec(force,velint,histvec,vderxy,vderxy2,funct,derxy,derxy2,
                   edeadng,aleconv,&presint,gradp,fac,visc,iel,hasext,
                   is_ale);
   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_int_res*/



/*!---------------------------------------------------------------------
\brief integration loop for integration of first derivatives 
       for projection

<pre>                                                         chfoe 01/05

This routine integrates the matrix and rhs to serve the L2-projection of 
first derivatives to the nodes. This is needed further to repare the 
Navier-Stokes equations residuum for linear elements which can not 
recover second derivatives.

see also: Jansen, KE, Collis SS, Whiting C, Shakib, F: 'A better consistency
          for low-oder stabilized finite element methods', Computer Methods
          174:153-170, 1999.

NOTE: the flag hasext is misused here in order to ensure assembly of all
      elements.

</pre>
\param  *ele	   ELEMENT	   (i)    actual element
\param  *hasext    INT             (o)    flag
\param **estif     DOUBLE          (o)    elemental stiffness matrix
\param  *force     DOUBLE	   (o)    elemental force vector to be filled
\param **xyze      DOUBLE          (-)    nodal coordinates
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **evelng    DOUBLE	   (i)    ele vel. at time n+g
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\return void                                                   

------------------------------------------------------------------------*/
void f2_int_stress_project(
	                   ELEMENT         *ele,
                           INT             *hasext,
                           DOUBLE         **estif,
                           DOUBLE          *force,
                           DOUBLE         **xyze,
                           DOUBLE          *funct,
                           DOUBLE         **deriv,
                           DOUBLE         **xjm,
                           DOUBLE         **derxy,
                           DOUBLE         **evelng,
                           DOUBLE         **vderxy
	                  )
{ 
INT       i,j;        /* couters                                        */
INT       iel;        /* number of nodes                                */
INT       intc;       /* "integration case" for tri for further infos
                          see f2_inpele.c and f2_intg.c                 */
INT       nir,nis;    /* number of integration nodesin r,s direction    */
INT       icode=2;    /* flag for eveluation of shape functions         */     
INT       lr, ls;     /* counter for integration                        */

DOUBLE    aux;
DOUBLE    fac;        /* total integration vactor                       */
DOUBLE    facr, facs; /* integration weights                            */
DOUBLE    det;        /* determinant of jacobian matrix                 */
DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
DIS_TYP   typ;	      /* element type                                   */

FLUID_DYNAMIC   *fdyn;
FLUID_DATA      *data;


#ifdef DEBUG 
dstrc_enter("f2_int_stress_project");
#endif

/*--------------------------------------------------- initialisation ---*/
iel    = ele->numnp;
typ    = ele->distyp;
fdyn   = alldyn[genprob.numff].fdyn;
data   = fdyn->data;

*hasext = 1; /* this is to ensure assembly */

/*------- get integraton data and check if elements are "higher order" */
switch (typ)
{
case quad4: case quad8: case quad9:  /* --> quad - element */
   icode   = 2;
   /* initialise integration */
   nir = ele->e.f2->nGP[0];
   nis = ele->e.f2->nGP[1];
break;
case tri6: /* --> tri - element */
   icode   = 2;
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
/*------------- get values of  shape functions and their derivatives ---*/
      switch(typ)
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
         e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
         e2   = data->qxg[ls][nis-1];
         facs = data->qwgt[ls][nis-1];
         f2_rec(funct,deriv,NULL,e1,e2,typ,icode);
      break;
      case tri3: case tri6:   /* --> tri - element */
         e1   = data->txgr[lr][intc];
         facr = data->twgt[lr][intc];
         e2   = data->txgs[lr][intc];
         facs = ONE;
         f2_tri(funct,deriv,NULL,e1,e2,typ,icode);
      break;
      default:
         dserror("typ unknown!");
      } /* end switch(typ) */
      
      /*------------------------------------ compute Jacobian matrix ---*/
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;

      /*----------------------------------- compute global derivates ---*/
      f2_gder(derxy,deriv,xjm,det,iel);

      /*-------------- get velocity derivatives at integration point ---*/
      f2_vder(vderxy,derxy,evelng,iel);
      
      /*------------------------------- integrate projection problem ---*/
      for (i=0; i<iel; i++)
      {
         for (j=0; j<iel; j++)
         {
            aux = funct[i] * funct[j] * fac;
            estif[i*3][j*3]     += aux;
            estif[i*3+1][j*3+1] += aux;
            estif[i*3+2][j*3+2] += aux;
         }
         /*  2*u_x,x  */
         force[i*3]   = vderxy[0][0] * funct[i] * fac;
         /*  2*u_y,y  */     
         force[i*3+1] = vderxy[1][1] * funct[i] * fac;
         /*u_x,y+u_y,x*/
         force[i*3+2] = 0.5 * (vderxy[0][1]+vderxy[1][0])*funct[i] * fac;
      }

   } /* end of loop over integration points ls*/
} /* end of loop over integration points lr */
 
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_int_stress*/



/*!---------------------------------------------------------------------
\brief routine to evaluate the stabilisation parameter within USFEM

<pre>                                                         chfoe 02/05

This routine evaluates the stabilisation parameters tau_M and tau_C
depending on the flags 'whichtau' and 'which_hk'. The respective values
are:

whichtau == 0 -> USFEM tau as reported by Barrenechea, Valentin and Franca

see also: Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite 
          element method for a generalized Stokes problem. Numerische 
          Mathematik, Vol. 92, pp. 652-677, 2002.
          http://www.lncc.br/~valentin/publication.htm
and:      Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized 
          Finite Element Method for the Advective-Reactive-Diffusive 
          Equation. Computer Methods in Applied Mechanics and Enginnering,
          Vol. 190, pp. 1785-1800, 2000.
          http://www.lncc.br/~valentin/publication.htm

  for this stabilisation parameter a number of different element length
  definitions is implemented.

  which_hk == 0  ->  hk = square root of area
  which_hk == 1  ->  hk = length in flow direction (-> Wall)
  which_hk == 2  ->  hk = approximate length in flow direction (Codina)
  which_hk == 3  ->  hk = second eigenvalue of left stretch tensor (Codina)

see also: Codina, R. and Soto, O.: Approximation of the incompressible
          Navier-Stokes equations using orthogonal subscale stabilisation 
          and pressure segregation on anisotropic finite element meshes. 
          Computer methods in Applied Mechanics and Engineering, 
          Vol 193, pp. 1403-1419, 2004.

further the parameter following Taylor and Hughes is implemented under 
whichtau == 1

see also: Taylor, Charles A. and Hughes, Thomas J.R. and Zarins, 
          Christopher K.: Finite element modeling of blood flow in 
          arteries. Computer Methods in Applied Mechanics and Engineering,
          Vol 158, pp. 155-196, 1998.
          (or better the dissertation of Whiting, Chapter 3)

This parameter uses an implicit definition for the element length which is
automatically used not regarding the value of which_hk.

</pre>
\param  *ele	   ELEMENT	   (i)    actual element
\param  *hasext    INT             (o)    flag
\param **estif     DOUBLE          (o)    elemental stiffness matrix
\param  *force     DOUBLE	   (o)    elemental force vector to be filled
\param **xyze      DOUBLE          (-)    nodal coordinates
\param  *funct     DOUBLE	   (-)    natural shape functions
\param **deriv     DOUBLE	   (-)	  deriv. of nat. shape funcs
\param **xjm	   DOUBLE	   (-)    jacobian matrix
\param **derxy     DOUBLE	   (-)	  global derivatives
\param **evelng    DOUBLE	   (i)    ele vel. at time n+g
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\return void                                                   

------------------------------------------------------------------------*/
void f2_get_tau(ELEMENT *ele,
                DOUBLE  **xjm,
                DOUBLE  **xyze,
                DOUBLE   *funct,
                DOUBLE    det, 
                DOUBLE   *velint, 
                DOUBLE    visc,
                INT       whichtau,
                INT       which_hk)
{

INT       i;              /* a counter                                  */
INT       iel;            /* number of element nodal points             */
DOUBLE    timefac;        /* factor from time integration (contains dt) */
DOUBLE    metric[3];      /* covariant metric tensor                    */
DOUBLE    gij[3];         /* contravariant metric tensor                */
DOUBLE    detJ, gijgij;
DOUBLE    advec, norm_p;
DOUBLE    trace, u0[2];
DOUBLE    mk,hk,h0;
DOUBLE    pe, re, xi1, xi2;
DOUBLE    fact;
DOUBLE    gcoor[2];       /* global coordinates of actual point         */
DOUBLE    velsquared;

/* nach Moeglichkeit zu entfernen:*/
DOUBLE    xyz[2][MAXNOD];
DOUBLE    area;

DIS_TYP        typ;

FLUID_DYNAMIC   *fdyn;

#ifdef DEBUG 
dstrc_enter("f2_get_tau");
#endif

/*--------------------------------------------------- initialisation ---*/
fdyn    = alldyn[genprob.numff].fdyn;
timefac = fdyn->thsl;

typ = ele->distyp;
iel = ele->numnp;

/*--------------------- decide the general way of calculation of tau ---*/
switch (whichtau)
{
case 0: /* Franca */
   /*--- get proper constant mk and h0 ---*/
   switch(typ)
   {
   case quad4:
      mk = 0.333333333333333333333;
      h0 = 2.0;
   break;
   case quad8:
      mk = 0.083333333333333333333;
      h0 = 2.0;
   break;
   case quad9:
      mk = 0.083333333333333333333;
      h0 = 2.0;
   break;
   default: dserror("element type not implemented!");
   }
   /*----------------- choose appropriate element length calculation ---*/
   switch (which_hk)
   {
   case 0: /* square root of area */
      /*--------------------- rewrite array of elemental coordinates ---*/
      for(i=0; i<iel; i++)
      {
         xyz[0][i] = xyze[0][i];
         xyz[1][i] = xyze[1][i];
      }
      /*-------------------------------- get area and element length ---*/
      area = area_lin_2d(ele,xyz);
      hk = sqrt(area);
   break;
   case 1: /* length in flow direction (Wall) */
      hk = 0.0;
      f2_gcoor(xyze,funct,iel,gcoor);
      f2_calstrlen(&hk,xyze,velint,ele,gcoor,typ);
   break;
   case 2: /* approximate length in flow direction */
      /*--------------------- get hk asuming length in parent domain ---*/
      velsquared = velint[0]*velint[0] + velint[1]*velint[1];
      
      if(velsquared == ZERO)
      {
         /*------------------ rewrite array of elemental coordinates ---*/
         for(i=0; i<iel; i++)
         {
            xyz[0][i] = xyze[0][i];
            xyz[1][i] = xyze[1][i];
         }
         /*----------------------------- get area and element length ---*/
         area = area_lin_2d(ele,xyz);
         hk = sqrt(area); 
      }
      else
      {
         /*-------------- transform velocity vector: u0 = J^{-1} * u ---*/
         u0[0] = 1.0/det * ( xjm[1][1] * velint[0] - xjm[0][1] * velint[1]);
         u0[1] = 1.0/det * (-xjm[1][0] * velint[0] + xjm[0][0] * velint[1]);
         hk = sqrt( velsquared / (u0[0]*u0[0] + u0[1]*u0[1]) ) * h0;
      }
   break;
   case 3: /* 'for anisotropic meshes' */
      /*-------------------------------- get half trace of (J * J^T) ---*/
      trace = (xjm[0][0]*xjm[0][0] + xjm[0][1]*xjm[0][1] 
             + xjm[1][0]*xjm[1][0] + xjm[1][1]*xjm[1][1]) * 0.5;
      /* a quick check for numerical rubbish */
      if (trace<det) /* this can be caused by numerical dust only */
      {
         hk = sqrt(trace) * 2.0;
      }
      else hk = sqrt(trace - sqrt(trace*trace - det*det)) * 2.0;
   break;
   default: dserror("Evaluation of element length unknown!");
   }

   /*---------------------------------------------------- get p-norm ---*/
   norm_p = sqrt(DSQR(velint[0]) + DSQR(velint[1]));

   pe = 4.0 * timefac * visc / (mk * DSQR(hk)); /* viscous : reactive forces */
   re = mk * norm_p * hk / (2.0 * visc);       /* advective : viscous forces */

   xi1 = DMAX(pe,1.0);
   xi2 = DMAX(re,1.0);

   fdyn->tau[0] = DSQR(hk) / (DSQR(hk) * xi1 + (4.0*timefac * visc/mk) * xi2);
   fdyn->tau[2] = norm_p * hk * 0.5 * xi2;
break;

case 1: /* Whiting tau */
   switch(typ)
   {
   case quad4:
      fact = 36.0;
   break;
   case quad8:
   case quad9:
      fact = 60.0;
   break;
   default: dserror("not implemented!");
   }
   /*----------------------------------- get covariant metric tensor ---*/
            /* 1st diagonal element */
   metric[0] = xjm[0][0]*xjm[0][0] + xjm[1][0]*xjm[1][0];
            /* 2nd diagonal element */
   metric[1] = xjm[0][1]*xjm[0][1] + xjm[1][1]*xjm[1][1];
            /* off diagonal element */
   metric[2] = xjm[0][0]*xjm[0][1] + xjm[1][0]*xjm[1][1];
   /*------------------------------- get contravariant metric tensor ---*/
   /*detJ = 1.0/(metric[0] * metric[1] - metric[2] * metric[2]);*/
   detJ = 1.0/( det * det ); 
   gij[0] =  metric[1] * detJ;   /* 1st diagonal element */
   gij[1] =  metric[0] * detJ;   /* 2nd diagonal element */
   gij[2] = -metric[2] * detJ;   /* off diagonal element */

   gijgij =gij[0]*gij[0] + gij[1]*gij[1] + 2.0 * gij[2]*gij[2];
   advec = velint[0]*gij[0]*velint[0] 
          +velint[1]*gij[1]*velint[1]
          +velint[0]*gij[2]*velint[1]*2.0;

   fdyn->tau[0] = 1.0 / sqrt(4.0 + timefac*timefac*(advec 
                             + fact * visc * visc * gijgij));
   /* tau_C = ( 8 * tau_M * tr(gij) )^{-1} */
   fdyn->tau[2] = 0.125 / (fdyn->tau[0] * (gij[0]+gij[1]));
break;

default: dserror(" Way of tau calculation unknown!");
} /* end switch (whichtau) */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
}



/*-----------------------------------------------------------------------*/
/*!
  \brief integration loop for the error calculation


  \param  *ele       ELEMENT   (i)    actual element
  \param **xyze      DOUBLE    (i)    nodal coordinates
  \param  *funct     DOUBLE    (i)    natural shape functions
  \param **deriv     DOUBLE    (i)    deriv. of nat. shape funcs
  \param **xjm       DOUBLE    (i)    jacobian matrix
  \param **evelng    DOUBLE    (i)    ele vel. at time n+g
  \param   visc      DOUBLE    (i)    viscosity
  \param  *epren     DOUBLE    (i)    ele pre
  \param  *container CONTAINER (i)    contains variables defined in container.h

  \return void

  \author mn
  \date   08/05

 */
/*-----------------------------------------------------------------------*/
void f2_int_kim_moin_err(
    ELEMENT         *ele,
    DOUBLE         **xyze,
    DOUBLE          *funct,
    DOUBLE         **deriv,
    DOUBLE         **xjm,
    DOUBLE         **evelng,
    DOUBLE           visc,
    DOUBLE          *epren,
    CONTAINER       *container
    )

{
  INT       i;          /* a couter                                       */
  INT       iel;        /* number of nodes                                */
  INT       intc;       /* "integration case" for tri for further infos
                           see f2_inpele.c and f2_intg.c                 */
  INT       is_ale;     /* ALE or Euler element flag                      */
  INT       nir,nis;    /* number of integration nodesin r,s direction    */
  INT       actmat;     /* material number of the element                 */
  INT       lr, ls;     /* counter for integration                        */

  DOUBLE    fac;        /* total integration vactor                       */
  DOUBLE    facr, facs; /* integration weights                            */
  DOUBLE    det;        /* determinant of jacobian matrix                 */
  DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
  DOUBLE    preint;     /* pressure at integration point                  */
  DOUBLE    velint[2];  /* velocity vector at integration point           */
  DOUBLE    xint[2];    /* coordinates at integration point               */
  DIS_TYP   typ;        /* element type                                   */
  DOUBLE    divuold;    /* velocity divergence at t=t^n                   */
  DOUBLE    presint;    /* pressure at integration point                  */
  DOUBLE    diffx, diffy, diffp;
  DOUBLE    solx, soly, solp;
  DOUBLE    a=2.0, t;
  FLUID_DYNAMIC   *fdyn;
  FLUID_DATA      *data;
  NODE     *actnode;


#ifdef DEBUG
  dstrc_enter("f2_int_kim_moin_err");
#endif


  /* initialisation */
  iel    = ele->numnp;
  typ    = ele->distyp;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  t = fdyn->acttime;

  is_ale = ele->e.f2->is_ale;

  /* get integraton data and check if elements are "higher order" */
  switch (typ)
  {
    case quad4: case quad8: case quad9:  /* --> quad - element */
      /* initialise integration */
      nir = ele->e.f2->nGP[0];
      nis = ele->e.f2->nGP[1];
      break;

    case tri6: /* --> tri - element */
    case tri3:
      /* initialise integration */
      nir  = ele->e.f2->nGP[0];
      intc = ele->e.f2->nGP[1];
      break;

    default:
      dserror("typ unknown!");

  }  /* switch (typ) */


  /*----------------------------------------------------------------------*
    |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (lr=0;lr<nir;lr++)
  {
    for (ls=0;ls<nis;ls++)
    {

      /* get values of  shape functions and their derivatives */
      switch(typ)
      {
        case quad4: case quad8: case quad9:   /* --> quad - element */
          e1   = data->qxg[lr][nir-1];
          facr = data->qwgt[lr][nir-1];
          e2   = data->qxg[ls][nis-1];
          facs = data->qwgt[ls][nis-1];
          f2_rec(funct,deriv,NULL,e1,e2,typ,2);
          break;

        case tri3: case tri6:   /* --> tri - element */
          e1   = data->txgr[lr][intc];
          facr = data->twgt[lr][intc];
          e2   = data->txgs[lr][intc];
          facs = ONE;
          f2_tri(funct,deriv,NULL,e1,e2,typ,2);
          break;

        default:
          dserror("typ unknown!");
      } /* end switch(typ) */

      /* compute Jacobian matrix */
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;


      /* get velocities at integration point */
      f2_veci(velint,funct,evelng,iel);


      /* get coordinates at integration point */
      f2_veci(xint,funct,xyze,iel);


      /* get pressure at integration point */
      presint = f2_scali(funct,epren,iel);


      if (fdyn->init == 9 )
      {
        /* calculate analytical solution */
        solx = - cos(a*PI*xint[0]) * sin(a*PI*xint[1])
          * exp(-2.0*a*a*PI*PI*t*visc);
        soly = sin(a*PI*xint[0]) * cos(a*PI*xint[1])
          * exp(-2.0*a*a*PI*PI*t*visc);
        solp = - 0.25 * ( cos(2.0*a*PI*xint[0])
            + cos(2.0*a*PI*xint[1]) )
          * exp(-4.0*a*a*PI*PI*t*visc);
      }
      else
        dserror(
            "Error calculation for fluid2 only possible for kim moin problem");


      /* calculate difference between anlytical and fe solution */
      diffx = (velint[0] - solx);
      diffy = (velint[1] - soly);
      diffp = (presint - solp);


      switch (container->error_norm)
      {
        case 0:
          /* infinity norm */
          container->vel_error = MAX( ABS(diffx), container->vel_error);
          container->vel_error = MAX( ABS(diffy), container->vel_error);
          container->pre_error = MAX( ABS(diffp), container->pre_error);

          container->vel_norm  = MAX( ABS(solx),  container->vel_norm);
          container->vel_norm  = MAX( ABS(soly),  container->vel_norm);
          container->pre_norm  = MAX( ABS(solp),  container->pre_norm);
          break;

        case 1:
          /* L1 norm */
          container->vel_error +=  ( ABS(diffx) + ABS(diffy) ) * fac ;
          container->pre_error +=  ( ABS(diffp) ) * fac ;

          container->vel_norm  +=  ( ABS(solx) + ABS(diffy) ) * fac ;
          container->pre_norm  +=  ( ABS(solp) ) * fac ;
          break;

        case 2:
          /* L2 norm */
          container->vel_error += (diffx*diffx + diffy*diffy) * fac;
          container->pre_error += diffp*diffp * fac;

          container->vel_norm  += (solx*solx + soly*soly) * fac;
          container->pre_norm  += solp*solp * fac;
          break;

        default:
          dserror("Error norm %2i not available!!",container->error_norm);
          break;
      }


    } /* end of loop over integration points ls*/
  } /* end of loop over integration points lr */

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f2_int_kim_moin_err */



#endif


/*! @} (documentation module close)*/


