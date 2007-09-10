/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2_is element

<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
#ifdef D_FLUID2_IS

#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2_is.h"

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

static FLUID_DYNAMIC *fdyn;

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
void f2is_int_usfem(
  ELEMENT         *ele,
  INT             *hasext,
  DOUBLE         **estif,
  DOUBLE          *eforce,
  DOUBLE         **xyze,
  DOUBLE          *funct,
  DOUBLE         **deriv,
  DOUBLE         **deriv2,
  DOUBLE          *pfunct,
  DOUBLE         **pderiv,
  DOUBLE         **pderiv2,
  DOUBLE         **xjm,
  DOUBLE         **derxy,
  DOUBLE         **derxy2,
  DOUBLE         **pderxy,
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
  DOUBLE           estress[3][MAXNOD_F2],
  INT              is_relax
  )
{
  INT       i;          /* a couter                                       */
  INT       iel;        /* number of nodes                                */
  INT       intc=0;       /* "integration case" for tri for further infos
			   see f2_inpele.c and f2_intg.c                 */
  INT       nir=0,nis=0;    /* number of integration nodesin r,s direction    */
  INT       ihoel=0;    /* flag for higher order elements                 */
  INT       icode=2;    /* flag for eveluation of shape functions         */
  INT       lr, ls;     /* counter for integration                        */
  INT       is_ale;
  DOUBLE    fac;        /* total integration vactor                       */
  DOUBLE    facr=0, facs=0; /* integration weights                            */
  DOUBLE    det;        /* determinant of jacobian matrix at time (n+1)   */
  DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
  DOUBLE    gradp[2];   /* pressure gradient at integration point         */
  DOUBLE    velint[2];  /* velocity vector at integration point           */
  DOUBLE    histvec[2]; /* history data at integration point              */
  DOUBLE    gridvelint[2]; /* grid velocity                               */
  DOUBLE    divuold;
  DOUBLE    press;
  DIS_TYP   typ;	/* element type                                   */

  FLUID_DYNAMIC   *fdyn;
  FLUID_DATA      *data;

  INT pdof=0;

#ifdef DEBUG
  dstrc_enter("f2is_int_usfem");
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
  case quad4:
    dserror("quad4 inf-sup not supported");
    break;
  case quad8: case quad9:  /* --> quad - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f2->nGP[0];
    nis = ele->e.f2->nGP[1];

    pdof = 4;
    break;
  case tri6: /* --> tri - element */
    icode   = 3;
    ihoel   = 1;
/* do NOT break at this point!!! */
  case tri3:
    dserror("currently not supported");

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
      /* get values of  shape functions and their derivatives ---*/
      switch(typ)
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
	e1   = data->qxg[lr][nir-1];
	facr = data->qwgt[lr][nir-1];
	e2   = data->qxg[ls][nis-1];
	facs = data->qwgt[ls][nis-1];
	f2_rec(funct,deriv,deriv2,e1,e2,typ,icode);
	f2_rec(pfunct,pderiv,pderiv2,e1,e2,quad4,icode);
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
      f2_gder(pderxy,pderiv,xjm,det,pdof);

      /*---------------- get velocities (n+1,i) at integration point ---*/
      f2_veci(velint,funct,evelng,iel);

      /*---------------- get history data (n,i) at integration point ---*/
      f2_veci(histvec,funct,evhist,iel);

      /*--------------------- get grid velocity at integration point ---*/
      if (is_ale)
      {
	f2_veci(gridvelint,funct,egridv,iel);
      }
      else
      {
        gridvelint[0] = 0;
        gridvelint[1] = 0;
      }

      /*-------- get velocity (n,i) derivatives at integration point ---*/
      f2_vder(vderxy,derxy,eveln,iel);
      divuold = vderxy[0][0] + vderxy[1][1];

      /*------ get velocity (n+1,i) derivatives at integration point ---*/
      f2_vder(vderxy,derxy,evelng,iel);

      /*--------------------------- compute second global derivative ---*/
      if (ihoel!=0)
      {
	f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
	f2_vder2(vderxy2,derxy2,evelng,iel);
      }

      /*------------------------------------- get pressure gradients ---*/
      gradp[0] = gradp[1] = 0.0;

      for (i=0; i<pdof; i++)
      {
	gradp[0] += pderxy[0][i] * epren[i];
	gradp[1] += pderxy[1][i] * epren[i];
      }

      press = f2_scali(pfunct,epren,pdof);

      /*-------------- perform integration for entire matrix and rhs ---*/
      f2is_calmat(estif,eforce,velint,histvec,gridvelint,press,vderxy,
		  vderxy2,gradp,funct,derxy,derxy2,pfunct,pderxy,edeadng,fac,
		  visc,iel,hasext,is_ale, is_relax);
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

/*------------------------------------------- assure assembly of rhs ---*/
  if (!is_relax)
    *hasext = 1;
  return;
}


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
\param **egridv    DOUBLE	   (i)    grid velocity
\param  *epren     DOUBLE	   (-)    ele pres. at time n
\param  *edeadng   DOUBLE	   (-)    ele dead load (selfweight) at n+1
\param **vderxy    DOUBLE	   (-)    global vel. derivatives
\param **vderxy2   DOUBLE	   (-)    2nd global vel. deriv.
\param   visc      DOUBLE          (i)    viscosity
\param **wa1	   DOUBLE	   (-)    working array
\param **wa2	   DOUBLE	   (-)    working array
\return void

------------------------------------------------------------------------*/
void f2is_int_res(
  ELEMENT         *ele,
  INT             *hasext,
  DOUBLE          *force,
  DOUBLE         **xyze,
  DOUBLE          *funct,
  DOUBLE         **deriv,
  DOUBLE         **deriv2,
  DOUBLE          *pfunct,
  DOUBLE         **pderiv,
  DOUBLE         **pderiv2,
  DOUBLE         **xjm,
  DOUBLE         **derxy,
  DOUBLE         **derxy2,
  DOUBLE         **pderxy,
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
  INT       intc=0;       /* "integration case" for tri for further infos
			   see f2_inpele.c and f2_intg.c                 */
  INT       is_ale;     /* ALE or Euler element flag                      */
  INT       nir=0,nis=0;    /* number of integration nodesin r,s direction    */
  INT       ihoel=0;    /* flag for higher order elements                 */
  INT       icode=2;    /* flag for eveluation of shape functions         */
  INT       lr, ls;     /* counter for integration                        */

  DOUBLE    fac;        /* total integration vactor                       */
  DOUBLE    facr=0, facs=0; /* integration weights                            */
  DOUBLE    det;        /* determinant of jacobian matrix                 */
  DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
  DOUBLE    velint[2];  /* velocity vector at integration point           */
  DOUBLE    histvec[2]; /* history data at integration point              */
  DOUBLE    gridvelint[2]; /* grid velocity                               */
  DIS_TYP   typ;	      /* element type                                   */
  DOUBLE    presint;    /* pressure at integration point                  */
  DOUBLE    gradp[2];   /* pressure gradient                              */
  FLUID_DYNAMIC   *fdyn;
  FLUID_DATA      *data;

  INT pdof=0;

#ifdef DEBUG
  dstrc_enter("f2is_int_usfem");
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
  case quad4:
    dserror("quad4 inf-sup not supported");
    break;
  case quad8: case quad9:  /* --> quad - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f2->nGP[0];
    nis = ele->e.f2->nGP[1];
    pdof = 4;
    break;
  case tri6: /* --> tri - element */
    icode   = 3;
    ihoel   = 1;
/* do NOT break at this point!!! */
  case tri3:
    dserror("currently not supported");

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
	f2_rec(pfunct,pderiv,pderiv2,e1,e2,quad4,icode);
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
      f2_gder(pderxy,pderiv,xjm,det,pdof);

      /*------------------------ get velocities at integration point ---*/
      f2_veci(velint,funct,evelng,iel);

      /*---------------------- get history data at integration point ---*/
      f2_veci(histvec,funct,evhist,iel);

      /*--------------------- get grid velocity at integration point ---*/
      if (is_ale)
      {
	f2_veci(gridvelint,funct,egridv,iel);
      }
      else
      {
        gridvelint[0] = 0;
        gridvelint[1] = 0;
      }

      /*-------------- get velocity derivatives at integration point ---*/
      f2_vder(vderxy,derxy,evelng,iel);

      /*--------------------------- compute second global derivative ---*/
      if (ihoel!=0)
      {
	f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,iel);
	f2_vder2(vderxy2,derxy2,evelng,iel);
      }

      /*------------------------------------- get pressure gradients ---*/
      gradp[0] = gradp[1] = 0.0;

      for (i=0; i<pdof; i++)
      {
	gradp[0] += pderxy[0][i] * epren[i];
	gradp[1] += pderxy[1][i] * epren[i];
      }
      /*-------------------------- get pressure at integration point ---*/
      presint = f2_scali(pfunct,epren,pdof);

      /*-------------- perform integration for entire matrix and rhs ---*/
      f2is_calresvec(force,velint,histvec,vderxy,vderxy2,funct,pfunct,derxy,derxy2,
		     edeadng,gridvelint,presint,gradp,fac,visc,iel,hasext,
		     is_ale);
    } /* end of loop over integration points ls*/
  } /* end of loop over integration points lr */

  for (i=0; i<MAXNOD*MAXDOFPERNODE; ++i)
  {
    force[i] /= fdyn->thsl;
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}


/*!---------------------------------------------------------------------
\brief evaluate fluid coefficient matrix

<pre>                                                        chfoe 04/04

In this routine the Gauss point contributions to the elemental coefficient
matrix of a stabilised fluid2 element are calculated. The procedure is
based on the Rothe method of first integrating in time. Hence the
resulting terms include coefficients containing time integration variables
such as theta or delta t which are represented by 'timefac'.

The routine was completed to contain ALE-terms also.         chfoe 11/04

The stabilisation is based on the residuum:

R_M = u + timefac u * grad u - timefac * 2 nu div epsilon(u)
    + timefac grad p - rhsint

R_C = div u

The corresponding weighting operators are
L_M = v + timefac u_old * grad v + timefac v * grad u_old
    - timefac * 2 nu alpha div epsilon (v) + timefac beta grad q

L_C = div v

where alpha = -1
      beta  = -1
are sign regulating factors and rhsint differs for different time
These factores are worked in now and cannot be changed any more.

integration schemes:

One-step-Theta:
rhsint = u_old + Theta dt f + (1-Theta) acc_old

BDF2:

generalised alpha:


The stabilisation by means of the momentum residuum R_M is of the unusual
type:
   Galerkin parts MINUS sum over elements (stabilising parts)
The stabilisation by means of the continuity equation R_C is done in the
usual way:
   Galerkin parts PLUS sum over elements (stabilising parts)

The calculation proceeds as follows.
1) obtain single (linearised) operators of R_M, R_C, L_M and L_C
2) build Galerkin terms from them
3) build stabilising terms from them
4) build Galerkin and stabilising terms of RHS

NOTE: u_old represents the last iteration value. (The most recent one
      we've got!)

NOTE: Galerkin and stabilisation matrices are calculated within one
      routine.

NOTE: In order to increase the performance plenty of terms are concentrated
      and worked into each other. A lengthy version of the file is available
      from the author.


Notational remarks:

                   /              \
                  | u_x,x   u_x,y |
vderxy = grad u = |               |
                  | u_y,x   u_y,y |
                  \               /

           /                         \
          | u_x,xx   u_x,yy   u_x,xy |
vderxy2 = |                          |
          | u_y,xx   u_y,yy   u_y,xy |
          \                          /

for further comments see comment lines within code.


</pre>
\param **estif      DOUBLE        (o)   ele stiffness matrix
\param  *eforce     DOUBLE        (o)   ele force vector
\param  *velint     DOUBLE        (i)   vel at INT point
\param  *histvec    DOUBLE        (i)   rhs at INT point
\param  *gridvint   DOUBLE        (i)   gridvel at INT point
\param **vderxy     DOUBLE        (i)   global vel derivatives
\param  *vderxy2    DOUBLE        (i)   2nd global vel derivatives
\param  *funct      DOUBLE        (i)   nat. shape funcs
\param **derxy      DOUBLE        (i)   global coord. deriv.
\param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
\param  *edeadng    DOUBLE        (i)   dead load at time n+1
\param   fac 	    DOUBLE        (i)   weighting factor
\param   visc       DOUBLE        (i)   fluid viscosity
\param   iel	    INT           (i)   number of nodes of act. ele
\param  *hasext     INT           (i)   flag, if element has volume load
\param   isale      INT           (i)   flag, if ALE or EULER
\return void

------------------------------------------------------------------------*/
void f2is_calmat( DOUBLE **estif,
		DOUBLE  *eforce,
		DOUBLE  *velint,
		DOUBLE   histvec[2],
		DOUBLE   gridvint[2],
		DOUBLE   press,
		DOUBLE **vderxy,
                DOUBLE **vderxy2,
                DOUBLE   gradp[2],
		DOUBLE  *funct,
		DOUBLE **derxy,
		DOUBLE **derxy2,
		  DOUBLE          *pfunct,
		  DOUBLE         **pderxy,
                DOUBLE  *edeadng,
		DOUBLE   fac,
		DOUBLE   visc,
		INT      iel,
                INT     *hasext,
                INT      isale,
                INT      is_relax
  )
{
  INT     i, j, ri, ci;
  DOUBLE  timefac;    /* One-step-Theta: timefac = theta*dt
			 BDF2:           timefac = 2/3 * dt               */
  DOUBLE  dt;         /* time step size*/
/*   DOUBLE  aux; */
/*   DOUBLE  auxmat[2][2]; */
  DOUBLE  tau_M, tau_C;             /* stabilisation parameter            */
  DOUBLE  tau_Mp;             /* stabilisation parameter            */
  DOUBLE  viscs2[2][2*MAXNOD]; /* viscous term incluiding 2nd derivatives */
#if 0
  DOUBLE  viscous[2][2][2*MAXNOD];   /* viscous term partially integrated */
#endif
  DOUBLE  conv_c[MAXNOD];    /* linearisation of convect, convective part */
  DOUBLE  conv_g[MAXNOD];          /* linearisation of convect, grid part */
  DOUBLE  conv_r[2][2*MAXNOD]; /* linearisation of convect, reactive part */
  DOUBLE  div[2*MAXNOD];             /* divergence of u or v              */
  DOUBLE  ugradv[MAXNOD][2*MAXNOD];  /* linearisation of u * grad v       */
  DOUBLE  conv_old[2]; /* convective term evalaluated with old velocities */
  DOUBLE  visc_old[2]; /* viscous term evaluated with old velocities      */
  DOUBLE  rhsint[2];   /* total right hand side terms at int.-point       */

  DOUBLE  time2nue, timetauM, timetauMp, ttimetauM, ttimetauMp, timefacfac;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("f2is_calmat");
#endif
/*========================== initialisation ============================*/
  fdyn = alldyn[genprob.numff].fdyn;

  tau_M  = fdyn->tau[0]*fac;
  tau_Mp = fdyn->tau[0]*fac;
  tau_C  = fdyn->tau[2]*fac;

  timefac = fdyn->thsl;
  dt      = fdyn->dta;

/* integration factors and koefficients of single terms */
  time2nue  = timefac * 2.0 * visc;
  timetauM   = timefac * tau_M;
  timetauMp  = timefac * tau_Mp;

  ttimetauM  = timefac * timetauM;
  ttimetauMp = timefac * timetauMp;
  timefacfac = timefac * fac;


/*------------------------- evaluate rhs vector at integration point ---*/
  if (!is_relax)
  {
    if (*hasext)
    {
      rhsint[0] = timefac * edeadng[0] + histvec[0];
      rhsint[1] = timefac * edeadng[1] + histvec[1];
    }
    else
    {
      rhsint[0] = histvec[0];
      rhsint[1] = histvec[1];
    }
  }
  else
  {
    rhsint[0]=0;
    rhsint[1]=0;
  }

/*----------------- get numerical representation of single operators ---*/

/* Convective term  u_old * grad u_old: */
  conv_old[0] = vderxy[0][0] * velint[0] + vderxy[0][1] * velint[1];
  conv_old[1] = vderxy[1][0] * velint[0] + vderxy[1][1] * velint[1];

/* Viscous term  div epsilon(u_old) */
  visc_old[0] = 0.5 * (2.0*vderxy2[0][0] + vderxy2[0][1] + vderxy2[1][2]);
  visc_old[1] = 0.5 * (2.0*vderxy2[1][1] + vderxy2[1][0] + vderxy2[0][2]);

  for (i=0; i<iel; i++) /* loop over nodes of element */
  {
    /* Reactive term  u:  funct */
    /* linearise convective term */

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
    conv_c[i] = derxy[0][i] * velint[0] + derxy[1][i] * velint[1] ;

    /*--- convective grid part u_G * grad (funct) -----------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
    if(isale)
    {
      conv_g[i] = - derxy[0][i] * gridvint[0] - derxy[1][i] * gridvint[1];
    }
    else
    {
      conv_g[i] = 0;
    }

    /*--- reactive part funct * grad (u_old) ----------------------------*/
    /* /                          \
       |  u_old_x,x   u_old_x,y   |
       |                          | * N   with  N .. form function matrix
       |  u_old_y,x   u_old_y,y   |
       \                         /                                       */
    conv_r[0][2*i]   = vderxy[0][0]*funct[i];
    conv_r[0][2*i+1] = vderxy[0][1]*funct[i];
    conv_r[1][2*i]   = vderxy[1][0]*funct[i];
    conv_r[1][2*i+1] = vderxy[1][1]*funct[i];

    /*--- viscous term  - grad * epsilon(u): ----------------------------*/
    /*   /                              \
	 1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
	 - - |                              |         N_y .. y-line of N
	 2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
	 \                             /                                 */
    viscs2[0][2*i]   = - 0.5 * ( 2.0 * derxy2[0][i] + derxy2[1][i] );
    viscs2[0][2*i+1] = - 0.5 * ( derxy2[2][i] );
    viscs2[1][2*i]   = - 0.5 * ( derxy2[2][i] );
    viscs2[1][2*i+1] = - 0.5 * ( derxy2[0][i] + 2.0 * derxy2[1][i] );

#if 0
    /*--- viscous term (after integr. by parts) -------------------------*/
    /*   /                            \
	 1 |  2 N_x,x    N_x,y + N_y,x  |    with N_x .. x-line of N
	 - |                            |         N_y .. y-line of N
	 2 |  N_y,x + N_x,y    2 N_y,y  |
	 \                           /                                   */
    viscous[0][0][2*i]   = derxy[0][i];
    viscous[0][0][2*i+1] = 0.0;                /* 1st index:             */
    viscous[0][1][2*i]   = 0.5 * derxy[1][i];  /*   line of epsilon      */
    viscous[0][1][2*i+1] = 0.5 * derxy[0][i];  /* 2nd index:             */
    viscous[1][0][2*i]   = 0.5 * derxy[1][i];  /*   column of epsilon    */
    viscous[1][0][2*i+1] = 0.5 * derxy[0][i];  /* 3rd index:             */
    viscous[1][1][2*i+1] = derxy[1][i];        /*   elemental vel dof    */
    viscous[1][1][2*i]   = 0.0;
#endif

    /* pressure gradient term derxy, funct without or with integration   *
     * by parts, respectively                                            */

    /*--- divergence u term ---------------------------------------------*/
    div[2*i]   = derxy[0][i];
    div[2*i+1] = derxy[1][i];

    /*--- ugradv-Term ---------------------------------------------------*/
    /*
      /                                                          \
      |  N1*N1,x  N1*N1,y  N2*N1,x  N2*N1,y  N3*N1,x ...       . |
      |                                                          |
      |  N1*N2,x  N1*N2,y  N2*N2,x  N2*N2,y  N3*N2,x ...       . |
      |                                                          |
      |  N1*N3,x  N1*N3,y  N2*N3,x  N2*N3,y  N3*N3,x ...       . |
      |                                           .              |
      |  . . .                                        .          |
      |                                                  Ni*Ni,y |
      \                                                          /       */
    /* remark: vgradu = ugradv^T */
    for (j=0; j<iel; j++)
    {
      ugradv[i][2*j]   = derxy[0][i] * funct[j];
      ugradv[i][2*j+1] = derxy[1][i] * funct[j];
    }
  }

/*--------------------------------- now build single stiffness terms ---*/

  /* An easy way to switch off terms. */

#define FLUID2_IS_TERM1
#define FLUID2_IS_TERM2
#define FLUID2_IS_TERM3
#define FLUID2_IS_TERM4
#define FLUID2_IS_TERM5
#define FLUID2_IS_TERM6
#define FLUID2_IS_TERM7
#define FLUID2_IS_TERM8
#define FLUID2_IS_TERM9
#define FLUID2_IS_TERM10
#define FLUID2_IS_TERM11
#define FLUID2_IS_TERM12
#define FLUID2_IS_TERM13
#define FLUID2_IS_TERM14
#define FLUID2_IS_TERM15
#define FLUID2_IS_TERM16
#define FLUID2_IS_TERM17


#define estif_(i,j)    estif[i][j]
#define eforce_(i)     eforce[i]
#define funct_(i)      funct[i]
#define funct_p_(i)    pfunct[i]
#define vderxy_(i,j)   vderxy[i][j]
#define conv_c_(j)     conv_c[j]
#define conv_g_(j)     conv_g[j]
#define conv_r_(i,j,k) conv_r[i][2*(k)+j]
#define conv_old_(j)   conv_old[j]
#define derxy_(i,j)    derxy[i][j]
#define pderxy_(i,j)   pderxy[i][j]
#define gridvint_(j)   gridvint[j]
#define velint_(j)     velint[j]
#define viscs2_(i,j,k) viscs2[i][2*(k)+j]
#define visc_old_(i)   visc_old[i]
#define rhsint_(i)     rhsint[i]
#define gradp_(j)      gradp[j]
#define ui             ci
#define vi             ri
#define nu_            visc
#define thsl           timefac

  /* This code is generated using MuPAD. Ask me for the MuPAD. u.kue */

#ifdef FLUID_INCREMENTAL

#ifdef QUASI_NEWTON
  if (!fdyn->qnewton || fdyn->itnum==1)
  {
#endif

    /* We keep two versions: with and without ale. The laster one is a
     * little faster. */

    if (isale)
    {
#include "f2is_stiff_ale.c"
    }
    else
    {
#include "f2is_stiff.c"
    }

#ifdef QUASI_NEWTON
  }
#endif

  if (isale)
  {
#include "f2is_rhs_incr_ale.c"
  }
  else
  {
#include "f2is_rhs_incr.c"
  }

#else  /* FLUID_INCREMENTAL */

#ifdef QUASI_NEWTON
  dserror("quasi newton not possible with non-incremental elements");
#endif

  if (isale)
  {
#include "f2is_stiff_ale.c"
#include "f2is_rhs_nonincr_ale.c"
  }
  else
  {
#include "f2is_stiff.c"
#include "f2is_rhs_nonincr.c"
  }

#endif	/* FLUID_INCREMENTAL */

#undef estif_
#undef eforce_
#undef conv_c_
#undef conv_g_
#undef conv_r_
#undef conv_old_
#undef derxy_
#undef pderxy_
#undef gridvint_
#undef velint_
#undef viscs2_
#undef gradp_
#undef ui
#undef vi
#undef funct_
#undef funct_p_
#undef vderxy_
#undef visc_old_
#undef rhsint_
#undef nu_
#undef thsl

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
}


/*!---------------------------------------------------------------------
\brief Gauss point contributions for integration of boundary forces

<pre>                                                         chfoe 03/05

This routine evaluates the Gauss point vaulues of the residual vector
of one element taking stabilisation effects into account. Only the
residual of the momentum equation R_M is calculated.

R_M = u + timefac u * grad u - timefac * 2 nu div epsilon(u)
    + timefac grad p - rhsint

The residual contains stabilisation of the type

Sum_over_k (R_M, tau L_M)_k with

L_M = v + timefac u_old * grad v + timefac v * grad u_old
    - timefac * 2 nu alpha div epsilon (v) + timefac beta grad q

where alpha = -1
      beta  = -1

timefac depends on the time integration scheme:

One-step theta:

timefac = theta * dt

BDF2:

timefac = 2/3 * dt

NOTE: this works perfectly only when the fluid is solved via usfem

</pre>
\param  *eforce	   DOUBLE    (o)    element force vector (residual)
\param  *velint    DOUBLE    (i)    (converged) vel. at int.-point
\param   histvec   DOUBLE    (i)    histroy data
\param **vderxy    DOUBLE    (i)    velocity gradient at int.-point
\param **vderxy2   DOUBLE    (i)    second vel. derivatives at int.-point
\param  *funct     DOUBLE    (i)    natural shape functions
\param **derxy     DOUBLE    (i)    shape function derivatives
\param **derxy2    DOUBLE    (i)    second shape funct. derivs
\param  *edeadng   DOUBLE    (i)    body forces
\param  *press     DOUBLE    (i)    pressure at Gauss point
\param   gradp[2]  DOUBLE    (i)    pressure gradient at GP
\param   fac       DOUBLE    (i)    integration factor
\param   visc      DOUBLE    (i)    fluid viscosity
\param   iel       INT       (i)    number of elemental nodes
\param  *hasext    INT       (i)    flag, if there is body force
\param   is_ale    INT       (i)    flag, if it's ale or Euler
\return void

------------------------------------------------------------------------*/
void f2is_calresvec(
  DOUBLE  *eforce,
  DOUBLE  *velint,
  DOUBLE   histvec[2],
  DOUBLE **vderxy,
  DOUBLE **vderxy2,
  DOUBLE  *funct,
  DOUBLE  *pfunct,
  DOUBLE **derxy,
  DOUBLE **derxy2,
  DOUBLE  *edeadng,
  DOUBLE   gridvint[2],
  DOUBLE   press,
  DOUBLE   gradp[2],
  DOUBLE   fac,
  DOUBLE   visc,
  INT      iel,
  INT     *hasext,
  INT      is_ale
  )
{
  INT     i, j, ri;
  DOUBLE  timefac;    /* One-step-Theta: timefac = theta*dt
			 BDF2:           timefac = 2/3 * dt               */
  DOUBLE  dt;         /* time step size*/
/*   DOUBLE  aux; */
/*   DOUBLE  auxmat[2][2]; */
  DOUBLE  tau_M, tau_C;             /* stabilisation parameter            */
  DOUBLE  tau_Mp;             /* stabilisation parameter            */
  DOUBLE  viscs2[2][2*MAXNOD]; /* viscous term incluiding 2nd derivatives */
  DOUBLE  viscous[2][2][2*MAXNOD];   /* viscous term partially integrated */
  DOUBLE  conv_c[MAXNOD];    /* linearisation of convect, convective part */
  DOUBLE  conv_g[MAXNOD];          /* linearisation of convect, grid part */
  DOUBLE  conv_r[2][2*MAXNOD]; /* linearisation of convect, reactive part */
  DOUBLE  div[2*MAXNOD];             /* divergence of u or v              */
  DOUBLE  ugradv[MAXNOD][2*MAXNOD];  /* linearisation of u * grad v       */
  DOUBLE  conv_old[2]; /* convective term evalaluated with old velocities */
  DOUBLE  visc_old[2]; /* viscous term evaluated with old velocities      */
  DOUBLE  rhsint[2];   /* total right hand side terms at int.-point       */

  DOUBLE  time2nue, timetauM, timetauMp, ttimetauM, ttimetauMp, timefacfac;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("f2is_calresvec");
#endif
/*========================== initialisation ============================*/
  fdyn = alldyn[genprob.numff].fdyn;

  tau_M  = fdyn->tau[0]*fac;
  tau_Mp = fdyn->tau[0]*fac;
  tau_C  = fdyn->tau[2]*fac;

  timefac = fdyn->thsl;
  dt      = fdyn->dta;

/* integration factors and koefficients of single terms */
  time2nue  = timefac * 2.0 * visc;
  timetauM   = timefac * tau_M;
  timetauMp  = timefac * tau_Mp;

  ttimetauM  = timefac * timetauM;
  ttimetauMp = timefac * timetauMp;
  timefacfac = timefac * fac;


/*------------------------- evaluate rhs vector at integration point ---*/
  if (*hasext)
  {
    rhsint[0] = timefac * edeadng[0] + histvec[0];
    rhsint[1] = timefac * edeadng[1] + histvec[1];
  }
  else
  {
    rhsint[0] = histvec[0];
    rhsint[1] = histvec[1];
  }

/*----------------- get numerical representation of single operators ---*/

/* Convective term  u_old * grad u_old: */
  conv_old[0] = vderxy[0][0] * velint[0] + vderxy[0][1] * velint[1];
  conv_old[1] = vderxy[1][0] * velint[0] + vderxy[1][1] * velint[1];

/* Viscous term  div epsilon(u_old) */
  visc_old[0] = 0.5 * (2.0*vderxy2[0][0] + vderxy2[0][1] + vderxy2[1][2]);
  visc_old[1] = 0.5 * (2.0*vderxy2[1][1] + vderxy2[1][0] + vderxy2[0][2]);

  for (i=0; i<iel; i++) /* loop over nodes of element */
  {
    /* Reactive term  u:  funct */
    /* linearise convective term */

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
    conv_c[i] = derxy[0][i] * velint[0] + derxy[1][i] * velint[1] ;

    /*--- convective grid part u_G * grad (funct) -----------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
    if(is_ale)
    {
      conv_g[i] = - derxy[0][i] * gridvint[0] - derxy[1][i] * gridvint[1];
    }
    else
    {
      conv_g[i] = 0;
    }

    /*--- reactive part funct * grad (u_old) ----------------------------*/
    /* /                          \
       |  u_old_x,x   u_old_x,y   |
       |                          | * N   with  N .. form function matrix
       |  u_old_y,x   u_old_y,y   |
       \                         /                                       */
    conv_r[0][2*i]   = vderxy[0][0]*funct[i];
    conv_r[0][2*i+1] = vderxy[0][1]*funct[i];
    conv_r[1][2*i]   = vderxy[1][0]*funct[i];
    conv_r[1][2*i+1] = vderxy[1][1]*funct[i];

    /*--- viscous term  - grad * epsilon(u): ----------------------------*/
    /*   /                              \
	 1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
	 - - |                              |         N_y .. y-line of N
	 2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
	 \                             /                                 */
    viscs2[0][2*i]   = - 0.5 * ( 2.0 * derxy2[0][i] + derxy2[1][i] );
    viscs2[0][2*i+1] = - 0.5 * ( derxy2[2][i] );
    viscs2[1][2*i]   = - 0.5 * ( derxy2[2][i] );
    viscs2[1][2*i+1] = - 0.5 * ( derxy2[0][i] + 2.0 * derxy2[1][i] );

    /*--- viscous term (after integr. by parts) -------------------------*/
    /*   /                            \
	 1 |  2 N_x,x    N_x,y + N_y,x  |    with N_x .. x-line of N
	 - |                            |         N_y .. y-line of N
	 2 |  N_y,x + N_x,y    2 N_y,y  |
	 \                           /                                   */
    viscous[0][0][2*i]   = derxy[0][i];
    viscous[0][0][2*i+1] = 0.0;                /* 1st index:             */
    viscous[0][1][2*i]   = 0.5 * derxy[1][i];  /*   line of epsilon      */
    viscous[0][1][2*i+1] = 0.5 * derxy[0][i];  /* 2nd index:             */
    viscous[1][0][2*i]   = 0.5 * derxy[1][i];  /*   column of epsilon    */
    viscous[1][0][2*i+1] = 0.5 * derxy[0][i];  /* 3rd index:             */
    viscous[1][1][2*i+1] = derxy[1][i];        /*   elemental vel dof    */
    viscous[1][1][2*i]   = 0.0;

    /* pressure gradient term derxy, funct without or with integration   *
     * by parts, respectively                                            */

    /*--- divergence u term ---------------------------------------------*/
    div[2*i]   = derxy[0][i];
    div[2*i+1] = derxy[1][i];

    /*--- ugradv-Term ---------------------------------------------------*/
    /*
      /                                                          \
      |  N1*N1,x  N1*N1,y  N2*N1,x  N2*N1,y  N3*N1,x ...       . |
      |                                                          |
      |  N1*N2,x  N1*N2,y  N2*N2,x  N2*N2,y  N3*N2,x ...       . |
      |                                                          |
      |  N1*N3,x  N1*N3,y  N2*N3,x  N2*N3,y  N3*N3,x ...       . |
      |                                           .              |
      |  . . .                                        .          |
      |                                                  Ni*Ni,y |
      \                                                          /       */
    /* remark: vgradu = ugradv^T */
    for (j=0; j<iel; j++)
    {
      ugradv[i][2*j]   = derxy[0][i] * funct[j];
      ugradv[i][2*j+1] = derxy[1][i] * funct[j];
    }

  }

/*--------------------------------- now build single stiffness terms ---*/

  /* An easy way to switch off terms. */

#define FLUID2_IS_TERM1
#define FLUID2_IS_TERM2
#define FLUID2_IS_TERM3
#define FLUID2_IS_TERM4
#define FLUID2_IS_TERM5
#define FLUID2_IS_TERM6
#define FLUID2_IS_TERM7
#define FLUID2_IS_TERM8
#define FLUID2_IS_TERM9
#define FLUID2_IS_TERM10
#define FLUID2_IS_TERM11
#define FLUID2_IS_TERM12
#define FLUID2_IS_TERM13
#define FLUID2_IS_TERM14
#define FLUID2_IS_TERM15
#define FLUID2_IS_TERM16
#define FLUID2_IS_TERM17


#define estif_(i,j)    estif[i][j]
#define eforce_(i)     eforce[i]
#define funct_(i)      funct[i]
#define funct_p_(i)    pfunct[i]
#define vderxy_(i,j)   vderxy[i][j]
#define conv_c_(j)     conv_c[j]
#define conv_g_(j)     conv_g[j]
#define conv_r_(i,j,k) conv_r[i][2*(k)+j]
#define conv_old_(j)   conv_old[j]
#define derxy_(i,j)    derxy[i][j]
#define pderxy_(i,j)   pderxy[i][j]
#define gridvint_(j)   gridvint[j]
#define velint_(j)     velint[j]
#define viscs2_(i,j,k) viscs2[i][2*(k)+j]
#define visc_old_(i)   visc_old[i]
#define rhsint_(i)     rhsint[i]
#define gradp_(j)      gradp[j]
#define ui             ci
#define vi             ri
#define nu_            visc
#define thsl           timefac

  /* This code is generated using MuPAD. Ask me for the MuPAD. u.kue */

  if (is_ale)
  {
#include "f2is_rhs_incr_ale.c"
  }
  else
  {
#include "f2is_rhs_incr.c"
  }

#undef estif_
#undef eforce_
#undef conv_c_
#undef conv_g_
#undef conv_r_
#undef conv_old_
#undef derxy_
#undef pderxy_
#undef gridvint_
#undef velint_
#undef viscs2_
#undef gradp_
#undef ui
#undef vi
#undef funct_
#undef funct_p_
#undef vderxy_
#undef visc_old_
#undef rhsint_
#undef nu_
#undef thsl

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif
#endif
