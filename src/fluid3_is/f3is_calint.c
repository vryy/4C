/*!----------------------------------------------------------------------
\file
\brief integration loop for one fluid2 element using USFEM

<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3
#include "../headers/standardtypes.h"
#include "fluid3_is.h"
#include "../ale3/ale3.h"
#include "../fluid3/fluid3_prototypes.h"

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
void f3is_int_usfem(
  ELEMENT         *ele,
  INT             *hasext,
  DOUBLE         **estif,
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
  DOUBLE         **wa2,
  INT              is_relax
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
  DIS_TYP   typ;        /* element type                                   */
  DOUBLE    press;

  FLUID_DYNAMIC   *fdyn;
  FLUID_DATA      *data;

  INT pdof=0;

#ifdef DEBUG
  dstrc_enter("f3is_int_usfem");
#endif

/*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  is_ale = ele->e.f3is->is_ale;
  intc = 0;
  nir  = 0;
  nis  = 0;
  nit  = 0;


  switch (typ)
  {
  case hex8:
    dserror("hex8 inf-sup not supported");
    break;
  case hex20: case hex27:  /* --> hex - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f3is->nGP[0];
    nis = ele->e.f3is->nGP[1];
    nit = ele->e.f3is->nGP[2];
    intc= 0;

    pdof = 8;
    break;
  case tet10: /* --> tet - element */
    icode   = 3;
    ihoel   = 1;
/* do NOT break at this point!!! */
  case tet4:    /* initialise integration */
    dserror("currently not supported");

    nir  = ele->e.f3is->nGP[0];
    nis  = 1;
    nit  = 1;
    intc = ele->e.f3is->nGP[1];
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
	  f3_hex(pfunct,pderiv,pderiv2,e1,e2,e3,hex8,icode);
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
	f3_gder(pderxy,pderiv,xjm,wa1,det,pdof);

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
	if(is_ale)
	  f3_veci(gridvelint,funct,egridv,iel);
	else
	{
	  gridvelint[0] = 0;
	  gridvelint[1] = 0;
	  gridvelint[2] = 0;
	}

	/*------------------------------------- get pressure gradients ---*/
	gradp[0] = gradp[1] = gradp[2] = 0.0;

	for (i=0; i<pdof; i++)
	{
	  gradp[0] += pderxy[0][i] * epren[i];
	  gradp[1] += pderxy[1][i] * epren[i];
	  gradp[2] += pderxy[2][i] * epren[i];
	}

	press = f3_scali(pfunct,epren,pdof);

	/*-------------- perform integration for entire matrix and rhs ---*/
	f3is_calmat(estif,force,velint,histvec,gridvelint,press,vderxy,
		    vderxy2,gradp,funct,derxy,derxy2,pfunct,edeadng,fac,
		    visc,iel,hasext,is_ale,is_relax);
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
\param **derxy      DOUBLE        (i)   global coord. deriv
\param **derxy2     DOUBLE        (i)   2nd global coord. deriv.
\param  *edeadng    DOUBLE        (i)   dead load at time n+1
\param   fac        DOUBLE        (i)   weighting factor
\param   visc       DOUBLE        (i)   fluid viscosity
\param   iel        INT           (i)   number of nodes of act. ele
\param  *hasext     INT           (i)   flag, if element has volume load
\param   isale      INT           (i)   flag, if ALE or EULER
\return void

------------------------------------------------------------------------*/
void f3is_calmat( DOUBLE **estif,
		  DOUBLE  *eforce,
		  DOUBLE  *velint,
		  DOUBLE   histvec[3],
		  DOUBLE   gridvint[3],
		  DOUBLE   press,
		  DOUBLE **vderxy,
		  DOUBLE **vderxy2,
		  DOUBLE   gradp[3],
		  DOUBLE  *funct,
		  DOUBLE **derxy,
		  DOUBLE **derxy2,
		  DOUBLE  *pfunct,
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
/*   DOUBLE  auxmat[3][3]; */
  DOUBLE  tau_M, tau_C;             /* stabilisation parameter            */
  DOUBLE  tau_Mp;                   /* stabilisation parameter            */
  DOUBLE  viscs2[3][3*MAXNOD_F3];   /* viscous term incluiding 2nd derivatives */
#if 0
  DOUBLE  viscous[3][3][3*MAXNOD_F3];/* viscous term partially integrated */
#endif
  DOUBLE  conv_c[MAXNOD_F3]; /* linearisation of convect, convective part */
  DOUBLE  conv_g[MAXNOD_F3];       /* linearisation of convect, grid part */
  DOUBLE  conv_r[3][3*MAXNOD_F3];/* linearisation of convect, reactive part */
  DOUBLE  div[3*MAXNOD_F3];          /* divergence of u or v              */
  DOUBLE  ugradv[MAXNOD_F3][3*MAXNOD_F3];/* linearisation of u * grad v   */
  DOUBLE  conv_old[3]; /* convective term evalaluated with old velocities */
  DOUBLE  conv_g_old[3];
  DOUBLE  visc_old[3]; /* viscous term evaluated with old velocities      */
  DOUBLE  rhsint[3];   /* total right hand side terms at int.-point       */
  DOUBLE  vconv_r[3][MAXNOD_F3];

  DOUBLE  time2nue, timetauM, timetauMp, ttimetauM, ttimetauMp, timefacfac;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("f3is_calmat");
#endif

/*========================== initialisation ============================*/
  fdyn = alldyn[genprob.numff].fdyn;

  tau_M  = fdyn->tau[0]*fac;
  tau_Mp = fdyn->tau[1]*fac;
  tau_C  = fdyn->tau[2]*fac;

  timefac = fdyn->thsl;
  dt      = fdyn->dta;

/* integration factors and koefficients of single terms */
  time2nue   = timefac * 2.0 * visc;
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
      rhsint[2] = timefac * edeadng[2] + histvec[2];
    }
    else
    {
      rhsint[0] = histvec[0];
      rhsint[1] = histvec[1];
      rhsint[2] = histvec[2];
    }
  }
  else
  {
    rhsint[0] = 0;
    rhsint[1] = 0;
    rhsint[2] = 0;
  }

/*----------------- get numerical representation of single operators ---*/

/* Convective term  u_old * grad u_old: */
  conv_old[0] = (vderxy[0][0] * velint[0] +
		 vderxy[0][1] * velint[1] +
		 vderxy[0][2] * velint[2]);
  conv_old[1] = (vderxy[1][0] * velint[0] +
		 vderxy[1][1] * velint[1] +
		 vderxy[1][2] * velint[2]);
  conv_old[2] = (vderxy[2][0] * velint[0] +
		 vderxy[2][1] * velint[1] +
		 vderxy[2][2] * velint[2]);

/* new for incremental formulation: */
/* Convective term  u_G_old * grad u_old: */
  conv_g_old[0] = (vderxy[0][0] * gridvint[0] +
		   vderxy[0][1] * gridvint[1] +
		   vderxy[0][2] * gridvint[2]);
  conv_g_old[1] = (vderxy[1][0] * gridvint[0] +
		   vderxy[1][1] * gridvint[1] +
		   vderxy[1][2] * gridvint[2]);
  conv_g_old[2] = (vderxy[2][0] * gridvint[0] +
		   vderxy[2][1] * gridvint[1] +
		   vderxy[2][2] * gridvint[2]);

/* Viscous term  div epsilon(u_old) */
  visc_old[0] = vderxy2[0][0] + 0.5 * (vderxy2[0][1] + vderxy2[1][3] + vderxy2[0][2] + vderxy2[2][4]);
  visc_old[1] = vderxy2[1][1] + 0.5 * (vderxy2[1][0] + vderxy2[0][3] + vderxy2[1][2] + vderxy2[2][5]);
  visc_old[2] = vderxy2[2][2] + 0.5 * (vderxy2[2][0] + vderxy2[0][4] + vderxy2[2][1] + vderxy2[1][5]);

  for (i=0; i<iel; i++) /* loop over nodes of element */
  {
    /* Reactive term  u:  funct */
    /* linearise convective term */

    /*--- convective part u_old * grad (funct) --------------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y + u_old_z * N,z
       with  N .. form function matrix                                   */
    conv_c[i] = derxy[0][i] * velint[0] + derxy[1][i] * velint[1]
      + derxy[2][i] * velint[2];

    /*--- convective grid part u_G * grad (funct) -----------------------*/
    /* u_old_x * N,x  +  u_old_y * N,y   with  N .. form function matrix */
    if(isale)
    {
      conv_g[i] = - derxy[0][i] * gridvint[0] - derxy[1][i] * gridvint[1]
	- derxy[2][i] * gridvint[2];
    }
    else
    {
      conv_g[i] = 0;
    }

    /*--- reactive part funct * grad (u_old) ----------------------------*/
    /* /                                     \
       |  u_old_x,x   u_old_x,y   u_old x,z  |
       |                                     |
       |  u_old_y,x   u_old_y,y   u_old_y,z  | * N
       |                                     |
       |  u_old_z,x   u_old_z,y   u_old_z,z  |
       \                                     /
       with  N .. form function matrix                                   */

    conv_r[0][3*i]   = vderxy[0][0]*funct[i];
    conv_r[0][3*i+1] = vderxy[0][1]*funct[i];
    conv_r[0][3*i+2] = vderxy[0][2]*funct[i];
    conv_r[1][3*i]   = vderxy[1][0]*funct[i];
    conv_r[1][3*i+1] = vderxy[1][1]*funct[i];
    conv_r[1][3*i+2] = vderxy[1][2]*funct[i];
    conv_r[2][3*i]   = vderxy[2][0]*funct[i];
    conv_r[2][3*i+1] = vderxy[2][1]*funct[i];
    conv_r[2][3*i+2] = vderxy[2][2]*funct[i];

    vconv_r[0][i] = conv_r[0][3*i]*velint[0] + conv_r[0][3*i+1]*velint[1] + conv_r[0][3*i+2]*velint[2];
    vconv_r[1][i] = conv_r[1][3*i]*velint[0] + conv_r[1][3*i+1]*velint[1] + conv_r[1][3*i+2]*velint[2];
    vconv_r[2][i] = conv_r[2][3*i]*velint[0] + conv_r[2][3*i+1]*velint[1] + conv_r[2][3*i+2]*velint[2];

    /*--- viscous term  - grad * epsilon(u): ----------------------------*/
    /*   /                                                \
	 |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
       1 |                                                |
     - - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
       2 |                                                |
	 |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
	 \                                                /

	 with N_x .. x-line of N
         N_y .. y-line of N                                             */

    viscs2[0][3*i]   = - 0.5 * (2.0 * derxy2[0][i] + derxy2[1][i] + derxy2[2][i]);
    viscs2[0][3*i+1] = - 0.5 *  derxy2[3][i];
    viscs2[0][3*i+2] = - 0.5 *  derxy2[4][i];
    viscs2[1][3*i]   = - 0.5 *  derxy2[3][i];
    viscs2[1][3*i+1] = - 0.5 * (derxy2[0][i] + 2.0 * derxy2[1][i] + derxy2[2][i]);
    viscs2[1][3*i+2] = - 0.5 *  derxy2[5][i];
    viscs2[2][3*i]   = - 0.5 *  derxy2[4][i];
    viscs2[2][3*i+1] = - 0.5 *  derxy2[5][i];
    viscs2[2][3*i+2] = - 0.5 * (derxy2[0][i] + derxy2[1][i] + 2.0 * derxy2[2][i]);

#if 0

    /*--- viscous term (after integr. by parts) -------------------------*/
    /*   /                                             \
	 |  2 N_x,x    N_x,y + N_y,x    N_x,z + N_z,x  |
	 1 |                                             |
	 - |  N_y,x + N_x,y   2 N_y,y     N_y,z + N_z,y  |
	 2 |                                             |
	 |  N_z,x + N_x,z   N_z,y + N_y,z    2 N_z,z   |
	 \                                             /
	 with N_x .. x-line of N
	 N_y .. y-line of N
	 N_z .. z-line of N                                              */

    viscous[0][0][3*i]   = derxy[0][i];
    viscous[0][0][3*i+1] = 0.0;
    viscous[0][0][3*i+2] = 0.0;                /* 1st index:             */
    viscous[0][1][3*i]   = 0.5 * derxy[1][i];  /*   line of epsilon      */
    viscous[0][1][3*i+1] = 0.5 * derxy[0][i];  /* 2nd index:             */
    viscous[0][1][3*i+2] = 0.0;                /*   column of epsilon    */
    viscous[0][2][3*i]   = 0.5 * derxy[2][i];  /* 3rd index:             */
    viscous[0][2][3*i+1] = 0.0;                /*   elemental vel dof    */
    viscous[0][2][3*i+2] = 0.5 * derxy[0][i];
    viscous[1][0][3*i]   = 0.5 * derxy[1][i];
    viscous[1][0][3*i+1] = 0.5 * derxy[0][i];
    viscous[1][0][3*i+2] = 0.0;
    viscous[1][1][3*i]   = 0.0;
    viscous[1][1][3*i+1] = derxy[1][i];
    viscous[1][1][3*i+2] = 0.0;
    viscous[1][2][3*i]   = 0.0;
    viscous[1][2][3*i+1] = 0.5 * derxy[2][i];
    viscous[1][2][3*i+2] = 0.5 * derxy[1][i];
    viscous[2][0][3*i]   = 0.5 * derxy[2][i];
    viscous[2][0][3*i+1] = 0.0;
    viscous[2][0][3*i+2] = 0.5 * derxy[0][i];
    viscous[2][1][3*i]   = 0.0;
    viscous[2][1][3*i+1] = 0.5 * derxy[2][i];
    viscous[2][1][3*i+2] = 0.5 * derxy[1][i];
    viscous[2][2][3*i]   = 0.0;
    viscous[2][2][3*i+1] = 0.0;
    viscous[2][2][3*i+2] = derxy[2][i];
#endif

    /* pressure gradient term derxy, funct without or with integration   *
     * by parts, respectively                                            */

    /*--- divergence u term ---------------------------------------------*/
    div[3*i]   = derxy[0][i];
    div[3*i+1] = derxy[1][i];
    div[3*i+2] = derxy[2][i];

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
      ugradv[i][3*j]   = derxy[0][i] * funct[j];
      ugradv[i][3*j+1] = derxy[1][i] * funct[j];
      ugradv[i][3*j+2] = derxy[2][i] * funct[j];
    }
  }

/*--------------------------------- now build single stiffness terms ---*/

#define FLUID3_IS_TERM1
#define FLUID3_IS_TERM2
#define FLUID3_IS_TERM3
#define FLUID3_IS_TERM4
#define FLUID3_IS_TERM5
#define FLUID3_IS_TERM6
#define FLUID3_IS_TERM7
#define FLUID3_IS_TERM8
#define FLUID3_IS_TERM9
#define FLUID3_IS_TERM10
#define FLUID3_IS_TERM11
#define FLUID3_IS_TERM12
#define FLUID3_IS_TERM13
#define FLUID3_IS_TERM14
#define FLUID3_IS_TERM15

#define estif_(i,j)    estif[i][j]
#define eforce_(i)     eforce[i]
#define funct_(i)      funct[i]
#define funct_p_(i)    pfunct[i]
#define vderxyz_(i,j)  vderxy[i][j]
#define conv_c_(j)     conv_c[j]
#define conv_g_(j)     conv_g[j]
#define conv_r_(i,j,k) conv_r[i][3*(k)+j]
#define vconv_r_(i,j)  vconv_r[i][j]
#define conv_old_(j)   conv_old[j]
#define conv_g_old_(j) conv_g_old[j]
#define derxyz_(i,j)   derxy[i][j]
#define gridvint_(j)   gridvint[j]
#define velint_(j)     velint[j]
#define viscs2_(i,j,k) viscs2[i][3*(k)+j]
#define visc_old_(i)   visc_old[i]
#define rhsint_(i)     rhsint[i]
#define gradp_(j)      gradp[j]
#define ui             ci
#define vi             ri
#define nu_            visc
#define thsl           timefac

#ifdef FLUID_INCREMENTAL

#ifdef QUASI_NEWTON
  if (!fdyn->qnewton || fdyn->itnum==1)
  {
#endif

    /* This code is generated using MuPAD. Ask me for the MuPAD. u.kue */

    /* We keep two versions: with and without ale. The laster one is a
     * little faster. (more than 10%) */

    if (isale)
    {
#include "f3is_stiff_ale.c"
    }
    else
    {
#include "f3is_stiff.c"
    }

#ifdef QUASI_NEWTON
  }
#endif

  if (isale)
  {
#include "f3is_rhs_incr_ale.c"
  }
  else
  {
#include "f3is_rhs_incr.c"
  }

#else  /* FLUID_INCREMENTAL */

#ifdef QUASI_NEWTON
  dserror("quasi newton not possible with non-incremental elements");
#endif

  if (isale)
  {
#include "f3is_stiff_ale.c"
#include "f3is_rhs_nonincr_ale.c"
  }
  else
  {
#include "f3is_stiff.c"
#include "f3is_rhs_nonincr.c"
  }

#endif  /* FLUID_INCREMENTAL */

#undef estif_
#undef eforce_
#undef conv_c_
#undef conv_g_
#undef conv_r_
#undef vconv_r_
#undef conv_old_
#undef conv_g_old_
#undef derxyz_
#undef gridvint_
#undef velint_
#undef viscs2_
#undef gradp_
#undef ui
#undef vi
#undef funct_
#undef funct_p_
#undef vderxyz_
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
