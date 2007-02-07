/*!----------------------------------------------------------------------
\file
\brief element control routine

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID3_PRO
*//*! @{ (documentation module open)*/
#ifdef D_FLUID3_PRO
#include "../headers/standardtypes.h"
#include "../fluid3/fluid3_prototypes.h"
#include "fluid3pro_prototypes.h"
#include "../fluid_full/fluid_prototypes.h"
#include "fluid3pro.h"
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
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;


/* Variables used by all element integration functions */

static ARRAY     evelnm_a;    /* element velocities at (n-1)               */
static DOUBLE  **evelnm;
static ARRAY     eveln_a;     /* element velocities at (n)                 */
static DOUBLE  **eveln;
static ARRAY     evelng_a;    /* element velocities at (n+gamma)           */
static DOUBLE  **evelng;
static ARRAY     evhist_a;    /* element velocities at (n+gamma)           */
static DOUBLE  **evhist;
static ARRAY     egridv_a;    /* element grid velocity                  */
static DOUBLE  **egridv;
static ARRAY     evnng_a;     /* element normal vector at n+1           */
static DOUBLE  **evnng;
static ARRAY     evnn_a;      /* element normal vector at n             */
static DOUBLE  **evnn;
static ARRAY     epren_a;     /* element pressures at (n+1)             */
static DOUBLE   *epren;
static ARRAY     eprenm_a;    /* element pressures at (n)               */
static DOUBLE   *eprenm;
static ARRAY     edeadn_a;    /* element dead load (selfweight)            */
static DOUBLE   *edeadng;
static ARRAY     edeadng_a;   /* element dead load (selfweight)            */
static DOUBLE   *edeadn;
static ARRAY     funct_a;     /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;     /* first natural derivatives                 */
static DOUBLE  **deriv;
static ARRAY     deriv2_a;    /* second natural derivatives                */
static DOUBLE  **deriv2;
static ARRAY     pfunct_a;    /* pressure shape functions                  */
static DOUBLE   *pfunct;
static ARRAY     pderiv_a;    /* pressure first natural derivatives        */
static DOUBLE  **pderiv;
static ARRAY     xyze_a;      /* actual element coordinates                */
static DOUBLE  **xyze;
static ARRAY     xjm_a;       /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     vderxy_a;    /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     pderxy_a;    /* pre -derivatives                          */
static DOUBLE  **pderxy;
static ARRAY     vderxy2_a;   /* vel - 2nd derivatives                     */
static DOUBLE  **vderxy2;
static ARRAY     vderxy_n_a;  /* vel - derivatives                         */
static DOUBLE  **vderxy_n;
static ARRAY     vderxy_nm_a; /* vel - derivatives                         */
static DOUBLE  **vderxy_nm;
static ARRAY     derxy_a;     /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     derxy2_a;    /* 2nd coordinate - derivatives              */
static DOUBLE  **derxy2;
static ARRAY     ekappan_a;   /* surface curvature at (n)                 */
static DOUBLE   *ekappan;
static ARRAY     ekappang_a;  /* surface curvature at (n+1)              */
static DOUBLE   *ekappang;
static ARRAY     ephin_a;     /* height function value at (n)            */
static DOUBLE   *ephin;
static ARRAY     ephing_a;    /* height function value at (n+1)          */
static DOUBLE   *ephing;
static ARRAY     iedgnod_a;
static INT      *iedgnod;
static ARRAY     w1_a;        /* working array of arbitrary chosen size    */
static DOUBLE  **wa1;         /* used in different element routines        */
static ARRAY     w2_a;        /* working array of arbitrary chosen size    */
static DOUBLE  **wa2;         /* used in different element routines        */
static DOUBLE  **estif;       /* pointer to global ele-stif                */
static DOUBLE  **emass;       /* pointer to galerkin ele-stif              */
static DOUBLE   *lmass;       /* pointer to galerkin lumped ele mass matrix */
static DOUBLE  **gradopr;     /* pointer to gradient operator               */
static DOUBLE   *eforce;      /* pointer to RHS                            */
static DOUBLE   *edforce;     /* pointer to RHS due to dirichl. conditions */
static DOUBLE   *gforce;
static FLUID_DYNAMIC   *fdyn;


/*----------------------------------------------------------------------*/
/*!
 \brief init the element calculating functions
 */
/*----------------------------------------------------------------------*/
void f3pro_calinit(
  ARRAY          *estif_global,
  ARRAY          *emass_global,
  ARRAY          *lmass_global,
  ARRAY          *gradopr_global,
  ARRAY          *eforce_global,
  ARRAY          *edforce_global,
  ARRAY          *gforce_global,
  ARRAY_POSITION *ipos
  )
{
#ifdef DEBUG
  dstrc_enter("f3pro_calinit");
#endif

  /* allocate working arrays and set pointers */

  evelnm    = amdef("evelnm"   ,&evelnm_a   ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
  eveln     = amdef("eveln"    ,&eveln_a    ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
  evelng    = amdef("evelng"   ,&evelng_a   ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
  evhist    = amdef("evhist"   ,&evhist_a   ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
  egridv    = amdef("egridv"   ,&egridv_a   ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
  evnng     = amdef("evnng"    ,&evnng_a    ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
  evnn      = amdef("evnn"     ,&evnn_a     ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
  epren     = amdef("epren"    ,&epren_a    ,MAXNOD_F3,1,"DV");
  eprenm    = amdef("eprenm"   ,&eprenm_a   ,MAXNOD_F3,1,"DV");
  edeadn    = amdef("edeadn"   ,&edeadn_a   ,3,1,"DV");
  edeadng   = amdef("edeadng"  ,&edeadng_a  ,3,1,"DV");
  funct     = amdef("funct"    ,&funct_a    ,MAXNOD_F3,1,"DV");
  deriv     = amdef("deriv"    ,&deriv_a    ,3,MAXNOD_F3,"DA");
  deriv2    = amdef("deriv2"   ,&deriv2_a   ,6,MAXNOD_F3,"DA");
  pfunct    = amdef("pfunct"   ,&pfunct_a   ,MAXNOD_F3,1,"DV");
  pderiv    = amdef("pderiv"   ,&pderiv_a   ,3,MAXNOD_F3,"DA");
  xjm       = amdef("xjm"      ,&xjm_a      ,3,3        ,"DA");
  xyze      = amdef("xyze"     ,&xyze_a     ,3,MAXNOD_F3,"DA");
  vderxy    = amdef("vderxy"   ,&vderxy_a   ,3,3,"DA");
  vderxy2   = amdef("vderxy2"  ,&vderxy2_a  ,3,6,"DA");
  vderxy_n  = amdef("vderxy_n" ,&vderxy_n_a ,3,3,"DA");
  vderxy_nm = amdef("vderxy_nm",&vderxy_nm_a,3,3,"DA");
  derxy     = amdef("derxy"    ,&derxy_a    ,3,MAXNOD_F3,"DA");
  derxy2    = amdef("derxy2"   ,&derxy2_a   ,6,MAXNOD_F3,"DA");
  pderxy    = amdef("pderxy"   ,&pderxy_a   ,3,MAXNOD_F3,"DA");
  ekappan   = amdef("ekappan"  ,&ekappan_a  ,MAXNOD_F3,1 ,"DV");
  ekappang  = amdef("ekappang" ,&ekappang_a ,MAXNOD_F3,1 ,"DV");
  ephin     = amdef("ephin"    ,&ephin_a    ,MAXNOD_F3,1 ,"DV");
  ephing    = amdef("ephing"   ,&ephing_a   ,MAXNOD_F3,1 ,"DV");
  iedgnod   = amdef("iedgnod"  ,&iedgnod_a  ,MAXNOD_F3,1 ,"IV");
  wa1       = amdef("wa1"      ,&w1_a       ,MAXDOFPERELE,MAXDOFPERELE,"DA");
  wa2       = amdef("wa2"      ,&w2_a       ,MAXDOFPERELE,MAXDOFPERELE,"DA");

  estif   = estif_global->a.da;
  emass   = emass_global->a.da;
  lmass   = lmass_global->a.dv;
  gradopr = gradopr_global->a.da;
  eforce  = eforce_global->a.dv;
  edforce = edforce_global->a.dv;
  gforce  = gforce_global->a.dv;

  fdyn = alldyn[genprob.numff].fdyn;

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
\brief control routine for element integration of fluid3pro

\param  *elev	         ELEMENT	(i)   actual velocity element
\param  *elep	         ELEMENT	(i)   actual pressure element
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *lmass_global    ARRAY	        (o)   lumped mass matrix
\param  *gradopr_global  ARRAY	        (o)   gradient operator
\param  *etforce_global  ARRAY	        (o)   element time force
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param  *ipos                           (i)   node array positions
\param  *hasdirich       int	        (o)   element flag

*/
/*----------------------------------------------------------------------*/
void f3pro_calele(
  ELEMENT        *ele,
  ARRAY          *estif_global,
  ARRAY          *emass_global,
  ARRAY          *lmass_global,
  ARRAY          *gradopr_global,
  ARRAY          *eforce_global,
  ARRAY          *edforce_global,
  ARRAY          *gforce_global,
  ARRAY_POSITION *ipos,
  INT            *hasdirich,
  INT            *hasext
  )
{
  DOUBLE visc;
  INT iel;

#ifdef DEBUG
  dstrc_enter("f3pro_calele");
#endif

#ifdef QUASI_NEWTON
    if (alldyn[genprob.numff].fdyn->qnewton && ele->e.f3pro->estif.fdim==0)
      amdef("estif",&ele->e.f3pro->estif,estif_global->fdim,estif_global->sdim,"DA");
#endif

  iel=ele->numnp;
  /*------------------------------------------------ initialise with ZERO */
  amzero(estif_global);
  amzero(emass_global);
  amzero(lmass_global);
  amzero(gradopr_global);
  amzero(eforce_global);
  amzero(edforce_global);
  amzero(gforce_global);
  *hasdirich=0;
  *hasext=0;

  switch(ele->e.f3pro->is_ale)
  {
  case 0:
    /*---------------------------------------------------- set element data */
    f3pro_calset(ele,xyze,evelnm,eveln,evelng,evhist,eprenm,epren,edeadng,ipos);

    /* We follow the Förster-element here. */

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f3_caltau(ele,xyze,funct,deriv,derxy,xjm,evelng,wa1,visc);

    /*-------------------------------- perform element integration ---*/
    f3pro_int_usfem(ele,hasext,estif,eforce,gforce,xyze,
                    funct,deriv,deriv2,pfunct,pderiv,xjm,
                    derxy,derxy2,pderxy,
                    evelng,eveln,evelnm,
                    evhist,NULL,eprenm,epren,edeadng,
                    vderxy,vderxy2,vderxy_n,vderxy_nm,visc,wa1,wa2,0);
    break;
  case 1:
    f3pro_calseta(ele,xyze,evelnm,eveln,evelng,evhist,egridv,eprenm,epren,edeadng,ipos);

    /* We follow the Förster-element here. */

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f3_caltau(ele,xyze,funct,deriv,derxy,xjm,evelng,wa1,visc);

    /*-------------------------------- perform element integration ---*/
    f3pro_int_usfem(ele,hasext,estif,eforce,gforce,xyze,
                    funct,deriv,deriv2,pfunct,pderiv,xjm,
                    derxy,derxy2,pderxy,
                    evelng,eveln,evelnm,
                    evhist,egridv,eprenm,epren,edeadng,
                    vderxy,vderxy2,vderxy_n,vderxy_nm,visc,wa1,wa2,0);
    break;
  default:
    dserror("parameter is_ale not 0 or 1!\n");
  }

#ifdef QUASI_NEWTON
  if (fdyn->qnewton)
  {
    if (fdyn->itnum==1)
    {
      amcopy(estif_global, &ele->e.f3pro->estif);
    }
    else
    {
      amcopy(&ele->e.f3pro->estif, estif_global);
    }
  }
#endif

#ifndef FLUID_INCREMENTAL
  /*------------------------------------------------ condensation of DBCs */
  /* estif is in xyz* so edforce is also in xyz* (but DBCs have to be
     tranformed before condensing the dofs                                */
  fluid_caldirich(ele,edforce,estif,hasdirich,ipos->velnp);
#endif

  /*----------------------------------------------------- local co-system */
  dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f3pro_calele */


/*!---------------------------------------------------------------------
\brief control routine for integration of element residual

<pre>                                                        chfoe 11/04

This routine controls the integration of the elemental residual which is
required to compute consistent nodal forces. These are also used to be
FSI coupling forces

</pre>

\param  *ele	         ELEMENT	(i)   actual element
\param  *eforce_global   ARRAY	        (o)   ele iteration force
\param  *ipos                           (i)   node array positions
\param  *hasdirich       INT	        (o)   element flag
\param  *hasext          INT	        (o)   element flag
\return void

------------------------------------------------------------------------*/
void f3pro_caleleres(
  ELEMENT        *ele,
  ARRAY          *eforce_global,
  ARRAY_POSITION *ipos,
  INT            *hasdirich,
  INT            *hasext
  )
{
  DOUBLE visc;

#ifdef DEBUG
  dstrc_enter("f2pro_caleleres");
#endif

#ifdef QUASI_NEWTON
  dserror("QUASI_NEWTON hack not supported!");
#endif

/*--------------------------------------------- initialise with ZERO ---*/
  amzero(eforce_global);
  *hasdirich=0;
  *hasext=0;

  switch (ele->e.f3pro->is_ale)
  {
  case 0:
    /*---------------------------------------------------- set element data */
    f3pro_calset(ele,xyze,evelnm,eveln,evelng,evhist,eprenm,epren,edeadng,ipos);

    /* We follow the Förster-element here. */

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f3_caltau(ele,xyze,funct,deriv,derxy,xjm,evelng,wa1,visc);

    /*-------------------------------- perform element integration ---*/
    f3pro_int_res(ele,hasext,eforce,xyze,funct,deriv,deriv2,pfunct,pderiv,xjm,derxy,
                  derxy2,pderxy,evelng,evhist,NULL,epren,edeadng,vderxy,
                  vderxy2,visc,wa1,wa2);
    break;
  case 1:
    f3pro_calseta(ele,xyze,evelnm,eveln,evelng,evhist,egridv,eprenm,epren,edeadng,ipos);

    /* We follow the Förster-element here. */

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f3_caltau(ele,xyze,funct,deriv,derxy,xjm,evelng,wa1,visc);

    /*----------------------------------- perform element integration ---*/
    f3pro_int_res(ele,hasext,eforce,xyze,funct,deriv,deriv2,pfunct,pderiv,xjm,derxy,
                  derxy2,pderxy,evelng,evhist,egridv,epren,edeadng,vderxy,
                  vderxy2,visc,wa1,wa2);
    break;
  default:
    dserror("parameter is_ale not 0 or 1!\n");
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*----------------------------------------------------------------------*/
/*!
  \brief weak form of pressure gradient

  Einfach den Druckgradienten berechnen. Der zu Grunde liegende Term
  ist  -p(x,y)*div_v  auf der linken Seite.
 */
/*----------------------------------------------------------------------*/
void f3pro_calgradp(
  ELEMENT* ele,
  ARRAY_POSITION *ipos)
{
  INT i;
  INT numpdof;
  INT       iel;		/* number of nodes                                */
  INT       intc=0;		/* "integration case" for tri for further infos
				   see f3_inpele.c and f3_intg.c                 */
  INT       nir=0,nis=0,nit=0;	/* number of integration nodesin r,s,t direction  */
  INT       ihoel=0;		/* flag for higher order elements                 */
  INT       icode=2;		/* flag for eveluation of shape functions         */
  INT       lr, ls, lt;		/* counter for integration                        */
  DOUBLE    fac;		/* total integration vactor                       */
  DOUBLE    facr=0, facs=0, fact=0;	/* integration weights                            */
  DOUBLE    det;		/* determinant of jacobian matrix at time (n+1)   */
  DOUBLE    e1,e2,e3;		/* natural coordinates of integr. point           */
  DIS_TYP   typ;		/* element type                                   */
  DISMODE   dm;

  FLUID_DATA      *data;

#ifdef DEBUG
  dstrc_enter("f3pro_calgradp");
#endif

  /* There is not that much to do, so we don't introduce subfunctions
   * here. */

  /* Initialize the fast way. We don't have the array at hand, so
   * amzero is out of reach.
   * Caution! The size here must match the size in calinit! */
  memset(gradopr[0],0,(MAXNOD*MAXDOFPERNODE)*(MAXNOD*MAXDOFPERNODE)*sizeof(DOUBLE));
  memset(emass[0]  ,0,(MAXNOD*MAXDOFPERNODE)*(MAXNOD*MAXDOFPERNODE)*sizeof(DOUBLE));
  /*memset(eforce,0,(MAXNOD*MAXDOFPERNODE)*sizeof(DOUBLE));*/

  /* set element coordinates */
  if (ele->e.f3pro->is_ale)
  {
    f3_alecoor(ele,xyze);
  }
  else
  {
    for(i=0;i<ele->numnp;i++)
    {
      xyze[0][i]=ele->node[i]->x[0];
      xyze[1][i]=ele->node[i]->x[1];
      xyze[2][i]=ele->node[i]->x[2];
    }
  }

  switch (ele->e.f3pro->dm)
  {
  case dm_q2pm1:
    numpdof = 4;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  case dm_q1q1:
  case dm_q2q2:
    numpdof = -1;
    break;
  case dm_q2q1:
    numpdof = -2;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f3pro->dm);
  }

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f3pro->dm;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  /*------- get integraton data and check if elements are "higher order" */
  switch (typ)
  {
  case hex8: case hex20: case hex27:  /* --> hex - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f3pro->nGP[0];
    nis = ele->e.f3pro->nGP[1];
    nit = ele->e.f3pro->nGP[2];
    intc= 0;
    break;
  case tet10: /* --> tet - element */
    icode   = 3;
    ihoel   = 1;
    /* do NOT break at this point!!! */
  case tet4:    /* initialise integration */
    nir  = ele->e.f3pro->nGP[0];
    nis  = 1;
    nit  = 1;
    intc = ele->e.f3pro->nGP[1];
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
	  if (numpdof>0)
	    f3pro_phex(pfunct, pderiv, e1, e2, e3, dm, &numpdof);
	  else if (numpdof==-2)
	    f3_hex(pfunct,pderiv,NULL,e1,e2,e3,hex8,2);
	  break;
	case tet4: case tet10:   /* --> tet - element */
	  e1   = data->txgr[lr][intc];
	  facr = data->twgt[lr][intc];
	  e2   = data->txgs[lr][intc];
	  facs = ONE;
	  e3   = data->txgt[lr][intc];
	  fact = ONE;
	  f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,icode);
	  dserror("illegal discretisation mode %d", dm);
	  break;
	default:
	  dserror("typ unknown!");
	}

	/*------------------------ compute Jacobian matrix at time n+1 ---*/
	f3_jaco(xyze,deriv,xjm,&det,ele,iel);
	fac = facr*facs*fact*det;

	/*----------------------------------- compute global derivates ---*/
	f3_gder(derxy,deriv,xjm,wa1,det,iel);

	/* loop over nodes of element */
	for (i=0; i<iel; i++)
	{
	  INT p;
	  if (numpdof==-1)
	  {
	    for (p=0; p<ele->e.f3pro->other->numnp; ++p)
	    {
	      gradopr[3*i  ][p] += fac*funct[p]*derxy[0][i] ;
	      gradopr[3*i+1][p] += fac*funct[p]*derxy[1][i] ;
	      gradopr[3*i+2][p] += fac*funct[p]*derxy[2][i] ;
	    }
	  }
	  else if (numpdof==-2)
	  {
	    for (p=0; p<ele->e.f3pro->other->numnp; ++p)
	    {
	      gradopr[3*i  ][p] += fac*pfunct[p]*derxy[0][i] ;
	      gradopr[3*i+1][p] += fac*pfunct[p]*derxy[1][i] ;
	      gradopr[3*i+2][p] += fac*pfunct[p]*derxy[2][i] ;
	    }
	  }
	  else
	  {
	    for (p=0; p<numpdof; ++p)
	    {
	      gradopr[3*i  ][p] += fac*pfunct[p]*derxy[0][i] ;
	      gradopr[3*i+1][p] += fac*pfunct[p]*derxy[1][i] ;
	      gradopr[3*i+2][p] += fac*pfunct[p]*derxy[2][i] ;
	    }
	  }
	}

	/* build mass matrix */
	for (i=0; i<iel; ++i)
	{
	  INT j;
	  for (j=0; j<iel; ++j)
	  {
	    DOUBLE aux = fac*funct[i]*funct[j];
	    emass[3*i  ][3*j  ] += aux ;
	    emass[3*i+1][3*j+1] += aux ;
	    emass[3*i+2][3*j+2] += aux ;
	  }
	}
      }
    }
  }

  /* lump mass matrix */
  /* A little explicit. It's supposed to be fast. */
  for (i=0; i<iel; ++i)
  {
    register INT j;
    register DOUBLE* row;
    register DOUBLE* dest;

    dest = &(lmass[3*i  ]);
    row = emass[3*i  ];
    *dest = 0;
    for (j=0; j<iel; ++j)
    {
      *dest += row[3*j  ];
    }

    /* all three entries are the same here */
    lmass[3*i+1] = lmass[3*i+2] = lmass[3*i];

#if 0
    dest = &(lmass[3*i+1]);
    row = emass[3*i+1];
    *dest = 0;
    for (j=0; j<iel; ++j)
    {
      *dest += row[3*j+1];
    }

    dest = &(lmass[3*i+2]);
    row = emass[3*i+2];
    *dest = 0;
    for (j=0; j<iel; ++j)
    {
      *dest += row[3*j+2];
    }
#endif
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*!
  \brief calculate rhs of pressure equation

  Maybe it would indeed be faster to do this multiplication globally.
 */
/*----------------------------------------------------------------------*/
void f3pro_calprhs(
  ELEMENT* ele,
  ARRAY_POSITION *ipos
  )
{
  INT i;
  INT numpdof;
  INT velnp;
  INT       iel;		/* number of nodes                                */
  INT       intc=0;		/* "integration case" for tri for further infos
				   see f3_inpele.c and f3_intg.c                 */
  INT       nir=0,nis=0,nit=0;	/* number of integration nodesin r,s,t direction  */
  INT       ihoel=0;		/* flag for higher order elements                 */
  INT       icode=2;		/* flag for eveluation of shape functions         */
  INT       lr, ls, lt;		/* counter for integration                        */
  DOUBLE    fac;		/* total integration vactor                       */
  DOUBLE    facr=0, facs=0, fact=0;	/* integration weights                            */
  DOUBLE    det;		/* determinant of jacobian matrix at time (n+1)   */
  DOUBLE    e1,e2,e3;		/* natural coordinates of integr. point           */
  DIS_TYP   typ;		/* element type                                   */
  DISMODE   dm;

  FLUID_DATA      *data;

#ifdef DEBUG
  dstrc_enter("f3pro_calgradp");
#endif

  /* There is not that much to do, so we don't introduce subfunctions
   * here. */

  /* Initialize the fast way. We don't have the array at hand, so
   * amzero is out of reach.
   * Caution! The size here must match the size in calinit! */
  memset(eforce,0,(MAXNOD*MAXDOFPERNODE)*sizeof(DOUBLE));

  /* set element coordinates */
  if (ele->e.f3pro->is_ale)
  {
    f3_alecoor(ele,xyze);
  }
  else
  {
    for(i=0;i<ele->numnp;i++)
    {
      xyze[0][i]=ele->node[i]->x[0];
      xyze[1][i]=ele->node[i]->x[1];
      xyze[2][i]=ele->node[i]->x[2];
    }
  }

  velnp = ipos->velnp;
  for(i=0;i<ele->numnp;i++)
  {
    NODE* actnode;
    actnode=ele->node[i];

    /*----------------------------------- set element velocities (n+gamma) */
    evelng[0][i]=actnode->sol_increment.a.da[velnp][0];
    evelng[1][i]=actnode->sol_increment.a.da[velnp][1];
    evelng[2][i]=actnode->sol_increment.a.da[velnp][2];
  }

  switch (ele->e.f3pro->dm)
  {
  case dm_q2pm1:
    numpdof = 4;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  case dm_q1q1:
  case dm_q2q2:
    numpdof = -1;
    break;
  case dm_q2q1:
    numpdof = -2;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f3pro->dm);
  }

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f3pro->dm;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  /*------- get integraton data and check if elements are "higher order" */
  switch (typ)
  {
  case hex8: case hex20: case hex27:  /* --> hex - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f3pro->nGP[0];
    nis = ele->e.f3pro->nGP[1];
    nit = ele->e.f3pro->nGP[2];
    intc= 0;
    break;
  case tet10: /* --> tet - element */
    icode   = 3;
    ihoel   = 1;
    /* do NOT break at this point!!! */
  case tet4:    /* initialise integration */
    nir  = ele->e.f3pro->nGP[0];
    nis  = 1;
    nit  = 1;
    intc = ele->e.f3pro->nGP[1];
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
      for (lt=0;lt<nit;lt++)
      {
	INT p;
	DOUBLE divu;

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
	  if (numpdof>-1)
	    f3pro_phex(pfunct, pderiv, e1, e2, e3, dm, &numpdof);
	  else if (numpdof==-2)
	    f3_hex(pfunct,pderiv,NULL,e1,e2,e3,hex8,2);
	  break;
	case tet4: case tet10:   /* --> tet - element */
	  e1   = data->txgr[lr][intc];
	  facr = data->twgt[lr][intc];
	  e2   = data->txgs[lr][intc];
	  facs = ONE;
	  e3   = data->txgt[lr][intc];
	  fact = ONE;
	  f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,icode);
	  dserror("illegal discretisation mode %d", dm);
	  break;
	default:
	  dserror("typ unknown!");
	}

	/*------------------------ compute Jacobian matrix at time n+1 ---*/
	f3_jaco(xyze,deriv,xjm,&det,ele,iel);
	fac = facr*facs*fact*det;

	/*----------------------------------- compute global derivates ---*/
	f3_gder(derxy,deriv,xjm,wa1,det,iel);

	/*------ get velocity (n+1,i) derivatives at integration point ---*/
	f3_vder(vderxy,derxy,evelng,iel);
	divu = vderxy[0][0] + vderxy[1][1] + vderxy[2][2];

	if (numpdof==-1)
	{
	  ELEMENT* pele;
	  pele = ele->e.f3pro->other;
	  for (p=0; p<pele->numnp; ++p)
	  {
	    eforce[p] -= fac*funct[p]*divu ;
	  }
	}
	else if (numpdof==-2)
	{
	  ELEMENT* pele;
	  pele = ele->e.f3pro->other;
	  for (p=0; p<pele->numnp; ++p)
	  {
	    eforce[p] -= fac*pfunct[p]*divu ;
	  }
	}
	else
	{
	  for (p=0; p<numpdof; ++p)
	  {
	    eforce[p] -= fac*pfunct[p]*divu ;
	  }
	}
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief calculate rhs
 */
/*----------------------------------------------------------------------*/
void f3pro_calvrhs(ELEMENT* ele, ARRAY_POSITION *ipos)
{
  INT i;
  INT numpdof;
  INT velnp;
  INT       iel;        /* number of nodes                                */
  INT       intc=0;       /* "integration case" for tri for further infos
                           see f2_inpele.c and f2_intg.c                 */
  INT       nir=0,nis=0,nit=0;/* number of integration nodesin r,s direction    */
  INT       ihoel=0;    /* flag for higher order elements                 */
  INT       icode=2;    /* flag for eveluation of shape functions         */
  INT       lr, ls, lt; /* counter for integration                        */
  DOUBLE    fac;        /* total integration vactor                       */
  DOUBLE    facr=0, facs=0, fact=0; /* integration weights                            */
  DOUBLE    det;        /* determinant of jacobian matrix at time (n+1)   */
  DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
  DIS_TYP   typ;        /* element type                                   */
  DISMODE   dm;

  FLUID_DATA      *data;
  DOUBLE gradp[3];
  DOUBLE velint[3];     /* velocity vector at integration point           */

#ifdef DEBUG
  dstrc_enter("f3pro_calvrhs");
#endif

  /* There is not that much to do, so we don't introduce subfunctions
   * here. */

  /* Initialize the fast way. We don't have the array at hand, so
   * amzero is out of reach.
   * Caution! The size here must match the size in calinit! */
  memset(eforce,0,(MAXNOD*MAXDOFPERNODE)*sizeof(DOUBLE));

  /* set element coordinates */
  if (ele->e.f3pro->is_ale)
  {
    f3_alecoor(ele,xyze);
  }
  else
  {
    for(i=0;i<ele->numnp;i++)
    {
      xyze[0][i]=ele->node[i]->x[0];
      xyze[1][i]=ele->node[i]->x[1];
      xyze[2][i]=ele->node[i]->x[2];
    }
  }

  velnp = ipos->velnp;
  for(i=0;i<ele->numnp;i++)
  {
    NODE* actnode;
    actnode=ele->node[i];

    /*----------------------------------- set element velocities (n+gamma) */
    evelng[0][i]=actnode->sol_increment.a.da[velnp][0];
    evelng[1][i]=actnode->sol_increment.a.da[velnp][1];
    evelng[2][i]=actnode->sol_increment.a.da[velnp][2];
  }

  switch (ele->e.f3pro->dm)
  {
  case dm_q2pm1:
    numpdof = 3;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  case dm_q1q1:
  case dm_q2q2:
    numpdof = -1;
    break;
  case dm_q2q1:
    numpdof = -2;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f3pro->dm);
  }

  /*---------------------------------------------- set pressures (n+1) ---*/
  if ((numpdof==-1) || (numpdof==-2))
  {
    ELEMENT* pele;
    pele = ele->e.f3pro->other;
    for (i=0; i<pele->numnp; ++i)
    {
      epren[i] = pele->node[i]->sol_increment.a.da[0][0];
    }
  }
  else
  {
    for (i=0; i<numpdof; ++i)
    {
      epren[i]   = ele->e.f3pro->phi[i];
    }
  }

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f3pro->dm;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  /*------- get integraton data and check if elements are "higher order" */
  switch (typ)
  {
  case hex8: case hex20: case hex27:  /* --> hex - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f3pro->nGP[0];
    nis = ele->e.f3pro->nGP[1];
    nit = ele->e.f3pro->nGP[2];
    intc= 0;
    break;
  case tet10: /* --> tet - element */
    icode   = 3;
    ihoel   = 1;
    /* do NOT break at this point!!! */
  case tet4:    /* initialise integration */
    nir  = ele->e.f3pro->nGP[0];
    nis  = 1;
    nit  = 1;
    intc = ele->e.f3pro->nGP[1];
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
	  if (numpdof>-1)
	    f3pro_phex(pfunct, pderiv, e1, e2, e3, dm, &numpdof);
	  else if (numpdof==-2)
	    f3_hex(pfunct,pderiv,NULL,e1,e2,e3,hex8,2);
	  break;
	case tet4: case tet10:   /* --> tet - element */
	  e1   = data->txgr[lr][intc];
	  facr = data->twgt[lr][intc];
	  e2   = data->txgs[lr][intc];
	  facs = ONE;
	  e3   = data->txgt[lr][intc];
	  fact = ONE;
	  f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,icode);
	  dserror("illegal discretisation mode %d", dm);
	  break;
	default:
	  dserror("typ unknown!");
	}

	/*------------------------ compute Jacobian matrix at time n+1 ---*/
	f3_jaco(xyze,deriv,xjm,&det,ele,iel);
	fac = facr*facs*fact*det;

	/*----------------------------------- compute global derivates ---*/
	f3_gder(derxy,deriv,xjm,wa1,det,iel);
        f3_gder(pderxy,pderiv,xjm,wa1,det,8);

	/*---------------- get velocities (n+1,i) at integration point ---*/
	f3_veci(velint,funct,evelng,iel);

	/*
	 * Now do the weak form of the pressure increment gradient. Used
	 * for the velocity update. */

	/*------------------------------------- get pressure gradients ---*/
	gradp[0] = gradp[1] = gradp[2] = 0.0;

	if (numpdof==-1)
	{
	  for (i=0; i<iel; i++)
	  {
	    gradp[0] += derxy[0][i] * epren[i];
	    gradp[1] += derxy[1][i] * epren[i];
	    gradp[2] += derxy[2][i] * epren[i];
	  }
	}
	else if (numpdof==-2)
	{
	  for (i=0; i<ele->e.f3pro->other->numnp; i++)
	  {
	    gradp[0] += pderxy[0][i] * epren[i];
	    gradp[1] += pderxy[1][i] * epren[i];
	    gradp[2] += pderxy[2][i] * epren[i];
	  }
	}
	else
	{
	  for (i=0; i<numpdof; i++)
	  {
	    gradp[0] += pderxy[0][i] * epren[i];
	    gradp[1] += pderxy[1][i] * epren[i];
	    gradp[2] += pderxy[2][i] * epren[i];
	  }
	}

	for (i=0; i<iel; i++)
	{
	  eforce[3*i  ] += -gradp[0]*fac*funct[i] ;
	  eforce[3*i+1] += -gradp[1]*fac*funct[i] ;
	  eforce[3*i+2] += -gradp[2]*fac*funct[i] ;
	}

#if 0
	for (i=0; i<iel; i++)
	{
	  eforce[3*i  ] += fac*funct[i]*velint[0] ;
	  eforce[3*i+1] += fac*funct[i]*velint[1] ;
	  eforce[3*i+2] += fac*funct[i]*velint[2] ;
	}
#endif
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief calculate the gradient of the pressure increment for the
  final velocity update
 */
/*----------------------------------------------------------------------*/
void f3pro_calvelupdate(
  ELEMENT* ele,
  ARRAY_POSITION *ipos
  )
{
  INT i;
  INT numpdof;
  INT       iel;		/* number of nodes                                */
  INT       intc=0;		/* "integration case" for tri for further infos
				   see f3_inpele.c and f3_intg.c                 */
  INT       nir=0,nis=0,nit=0;	/* number of integration nodesin r,s,t direction  */
  INT       ihoel=0;		/* flag for higher order elements                 */
  INT       icode=2;		/* flag for eveluation of shape functions         */
  INT       lr, ls, lt;		/* counter for integration                        */
  DOUBLE    fac;		/* total integration vactor                       */
  DOUBLE    facr=0, facs=0, fact=0;	/* integration weights                            */
  DOUBLE    det;		/* determinant of jacobian matrix at time (n+1)   */
  DOUBLE    e1,e2,e3;		/* natural coordinates of integr. point           */
  DIS_TYP   typ;		/* element type                                   */
  DISMODE   dm;

  FLUID_DATA      *data;

#ifdef DEBUG
  dstrc_enter("f3pro_calvelupdate");
#endif

  /* There is not that much to do, so we don't introduce subfunctions
   * here. */

  /* Initialize the fast way. We don't have the array at hand, so
   * amzero is out of reach.
   * Caution! The size here must match the size in calinit! */
  memset(eforce,0,(MAXNOD*MAXDOFPERNODE)*sizeof(DOUBLE));

  /* set element coordinates */
  if (ele->e.f3pro->is_ale)
  {
    f3_alecoor(ele,xyze);
  }
  else
  {
    for(i=0;i<ele->numnp;i++)
    {
      xyze[0][i]=ele->node[i]->x[0];
      xyze[1][i]=ele->node[i]->x[1];
      xyze[2][i]=ele->node[i]->x[2];
    }
  }

  switch (ele->e.f3pro->dm)
  {
  case dm_q2pm1:
    numpdof = 4;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  case dm_q1q1:
  case dm_q2q2:
    numpdof = -1;
    break;
  case dm_q2q1:
    numpdof = -2;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f3pro->dm);
  }

  /*---------------------------------------------- set pressures (n+1) ---*/
  if ((numpdof==-1) || (numpdof==-2))
  {
    ELEMENT* pele;
    pele = ele->e.f3pro->other;
    for (i=0; i<pele->numnp; ++i)
    {
      epren[i] = pele->node[i]->sol_increment.a.da[0][0];
    }
  }
  else
  {
    for (i=0; i<numpdof; ++i)
    {
      epren[i]   = ele->e.f3pro->phi[i];
    }
  }

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f3pro->dm;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  /*------- get integraton data and check if elements are "higher order" */
  switch (typ)
  {
  case hex8: case hex20: case hex27:  /* --> hex - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f3pro->nGP[0];
    nis = ele->e.f3pro->nGP[1];
    nit = ele->e.f3pro->nGP[2];
    intc= 0;
    break;
  case tet10: /* --> tet - element */
    icode   = 3;
    ihoel   = 1;
    /* do NOT break at this point!!! */
  case tet4:    /* initialise integration */
    nir  = ele->e.f3pro->nGP[0];
    nis  = 1;
    nit  = 1;
    intc = ele->e.f3pro->nGP[1];
    break;
  default:
    dserror("typ unknown!");
  }


  /*----------------------------------------------------------------------*
   |               start loop over integration points                     |
   *----------------------------------------------------------------------*/
  for (lr=0;lr<nir;lr++)
  {
    for (ls=0;ls<nis;ls++)
    {
      for (lt=0;lt<nit;lt++)
      {
	DOUBLE ipress;

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
	  if (numpdof>-1)
	    f3pro_phex(pfunct, pderiv, e1, e2, e3, dm, &numpdof);
	  else if (numpdof==-2)
	    f3_hex(pfunct,pderiv,NULL,e1,e2,e3,hex8,2);
	  break;
	case tet4: case tet10:   /* --> tet - element */
	  e1   = data->txgr[lr][intc];
	  facr = data->twgt[lr][intc];
	  e2   = data->txgs[lr][intc];
	  facs = ONE;
	  e3   = data->txgt[lr][intc];
	  fact = ONE;
	  f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,icode);
	  dserror("illegal discretisation mode %d", dm);
	  break;
	default:
	  dserror("typ unknown!");
	}

	/*------------------------ compute Jacobian matrix at time n+1 ---*/
	f3_jaco(xyze,deriv,xjm,&det,ele,iel);
	fac = facr*facs*fact*det;

	/*----------------------------------- compute global derivates ---*/
	f3_gder(derxy,deriv,xjm,wa1,det,iel);

	/*
	 * Now do the weak form of the pressure increment gradient. Used
	 * for the velocity update. */

	ipress = 0;
	if (numpdof==-1)
	{
	  for (i=0; i<ele->e.f3pro->other->numnp; ++i)
	  {
	    ipress += funct[i] * epren[i];
	  }
	}
	else if (numpdof==-2)
	{
	  for (i=0; i<ele->e.f3pro->other->numnp; ++i)
	  {
	    ipress += pfunct[i] * epren[i];
	  }
	}
	else
	{
	  for (i=0; i<numpdof; ++i)
	    ipress += pfunct[i] * epren[i];
	}

	/* loop over nodes of element */
	for (i=0; i<iel; i++)
	{
	  eforce[3*i  ] += -ipress*fac*derxy[0][i] ;
	  eforce[3*i+1] += -ipress*fac*derxy[1][i] ;
	  eforce[3*i+2] += -ipress*fac*derxy[2][i] ;
	}
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief add this elements share of node pressure

  \param ele       (i) the element
  \param k         (i) the node those pressure is calculated
  \param pressure  (o) the pressure value
 */
/*----------------------------------------------------------------------*/
void f3pro_addnodepressure(ELEMENT* ele, INT k, DOUBLE* pressure)
{
  INT i;
  INT numpdof;
  INT iel;
  DIS_TYP   typ;
  DISMODE   dm;

  /* edge positions in element coordinates */
  static INT edges[] = {
    1,   1,   1,
    1,  -1,   1,
    -1,  -1,   1,
    -1,   1,   1,

    1,   1,  -1,
    1,  -1,  -1,
    -1,  -1,  -1,
    -1,   1,  -1,

    1,   0,   1,
    0,  -1,   1,
    -1,   0,   1,
    0,   1,   1,

    1,   1,   0,
    1,  -1,   0,
    -1,  -1,   0,
    -1,   1,   0,

    1,   0,  -1,
    0,  -1,  -1,
    -1,   0,  -1,
    0,   1,  -1,

    0,   0,   1,

    1,   0,   0,
    0,  -1,   0,
    -1,   0,   0,
    0,   1,   0,

    0,   0,  -1,

    0,   0,   0
  };

#ifdef DEBUG
  dstrc_enter("f2pro_addnodepressure");
#endif

  switch (ele->e.f3pro->dm)
  {
  case dm_q2pm1:
    numpdof = 4;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f3pro->dm);
  }

  for (i=0; i<numpdof; ++i)
  {
    epren[i]   = ele->e.f3pro->press[i];
  }

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f3pro->dm;

  /* This used to be a loop the binary io demands a node wise
   * calculation. */
  /*for (k=0; k<iel; ++k)*/
  {
    DOUBLE press;
    INT lr;
    INT ls;
    INT lt;

    lr = edges[3*k];
    ls = edges[3*k+1];
    lt = edges[3*k+2];

    /*------------- get values of  shape functions and their derivatives ---*/
    switch(typ)
    {
    case hex8: case hex20: case hex27:   /* --> hex - element */
      f3pro_phex(pfunct, pderiv, lr, ls, lt, dm, &numpdof);
      break;
    case tet4: case tet10:   /* --> tet - element */
      dserror("illegal discretisation mode %d", dm);
      break;
    default:
      dserror("typ unknown!");
    }

    press = 0;
    for (i=0; i<numpdof; ++i)
      press += pfunct[i] * epren[i];

    *pressure = press/ele->node[k]->numele;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief calculate laplacian pressure matrix

  \param ele       (i) the element
 */
/*----------------------------------------------------------------------*/
void f3pro_calpress(ELEMENT* ele,
		    ARRAY* estif_global,
		    ARRAY* emass_global,
		    ARRAY_POSITION* ipos)
{
  INT i;
  INT numpdof;
  INT       iel;        /* number of nodes                                */
  INT       intc=0;       /* "integration case" for tri for further infos
                           see f2_inpele.c and f2_intg.c                 */
  INT       nir=0,nis=0,nit=0;/* number of integration nodesin r,s direction    */
  INT       ihoel=0;    /* flag for higher order elements                 */
  INT       icode=2;    /* flag for eveluation of shape functions         */
  INT       lr,ls,lt;   /* counter for integration                        */
  DOUBLE    fac;        /* total integration vactor                       */
  DOUBLE    facr=0, facs=0, fact=0; /* integration weights                            */
  DOUBLE    det;        /* determinant of jacobian matrix at time (n+1)   */
  DOUBLE    e1,e2,e3;   /* natural coordinates of integr. point           */
  DIS_TYP   typ;        /* element type                                   */
  DISMODE   dm;

  FLUID_DATA      *data;

#ifdef DEBUG
  dstrc_enter("f3pro_calpress");
#endif

  amzero(estif_global);
  amzero(emass_global);

  /* set element coordinates */
  if (ele->e.f3pro->is_ale)
  {
    f3_alecoor(ele,xyze);
  }
  else
  {
    for(i=0;i<ele->numnp;i++)
    {
      xyze[0][i]=ele->node[i]->x[0];
      xyze[1][i]=ele->node[i]->x[1];
      xyze[2][i]=ele->node[i]->x[2];
    }
  }

  switch (ele->e.f3pro->dm)
  {
  case dm_q2pm1:
    numpdof = 3;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  case dm_q1q1:
  case dm_q2q2:
    numpdof = -1;
    break;
  case dm_q2q1:
    numpdof = -2;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f3pro->dm);
  }

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f3pro->dm;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  /*------- get integraton data and check if elements are "higher order" */
  switch (typ)
  {
  case hex8: case hex20: case hex27:  /* --> hex - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f3pro->nGP[0];
    nis = ele->e.f3pro->nGP[1];
    nit = ele->e.f3pro->nGP[2];
    intc= 0;
    break;
  case tet10: /* --> tet - element */
    icode   = 3;
    ihoel   = 1;
    /* do NOT break at this point!!! */
  case tet4:    /* initialise integration */
    nir  = ele->e.f3pro->nGP[0];
    nis  = 1;
    nit  = 1;
    intc = ele->e.f3pro->nGP[1];
    break;
  default:
    dserror("typ unknown!");
  }

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
	  if (numpdof>-1)
	    f3pro_phex(pfunct, pderiv, e1, e2, e3, dm, &numpdof);
	  else if (numpdof==-2)
	    f3_hex(pfunct,pderiv,NULL,e1,e2,e3,hex8,2);
	  break;
	case tet4: case tet10:   /* --> tet - element */
	  e1   = data->txgr[lr][intc];
	  facr = data->twgt[lr][intc];
	  e2   = data->txgs[lr][intc];
	  facs = ONE;
	  e3   = data->txgt[lr][intc];
	  fact = ONE;
	  f3_tet(funct,deriv,deriv2,e1,e2,e3,typ,icode);
	  dserror("illegal discretisation mode %d", dm);
	  break;
	default:
	  dserror("typ unknown!");
	}

	/*------------------------ compute Jacobian matrix at time n+1 ---*/
	f3_jaco(xyze,deriv,xjm,&det,ele,iel);
	fac = facr*facs*fact*det;

	/* build mass matrix */
	for (i=0; i<iel; ++i)
	{
	  INT j;
	  for (j=0; j<iel; ++j)
	  {
	    DOUBLE aux = fac*funct[i]*funct[j];
	    emass[3*i  ][3*j  ] += aux ;
	    emass[3*i+1][3*j+1] += aux ;
	    emass[3*i+2][3*j+2] += aux ;
	  }
	}

	if (numpdof == -2)
	{
	  /*----------------------------------- compute global derivates ---*/
	  f3_gder(derxy,pderiv,xjm,wa1,det,8);

	  for (i=0; i<8; i++)
	  {
	    INT j;
	    for (j=0; j<8; j++)
	    {
	      estif[i][j] += fac*(derxy[0][i]*derxy[0][j] +
				  derxy[1][i]*derxy[1][j] +
				  derxy[2][i]*derxy[2][j]) ;
	    }
	  }
	}
	else
	{
	  dserror("not supported");
	}
      }
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif
/*! @} (documentation module close)*/
