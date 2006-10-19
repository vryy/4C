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
\addtogroup FLUID2_PRO
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2_PRO
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2pro_prototypes.h"
#include "../fluid_full/fluid_prototypes.h"
#include "fluid2pro.h"
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
static ARRAY     epren_a;     /* element pressures at (n)	                */
static DOUBLE   *epren;
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
static ARRAY     derxy_a;     /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     derxy2_a;    /* 2nd coordinate - derivatives              */
static DOUBLE  **derxy2;
static ARRAY     sigmaint_a;  /* fluid stresses at integration point     */
static double  **sigmaint;
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
void f2pro_calinit(
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
  dstrc_enter("f2pro_calinit");
#endif

  /* allocate working arrays and set pointers */

  eveln     = amdef("eveln"    ,&eveln_a    ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
  evelng    = amdef("evelng"   ,&evelng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
  evhist    = amdef("evhist"   ,&evhist_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
  egridv    = amdef("egridv"   ,&egridv_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
  evnng     = amdef("evnng"    ,&evnng_a    ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
  evnn      = amdef("evnn"     ,&evnn_a     ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
  epren     = amdef("epren"    ,&epren_a    ,MAXNOD_F2,1,"DV");
  edeadn    = amdef("edeadn"   ,&edeadn_a   ,2,1,"DV");
  edeadng   = amdef("edeadng"  ,&edeadng_a  ,2,1,"DV");
  funct     = amdef("funct"    ,&funct_a    ,MAXNOD_F2,1,"DV");
  deriv     = amdef("deriv"    ,&deriv_a    ,2,MAXNOD_F2,"DA");
  deriv2    = amdef("deriv2"   ,&deriv2_a   ,3,MAXNOD_F2,"DA");
  pfunct    = amdef("pfunct"   ,&pfunct_a   ,MAXNOD_F2,1,"DV");
  pderiv    = amdef("pderiv"   ,&pderiv_a   ,2,MAXNOD_F2,"DA");
  xjm       = amdef("xjm"      ,&xjm_a      ,2,2        ,"DA");
  xyze      = amdef("xyze"     ,&xyze_a     ,2,MAXNOD_F2,"DA");
  vderxy    = amdef("vderxy"   ,&vderxy_a   ,2,2,"DA");
  vderxy2   = amdef("vderxy2"  ,&vderxy2_a  ,2,3,"DA");
  derxy     = amdef("derxy"    ,&derxy_a    ,2,MAXNOD_F2,"DA");
  derxy2    = amdef("derxy2"   ,&derxy2_a   ,3,MAXNOD_F2,"DA");
  pderxy    = amdef("pderxy"   ,&pderxy_a   ,2,MAXNOD_F2,"DA");
  sigmaint  = amdef("sigmaint" ,&sigmaint_a ,3,MAXGAUSS ,"DA");
  ekappan   = amdef("ekappan"  ,&ekappan_a  ,MAXNOD_F2,1 ,"DV");
  ekappang  = amdef("ekappang" ,&ekappang_a ,MAXNOD_F2,1 ,"DV");
  ephin     = amdef("ephin"    ,&ephin_a    ,MAXNOD_F2,1 ,"DV");
  ephing    = amdef("ephing"   ,&ephing_a   ,MAXNOD_F2,1 ,"DV");
  iedgnod   = amdef("iedgnod"  ,&iedgnod_a  ,MAXNOD_F2,1 ,"IV");
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
\brief control routine for element integration of fluid2pro

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
void f2pro_calele(
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
  dstrc_enter("f2pro_calele");
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

  /*---------------------------------------------------- set element data */
  f2pro_calset(ele,xyze,eveln,evelng,evhist,epren,ipos);

  /* We follow the förster-element here. */

  /*---------------------------------------------- get viscosity ---*/
  visc = mat[ele->mat-1].m.fluid->viscosity;

  /*--------------------------------------------- stab-parameter ---*/
  f2_caltau(ele,xyze,funct,deriv,xjm,evelng,visc);

  /*-------------------------------- perform element integration ---*/
  f2pro_int_usfem(ele,hasext,estif,eforce,gforce,xyze,
                  funct,deriv,deriv2,pfunct,pderiv,xjm,
                  derxy,derxy2,pderxy,
                  evelng,eveln,
                  evhist,NULL,epren,edeadng,
                  vderxy,vderxy2,visc,wa1,wa2);

  /*------------------------------------------------ condensation of DBCs */
  /* estif is in xyz* so edforce is also in xyz* (but DBCs have to be
     tranformed before condensing the dofs                                */
  fluid_caldirich(ele,edforce,estif,hasdirich,ipos->velnp);

  /*----------------------------------------------------- local co-system */
  dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of f2pro_calele */



/*----------------------------------------------------------------------*/
/*!
  \brief weak form of pressure gradient

  Einfach den Druckgradienten berechnen. Der zu Grunde liegende Term
  ist  -p(x,y)*div_v  auf der linken Seite.
 */
/*----------------------------------------------------------------------*/
void f2pro_calgradp(
  ELEMENT* ele,
  ARRAY_POSITION *ipos)
{
  INT i;
  INT numpdof;
  INT       iel;        /* number of nodes                                */
  INT       intc;       /* "integration case" for tri for further infos
                           see f2_inpele.c and f2_intg.c                 */
  INT       nir,nis;    /* number of integration nodesin r,s direction    */
  INT       ihoel=0;    /* flag for higher order elements                 */
  INT       icode=2;    /* flag for eveluation of shape functions         */
  INT       lr, ls;     /* counter for integration                        */
  DOUBLE    fac;        /* total integration vactor                       */
  DOUBLE    facr, facs; /* integration weights                            */
  DOUBLE    det;        /* determinant of jacobian matrix at time (n+1)   */
  DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
  DIS_TYP   typ;	      /* element type                                   */
  DISMODE   dm;

  FLUID_DATA      *data;

#ifdef DEBUG
  dstrc_enter("f2pro_calgradp");
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
  for(i=0;i<ele->numnp;i++)
  {
    xyze[0][i]=ele->node[i]->x[0];
    xyze[1][i]=ele->node[i]->x[1];
  }

  switch (ele->e.f2pro->dm)
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
    dserror("unsupported discretization mode %d", ele->e.f2pro->dm);
  }

#if 0
  /*---------------------------------------------- set pressures (n+1) ---*/
  if (numpdof==-1)
  {
    ELEMENT* pele;
    pele = ele->e.f2pro->other;
    for (i=0; i<pele->numnp; ++i)
    {
      epren[i] = pele->node[i]->sol_increment.a.da[1][0];
    }
  }
  else
  {
    for (i=0; i<numpdof; ++i)
    {
      epren[i]   = ele->e.f2pro->press[i];
    }
  }
#endif

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f2pro->dm;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  /*------- get integraton data and check if elements are "higher order" */
  switch (typ)
  {
  case quad4: case quad8: case quad9:  /* --> quad - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f2pro->nGP[0];
    nis = ele->e.f2pro->nGP[1];
    break;
  case tri6: /* --> tri - element */
    icode   = 3;
    ihoel   = 1;
    /* do NOT break at this point!!! */
  case tri3:
    /* initialise integration */
    nir  = ele->e.f2pro->nGP[0];
    nis  = 1;
    intc = ele->e.f2pro->nGP[1];
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
	if (numpdof>0)
	  f2pro_prec(pfunct, pderiv, e1, e2, dm, &numpdof);
	else if (numpdof==-2)
	  f2_rec(pfunct,pderiv,NULL,e1,e2,quad4,2);
        break;
      case tri3: case tri6:   /* --> tri - element */
        e1   = data->txgr[lr][intc];
        facr = data->twgt[lr][intc];
        e2   = data->txgs[lr][intc];
        facs = ONE;
        f2_tri(funct,deriv,deriv2,e1,e2,typ,icode);
        dserror("illegal discretisation mode %d", dm);
        break;
      default:
        dserror("typ unknown!");
      } /* end switch(typ) */

      /*------------------------ compute Jacobian matrix at time n+1 ---*/
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr * facs * det;

      /*----------------------------------- compute global derivates ---*/
      f2_gder(derxy,deriv,xjm,det,iel);

      /* loop over nodes of element */
      for (i=0; i<iel; i++)
      {
	INT p;
	if (numpdof==-1)
	{
	  /*
	   * The number of nodes of the pressure element is supposed
	   * to be equal to the number of my nodes. The shape
	   * functions are supposed to be the same. */
	  for (p=0; p<ele->e.f2pro->other->numnp; ++p)
	  {
	    gradopr[2*i  ][p] += fac*funct[p]*derxy[0][i] ;
	    gradopr[2*i+1][p] += fac*funct[p]*derxy[1][i] ;
	  }
	}
	else if (numpdof==-2)
	{
	  for (p=0; p<ele->e.f2pro->other->numnp; ++p)
	  {
	    gradopr[2*i  ][p] += fac*pfunct[p]*derxy[0][i] ;
	    gradopr[2*i+1][p] += fac*pfunct[p]*derxy[1][i] ;
	  }
	}
	else
	{
	  for (p=0; p<numpdof; ++p)
	  {
	    gradopr[2*i  ][p] += fac*pfunct[p]*derxy[0][i] ;
	    gradopr[2*i+1][p] += fac*pfunct[p]*derxy[1][i] ;
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
          emass[2*i  ][2*j  ] += aux ;
          emass[2*i+1][2*j+1] += aux ;
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

    dest = &(lmass[2*i  ]);
    row = emass[2*i  ];
    *dest = 0;
    for (j=0; j<iel; ++j)
    {
      *dest += row[2*j  ];
    }

    dest = &(lmass[2*i+1]);
    row = emass[2*i+1];
    *dest = 0;
    for (j=0; j<iel; ++j)
    {
      *dest += row[2*j+1];
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*!
  \brief calculate rhs of pressure equation
 */
/*----------------------------------------------------------------------*/
void f2pro_calprhs(
  ELEMENT* ele,
  ARRAY_POSITION *ipos
  )
{
  INT i;
  INT numpdof;
  INT velnp;
  INT       iel;        /* number of nodes                                */
  INT       intc;       /* "integration case" for tri for further infos
                           see f2_inpele.c and f2_intg.c                 */
  INT       nir,nis;    /* number of integration nodesin r,s direction    */
  INT       ihoel=0;    /* flag for higher order elements                 */
  INT       icode=2;    /* flag for eveluation of shape functions         */
  INT       lr, ls;     /* counter for integration                        */
  DOUBLE    fac;        /* total integration vactor                       */
  DOUBLE    facr, facs; /* integration weights                            */
  DOUBLE    det;        /* determinant of jacobian matrix at time (n+1)   */
  DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
  DIS_TYP   typ;        /* element type                                   */
  DISMODE   dm;

  FLUID_DATA      *data;

#ifdef DEBUG
  dstrc_enter("f2pro_calgradp");
#endif

  /* There is not that much to do, so we don't introduce subfunctions
   * here. */

  /* Initialize the fast way. We don't have the array at hand, so
   * amzero is out of reach.
   * Caution! The size here must match the size in calinit! */
  memset(eforce,0,(MAXNOD*MAXDOFPERNODE)*sizeof(DOUBLE));

  /* set element coordinates */
  for(i=0;i<ele->numnp;i++)
  {
    xyze[0][i]=ele->node[i]->x[0];
    xyze[1][i]=ele->node[i]->x[1];
  }

  velnp = ipos->velnp;
  for(i=0;i<ele->numnp;i++)
  {
    NODE* actnode;
    actnode=ele->node[i];

    /*----------------------------------- set element velocities (n+gamma) */
    evelng[0][i]=actnode->sol_increment.a.da[velnp][0];
    evelng[1][i]=actnode->sol_increment.a.da[velnp][1];
  }

  switch (ele->e.f2pro->dm)
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
    dserror("unsupported discretization mode %d", ele->e.f2pro->dm);
  }

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f2pro->dm;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  /*------- get integraton data and check if elements are "higher order" */
  switch (typ)
  {
  case quad4: case quad8: case quad9:  /* --> quad - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f2pro->nGP[0];
    nis = ele->e.f2pro->nGP[1];
    break;
  case tri6: /* --> tri - element */
    icode   = 3;
    ihoel   = 1;
    /* do NOT break at this point!!! */
  case tri3:
    /* initialise integration */
    nir  = ele->e.f2pro->nGP[0];
    nis  = 1;
    intc = ele->e.f2pro->nGP[1];
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
      INT p;
      DOUBLE divu;

      /*------------- get values of  shape functions and their derivatives ---*/
      switch(typ)
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
        e1   = data->qxg[lr][nir-1];
        facr = data->qwgt[lr][nir-1];
        e2   = data->qxg[ls][nis-1];
        facs = data->qwgt[ls][nis-1];
        f2_rec(funct,deriv,deriv2,e1,e2,typ,icode);
	if (numpdof>-1)
	  f2pro_prec(pfunct, pderiv, e1, e2, dm, &numpdof);
	else if (numpdof==-2)
	  f2_rec(pfunct,pderiv,NULL,e1,e2,quad4,2);
        break;
      case tri3: case tri6:   /* --> tri - element */
        e1   = data->txgr[lr][intc];
        facr = data->twgt[lr][intc];
        e2   = data->txgs[lr][intc];
        facs = ONE;
        f2_tri(funct,deriv,deriv2,e1,e2,typ,icode);
        dserror("illegal discretisation mode %d", dm);
        break;
      default:
        dserror("typ unknown!");
      } /* end switch(typ) */

      /*------------------------ compute Jacobian matrix at time n+1 ---*/
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr * facs * det;

      /*----------------------------------- compute global derivates ---*/
      f2_gder(derxy,deriv,xjm,det,iel);

      /*------ get velocity (n+1,i) derivatives at integration point ---*/
      f2_vder(vderxy,derxy,evelng,iel);
      divu = vderxy[0][0] + vderxy[1][1];

      if (numpdof==-1)
      {
	ELEMENT* pele;
	pele = ele->e.f2pro->other;
	for (p=0; p<pele->numnp; ++p)
	{
	  eforce[p] -= fac*funct[p]*divu ;
	}
      }
      else if (numpdof==-2)
      {
	ELEMENT* pele;
	pele = ele->e.f2pro->other;
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
void f2pro_calvelupdate(
  ELEMENT* ele,
  ARRAY_POSITION *ipos
  )
{
  INT i;
  INT numpdof;
  INT       iel;        /* number of nodes                                */
  INT       intc;       /* "integration case" for tri for further infos
                           see f2_inpele.c and f2_intg.c                 */
  INT       nir,nis;    /* number of integration nodesin r,s direction    */
  INT       ihoel=0;    /* flag for higher order elements                 */
  INT       icode=2;    /* flag for eveluation of shape functions         */
  INT       lr, ls;     /* counter for integration                        */
  DOUBLE    fac;        /* total integration vactor                       */
  DOUBLE    facr, facs; /* integration weights                            */
  DOUBLE    det;        /* determinant of jacobian matrix at time (n+1)   */
  DOUBLE    e1,e2;      /* natural coordinates of integr. point           */
  DIS_TYP   typ;        /* element type                                   */
  DISMODE   dm;

  FLUID_DATA      *data;

#ifdef DEBUG
  dstrc_enter("f2pro_calvelupdate");
#endif

  /* There is not that much to do, so we don't introduce subfunctions
   * here. */

  /* Initialize the fast way. We don't have the array at hand, so
   * amzero is out of reach.
   * Caution! The size here must match the size in calinit! */
  memset(eforce,0,(MAXNOD*MAXDOFPERNODE)*sizeof(DOUBLE));

  /* set element coordinates */
  for(i=0;i<ele->numnp;i++)
  {
    xyze[0][i]=ele->node[i]->x[0];
    xyze[1][i]=ele->node[i]->x[1];
  }

  switch (ele->e.f2pro->dm)
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
    dserror("unsupported discretization mode %d", ele->e.f2pro->dm);
  }

  /*---------------------------------------------- set pressures (n+1) ---*/
  if ((numpdof==-1) || (numpdof==-2))
  {
    ELEMENT* pele;
    pele = ele->e.f2pro->other;
    for (i=0; i<pele->numnp; ++i)
    {
      epren[i] = pele->node[i]->sol_increment.a.da[0][0];
    }
  }
  else
  {
    for (i=0; i<numpdof; ++i)
    {
      epren[i]   = ele->e.f2pro->phi[i];
    }
  }

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f2pro->dm;
  fdyn   = alldyn[genprob.numff].fdyn;
  data   = fdyn->data;

  /*------- get integraton data and check if elements are "higher order" */
  switch (typ)
  {
  case quad4: case quad8: case quad9:  /* --> quad - element */
    icode   = 3;
    ihoel   = 1;
    /* initialise integration */
    nir = ele->e.f2pro->nGP[0];
    nis = ele->e.f2pro->nGP[1];
    break;
  case tri6: /* --> tri - element */
    icode   = 3;
    ihoel   = 1;
    /* do NOT break at this point!!! */
  case tri3:
    /* initialise integration */
    nir  = ele->e.f2pro->nGP[0];
    nis  = 1;
    intc = ele->e.f2pro->nGP[1];
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
      DOUBLE ipress;

      /*------------- get values of  shape functions and their derivatives ---*/
      switch(typ)
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
        e1   = data->qxg[lr][nir-1];
        facr = data->qwgt[lr][nir-1];
        e2   = data->qxg[ls][nis-1];
        facs = data->qwgt[ls][nis-1];
        f2_rec(funct,deriv,deriv2,e1,e2,typ,icode);
	if (numpdof>-1)
	  f2pro_prec(pfunct, pderiv, e1, e2, dm, &numpdof);
	else if (numpdof==-2)
	  f2_rec(pfunct,pderiv,NULL,e1,e2,quad4,2);
        break;
      case tri3: case tri6:   /* --> tri - element */
        e1   = data->txgr[lr][intc];
        facr = data->twgt[lr][intc];
        e2   = data->txgs[lr][intc];
        facs = ONE;
        f2_tri(funct,deriv,deriv2,e1,e2,typ,icode);
        dserror("illegal discretisation mode %d", dm);
        break;
      default:
        dserror("typ unknown!");
      } /* end switch(typ) */

      /*------------------------ compute Jacobian matrix at time n+1 ---*/
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr * facs * det;

      /*----------------------------------- compute global derivates ---*/
      f2_gder(derxy,deriv,xjm,det,iel);

      /*
       * Now do the weak form of the pressure increment gradient. Used
       * for the velocity update. */

      ipress = 0;
      if (numpdof==-1)
      {
	for (i=0; i<ele->e.f2pro->other->numnp; ++i)
	{
	  ipress += funct[i] * epren[i];
	}
      }
      else if (numpdof==-2)
      {
	for (i=0; i<ele->e.f2pro->other->numnp; ++i)
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
	eforce[2*i  ] += -ipress*fac*derxy[0][i] ;
	eforce[2*i+1] += -ipress*fac*derxy[1][i] ;
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
void f2pro_addnodepressure(ELEMENT* ele, INT k, DOUBLE* pressure)
{
  INT i;
  INT numpdof;
  INT iel;
  DIS_TYP   typ;
  DISMODE   dm;

  /* edge positions in element coordinates */
  static INT edges[] = {
    -1, -1,
    1, -1,
    1, 1,
    -1, 1,
    0, -1,
    1, 0,
    0, 1,
    -1, 0,
    0, 0
  };

#ifdef DEBUG
  dstrc_enter("f2pro_addnodepressure");
#endif

  switch (ele->e.f2pro->dm)
  {
  case dm_q2pm1:
    numpdof = 3;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f2pro->dm);
  }

  for (i=0; i<numpdof; ++i)
  {
    epren[i]   = ele->e.f2pro->press[i];
  }

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f2pro->dm;

  /* This used to be a loop the binary io demands a node wise
   * calculation. */
  /*for (k=0; k<iel; ++k)*/
  {
    DOUBLE press;
    INT lr;
    INT ls;

    lr = edges[2*k];
    ls = edges[2*k+1];

    switch(typ)
    {
    case quad4: case quad8: case quad9:   /* --> quad - element */
      f2pro_prec(pfunct, pderiv, lr, ls, dm, &numpdof);
      break;
    case tri3: case tri6:   /* --> tri - element */
      f2_tri(funct,deriv,deriv2,lr,ls,typ,0);
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

#if defined(DEBUG) && !defined(PARALLEL)

void f2pro_debugoutpressure(ELEMENT* ele, FILE* f)
{
  INT i;
  INT k;
  INT numpdof;
  INT iel;
  DIS_TYP   typ;
  DISMODE   dm;

  INT edges[] = {
    -1, -1,
    1, -1,
    1, 1,
    -1, 1,
    -1, -1
  };

  for(i=0;i<ele->numnp;i++)
  {
    xyze[0][i]=ele->node[i]->x[0];
    xyze[1][i]=ele->node[i]->x[1];
  }

  switch (ele->e.f2pro->dm)
  {
  case dm_q2pm1:
    numpdof = 3;
    break;
  case dm_q1p0:
    numpdof = 1;
    break;
  default:
    dserror("unsupported discretization mode %d", ele->e.f2pro->dm);
  }

  for (i=0; i<numpdof; ++i)
  {
    epren[i]   = ele->e.f2pro->press[i];
  }

  /*--------------------------------------------------- initialisation ---*/
  iel    = ele->numnp;
  typ    = ele->distyp;
  dm     = ele->e.f2pro->dm;

  /* Loop the corners of a quad. We do the start node twice to get a
   * closed line in gnuplot. */
  for (k=0; k<10; k+=2)
  {
    INT lr;
    INT ls;
    DOUBLE x, y;
    DOUBLE press;

    lr = edges[k];
    ls = edges[k+1];

    switch(typ)
    {
    case quad4: case quad8: case quad9:   /* --> quad - element */
      f2_rec(funct,deriv,deriv2,lr,ls,typ,0);
      f2pro_prec(pfunct, pderiv, lr, ls, dm, &numpdof);
      break;
    case tri3: case tri6:   /* --> tri - element */
      f2_tri(funct,deriv,deriv2,lr,ls,typ,0);
      dserror("illegal discretisation mode %d", dm);
      break;
    default:
      dserror("typ unknown!");
    }

    /*
     * That is ridiculous but it frees us from matching the element
     * edges to the node coordinates. */
    x = 0;
    y = 0;
    for (i=0; i<iel; ++i)
    {
      x += funct[i] * xyze[0][i];
      y += funct[i] * xyze[1][i];
    }

    press = 0;
    for (i=0; i<numpdof; ++i)
      press += pfunct[i] * epren[i];

    fprintf(f, "%f %f %f\n", x, y, press);
  }
  fprintf(f, "\n\n");
}

#endif

#endif
/*! @} (documentation module close)*/
