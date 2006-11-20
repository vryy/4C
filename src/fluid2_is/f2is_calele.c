/*!----------------------------------------------------------------------
\file
\brief element control routine

<pre>
Maintainer: Ulrich Küttler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

------------------------------------------------------------------------*/
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "fluid2_is.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../fluid2/fluid2_prototypes.h"
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
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*/
static ARRAY     eveln_a;  /* element velocities at (n)                 */
static DOUBLE  **eveln;
static ARRAY     evelng_a; /* element velocities at (n+gamma)           */
static DOUBLE  **evelng;
static ARRAY     evhist_a; /* element velocities at (n+gamma)           */
static DOUBLE  **evhist;
static ARRAY     ealecovn_a;  /* element ale-convective velocities      */
static DOUBLE  **ealecovn;    /* at (n)                                 */
static ARRAY     ealecovng_a; /* element ale-convective velocities      */
static DOUBLE  **ealecovng;   /* at (n+gamma)                           */
static ARRAY     egridv_a;    /* element grid velocity                  */
static DOUBLE  **egridv;
static ARRAY     evnng_a;     /* element normal vector at n+1           */
static DOUBLE  **evnng;
static ARRAY     evnn_a;      /* element normal vector at n             */
static DOUBLE  **evnn;
static ARRAY     epren_a;  /* element pressures at (n)	                */
static DOUBLE   *epren;
static ARRAY     edeadn_a; /* element dead load (selfweight)            */
static DOUBLE   *edeadng;
static ARRAY     edeadng_a;/* element dead load (selfweight)            */
static DOUBLE   *edeadn;
static ARRAY     funct_a;  /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;  /* first natural derivatives                 */
static DOUBLE  **deriv;
static ARRAY     deriv2_a; /* second natural derivatives                */
static DOUBLE  **deriv2;
static ARRAY     pfunct_a;  /* shape functions                           */
static DOUBLE   *pfunct;
static ARRAY     pderiv_a;  /* first natural derivatives                 */
static DOUBLE  **pderiv;
static ARRAY     pderiv2_a; /* second natural derivatives                */
static DOUBLE  **pderiv2;
static ARRAY     xyze_a;   /* actual element coordinates                */
static DOUBLE  **xyze;
static ARRAY     xjm_a;    /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     vderxy_a; /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     pderxy_a; /* pre -derivatives                          */
static DOUBLE  **pderxy;
static ARRAY     vderxy2_a;/* vel - 2nd derivatives                     */
static DOUBLE  **vderxy2;
static ARRAY     derxy_a;  /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     derxy2_a; /* 2nd coordinate - derivatives              */
static DOUBLE  **derxy2;
static ARRAY     sigmaint_a; /* fluid stresses at integration point     */
static double  **sigmaint;
static ARRAY     ekappan_a; /* surface curvature at (n)                 */
static DOUBLE   *ekappan;
static ARRAY     ekappang_a; /* surface curvature at (n+1)              */
static DOUBLE   *ekappang;
static ARRAY     ephin_a;    /* height function value at (n)            */
static DOUBLE   *ephin;
static ARRAY     ephing_a;   /* height function value at (n+1)          */
static DOUBLE   *ephing;
static ARRAY     iedgnod_a;
static INT      *iedgnod;
static ARRAY     w1_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa1;      /* used in different element routines        */
static ARRAY     w2_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa2;      /* used in different element routines        */
static DOUBLE  **estif;    /* pointer to global ele-stif                */
static DOUBLE  **emass;    /* pointer to galerkin ele-stif              */
static DOUBLE   *eforce;   /* pointer to RHS                            */
static DOUBLE   *edforce;  /* pointer to RHS due to dirichl. conditions */

static ARRAY     eddy_a;       /* element turbulent ken. energy  at (n+gamma)*/
static DOUBLE   *eddy;

static DOUBLE    visc;

static FLUID_DYNAMIC   *fdyn;
/*!---------------------------------------------------------------------
\brief control routine for element integration of fluid2

<pre>                                                         genk 03/02

This routine controls the element evaluation:
-actual vel. and pres. variables are set
-stabilisation parameters are calculated
-element integration is performed --> element stiffness matrix and
                                  --> element load vectors
-stiffness matrix and load vectors are permuted for assembling
-element load vector due to dirichlet conditions is calculated

</pre>
\param  *ele	         ELEMENT	(i)   actual element
\param  *eleke	         ELEMENT	(i)   element for turbulence-model
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *eforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param  *hasdirich       INT	        (o)   element flag
\param  *hasext          INT	        (o)   element flag
\param   imyrank         INT            (i)   proc number
\param  *ipos                           (i)   node array positions
\param   is_relax        INT            (i)   flag
\param   init	         INT	        (i)   init flag
\return void

------------------------------------------------------------------------*/
void f2is_calele(
  ELEMENT        *ele,
  ELEMENT        *eleke,
  ARRAY          *estif_global,
  ARRAY          *emass_global,
  ARRAY          *eforce_global,
  ARRAY          *edforce_global,
  INT            *hasdirich,
  INT            *hasext,
  INT             imyrank,
  ARRAY_POSITION *ipos,
  INT             is_relax,
  INT             init
  )
{
  INT		readfrom;	/* where to read dbc from 		*/

  DOUBLE estress[3][MAXNOD_F2];   /* ele stresses reprojected for lin eles*/

#ifdef DEBUG
  dstrc_enter("f2is_calele");
#endif

  if (init==1) /* allocate working arrays and set pointers */
  {
    eveln     = amdef("eveln"    ,&eveln_a    ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
    evelng    = amdef("evelng"   ,&evelng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
    evhist    = amdef("evhist"   ,&evhist_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
    ealecovn  = amdef("ealecovn" ,&ealecovn_a ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
    ealecovng = amdef("ealecovng",&ealecovng_a,NUM_F2_VELDOF,MAXNOD_F2,"DA");
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
    pderiv2   = amdef("pderiv2"  ,&pderiv2_a  ,3,MAXNOD_F2,"DA");
    xjm       = amdef("xjm"      ,&xjm_a      ,2,2        ,"DA");
    xyze      = amdef("xyze"     ,&xyze_a     ,2,MAXNOD_F2,"DA");
    vderxy    = amdef("vderxy"   ,&vderxy_a   ,2,2,"DA");
    pderxy    = amdef("pderxy"   ,&pderxy_a   ,2,MAXNOD_F2,"DA");
    vderxy2   = amdef("vderxy2"  ,&vderxy2_a  ,2,3,"DA");
    derxy     = amdef("derxy"    ,&derxy_a    ,2,MAXNOD_F2,"DA");
    derxy2    = amdef("derxy2"   ,&derxy2_a   ,3,MAXNOD_F2,"DA");
    sigmaint  = amdef("sigmaint" ,&sigmaint_a ,3,MAXGAUSS ,"DA");
    ekappan   = amdef("ekappan"  ,&ekappan_a  ,MAXNOD_F2,1 ,"DV");
    ekappang  = amdef("ekappang" ,&ekappang_a ,MAXNOD_F2,1 ,"DV");
    ephin     = amdef("ephin"    ,&ephin_a    ,MAXNOD_F2,1 ,"DV");
    ephing    = amdef("ephing"   ,&ephing_a   ,MAXNOD_F2,1 ,"DV");
    iedgnod   = amdef("iedgnod"  ,&iedgnod_a  ,MAXNOD_F2,1 ,"IV");
    eddy      = amdef("eddy"     ,&eddy_a  ,MAXNOD_F2,1,"DV");
    wa1       = amdef("wa1"      ,&w1_a       ,MAXDOFPERELE,MAXDOFPERELE,"DA");
    wa2       = amdef("wa2"      ,&w2_a       ,MAXDOFPERELE,MAXDOFPERELE,"DA");
/*                                               \- size is arbitrary chosen!  */
    estif   = estif_global->a.da;
    emass   = emass_global->a.da;
    eforce  = eforce_global->a.dv;
    edforce = edforce_global->a.dv;

    fdyn = alldyn[genprob.numff].fdyn;
    goto end;
  } /* endif (init==1) */

#ifdef QUASI_NEWTON
  if (fdyn->qnewton && ele->e.f2is->estif.fdim==0)
    amdef("estif",&ele->e.f2is->estif,estif_global->fdim,estif_global->sdim,"DA");
#endif

/*------------------------------------------------ initialise with ZERO */
#if 0
  amzero(estif_global);
  amzero(emass_global);
  amzero(eforce_global);
  amzero(edforce_global);
#else
  memset(estif_global->a.da[0],0,estif_global->fdim*estif_global->sdim*sizeof(DOUBLE));
  memset(emass_global->a.da[0],0,emass_global->fdim*emass_global->sdim*sizeof(DOUBLE));
  memset(eforce_global->a.dv,0,eforce_global->fdim*sizeof(DOUBLE));
  memset(edforce_global->a.dv,0,edforce_global->fdim*sizeof(DOUBLE));
#endif
  *hasdirich=0;
  *hasext=0;

  switch(ele->e.f2is->is_ale)
  {
  case 0:
    /* set element data */
    f2is_calset(ele,xyze,eveln,evelng,evhist,epren,edeadn,edeadng,ipos,hasext);

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f2_caltau(ele,xyze,funct,deriv,xjm,evelng,visc);

    /*-------------------------------- perform element integration ---*/
    f2is_int_usfem(ele,hasext,estif,eforce,xyze,
		   funct,deriv,deriv2,pfunct,pderiv,pderiv2,
		   xjm,derxy,derxy2,pderxy,evelng,eveln,
		   evhist,NULL,epren,edeadng,
		   vderxy,vderxy2,visc,wa1,wa2,estress, is_relax);
    break;
  case 1:
    /*---------------------------------------------- set element data ---*/
    f2is_calseta(ele,xyze,eveln,evelng,evhist,ealecovn,
	       ealecovng,egridv,epren,edeadn,edeadng,ekappan,ekappang,
	       ephin,ephing,evnng,evnn,ipos,hasext,is_relax);
    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f2_caltau(ele,xyze,funct,deriv,xjm,ealecovng,visc);

    /*-------------------------------- perform element integration ---*/
    f2is_int_usfem(ele,hasext,estif,eforce,xyze,
		   funct,deriv,deriv2,pfunct,pderiv,pderiv2,
		   xjm,derxy,derxy2,pderxy,evelng,eveln,
		   evhist,egridv,epren,edeadng,
		   vderxy,vderxy2,visc,wa1,wa2,estress, is_relax);
    break;
  default:
    dserror("parameter is_ale not 0 or 1!\n");
  }

/* look for neumann bc */
  f2_calneumann(ele, imyrank, eforce, xyze, funct, deriv, xjm, edeadn, edeadng);

/*-------------------------------------------- local co-ordinate system */
  if (ele->locsys==locsys_yes)
    locsys_trans(ele,estif,NULL,NULL,eforce);


/*------------------------------- calculate element load vector edforce */
  if (is_relax)			/* calculation for relaxation parameter	*/
    readfrom = ipos->relax;
  else				/* standard case			*/
    readfrom = ipos->velnp;

#ifdef QUASI_NEWTON
  if (fdyn->qnewton)
  {
    if (fdyn->itnum==1)
    {
      amcopy(estif_global, &ele->e.f2is->estif);
    }
    else
    {
      amcopy(&ele->e.f2is->estif, estif_global);
    }
  }
#endif

#ifndef FLUID_INCREMENTAL
/*------------------------------------------------ condensation of DBCs */
/* estif is in xyz* so edforce is also in xyz* (but DBCs have to be
   tranformed before condensing the dofs                                */
  fluid_caldirich(ele,edforce,estif,hasdirich,readfrom);
#endif

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}


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
void f2is_caleleres(
  ELEMENT        *ele,
  ARRAY          *eforce_global,
  ARRAY_POSITION *ipos,
  INT            *hasdirich,
  INT            *hasext
  )
{
  INT             i;
  NODE           *actnode;

#ifdef DEBUG
  dstrc_enter("f2is_caleleres");
#endif

/*--------------------------------------------- initialise with ZERO ---*/
  amzero(eforce_global);
  *hasdirich=0;
  *hasext=0;

  switch(ele->e.f2->is_ale)
  {
  case 0:
/*---------------------------------------------------- set element data */
    f2is_calset(ele,xyze,eveln,evelng,evhist,epren,edeadn,edeadng,ipos,hasext);

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f2_caltau(ele,xyze,funct,deriv,xjm,evelng,visc);

    /*-------------------------------- perform element integration ---*/
    f2is_int_res(ele,hasext,eforce,xyze,
		 funct,deriv,deriv2,pfunct,pderiv,pderiv2,
		 xjm,derxy,
		 derxy2,pderxy,evelng,evhist,NULL,epren,edeadng,vderxy,
		 vderxy2,visc,wa1,wa2);
    break;
  case 1:
    /*---------------------------------------------- set element data ---*/
    f2is_calseta(ele,xyze,eveln,evelng,evhist,ealecovn,
		 ealecovng,egridv,epren,edeadn,edeadng,ekappan,ekappang,
		 ephin,ephing,evnng,evnn,ipos,hasext,0);

    /*------------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*------------------------------------------------ stab-parameter ---*/
    f2_caltau(ele,xyze,funct,deriv,xjm,ealecovng,visc);

    /*----------------------------------- perform element integration ---*/
    f2is_int_res(ele,hasext,eforce,xyze,
		 funct,deriv,deriv2,pfunct,pderiv,pderiv2,
		 xjm,derxy,
		 derxy2,pderxy,evelng,evhist,egridv,epren,edeadng,
		 vderxy,vderxy2,visc,wa1,wa2);
    break;
  default:
    dserror("parameter is_ale not 0 or 1!\n");
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}

#endif
/*! @} (documentation module close)*/

