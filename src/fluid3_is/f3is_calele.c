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
#ifdef D_FLUID3_IS
#include "../headers/standardtypes.h"
#include "fluid3_is.h"
#include "../fluid3/fluid3_prototypes.h"
#include "../fluid_full/fluid_prototypes.h"

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
/*!----------------------------------------------------------------------
\brief positions of physical values in node arrays
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
static ARRAY     ehist_a;  /* element history data                      */
static DOUBLE  **ehist;
static ARRAY     eveln_a;  /* element velocities at (n)	      	        */
static DOUBLE  **eveln;
static ARRAY     evelng_a; /* element velocities at (n+gamma)		*/
static DOUBLE  **evelng;
static ARRAY     ealecovn_a;  /* element ale-convective velocities      */
static DOUBLE  **ealecovn;    /* at (n)                                 */
static ARRAY     ealecovng_a; /* element ale-convective velocities      */
static DOUBLE  **ealecovng;   /* at (n+gamma)                           */
static ARRAY     egridv_a; /* element grid velocity                     */
static DOUBLE  **egridv;
static ARRAY     epren_a;  /* element pressures at (n)  		*/
static DOUBLE   *epren;
static DOUBLE   *edeadng;
static ARRAY     edeadng_a;/* element dead load (selfweight)            */
static ARRAY     funct_a;  /* shape functions				*/
static DOUBLE   *funct;
static ARRAY     deriv_a;  /* first natural derivatives 		*/
static DOUBLE  **deriv;
static ARRAY     deriv2_a; /* second natural derivatives		*/
static DOUBLE  **deriv2;
static ARRAY     pfunct_a;  /* shape functions                           */
static DOUBLE   *pfunct;
static ARRAY     pderiv_a;  /* first natural derivatives                 */
static DOUBLE  **pderiv;
static ARRAY     pderiv2_a; /* second natural derivatives                */
static DOUBLE  **pderiv2;
static ARRAY     xyze_a;
static DOUBLE  **xyze;
static ARRAY     xjm_a;    /* jocobian matrix				*/
static DOUBLE  **xjm;
static ARRAY     vderxy_a; /* vel - derivatives 			*/
static DOUBLE  **vderxy;
static ARRAY     pderxy_a; /* pre -derivatives  			*/
static DOUBLE  **pderxy;
static ARRAY     vderxy2_a;/* vel - 2nd derivatives			*/
static DOUBLE  **vderxy2;
static ARRAY     derxy_a;  /* coordinate - derivatives  		*/
static DOUBLE  **derxy;
static ARRAY     derxy2_a; /* 2nd coordinate - derivatives		*/
static DOUBLE  **derxy2;
static ARRAY     sigmaint_a; /* fluid stresses at integration point     */
static DOUBLE  **sigmaint;
static ARRAY     ephin_a;    /* height function value at (n)            */
static DOUBLE   *ephin;
static ARRAY     ephing_a;   /* height function value at (n+1)          */
static DOUBLE   *ephing;
static ARRAY     iedgnod_a;
static INT      *iedgnod;
static ARRAY     w1_a;     /* working array of arbitrary chosen size	*/
static DOUBLE  **wa1;      /* used in different element routines	*/
static ARRAY     w2_a;     /* working array of arbitrary chosen size	*/
static DOUBLE  **wa2;      /* used in different element routines	*/
static DOUBLE  **estif;    /* pointer to global ele-stif		*/
static DOUBLE  **emass;    /* pointer to galerkin ele-stif		*/
static DOUBLE   *eforce;   /* pointer to RHS                            */
static DOUBLE   *edforce;  /* pointer to RHS due to dirichl. conditions */

static DOUBLE    visc;

static FLUID_DYNAMIC   *fdyn;
/*!---------------------------------------------------------------------
\brief control routine for element integration of fluid3

<pre>                                                         genk 05/02

This routine controls the element evaluation:
-actual vel. and pres. variables are set
-stabilisation parameters are calculated
-element integration is performed --> element stiffness matrix and
                                  --> element load vectors
-stiffness matrix and load vectors are permuted for assembling
-element load vector due to dirichlet conditions is calculated

</pre>
\param  *ele	         ELEMENT	(i)   actual element
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *eforce_global   ARRAY	        (o)   ele force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param  *ipos                           (i)   node array positions
\param  *hasdirich       INT	        (o)   element flag
\param  *hasext          INT	        (o)   element flag
\param   init	         INT	        (i)   init flag
\return void

------------------------------------------------------------------------*/
void f3is_calele(
  ELEMENT        *ele,
  ARRAY          *estif_global,
  ARRAY          *emass_global,
  ARRAY          *eforce_global,
  ARRAY          *edforce_global,
  ARRAY_POSITION *ipos,
  INT            *hasdirich,
  INT            *hasext,
  INT             is_relax,
  INT             init
  )
{
  INT		readfrom;	/* where to read dbc from 		*/

#ifdef DEBUG
  dstrc_enter("f3_calele");
#endif

  if (init==1) /* allocate working arrays and set pointers */
  {
    ehist     = amdef("ehist"  ,&ehist_a  ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
    eveln     = amdef("eveln"  ,&eveln_a  ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
    evelng    = amdef("evelng" ,&evelng_a ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
    ealecovn  = amdef("ealecovn" ,&ealecovn_a ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
    ealecovng = amdef("ealecovng",&ealecovng_a,NUM_F3_VELDOF,MAXNOD_F3,"DA");
    egridv    = amdef("egridv"   ,&egridv_a   ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
    epren     = amdef("epren"  ,&epren_a  ,MAXNOD_F3,1,"DV");
    edeadng   = amdef("edeadng",&edeadng_a,3,1,"DV");
    funct     = amdef("funct"  ,&funct_a  ,MAXNOD_F3,1,"DV");
    deriv     = amdef("deriv"  ,&deriv_a  ,3,MAXNOD_F3,"DA");
    deriv2    = amdef("deriv2" ,&deriv2_a ,6,MAXNOD_F3,"DA");
    pfunct    = amdef("pfunct"   ,&pfunct_a   ,MAXNOD_F3,1,"DV");
    pderiv    = amdef("pderiv"   ,&pderiv_a   ,3,MAXNOD_F3,"DA");
    pderiv2   = amdef("pderiv2"  ,&pderiv2_a  ,6,MAXNOD_F3,"DA");
    xjm       = amdef("xjm"    ,&xjm_a    ,3,3        ,"DA");
    xyze      = amdef("xyze"   ,&xyze_a   ,3,MAXNOD_F3,"DA");
    vderxy    = amdef("vderxy" ,&vderxy_a ,3,3,"DA");
    pderxy    = amdef("pderxy" ,&pderxy_a ,3,MAXNOD_F3,"DA");
    vderxy2   = amdef("vderxy2",&vderxy2_a,3,6,"DA");
    derxy     = amdef("derxy"  ,&derxy_a  ,3,MAXNOD_F3,"DA");
    derxy2    = amdef("derxy2" ,&derxy2_a ,6,MAXNOD_F3,"DA");
    sigmaint  = amdef("sigmaint" ,&sigmaint_a ,MAXGAUSS,6,"DA");
    ephin     = amdef("ephin"    ,&ephin_a    ,MAXNOD_F3,1 ,"DV");
    ephing    = amdef("ephing"   ,&ephing_a   ,MAXNOD_F3,1 ,"DV");
    iedgnod   = amdef("iedgnod"  ,&iedgnod_a  ,MAXNOD_F3,1 ,"IV");
    wa1       = amdef("wa1"    ,&w1_a      ,300,300        ,"DA");
    wa2       = amdef("wa2"    ,&w2_a      ,300,300        ,"DA");
/*                                        \- size is chosen arbitrarily! */
    estif   = estif_global->a.da;
    emass   = emass_global->a.da;
    eforce  = eforce_global->a.dv;
    edforce = edforce_global->a.dv;

    fdyn    = alldyn[genprob.numff].fdyn;
    goto end;
  } /* endif (init==1) */

#ifdef QUASI_NEWTON
  if (fdyn->qnewton && ele->e.f3->estif.fdim==0)
    amdef("estif",&ele->e.f3->estif,estif_global->fdim,estif_global->sdim,"DA");
#endif

/*------------------------------------------------ initialise with ZERO */
#if 1
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

  switch(ele->e.f3->is_ale)
  {
  case 0:
/*---------------------------------------------------- set element data */
    f3is_calset(ele,xyze,ehist,evelng,epren,edeadng,ipos,hasext);

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f3_caltau(ele,xyze,funct,deriv,derxy,xjm,evelng,wa1,visc);

    /*-------------------------------- perform element integration ---*/
    f3is_int_usfem(ele,hasext,estif,eforce,xyze,
                   funct,deriv,deriv2,pfunct,pderiv,pderiv2,
		   xjm,derxy,derxy2,pderxy,evelng,
                   ehist,NULL,epren,edeadng,
                   vderxy,vderxy2,visc,wa1,wa2,is_relax);
    break;
  case 1:
    /* set element data */
    f3is_calseta(ele,xyze,ehist,evelng,
		 ealecovng,egridv,epren,edeadng,ipos,hasext,is_relax);

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f3_caltau(ele,xyze,funct,deriv,derxy,xjm,ealecovng,wa1,visc);

    /*-------------------------------- perform element integration ---*/
    f3is_int_usfem(ele,hasext,estif,eforce,xyze,
                   funct,deriv,deriv2,pfunct,pderiv,pderiv2,
		   xjm,derxy,derxy2,pderxy,evelng,
                   ehist,egridv,epren,edeadng,
                   vderxy,vderxy2,visc,wa1,wa2,is_relax);
    break;
  default:
    dserror("parameter is_ale not 0 or 1!\n");
  }  /* end switch */

  /* local co-ordinate system */
  if(ele->locsys==locsys_yes)
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
      amcopy(estif_global, &ele->e.f3->estif);
    }
    else
    {
      amcopy(&ele->e.f3->estif, estif_global);
    }
  }
#endif

/*------------------------------------------------ condensation of DBCs */
/* estif is in xyz* so edforce is also in xyz* (but DBCs have to be
   tranformed before condensing the dofs                                */
#ifdef FLUID_INCREMENTAL
  /* with incremental fluid we want dirichlet forces only during
   * steepest descent relaxation factor calculation */
  if (is_relax)
#endif
  fluid_caldirich(ele,edforce,estif,hasdirich,readfrom);

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
void f3is_caleleres(
  ELEMENT        *ele,
  ARRAY          *eforce_global,
  ARRAY_POSITION *ipos,
  INT            *hasdirich,
  INT            *hasext
  )
{
#ifdef DEBUG
  dstrc_enter("f3is_caleleres");
#endif

/*--------------------------------------------- initialise with ZERO ---*/
  amzero(eforce_global);
  *hasdirich=0;
  *hasext=0;

  switch(ele->e.f3->is_ale)
  {
  case 0:
/*---------------------------------------------------- set element data */
    f3is_calset(ele,xyze,ehist,evelng,epren,edeadng,ipos,hasext);

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f3_caltau(ele,xyze,funct,deriv,derxy,xjm,evelng,wa1,visc);

    /*-------------------------------- perform element integration ---*/
    f3is_int_res(ele,hasext,eforce,xyze,
		 funct,deriv,deriv2,pfunct,pderiv,pderiv2,
		 xjm,derxy,
		 derxy2,pderxy,evelng,ehist,NULL,epren,edeadng,vderxy,
		 vderxy2,visc,wa1,wa2);
    break;
  case 1:
    /* set element data */
    f3is_calseta(ele,xyze,ehist,evelng,
		 ealecovng,egridv,epren,edeadng,ipos,hasext,0);

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f3_caltau(ele,xyze,funct,deriv,derxy,xjm,ealecovng,wa1,visc);

    /*----------------------------------- perform element integration ---*/
    f3is_int_res(ele,hasext,eforce,xyze,
		 funct,deriv,deriv2,pfunct,pderiv,pderiv2,
		 xjm,derxy,
		 derxy2,pderxy,evelng,ehist,egridv,epren,edeadng,
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


/*----------------------------------------------------------------------*/
/*!
  \brief calculate reaction forces for SD relaxation

  We just calculated a linear fluid solution at the current state
  without any rhs and with the residuum prescribes at the fsi
  interface. Now we need to know the fluid reaction forces. We simply
  recalculate the element matrices at the interface and multiply with
  the known solution. This way we get consistent nodal forces.

  Note: Only the dofs belonging to the interface are calculated.

  \param ele           (i) the element
  \param estif_global  (-) global stiffness matrix
  \param eforce_global (o) consistent nodal forces at the interface
  \param ipos          (i) fluid field array positions
  \param hasdirich     (-) dirichlet flag
  \param hasext        (-) ext flag

  \author u.kue
  \date 01/07
 */
/*----------------------------------------------------------------------*/
void f3is_caleleres_relax(ELEMENT        *ele,
			  ARRAY          *estif_global,
			  ARRAY          *eforce_global,
			  ARRAY_POSITION *ipos,
			  INT            *hasdirich,
			  INT            *hasext)
{
  INT is_relax = 1;
#ifdef DEBUG
  dstrc_enter("f3is_caleleres_relax");
#endif

#ifdef QUASI_NEWTON
  dserror("quasi newton hack not supported with SD");
#endif

  /*------------------------------------------------ initialise with ZERO */
  amzero(estif_global);
  amzero(eforce_global);
  *hasdirich=0;
  *hasext=0;

  memset(emass[0],0,estif_global->fdim*estif_global->sdim*sizeof(DOUBLE));

  /* The point here is to calculate the element matrix and to apply
   * the (independent) solution afterwards. */

  switch(ele->e.f3->is_ale)
  {
  case 0:
/*---------------------------------------------------- set element data */
    f3is_calset(ele,xyze,ehist,evelng,epren,edeadng,ipos,hasext);

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f3_caltau(ele,xyze,funct,deriv,derxy,xjm,evelng,wa1,visc);

    /*-------------------------------- perform element integration ---*/
    f3is_int_usfem(ele,hasext,estif,eforce,xyze,
                   funct,deriv,deriv2,pfunct,pderiv,pderiv2,
		   xjm,derxy,derxy2,pderxy,evelng,
                   ehist,NULL,epren,edeadng,
                   vderxy,vderxy2,visc,wa1,wa2,is_relax);
    break;
  case 1:
    /* set element data */
    f3is_calseta(ele,xyze,ehist,evelng,
		 ealecovng,egridv,epren,edeadng,ipos,hasext,is_relax);

    /*---------------------------------------------- get viscosity ---*/
    visc = mat[ele->mat-1].m.fluid->viscosity;

    /*--------------------------------------------- stab-parameter ---*/
    f3_caltau(ele,xyze,funct,deriv,derxy,xjm,ealecovng,wa1,visc);

    /*-------------------------------- perform element integration ---*/
    f3is_int_usfem(ele,hasext,estif,eforce,xyze,
                   funct,deriv,deriv2,pfunct,pderiv,pderiv2,
		   xjm,derxy,derxy2,pderxy,evelng,
                   ehist,egridv,epren,edeadng,
                   vderxy,vderxy2,visc,wa1,wa2,is_relax);
    break;
  default:
    dserror("parameter is_ale not 0 or 1!\n");
  }

  /* Use stiffness matrix to calculate reaction forces. */
  fluid_reaction_forces(ele, fdyn,
			estif_global->a.da,
			eforce_global->a.dv,
			ipos->relax);

#ifdef DEBUG
  dstrc_exit();
#endif
}

#endif
