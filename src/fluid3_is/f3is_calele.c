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
#ifdef D_FLUID3
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
		 ealecovng,egridv,epren,edeadng,ipos,hasext);

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

#ifndef FLUID_INCREMENTAL
/*------------------------------------------------ condensation of DBCs */
/* estif is in xyz* so edforce is also in xyz* (but DBCs have to be
   tranformed before condensing the dofs                                */
  fluid_caldirich(ele,edforce,estif,hasdirich,ipos->velnp);
#endif

  /* amprint(stdout,estif_global,27*3+8,27*3+8); */
#if 0
  {
    static FILE* f;
    INT i;
    FILE* bin;
    char buf[50];

    sprintf(buf,"erhs.%d.data",ele->Id);
    bin = fopen(buf,"w");
    fwrite(eforce,sizeof(DOUBLE),4*8+3*19,bin);
    fclose(bin);

    sprintf(buf,"edrhs.%d.data",ele->Id);
    bin = fopen(buf,"w");
    fwrite(edforce,sizeof(DOUBLE),4*8+3*19,bin);
    fclose(bin);

    sprintf(buf,"estiff.%d.data",ele->Id);
    bin = fopen(buf,"w");
    for (i=0; i<4*8+3*19; ++i)
    {
      fwrite(estif[i],sizeof(DOUBLE),4*8+3*19,bin);
    }
    fclose(bin);

    if (f==NULL)
    {
      f = fopen("elements.py","w");
      fprintf(f,"from Matrix import Matrix\n\n"
	      "element = []\n"
	      "rhs = []\n"
	      "vn = []\n\n");
    }

    if (ele->Id < 10)
    {
      fprintf(f,"element.append(([\n");
      for (i=0; i<4*8+3*19; ++i)
      {
	INT j;
	fprintf(f,"[");
	for (j=0; j<4*8+3*19; ++j)
	{
	  fprintf(f,"%.20e,",estif[i][j]);
	}
	fprintf(f,"],\n");
      }
      fprintf(f,"]))\n\n");

      fprintf(f,"rhs.append(([\n");
      for (i=0; i<4*8+3*19; ++i)
      {
	fprintf(f,"%.20e,",eforce[i]);
      }
      fprintf(f,"\n]))\n\n");

      fprintf(f,"vn.append(([\n");
      for (i=0; i<ele->numnp; ++i)
      {
	NODE* node = ele->node[i];
	fprintf(f,"%.20e,",node->sol_increment.a.da[ipos->velnp][0]);
	fprintf(f,"%.20e,",node->sol_increment.a.da[ipos->velnp][1]);
	fprintf(f,"%.20e,",node->sol_increment.a.da[ipos->velnp][2]);
	if (node->sol_increment.sdim > 3)
	  fprintf(f,"%.20e,",node->sol_increment.a.da[ipos->velnp][3]);
      }
      fprintf(f,"\n]))\n\n");

      fflush(f);
    }
  }
#endif

end:

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}

#endif
