/*!----------------------------------------------------------------------
\file
\brief multilevel element control routine for fluid2

<pre>
Maintainer: Volker Gravemeier
            vgravem@stanford.edu


</pre>

------------------------------------------------------------------------*/
#ifdef FLUID2_ML
#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../fluid2/fluid2_prototypes.h"
#include "fluid2ml_prototypes.h"
#include "../fluid2/fluid2.h"
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

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------- general static arrays */
static ARRAY     eveln_a;  /* element velocities at (n)                 */
static DOUBLE  **eveln;
static ARRAY     evel_a;   /* element velocities at (n+1)               */
static DOUBLE  **evel;
static ARRAY     epren_a;  /* element pressures at (n)	                */
static DOUBLE   *epren;
static ARRAY     epre_a;   /* element pressures at (n+1)                */
static DOUBLE   *epre;
static ARRAY     edeadn_a; /* element dead load at (n)                  */
static DOUBLE   *edeadn;
static ARRAY     edead_a;  /* element dead load at (n+1)                */
static DOUBLE   *edead;
static ARRAY     funct_a;  /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;  /* local shape function derivatives          */
static DOUBLE  **deriv;
static ARRAY     deriv2_a; /* 2nd local shape function derivatives      */
static DOUBLE  **deriv2;
static ARRAY     xjm_a;    /* jacobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     velint_a; /* velocities at integration point           */
static DOUBLE   *velint;
static ARRAY     vel2int_a;/* velocities at integration point           */
static DOUBLE   *vel2int;
static ARRAY	 covint_a; /* convective velocities at integr. point    */
static DOUBLE   *covint;
static ARRAY     vderxy_a; /* velocity derivatives                      */
static DOUBLE  **vderxy;
static ARRAY     pderxy_a; /* pressure derivatives                      */
static DOUBLE   *pderxy;
static ARRAY     vderxy2_a;/* 2nd velocity derivatives                  */
static DOUBLE  **vderxy2;
static ARRAY     derxy_a;  /* global shape function derivatives         */
static DOUBLE  **derxy;
static ARRAY     derxy2_a; /* 2nd global shape function derivatives     */
static DOUBLE  **derxy2;
static ARRAY     xyze_a;   /* (large-scale) element coordinates         */
static DOUBLE  **xyze;
static ARRAY     vderint_a;/* velocity derivatives at integration point */
static DOUBLE  **vderint;
static ARRAY     w1_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa1;      /* used in different element routines        */
static ARRAY     w2_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa2;      /* used in different element routines        */
static DOUBLE  **estif;    /* pointer to element stiffness matrix       */
static DOUBLE  **emass;    /* pointer to element mass matrix            */
static DOUBLE   *etforce;  /* pointer to element time rhs               */
static DOUBLE   *eiforce;  /* pointer to element iteration rhs          */
static DOUBLE   *edforce;  /* pointer to element rhs due to dirichl. bc */

/*----------------------------------- static arrays for multi-level FEM */
static ARRAY     smlme_a;  /* submesh element location matrix           */
static INT      *smlme;
static ARRAY     smitope_a;/* submesh element topology array            */
static INT      *smitope;
static ARRAY     smxyze_a; /* submesh element coordinates               */
static DOUBLE  **smxyze;
static ARRAY     smxyzep_a;/* sm element coordinates on parent domain   */
static DOUBLE  **smxyzep;
static ARRAY     smfunct_a;/* submesh shape functions                   */
static DOUBLE   *smfunct;
static ARRAY     smderiv_a;/* submesh local shape function derivatives  */
static DOUBLE  **smderiv;
static ARRAY     smderxy_a;/* submesh global shape function derivatives */
static DOUBLE  **smderxy;
static ARRAY     smderiv2_a;/* sm 2nd local shape function derivatives  */
static DOUBLE  **smderiv2;
static ARRAY     smderxy2_a;/* sm 2nd global shape function derivatives */
static DOUBLE  **smderxy2;
static ARRAY     smxjm_a;  /* submesh jacobian matrix                   */
static DOUBLE  **smxjm;
static ARRAY     smestif_a;/* submesh element stiffness matrix          */
static DOUBLE  **smestif;
static ARRAY     smemass_a;/* submesh element mass matrix               */
static DOUBLE  **smemass;
static ARRAY     smevfor_a;/* submesh element vmm rhs                   */
static DOUBLE  **smevfor;
static ARRAY     smetfor_a;/* submesh element time rhs                  */
static DOUBLE  **smetfor;

static ARRAY     evbub_a;  /* submesh element velocity bubble functions */
static DOUBLE  **evbub;
static ARRAY     epbub_a;  /* submesh element pressure bubble functions */
static DOUBLE  **epbub;
static ARRAY     efbub_a;  /* submesh element rhs bubble functions      */
static DOUBLE  **efbub;
static ARRAY     evbubn_a; /* sm elem. velocity bubble functions at (n) */
static DOUBLE  **evbubn;
static ARRAY     epbubn_a; /* sm elem. pressure bubble functions at (n) */
static DOUBLE  **epbubn;
static ARRAY     efbubn_a; /* sm elem. rhs bubble functions at (n)      */
static DOUBLE  **efbubn;
static ARRAY     vbubint_a;/* velocity bubble functions at int. point   */
static DOUBLE   *vbubint;
static ARRAY     vbubderxy_a;/* vel. bubble funct. deriv. at int. p.    */
static DOUBLE  **vbubderxy;
static ARRAY     vbubderxy2_a;/* 2nd vel. bub. funct. deriv. at int. p. */
static DOUBLE  **vbubderxy2;
static ARRAY     pbubint_a;/* pressure bubble functions at int. point   */
static DOUBLE  **pbubint;
static ARRAY4D   pbubderxy_a;/* pre. bubble funct. deriv. at int. p.    */
static DOUBLE ***pbubderxy;
static ARRAY4D   pbubderxy2_a;/* 2nd pre. bub. funct. deriv. at int. p. */
static DOUBLE ***pbubderxy2;
static ARRAY     vbubintn_a;/* vel. bubble funct. at int. point at (n)  */
static DOUBLE   *vbubintn;
static ARRAY     vbubderxyn_a;/* vel. bub. fun. der. at int. p. at (n)  */
static DOUBLE  **vbubderxyn;
static ARRAY     vbubderxy2n_a;/* 2nd vel. bub. fun. der. at i.p. at (n)*/
static DOUBLE  **vbubderxy2n;
static ARRAY     pbubintn_a;/* pre. bubble funct. at int. point at (n)  */
static DOUBLE  **pbubintn;
static ARRAY4D   pbubderxyn_a;/* pre. bub. fun. der. at int. p. at (n)  */
static DOUBLE ***pbubderxyn;
static ARRAY4D   pbubderxy2n_a;/* 2nd pre. bub. fun. der. at i.p. at (n)*/
static DOUBLE ***pbubderxy2n;
static ARRAY     velintn_a; /* velocities at integration point at (n)   */
static DOUBLE   *velintn;
static ARRAY     velintnt_a;/* 'temporal' velocities at int. p. at (n)  */
static DOUBLE   *velintnt;
static ARRAY     velintnc_a;/* 'convective' velocit. at int. p. at (n)  */
static DOUBLE   *velintnc;
static ARRAY     covintn_a; /* convective velocities at int. p. at (n)  */
static DOUBLE   *covintn;
static ARRAY     vderxyn_a; /* velocity derivatives at int. p. at (n)   */
static DOUBLE  **vderxyn;
static ARRAY     vderxync_a;/* 'convective' vel. deriv. at i. p. at (n) */
static DOUBLE  **vderxync;
static ARRAY     vderxynv_a;/* 'viscous' vel. deriv. at i. p. at (n)    */
static DOUBLE  **vderxynv;
static ARRAY     vderxy2n_a;/* 2nd vel. derivatives at int. p. at (n)   */
static DOUBLE  **vderxy2n;
static ARRAY     vderxy2nv_a;/* 2nd 'viscous' vel. der. at i. p. at (n) */
static DOUBLE  **vderxy2nv;
static ARRAY     pderxyn_a; /* pressure derivatives at int. p. at (n)   */
static DOUBLE   *pderxyn;
static ARRAY     smvelint_a;/* small-scale velocities at int. point     */
static DOUBLE   *smvelint;
static ARRAY     smvderxy_a;/* small-scale velocity deriv. at int. p.   */
static DOUBLE  **smvderxy;
static ARRAY     smpreint_a;/* small-scale 'pressures' at int. point    */
static DOUBLE   *smpreint;
static ARRAY     smpderxy_a;/* small-scale 'pressure' deriv. at int. p. */
static DOUBLE  **smpderxy;
static ARRAY     smvelintn_a;/* small-scale velocities at i. p. at (n)  */
static DOUBLE   *smvelintn;
static ARRAY     smvderxyn_a;/* small-scale vel. deriv. at i. p. at (n) */
static DOUBLE  **smvderxyn;
static ARRAY     smvderxy2n_a;/* 2nd s-s vel. deriv. at i. p. at (n)    */
static DOUBLE  **smvderxy2n;
static ARRAY     smpreintn_a;/* small-scale 'pressures' at i. p. at (n) */
static DOUBLE   *smpreintn;
static ARRAY     smpderxyn_a;/* small-scale 'pre.' der. at i. p. at (n) */
static DOUBLE  **smpderxyn;
static ARRAY     smpderxy2n_a;/* 2nd s-s 'pre.' der. at i. p. at (n)    */
static DOUBLE  **smpderxy2n;
static ARRAY     smfint_a;  /* small-scale 'rhs' at int. point          */
static DOUBLE   *smfint;
static ARRAY     smfderxy_a;/* small-scale 'rhs' deriv. at int. point   */
static DOUBLE  **smfderxy;
static ARRAY     smfintn_a; /* small-scale 'rhs' at int. point at (n)   */
static DOUBLE   *smfintn;
static ARRAY     smfderxyn_a;/* small-scale 'rhs' der. at i. p. at (n)  */
static DOUBLE  **smfderxyn;
static ARRAY     smfderxy2n_a;/* 2nd s-s 'rhs' der. at i. p. at (n)     */
static DOUBLE  **smfderxy2n;

static ARRAY     smiediff_a;/* sm element lhs integral (diffusive part) */
static DOUBLE  **smiediff;
static ARRAY     smierhs_a; /* sm element rhs integral                  */
static DOUBLE   *smierhs;

static ARRAY     sslme_a;  /* sub-submesh element location matrix       */
static INT      *sslme;
static ARRAY     ssitope_a;/* sub-submesh element topology array        */
static INT      *ssitope;
static ARRAY     ssxyze_a; /* sub-submesh element coordinates           */
static DOUBLE  **ssxyze;
static ARRAY     ssxyzep_a;/* ssm element coordinates on parent domain  */
static DOUBLE  **ssxyzep;
static ARRAY     ssfunct_a;/* sub-submesh shape functions               */
static DOUBLE   *ssfunct;
static ARRAY     ssderiv_a;/* sub-submesh local shape function deriv.   */
static DOUBLE  **ssderiv;
static ARRAY     ssderxy_a;/* sub-submesh global shape function deriv.  */
static DOUBLE  **ssderxy;
static ARRAY     ssderiv2_a;/* ssm 2nd local shape function deriv.      */
static DOUBLE  **ssderiv2;
static ARRAY     ssderxy2_a;/* ssm 2nd global shape function deriv.     */
static DOUBLE  **ssderxy2;
static ARRAY     ssxjm_a;  /* sub-submesh jacobian matrix               */
static DOUBLE  **ssxjm;
static ARRAY     ssestif_a;/* sub-submesh element stiffness matrix      */
static DOUBLE  **ssestif;
static ARRAY     ssenfor_a;/* sub-submesh element normalized rhs        */
static DOUBLE   *ssenfor;

static ARRAY     ebub_a;   /* ssm normalized element bubble functions   */
static DOUBLE   *ebub;

static DOUBLE  **smmat;    /* pointer to global submesh matrix          */
static DOUBLE  **smrhs;    /* pointer to global submesh rhs             */
static DOUBLE  **ssmat;    /* pointer to global sub-submesh matrix  	*/
static DOUBLE   *ssrhs;    /* pointer to global sub-submesh rhs		*/

/*!---------------------------------------------------------------------
\brief control routine for large-scale element integration of fluid2

<pre>                                                       gravem 07/03

This routine controls the element evaluation of the large-scale element:
-actual vel. and pres. variables are set
-for 3-level: control routine for dynamic subgrid viscosity is called
-the control routine for the small-scale solution is called
-the small-scale problem is solved
-the small-scale bubble functions are integrated for this l-s element
-for additional stabilization: stabilization parameters are calculated
-element integration is performed --> element stiffness matrix and
                                  --> element load vectors
-stiffness matrix and load vectors are permuted for assembling
-element load vector due to dirichlet conditions is calculated

</pre>
\param  *data	         FLUID_DATA     (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ssmesh	         FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual element
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *etforce_global  ARRAY	        (o)   element time force
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param  *hasdirich       int	        (o)   element flag
\param  *hasext          int	        (o)   element flag
\param   init	         int	        (i)   init flag
\return void

------------------------------------------------------------------------*/
void f2_lsele(FLUID_DATA     *data,
              FLUID_DYN_ML   *mlvar,
              FLUID_ML_SMESH *submesh,
              FLUID_ML_SMESH *ssmesh,
	      ELEMENT	     *ele,
              ARRAY	     *estif_global,
              ARRAY	     *emass_global,
	      ARRAY	     *etforce_global,
	      ARRAY	     *eiforce_global,
	      ARRAY	     *edforce_global,
	      INT	     *hasdirich,
              INT	     *hasext,
	      INT	      init)
{
INT              hasdead;
INT              info=0;
INT              infrhs=0;
INT              i,j;
DOUBLE           matvec[7000],rhsvec[2000];

#ifdef DEBUG
dstrc_enter("f2_lsele");
#endif

if (init==1) /* allocate working arrays and set pointers */
{
   fdyn = alldyn[genprob.numff].fdyn;

   eveln   = amdef("eveln"  ,&eveln_a  ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   evel    = amdef("evel"   ,&evel_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   epren   = amdef("epren"  ,&epren_a  ,MAXNOD_F2,1,"DV");
   epre    = amdef("epre"   ,&epre_a   ,MAXNOD_F2,1,"DV");
   edeadn  = amdef("edeadn" ,&edeadn_a ,2,1,"DV");
   edead   = amdef("edead"  ,&edead_a  ,2,1,"DV");
   funct   = amdef("funct"  ,&funct_a  ,MAXNOD_F2,1,"DV");
   deriv   = amdef("deriv"  ,&deriv_a  ,2,MAXNOD_F2,"DA");
   deriv2  = amdef("deriv2" ,&deriv2_a ,3,MAXNOD_F2,"DA");
   xjm     = amdef("xjm"    ,&xjm_a    ,2,2        ,"DA");
   velint  = amdef("velint" ,&velint_a ,NUM_F2_VELDOF,1,"DV");
   vel2int = amdef("vel2int",&vel2int_a,NUM_F2_VELDOF,1,"DV");
   covint  = amdef("covint" ,&covint_a ,NUM_F2_VELDOF,1,"DV");
   vderxy  = amdef("vderxy" ,&vderxy_a ,2,2,"DA");
   pderxy  = amdef("pderxy" ,&pderxy_a ,2,1,"DV");
   vderxy2 = amdef("vderxy2",&vderxy2_a,2,3,"DA");
   derxy   = amdef("derxy"  ,&derxy_a  ,2,MAXNOD_F2,"DA");
   derxy2  = amdef("derxy2" ,&derxy2_a ,3,MAXNOD_F2,"DA");
   xyze    = amdef("xyze"   ,&xyze_a   ,2,MAXNOD_F2,"DA");
   vderint = amdef("vderint",&vderint_a,4,MAXGAUSS ,"DA");
   wa1     = amdef("wa1"    ,&w1_a     ,30,30,"DA");
   wa2     = amdef("wa2"    ,&w2_a     ,30,30,"DA");
/*                                        \- size is arbitrary chosen!  */
   estif   = estif_global->a.da;
   emass   = emass_global->a.da;
   eiforce = eiforce_global->a.dv;
   etforce = etforce_global->a.dv;
   edforce = edforce_global->a.dv;

/*------------------ allocation for submesh of this large-scale element */
   f2_smele(data,mlvar,submesh,ele,init);
/*-------------- allocation for sub-submesh of this large-scale element */
   if (mlvar->smsgvi>2) f2_ssele(data,mlvar,ssmesh,ele,init);
   goto end;
} /* endif (init==1) */

/*------------------------------------------------ initialize with ZERO */
amzero(estif_global);
amzero(emass_global);
amzero(eiforce_global);
amzero(etforce_global);
amzero(edforce_global);
*hasdirich=0;
*hasext=0;

/*---------------------------- set element data for large-scale element */
f2_lsset(ele,eveln,evel,epren,epre,edeadn,edead,hasext);

/*------------------------------- dynamic subgrid viscosity calculation */
if (mlvar->smsgvi>2) f2_dynsgv(data,mlvar,submesh,ssmesh,ele);

/*------------------ initialize submesh global matrix and rhs with ZERO */
amzero(&(submesh->mat));
amzero(&(submesh->rhs));
/*------------ small-scale solution on submesh: element control routine */
f2_smele(data,mlvar,submesh,ele,init);

/*fluid_prgmr(smmat,smrhs,submesh->numeq,mlvar->nelbub); */
/*--------------------- change order in matrix and rhs array for solver */
for (i=0; i<submesh->numeq; i++)
{
  for (j=0; j<submesh->numeq; j++)
  {
    matvec[info] = submesh->mat.a.da[j][i];
    info++;
  }
}
for (i=0; i<mlvar->nelbub; i++)
{
  for (j=0; j<submesh->numeq; j++)
  {
    rhsvec[infrhs] = submesh->rhs.a.da[j][i];
    infrhs++;
  }
}
/*-------------------------------------------- solve small-scale system */
dgesv(&(submesh->numeq),&(mlvar->nelbub),&(matvec[0]),&(submesh->numeq),
      submesh->ipiv.a.iv,&(rhsvec[0]),&(submesh->numeq),&info);
/*---------------------------------------- error in solution procedure? */
if (info<0)      dserror("illegal value in small-scale solution");
else if (info>0) dserror("small-scale solution could not be computed");

/*----------------------------------- re-change order in solution array */
infrhs=0;
for (i=0; i<mlvar->nelbub; i++)
{
for (j=0; j<submesh->numeq; j++)
{
   submesh->rhs.a.da[j][i] = rhsvec[infrhs];
   infrhs++;
}
}
/*fluid_prgmr(smmat,smrhs,submesh->numeq,mlvar->nelbub);*/
/*f2_cgsbub(submesh,smrhs,mlvar->nelbub);*/

/*-------------------------- copy small-scale solution to element array */
f2_smcopy(smrhs,ele,submesh->numeq,mlvar->nelbub);

/*----------------- calculate small-scale part of element matrices and
                               element force vectors (bubble functions) */
f2_bubele(data,mlvar,submesh,ele);

/*- calculate charact. l-s element length and stab. param. if necessary */
dsassert(ele->e.f2->stab_type == stab_gls, /* check for proper stabilisation */
     "wrong stabilisation in mulitlevel case");

f2_mlcalelesize(ele,data,funct,deriv,deriv2,derxy,xjm,evel,velint,
                  vderxy,wa1);
/*------------------ calculate large-scale part of element matrices and
                                                  element force vectors */
f2_lsint(data,ele,mlvar,hasext,estif,emass,eiforce,etforce,funct,
         deriv,deriv2,xjm,derxy,derxy2,evel,eveln,epren,edeadn,edead,
	 velint,velintn,covint,covintn,vderxy,vderxyn,wa1,wa2);

/*----------------- add emass and estif to estif and permute the matrix */
f2_mlpermestif(estif,emass,wa1,ele->numnp);

/*--------------------------------- permute element load vector etforce */
if (fdyn->nif!=0) f2_permeforce(etforce,wa1,ele->numnp);

/*--------------------------------- permute element load vector eiforce */
if (fdyn->nii+(*hasext)!=0) f2_permeforce(eiforce,wa1,ele->numnp);

/*------------------------------- calculate element load vector edforce */
fluid_mlcaldirich(ele,edforce,estif,hasdirich);

/*----------------------------------------------------- local co-system */
dsassert(ele->locsys==locsys_no,"locsys not implemented for this element!\n");

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_lsele */

/*!---------------------------------------------------------------------
\brief control routine for submesh element integration of fluid2

<pre>                                                       gravem 07/03

This routine controls the element evaluation of the submesh element:
-element data is set
-stabilization parameter or subgrid viscosity is calculated if necessary
-element integration is performed --> element stiffness matrix and
                                  --> element load vectors

</pre>
\param  *data	         FLUID_DATA     (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual large-scale element
\param   init	         int	        (i)   init flag
\return void

------------------------------------------------------------------------*/
void f2_smele(FLUID_DATA     *data,
              FLUID_DYN_ML   *mlvar,
              FLUID_ML_SMESH *submesh,
	      ELEMENT        *ele,
              INT             init)
{
INT              iele;    /* submesh element counter                   */

#ifdef DEBUG
dstrc_enter("f2_smele");
#endif

if (init==1) /* allocate working arrays and set pointers */
{
   fdyn = alldyn[genprob.numff].fdyn;

   smlme    = amdef("smlme"   ,&smlme_a   ,submesh->numen,1,"IV");
   smitope  = amdef("smitope" ,&smitope_a ,submesh->numen,1,"IV");
   smxyze   = amdef("smxyze"  ,&smxyze_a  ,2,submesh->numen,"DA");
   smxyzep  = amdef("smxyzep" ,&smxyzep_a ,2,submesh->numen,"DA");
   smfunct  = amdef("smfunct" ,&smfunct_a ,submesh->numen,1,"DV");
   smderiv  = amdef("smderiv" ,&smderiv_a ,2,submesh->numen,"DA");
   smderxy  = amdef("smderxy" ,&smderxy_a ,2,submesh->numen,"DA");
   smderiv2 = amdef("smderiv2",&smderiv2_a,3,submesh->numen,"DA");
   smderxy2 = amdef("smderxy2",&smderxy2_a,3,submesh->numen,"DA");
   smxjm    = amdef("smxjm"   ,&smxjm_a   ,2,2             ,"DA");

   smestif = amdef("smestif",&smestif_a,submesh->numen,submesh->numen,"DA");
   smemass = amdef("smemass",&smemass_a,submesh->numen,submesh->numen,"DA");
   smevfor = amdef("smevfor",&smevfor_a,submesh->numen,mlvar->nelbub ,"DA");
   smetfor = amdef("smetfor",&smetfor_a,submesh->numen,mlvar->nelbub ,"DA");

   evbub  = amdef("evbub" ,&evbub_a ,submesh->numen,mlvar->nvbub,"DA");
   epbub  = amdef("epbub" ,&epbub_a ,submesh->numen,mlvar->npbub,"DA");
   efbub  = amdef("efbub" ,&efbub_a ,submesh->numen,2		,"DA");
   evbubn = amdef("evbubn",&evbubn_a,submesh->numen,mlvar->nvbub,"DA");
   epbubn = amdef("epbubn",&epbubn_a,submesh->numen,mlvar->npbub,"DA");
   efbubn = amdef("efbubn",&efbubn_a,submesh->numen,2		,"DA");

   vbubint    = amdef("vbubint"    ,&vbubint_a    ,MAXNOD_F2,1,"DV");
   vbubderxy  = amdef("vbubderxy"  ,&vbubderxy_a  ,2,MAXNOD_F2,"DA");
   vbubderxy2 = amdef("vbubderxy2" ,&vbubderxy2_a ,3,MAXNOD_F2,"DA");

   pbubint    = amdef ("pbubint"    ,&pbubint_a    ,2,MAXNOD_F2,"DA");
   pbubderxy  = am4def("pbubderxy"  ,&pbubderxy_a  ,2,2,MAXNOD_F2,0,"D3");
   pbubderxy2 = am4def("pbubderxy2" ,&pbubderxy2_a ,3,2,MAXNOD_F2,0,"D3");

   vbubintn    = amdef("vbubintn"    ,&vbubintn_a    ,MAXNOD_F2,1,"DV");
   vbubderxyn  = amdef("vbubderxyn"  ,&vbubderxyn_a  ,2,MAXNOD_F2,"DA");
   vbubderxy2n = amdef("vbubderxy2n" ,&vbubderxy2n_a ,3,MAXNOD_F2,"DA");

   pbubintn    = amdef("pbubintn"    ,&pbubintn_a    ,2,MAXNOD_F2,"DA");
   pbubderxyn  = am4def("pbubderxyn"  ,&pbubderxyn_a  ,2,2,MAXNOD_F2,0,"D3");
   pbubderxy2n = am4def("pbubderxy2n" ,&pbubderxy2n_a ,3,2,MAXNOD_F2,0,"D3");

   velintn  = amdef("velintn"  ,&velintn_a  ,2,1,"DV");
   velintnt = amdef("velintnt" ,&velintnt_a ,2,1,"DV");
   velintnc = amdef("velintnc" ,&velintnc_a ,2,1,"DV");
   covintn  = amdef("covintn"  ,&covintn_a  ,2,1,"DV");

   vderxyn   = amdef("vderxyn"  ,&vderxyn_a  ,2,2,"DA");
   vderxync  = amdef("vderxync" ,&vderxync_a ,2,2,"DA");
   vderxynv  = amdef("vderxynv" ,&vderxynv_a ,2,2,"DA");
   vderxy2n  = amdef("vderxy2n" ,&vderxy2n_a ,2,3,"DA");
   vderxy2nv = amdef("vderxy2nv",&vderxy2nv_a,2,3,"DA");
   pderxyn   = amdef("pderxyn"  ,&pderxyn_a  ,2,1,"DV");

   smvelint  = amdef("smvelint" ,&smvelint_a ,2,1,"DV");
   smvderxy  = amdef("smvderxy" ,&smvderxy_a ,2,2,"DA");
   smpreint  = amdef("smpreint" ,&smpreint_a ,2,1,"DV");
   smpderxy  = amdef("smpderxy" ,&smpderxy_a ,2,2,"DA");

   smvelintn  = amdef("smvelintn" ,&smvelintn_a ,2,1,"DV");
   smvderxyn  = amdef("smvderxyn" ,&smvderxyn_a ,2,2,"DA");
   smvderxy2n = amdef("smvderxy2n",&smvderxy2n_a,2,3,"DA");
   smpreintn  = amdef("smpreintn" ,&smpreintn_a ,2,1,"DV");
   smpderxyn  = amdef("smpderxyn" ,&smpderxyn_a ,2,2,"DA");
   smpderxy2n = amdef("smpderxy2n",&smpderxy2n_a,2,3,"DA");

   smfint   = amdef("smfint"  ,&smfint_a  ,2,1,"DV");
   smfderxy = amdef("smfderxy",&smfderxy_a,2,2,"DA");

   smfintn    = amdef("smfintn"   ,&smfintn_a   ,2,1,"DV");
   smfderxyn  = amdef("smfderxyn" ,&smfderxyn_a ,2,2,"DA");
   smfderxy2n = amdef("smfderxy2n",&smfderxy2n_a,2,3,"DA");

   smiediff = amdef("smiediff",&smiediff_a,submesh->numen,submesh->numen,"DA");
   smierhs  = amdef("smierhs" ,&smierhs_a ,submesh->numen,1             ,"DV");

   smmat  = submesh->mat.a.da;
   smrhs  = submesh->rhs.a.da;
   goto end;
} /* endif (init==1) */

for (iele=0; iele<submesh->numele; iele++)/* loop over submesh elements */
{
/*---------------------------------------------------- set element data */
  f2_smset(submesh,ele,smlme,smitope,smxyze,smxyzep,iele,0);

/*----------------- set iterative bubble functions at current time step */
  if (mlvar->convel==0)
    f2_bubset(mlvar,submesh,ele,smlme,evbub,epbub,efbub,0);
/*------------------------------- set bubble functions at time step (n) */
  if (fdyn->nis==0 && mlvar->quastabub==0)
    f2_bubset(mlvar,submesh,ele,smlme,evbubn,epbubn,efbubn,1);

/*--- calculate charact. element length and stab. param. / subgr. visc. */
  if (mlvar->smstabi>0 || mlvar->smsgvi==1 || mlvar->smsgvi==2)
    f2_smelesize(ele,data,mlvar,funct,deriv,deriv2,smfunct,smderiv,
                 smderiv2,derxy,xjm,evel,velint,vderxy,smxyze,smxyzep,wa1);

/*------------ initialize submesh global matrices and vectors with ZERO */
  amzero(&smestif_a);
  amzero(&smemass_a);
  amzero(&smevfor_a);
  amzero(&smetfor_a);

/*-------- calculate submesh element matrices and element force vectors */
  f2_smint(data,ele,mlvar,submesh,smestif,smemass,smevfor,smetfor,
	   smxyze,smxyzep,funct,deriv,deriv2,xjm,derxy,derxy2,smfunct,
	   smderiv,smderiv2,smxjm,smderxy,smderxy2,eveln,evel,epren,epre,
	   evbub,epbub,efbub,evbubn,epbubn,efbubn,vbubint,vbubderxy,
	   vbubderxy2,pbubint,pbubderxy,pbubderxy2,vbubintn,vbubderxyn,
	   vbubderxy2n,pbubintn,pbubderxyn,pbubderxy2n,velint,velintn,
	   velintnt,velintnc,vderxy,vderxyn,vderxync,vderxynv,vderxy2,
	   vderxy2n,vderxy2nv,pderxyn,smvelint,smvderxy,smpreint,smpderxy,
	   smvelintn,smvderxyn,smvderxy2n,smpreintn,smpderxyn,smpderxy2n,
	   smfint,smfderxy,smfintn,smfderxyn,smfderxy2n,wa1,wa2);

/*------------- add element stiffness matrix to global stiffness matrix */
  if (mlvar->quastabub!=0)
    fluid_add_smat(smmat,smestif,submesh->numen,smlme,fdyn->dta);
  else fluid_add_smat(smmat,smestif,submesh->numen,smlme,fdyn->thsl);

/*----------------------- add element mass matrix to global mass matrix */
  if (fdyn->nis==0 && mlvar->quastabub==0)
    fluid_add_smat(smmat,smemass,submesh->numen,smlme,ONE);

/*------------------------------- add element vmm rhs to global vmm rhs */
  fluid_add_smrhs(smrhs,smevfor,mlvar->nelbub,submesh->numen,smlme);

/*----------------------------- add element time rhs to global time rhs */
  if (fdyn->nis==0)
    fluid_add_smrhs(smrhs,smetfor,mlvar->nelbub,submesh->numen,smlme);
}

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_smele */

/*!---------------------------------------------------------------------
\brief control routine for submesh bubble function integration of fluid2

<pre>                                                       gravem 07/03

This routine controls the element integration of bubble functions on
the submesh element:
-element data is set
-stabilization parameter or subgrid viscosity is calculated if necessary
-element integration is performed --> element stiffness matrix and
                                  --> element load vectors

</pre>
\param  *data	         FLUID_DATA     (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual large-scale element
\return void

------------------------------------------------------------------------*/
void f2_bubele(FLUID_DATA     *data,
               FLUID_DYN_ML   *mlvar,
               FLUID_ML_SMESH *submesh,
	       ELEMENT        *ele)
{
INT              iele;    /* submesh element counter                   */

#ifdef DEBUG
dstrc_enter("f2_bubele");
#endif

for (iele=0; iele<submesh->numele; iele++)/* loop over submesh elements */
{
/*---------------------------------------------------- set element data */
  f2_smset(submesh,ele,smlme,smitope,smxyze,smxyzep,iele,0);

/*----------------- set iterative bubble functions at current time step */
  f2_bubset(mlvar,submesh,ele,smlme,evbub,epbub,efbub,0);

/*---------- calculate submesh bubble matrices and bubble force vectors */
  f2_bubint(data,ele,mlvar,submesh,estif,emass,eiforce,smxyze,
            smxyzep,funct,deriv,deriv2,xjm,derxy,smfunct,smderiv,smderiv2,
	    smxjm,smderxy,evel,epre,evbub,epbub,efbub,vbubint,vbubderxy,
	    pbubint,pbubderxy,covint,velint,vderxy,smvelint,smvderxy,
	    smpreint,smpderxy,smfint,smfderxy,wa1,wa2);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_bubele */

/*!---------------------------------------------------------------------
\brief control routine for dynamic calc. of subgrid viscosity for fluid2

<pre>                                                       gravem 07/03

This routine controls the dynamic calculation of the subgrid viscosity
on the large-scale element:
-element data is set
-stabilization parameter or subgrid viscosity is calculated if necessary
-element integration is performed --> element stiffness matrix and
                                  --> element load vectors

</pre>
\param  *data	         FLUID_DATA     (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ssmesh	         FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual large-scale element
\return void

------------------------------------------------------------------------*/
void f2_dynsgv(FLUID_DATA     *data,
               FLUID_DYN_ML   *mlvar,
               FLUID_ML_SMESH *submesh,
               FLUID_ML_SMESH *ssmesh,
	       ELEMENT	      *ele)
{
INT              i,j,ite;    /* just some counters                      */
INT              info=0;     /* info variable for solver                */
INT              itemax=25;  /* maximum number of iterations            */
INT              actmat;     /* material number of the element          */
INT              nrhs=1;     /* number of rhs                           */
DOUBLE           visc;       /* physical viscosity                      */
DOUBLE           smidiff=0.0;/* sm global lhs integral (diffusive part) */
DOUBLE           smirhs=0.0; /* sm global rhs integral                  */
DOUBLE           sgtol;      /* tolerance for subgr. visc. iteration    */
DOUBLE           sgdiff;     /* iterative difference in subgr. visc.    */
DOUBLE           sgvisc;     /* subgrid viscosity                       */
DOUBLE           ssinbu;     /* ssm integral of normalized bubble fun.  */
DOUBLE           matvec[7000],rhsvec[2000];

#ifdef DEBUG
dstrc_enter("f2_dynsgv");
#endif

/*------------------------------------------------- set some parameters */
actmat = ele->mat-1;
visc   = mat[actmat].m.fluid->viscosity;
sgtol  = visc/TEN;
mlvar->smsgvisc = ZERO;

/*-- calculate global lhs (diffusive part) and rhs integrals on submesh */
f2_smintele(data,submesh,ele,&smidiff,&smirhs);

/*-----------------------------------------------------------------------
   iterative loop for dynamic subgrid viscosity calculation
----------------------------------------------------------------------- */
for (ite=0;ite<itemax;ite++)
{
/*------------------------------------- set iterative subgrid viscosity */
  sgvisc = mlvar->smsgvisc;
/* initialize integral of normalized bubble fun. on sub-submesh to ZERO */
  ssinbu = ZERO;

/*---------------- initialize sub-submesh global matrix and rhs to ZERO */
  amzero(&(ssmesh->mat));
  amzero(&(ssmesh->rhs));
/*-------------------- solution on sub-submesh: element control routine */
  f2_ssele(data,mlvar,ssmesh,ele,0);

/*----------------------------- change order in matrix array for solver */
  for (i=0; i<ssmesh->numeq; i++)
  {
    for (j=0; j<ssmesh->numeq; j++)
    {
      matvec[info] = ssmesh->mat.a.da[j][i];
      info++;
    }
  }
/*-------------------------------------------------------- solve system */
  dgesv(&(ssmesh->numeq),&nrhs,&(matvec[0]),&(ssmesh->numeq),ssmesh->ipiv.a.iv,
        ssrhs,&(ssmesh->numeq),&info);
/*---------------------------------------- error in solution procedure? */
  if (info<0)      dserror("illegal value in small-scale solution");
  else if (info>0) dserror("small-scale solution could not be computed");

/*----- calculate integral of normalized bubble function on sub-submesh */
  f2_ssintele(data,ssmesh,ele,&ssinbu);

/*------------------------ calculate updated value of subgrid viscosity */
  mlvar->smsgvisc = smirhs*smirhs/smidiff/ssinbu - visc;
/*--- if subgrid viscosity < 0: subgrid viscosity = 0 -> stop iteration */
  if (mlvar->smsgvisc<ZERO) goto end;
/*--------------------------- iterative difference in subgrid viscosity */
  sgdiff = FABS(mlvar->smsgvisc-sgvisc);
/*------------------------------------------- check iteration tolerance */
  if (sgdiff<sgtol) goto end;
/*---------------------------------- check maximum number of iterations */
  else if (ite==itemax)
  {
    printf("|            |                   |              |             | \n");
    printf("| >>>>>> dynamic subgrid visc. not converged in itemax steps! | \n");
    printf("|            |                   |              |             | \n");
  }
}

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_dynsgv */

/*!---------------------------------------------------------------------
\brief submesh global lhs and rhs element integration for fluid2

<pre>                                                       gravem 07/03

This routine controls the elementwise integration of the global lhs
(only the diffusive part) and rhs on the submesh.

</pre>
\param  *data	         FLUID_DATA     (i)
\param  *submesh	 FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)  actual large-scale element
\param  *smidiff	 DOUBLE	        (o)  submesh global lhs integral
\param  *smirhs  	 DOUBLE	        (o)  submesh global rhs integral
\return void

------------------------------------------------------------------------*/
void f2_smintele(FLUID_DATA     *data,
                 FLUID_ML_SMESH *submesh,
	         ELEMENT        *ele,
	         DOUBLE         *smidiff,
	         DOUBLE         *smirhs)
{
INT              iele;    /* submesh element counter                   */

#ifdef DEBUG
dstrc_enter("f2_smintele");
#endif

for (iele=0; iele<submesh->numele; iele++)/* loop over submesh elements */
{
/*---------------------------------------------------- set element data */
  f2_smset(submesh,ele,smlme,smitope,smxyze,smxyzep,iele,0);

/*------------------- initialize element lhs and rhs integral with ZERO */
  amzero(&smiediff_a);
  amzero(&smierhs_a);
/*---------------------- calculate submesh element lhs and rhs integral */
  f2_smint2(data,ele,submesh,smiediff,smierhs,smxyze,smfunct,smderiv,
            smderiv2,smxjm,smderxy);

/*--------------------- add element lhs integral to global lhs integral */
  fluid_add_intlhs(smidiff,smiediff,submesh->numen,smlme);

/*--------------------- add element rhs integral to global rhs integral */
  fluid_add_intrhs(smirhs,smierhs,submesh->numen,smlme);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_smintele */

/*!---------------------------------------------------------------------
\brief control routine for sub-submesh element integration for fluid2

<pre>                                                       gravem 07/03

This routine controls the element evaluation of the sub-submesh element:
-element data is set
-element integration is performed --> element stiffness matrix and
                                  --> element load vectors

</pre>
\param  *data	         FLUID_DATA     (i)
\param  *mlvar	         FLUID_DYN_ML   (i)
\param  *ssmesh	         FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)   actual large-scale element
\param   init	         int	        (i)   init flag
\return void

------------------------------------------------------------------------*/
void f2_ssele(FLUID_DATA      *data,
              FLUID_DYN_ML    *mlvar,
              FLUID_ML_SMESH  *ssmesh,
	      ELEMENT         *ele,
	      INT              init)
{
INT              iele;    /* sub-submesh element counter               */

#ifdef DEBUG
dstrc_enter("f2_ssele");
#endif

if (init==1) /* allocate working arrays and set pointers */
{
   sslme    = amdef("sslme"   ,&sslme_a   ,ssmesh->numen,1,"IV");
   ssitope  = amdef("ssitope" ,&ssitope_a ,ssmesh->numen,1,"IV");
   ssxyze   = amdef("ssxyze"  ,&ssxyze_a  ,2,ssmesh->numen,"DA");
   ssxyzep  = amdef("ssxyzep" ,&ssxyzep_a ,2,ssmesh->numen,"DA");
   ssfunct  = amdef("ssfunct" ,&ssfunct_a ,ssmesh->numen,1,"DV");
   ssderiv  = amdef("ssderiv" ,&ssderiv_a ,2,ssmesh->numen,"DA");
   ssderxy  = amdef("ssderxy" ,&ssderxy_a ,2,ssmesh->numen,"DA");
   ssderiv2 = amdef("ssderiv2",&ssderiv2_a,3,ssmesh->numen,"DA");
   ssderxy2 = amdef("ssderxy2",&ssderxy2_a,3,ssmesh->numen,"DA");
   ssxjm    = amdef("ssxjm"   ,&ssxjm_a   ,2,2             ,"DA");

   ssestif = amdef("ssestif",&ssestif_a,ssmesh->numen,ssmesh->numen,"DA");
   ssenfor = amdef("ssenfor",&ssenfor_a,ssmesh->numen,1            ,"DV");

   ebub    = amdef("ebub"   ,&ebub_a   ,ssmesh->numen,1,"DV");

   ssmat  = ssmesh->mat.a.da;
   ssrhs  = ssmesh->rhs.a.dv;
   goto end;
} /* endif (init==1) */

for (iele=0; iele<ssmesh->numele; iele++)/* loop over sub-submesh elem. */
{
/*---------------------------------------------------- set element data */
  f2_smset(ssmesh,ele,sslme,ssitope,ssxyze,ssxyzep,iele,1);

/*----------- initialize sub-submesh global matrix and vector with ZERO */
  amzero(&ssestif_a);
  amzero(&ssenfor_a);

/*------- calculate sub-submesh element matrix and element force vector */
  f2_ssint(data,ele,mlvar,ssmesh,ssestif,ssenfor,ssxyze,ssxyzep,
           funct,deriv,deriv2,xjm,derxy,ssfunct,ssderiv,ssderiv2,ssxjm,
	   ssderxy,evel,velint,vderxy);

/*--------------------------------- add element matrix to global matrix */
  fluid_add_smat(ssmat,ssestif,ssmesh->numen,sslme,ONE);

/*--------------------------------------- add element rhs to global rhs */
  fluid_add_ssrhs(ssrhs,ssenfor,ssmesh->numen,sslme);
}

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_ssele */

/*!---------------------------------------------------------------------
\brief sub-submesh element integration of normalized bubble for fluid2

<pre>                                                       gravem 07/03

This routine controls the elementwise integration of the normalized
bubble function on the sub-submesh.

</pre>
\param  *data	         FLUID_DATA     (i)
\param  *ssmesh  	 FLUID_ML_SMESH (i)
\param  *ele	         ELEMENT	(i)  actual large-scale element
\param  *ssinbu    	 DOUBLE	        (o)  sub-submesh bubble integral
\return void

------------------------------------------------------------------------*/
void f2_ssintele(FLUID_DATA     *data,
                 FLUID_ML_SMESH *ssmesh,
	         ELEMENT        *ele,
		 DOUBLE         *ssinbu)
{
INT              iele;    /* sub-submesh element counter               */

#ifdef DEBUG
dstrc_enter("f2_ssintele");
#endif

for (iele=0; iele<ssmesh->numele; iele++)/* loop over sub-submesh elem. */
{
/*------------------------------ set element normalized bubble function */
  f2_ssset(ssmesh,ele,sslme,ssitope,ssxyze,ebub,iele);

/* calculate sub-submesh element integral of normalized bubble function */
  f2_inbu(data,ele,ssmesh,ssinbu,ebub,ssxyze,ssfunct,ssderiv,ssderiv2,
          ssxjm);
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_ssintele */

#endif
