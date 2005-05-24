/*!----------------------------------------------------------------------
\file
\brief element control routine

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID3
#include "../headers/standardtypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
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
static ARRAY     edeadn_a; /* element dead load (selfweight)            */
static DOUBLE   *edeadng;
static ARRAY     edeadng_a;/* element dead load (selfweight)            */
static DOUBLE   *edeadn;
static ARRAY     funct_a;  /* shape functions				*/
static DOUBLE   *funct;
static ARRAY     deriv_a;  /* first natural derivatives 		*/
static DOUBLE  **deriv;
static ARRAY     deriv2_a; /* second natural derivatives		*/
static DOUBLE  **deriv2;
static ARRAY     xyze_a;
static DOUBLE  **xyze;
static ARRAY     xjm_a;    /* jocobian matrix				*/
static DOUBLE  **xjm;
static ARRAY     vderxy_a; /* vel - derivatives 			*/
static DOUBLE  **vderxy;
static ARRAY     pderxy_a; /* pre -derivatives  			*/
static DOUBLE   *pderxy;
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
void f3_calele(
	        ELEMENT        *ele,
                ARRAY          *estif_global,
                ARRAY          *emass_global,
	        ARRAY          *eforce_global,
		ARRAY          *edforce_global,
                ARRAY_POSITION *ipos,
		INT            *hasdirich,
		INT            *hasext,
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
   edeadn    = amdef("edeadn" ,&edeadn_a ,3,1,"DV");
   edeadng   = amdef("edeadng",&edeadng_a,3,1,"DV");
   funct     = amdef("funct"  ,&funct_a  ,MAXNOD_F3,1,"DV");
   deriv     = amdef("deriv"  ,&deriv_a  ,3,MAXNOD_F3,"DA");
   deriv2    = amdef("deriv2" ,&deriv2_a ,6,MAXNOD_F3,"DA");
   xjm       = amdef("xjm"    ,&xjm_a    ,3,3        ,"DA");
   xyze      = amdef("xyze"   ,&xyze_a   ,3,MAXNOD_F3,"DA");
   vderxy    = amdef("vderxy" ,&vderxy_a ,3,3,"DA");
   pderxy    = amdef("pderxy" ,&pderxy_a ,3,1,"DV");
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

/*------------------------------------------------ initialise with ZERO */
amzero(estif_global);
amzero(emass_global);
amzero(eforce_global);
amzero(edforce_global);
*hasdirich=0;
*hasext=0;
switch(ele->e.f3->is_ale)
{
case 0:
/*---------------------------------------------------- set element data */
   f3_calset(ele,xyze,ehist,evelng,epren,edeadn,edeadng,ipos,hasext);

/*------------------------- calculate element size and stab-parameter: */
   f3_calelesize(ele,xyze,funct,deriv,deriv2,derxy,xjm,evelng,wa1,0);

/*------------------------------- calculate element stiffness matrices */
/*                                           and element force vectors */
   f3_calint(ele,estif,emass,eforce,
             xyze,funct,deriv,deriv2,xjm,derxy,derxy2,
             evelng,vderxy,wa1,wa2);
break;
case 1:
/*---------------------------------------------------- set element data */
   f3_calseta(ele,xyze,ehist,evelng,ealecovn,
              ealecovng,egridv,epren,edeadn,edeadng,ipos,hasext);
/*------------------------- calculate element size and stab-parameter: */
   f3_calelesize(ele,xyze,funct,deriv,deriv2,derxy,xjm,evelng,wa1,0);
/*------------------------------- calculate element stiffness matrices */
/*                                           and element force vectors */
   f3_calinta(ele,estif,emass,eforce,
              xyze,funct,deriv,deriv2,xjm,derxy,derxy2,
              evelng,ealecovng,egridv,vderxy,wa1,wa2);
break;
default:
   dserror("parameter is_ale not 0 or 1!\n");
}

#ifdef PERF
  perf_begin(21);
#endif

switch(ele->e.f3->fs_on)
{
case 0: case 1: case 3: /* no or explict free surface */
   /*-------------- add emass and estif to estif and permute the matrix */
   f3_permestif(estif,emass,wa1,ele->numnp);
   /*------------------------------- permute element load vector eforce */
   if (fdyn->nii+(*hasext)!=0)
      f3_permeforce(eforce,wa1,ele->numnp);
break;
case 2: case 5: /* partitioned implict free surface */
   dsassert(ele->e.f3->is_ale!=0,"element at free surface has to be ALE!\n");
   /*-------------- add emass and estif to estif and permute the matrix */
   f3_permestif_ifs(estif,emass,wa1,ele);
   /*------------------------------- permute element load vector eforce */
   if (fdyn->nii+(*hasext)!=0)
      f3_permeforce_ifs(eforce,wa1,ele);
break;
default:
   dserror("parameter fs_on out of range!\n");
}

#ifdef PERF
  perf_end(21);
#endif

/*-------------------------------------------- local co-ordinate system */
if(ele->locsys==locsys_yes)
   locsys_trans(ele,estif,NULL,NULL,eforce);

/*------------------------------------------------ condensation of DBCs */
/* estif is in xyz* so edforce is also in xyz* (but DBCs have to be
   tranformed before condensing the dofs                                */
fluid_caldirich(ele,edforce,estif,hasdirich,ipos->velnp);

/*------------------------------------------ calculate emass * ehist ---*/
f3_massrhs(ele,emass,ehist,edeadng,eforce,hasext);

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f3_calele */


/*!---------------------------------------------------------------------
\brief control routine for stress calculation

<pre>                                                         genk 03/04

</pre>

\param     str       FLUID_STRESS   (i)    flag for stress calculation
\param     viscstr   INT            (i)    viscose stresses yes/no?
\param    *ele       ELEMENt        (i)    actual element
\param   *ipos                      (i)    node array positions
\return void

------------------------------------------------------------------------*/
void f3_stress(FLUID_STRESS  str,
               INT           viscstr,
	       ELEMENT      *ele,
               ARRAY_POSITION *ipos)
{
#ifdef D_FSI
INT       i;
INT       coupled;      /* flag for fsi interface element */
INT       iel;          /* number of nodes per element */
GNODE    *actgnode;     /* actual gnode */
#endif

#ifdef DEBUG
dstrc_enter("f3_stress");
#endif

switch(str)
{
case str_none: /* do nothing */
break;
#ifdef D_FSI
case str_fsicoupling:
/* check if fluid element is coupled to struct element */
   iel=ele->numnp;
   coupled=0;
   for (i=0;i<iel;i++)
   {
      actgnode = ele->node[i]->gnode;
      /* check if there is a coupled struct node */
      if (actgnode->mfcpnode[genprob.numsf]==NULL) continue;
      coupled=1;
      break;
   }
   if (coupled==1)
   f3_calelestress(viscstr,ele,eveln,epren,funct,
                   deriv,derxy,vderxy,xjm,wa1,xyze,sigmaint,ipos);
break;
#endif
case str_liftdrag:
case str_all:
   f3_calelestress(viscstr,ele,eveln,epren,funct,
                  deriv,derxy,vderxy,xjm,wa1,xyze,sigmaint,ipos);
break;
default:
  dserror("stress calculation not possible!\n");
}

  /*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of f3_stress */

/*!---------------------------------------------------------------------
\brief control routine for heightfunction

<pre>                                                       genk 04/04

evaluation of heightfunction

</pre>
\param   *ele              ELEMENT           the actual element
\param   *estif_global     DOUBLE            element stiffness matrix
\param   *eforce_global    DOUBLE            ele RHS
\param   *container        CONTAINER
\param   *ipos                           (i) node array positions
\return void

------------------------------------------------------------------------*/
void f3_heightfunc(
                   ELEMENT              *ele,
                   ARRAY                *estif_global,
		   ARRAY                *eforce_global,
		   CONTAINER            *container,
                   ARRAY_POSITION       *ipos
		   )
{
#ifdef D_FSI
INT       i,surf;
INT       nir ,nil;
INT       ngsurf;
INT       ngnode;
INT       foundsurf;
DOUBLE    velint[3],vel2int[3];
DIS_TYP   typ;
NODE     *actfnode;
GSURF    *gsurf[6];
FLUID_FREESURF_CONDITION *surffs[6];

#ifdef DEBUG
dstrc_enter("f3_heightfunc");
#endif


fdyn   = alldyn[genprob.numff].fdyn;

amzero(estif_global);
amzero(eforce_global);

/*--------------------------------------------- set element coordinates */
f3_alecoor(ele,xyze);

/*-------------------------------------------------- set element values */
for(i=0;i<ele->numnp;i++) /* loop nodes of element */
{
   actfnode=ele->node[i];
   if(actfnode->hfdof==NULL) continue;
   /*---------------------------------------- set element values at n+1 */
   evelng[0][i]   =actfnode->sol_increment.a.da[ipos->velnp][0];
   evelng[1][i]   =actfnode->sol_increment.a.da[ipos->velnp][1];
   evelng[2][i]   =actfnode->sol_increment.a.da[ipos->velnp][2];
   ephing[i]      =actfnode->xfs[2];
   /*------------------------------------------ set element values at n */
   eveln[0][i]    =actfnode->sol_increment.a.da[ipos->veln][0];
   eveln[1][i]    =actfnode->sol_increment.a.da[ipos->veln][1];
   eveln[2][i]    =actfnode->sol_increment.a.da[ipos->veln][2];
   ephin[i]       =actfnode->sol_increment.a.da[ipos->veln][4];
} /* end of loop over nodes of element */

typ  = ele->distyp;
nir  = ele->e.f3->nGP[0];
foundsurf=0;
/*---------------------------------- number of surfaces to this element */
ngsurf=ele->g.gvol->ngsurf;

/*------- loop over surfaces, check for freesurface conditions on surfs */
for (i=0; i<ngsurf; i++)
{
   gsurf[i] = ele->g.gvol->gsurf[i];
   surffs[i] = gsurf[i]->freesurf;
   if(surffs[i]==NULL) continue;
   foundsurf++;
}

if (foundsurf!=1)
   dserror("no or too many element edges at free surface!\n");

/*------------------------------------------ set number of gauss points */
nil = IMAX(nir,2);

/*------------------------------------- loop over surfs at free surface */
for (surf=0; surf<ngsurf; surf++)
{
   if (surffs[surf]==NULL) continue;
   /*------------------------------------ check number of nodes on surf */
   ngnode = gsurf[surf]->ngnode;
   /*--------------------------------------------------- distyp of edge */
   switch (typ)
   {
   case hex8: typ=quad4; break;
   default: dserror("distyp not allowed for implicit free surface!\n");
   }
   /*--------------------------------------------------- get edge nodes */
   f3_iedg(iedgnod,ele,surf);
   /*--------------------------------- integration loop on actual gline */
   f3_calint_hfsep(ele,funct,deriv,wa1,wa2,xyze,ngnode,nil,
                   iedgnod,velint,vel2int,evelng,eveln,ephing,ephin,derxy,typ,
                   estif,eforce);
}

/*------------------------------------------- copy iedgnod to container */
container->ngnode=ngnode;
container->iedgnod=iedgnod;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
#else
dserror("FSI-functions not compiled!\n");
#endif
return;
} /* end of f3_heightfunc */

/*!---------------------------------------------------------------------
\brief control routine for stabilisation calculation

<pre>                                                         genk 05/04

evaluation of stabilisation parameter at the end of a time step

</pre>
\param   *ele     ELEMENT        the acutal element
\param   *ipos                           (i)   node array positions
\return void

------------------------------------------------------------------------*/
void f3_calstab(ELEMENT *ele, ARRAY_POSITION *ipos)
{
INT      i;
NODE    *actfnode;

#ifdef DEBUG
dstrc_enter("f3_calstab");
#endif

fdyn   = alldyn[genprob.numff].fdyn;

/*--------------------------------------------- get actual co-ordinates */
if (ele->e.f3->is_ale==0)
for(i=0;i<ele->numnp;i++)
{
   xyze[0][i]=ele->node[i]->x[0];
   xyze[1][i]=ele->node[i]->x[1];
   xyze[2][i]=ele->node[i]->x[2];
}
else
   f3_alecoor(ele,xyze);

/*------------------------------------------------ get actual velocity */
for(i=0;i<ele->numnp;i++) /* loop nodes of element */
{
   actfnode=ele->node[i];
   evelng[0][i]   =actfnode->sol_increment.a.da[ipos->velnp][0];
   evelng[1][i]   =actfnode->sol_increment.a.da[ipos->velnp][1];
   evelng[2][i]   =actfnode->sol_increment.a.da[ipos->velnp][2];
}

/*-------------------------- calculate element size and stab-parameter: */
f3_calelesize(ele,xyze,funct,deriv,deriv2,derxy,xjm,evelng,wa1,1);


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f2_calstab */



#endif
