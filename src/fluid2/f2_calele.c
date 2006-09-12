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
/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2
#include "../headers/standardtypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
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
static ARRAY     xyze_a;   /* actual element coordinates                */
static DOUBLE  **xyze;
static ARRAY     xjm_a;    /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     vderxy_a; /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     pderxy_a; /* pre -derivatives                          */
static DOUBLE   *pderxy;
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
void f2_calele(
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
dstrc_enter("f2_calele");
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
   xjm       = amdef("xjm"      ,&xjm_a      ,2,2        ,"DA");
   xyze      = amdef("xyze"     ,&xyze_a     ,2,MAXNOD_F2,"DA");
   vderxy    = amdef("vderxy"   ,&vderxy_a   ,2,2,"DA");
   pderxy    = amdef("pderxy"   ,&pderxy_a   ,2,1,"DV");
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

/*------------------------------------------------ initialise with ZERO */
amzero(estif_global);
amzero(emass_global);
amzero(eforce_global);
amzero(edforce_global);
*hasdirich=0;
*hasext=0;

switch(ele->e.f2->is_ale)
{
case 0:
/*---------------------------------------------------- set element data */
   f2_calset(ele,xyze,eveln,evelng,evhist,epren,edeadn,edeadng,ipos,hasext);

   switch (ele->e.f2->stab_type)
   {
   case stab_gls:
/*-------------------------- calculate element size and stab-parameter: */
      f2_calelesize(ele,eleke,xyze,funct,deriv,deriv2,xjm,derxy,
                    vderxy,evelng,ephin,ephing,eddy,&visc,ipos,0);
/*-------------------------------- calculate element stiffness matrices */
/*                                            and element force vectors */
      f2_calint(ele,estif,emass,eforce,
                xyze,funct,deriv,deriv2,xjm,derxy,derxy2,
                evelng,vderxy,wa1,wa2,visc);
   break;
   case stab_usfem:
      /*---------------------------------------------- get viscosity ---*/
      visc = mat[ele->mat-1].m.fluid->viscosity;

      /*--------------------------------------------- stab-parameter ---*/
      f2_caltau(ele,xyze,funct,deriv,xjm,evelng,visc);

      /*-------------------------------- perform element integration ---*/
      f2_int_usfem(ele,hasext,estif,eforce,xyze,
                   funct,deriv,deriv2,xjm,derxy,derxy2,evelng,eveln,
                   evhist,NULL,epren,edeadng,
                   vderxy,vderxy2,visc,wa1,wa2,estress, is_relax);
   break;
   default: dserror("unknown stabilisation type");
   }
break;
case 1:
   /*---------------------------------------------- set element data ---*/
   f2_calseta(ele,xyze,eveln,evelng,evhist,ealecovn,
              ealecovng,egridv,epren,edeadn,edeadng,ekappan,ekappang,
	      ephin,ephing,evnng,evnn,ipos,hasext,is_relax);
   switch (ele->e.f2->stab_type)
   {
   case stab_gls:
      /*-------------------- calculate element size and stab-parameter: */
      f2_calelesize(ele,eleke,xyze,funct,deriv,deriv2,xjm,derxy,
                    vderxy,ealecovng,ephin,ephing,eddy,&visc,ipos,0);
      /*-------------------------- calculate element stiffness matrices */
      /*                                         and element force vectors */
      f2_calinta(ele,imyrank,
                 estif,emass,eforce,
	         xyze,funct,deriv,deriv2,xjm,derxy,derxy2,
	         eveln,evelng,ealecovng,egridv,vderxy,
	         ekappan,ekappang,ephin,ephing,evnng,wa1,wa2);
   break;
   case stab_usfem:
      {
      /*---------------------------------------------- get viscosity ---*/
      visc = mat[ele->mat-1].m.fluid->viscosity;

      /*--------------------------------------------- stab-parameter ---*/
      f2_caltau(ele,xyze,funct,deriv,xjm,ealecovng,visc);

      /*--------------------------------- care for stress projection ---*/
      if (fdyn->stresspro)
      {
      INT i;
      NODE *actnode;
      for(i=0;i<ele->numnp;i++) /* loop nodes of element */
         {
            actnode=ele->node[i];
            /*---------------- set interpolated velocity derivatives ---*/
            estress[0][i]=actnode->sol_increment.a.da[ipos->stresspro][0];
            estress[1][i]=actnode->sol_increment.a.da[ipos->stresspro][1];
            estress[2][i]=actnode->sol_increment.a.da[ipos->stresspro][2];
         }
      }

      /*-------------------------------- perform element integration ---*/
      f2_int_usfem(ele,hasext,estif,eforce,xyze,
                   funct,deriv,deriv2,xjm,derxy,derxy2,evelng,eveln,
                   evhist,egridv,epren,edeadng,
                   vderxy,vderxy2,visc,wa1,wa2,estress, is_relax);
      break;
      }
    default: dserror("unknown stabilisation type");
   }
break;
default:
   dserror("parameter is_ale not 0 or 1!\n");
} /* end switch */

if (ele->e.f2->stab_type != stab_usfem)
{
#ifdef PERF
  perf_begin(21);
#endif

switch(ele->e.f2->fs_on)
{
case 0: case 1: case 3: /* no or explict free surface */
   /*---------------------- add up emassto estif and permute the matrix */
   f2_permestif(estif,emass,wa1,ele);
   /*------------------------------ permute element load vector eforce */
   if (fdyn->nii+(*hasext)!=0)
      f2_permeforce(eforce,wa1,ele->numnp);
break;
case 2: case 5: case 6:/* partitioned implict free surface */
   dsassert(ele->e.f2->is_ale!=0,"element at free surface has to be ALE!\n");
   /*--------------------- add up emass to estif and permute the matrix */
   f2_permestif_ifs(estif,emass,wa1,ele);
   /*------------------------------ permute element load vector eforce */
   if (fdyn->nii+(*hasext)!=0)
      f2_permeforce_ifs(eforce,wa1,ele);
break;
default:
   dserror("parameter fs_on out of range!\n");
}

#ifdef PERF
  perf_end(21);
#endif
/*---------------------------------- calculate rhs: emass * vel_hist ---*/
if (!is_relax)			/* calculation for relaxation parameter	*/
f2_massrhs(ele,emass,evhist,edeadng,eforce,hasext);
}

/* look for neumann bc */
f2_calneumann(ele, imyrank, eforce, xyze, funct, deriv, xjm, edeadn, edeadng);

/*-------------------------------------------- local co-ordinate system */
if(ele->locsys==locsys_yes)
   locsys_trans(ele,estif,NULL,NULL,eforce);


/*------------------------------- calculate element load vector edforce */
if (is_relax)			/* calculation for relaxation parameter	*/
   readfrom = ipos->relax;
else				/* standard case			*/
   readfrom = ipos->velnp;

/*------------------------------------------------ condensation of DBCs */
/* estif is in xyz* so edforce is also in xyz* (but DBCs have to be
   tranformed before condensing the dofs                                */
fluid_caldirich(ele,edforce,estif,hasdirich,readfrom);

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calele */

/*!---------------------------------------------------------------------
\brief control routine for stress calculation

<pre>                                                         genk 10/02

</pre>

\param     str       FLUID_STRESS   (i)    flag for stress calculation
\param     viscstr   INT            (i)    viscose stresses yes/no?
\param    *ele       ELEMENt        (i)    actual element
\param    *ipos                     (i)    node array positions
\param     is_relax  INT            (i)    flag
\return void

------------------------------------------------------------------------*/
void f2_stress(FLUID_STRESS  str,
               INT           viscstr,
	       ELEMENT      *ele,
               ARRAY_POSITION *ipos,
	       INT           is_relax  )
{
INT       i;        /* simply a counter                                 */
#ifdef D_FSI
INT       coupled;  /* flag for fsi interface element                   */
INT       iel;      /* number of nodes per element                      */
GNODE    *actgnode; /* actual gnode                                     */
#endif
INT       ldflag;
GSURF    *actgsurf;
GLINE    *actgline;
DLINE    *actdline;

#ifdef DEBUG
dstrc_enter("f2_stress");
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
      /* this approach does not work with a nonconforming discretization
         of the interface, thus it is replaced by the second one */
      /*
      if (actgnode->mfcpnode[genprob.numsf]==NULL) continue;
      */
      if (actgnode->fsicouple == NULL) continue;
      coupled=1;
      break;
   }
   if (coupled==1)
   f2_calelestress(viscstr,ele,evelng,epren,funct,
                   deriv,derxy,vderxy,xjm,xyze,sigmaint,ipos,is_relax);
break;
#endif
case str_liftdrag:
   /* check if element is on liftdrag-dline */
   actgsurf=ele->g.gsurf;
   ldflag=0;
   for (i=0;i<actgsurf->ngline;i++)
   {
      actgline=actgsurf->gline[i];
      actdline=actgline->dline;
      if (actdline==NULL) continue;
      if (actdline->liftdrag==0) continue;
      ldflag++;
      break;
   }
   if (ldflag>0)
   /* calculate element stresses for selected elements */
   f2_calelestress(viscstr,ele,evelng,epren,funct,
                   deriv,derxy,vderxy,xjm,xyze,sigmaint,ipos,0);
break;
case str_all:
   dserror(
       "stress computation for all elements not implemented yet\n");
break;
default:
   dserror("stress calculation not possible!\n");
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_stress */

/*!---------------------------------------------------------------------
\brief calculate shearstress for TURBULENCE MODEL

<pre>                                                         he  12/02


</pre>
\param   *ele      ELEMENT	     (i)    actual element
\param   *ipos                       (i)    node array positions
\return void

------------------------------------------------------------------------*/
void f2_shearstress(
	            ELEMENT    *ele,
                    ARRAY_POSITION *ipos
                   )
{
DOUBLE           det;
INT              actmat;       /* material number of the element */
INT              iel;
DOUBLE           r,s;
INT              icode,ihoel;
INT              node,i;
DIS_TYP          typ;	      /* element type                    */

GNODE *actgnode;
NODE  *actnode;               /* actual node                     */

#ifdef DEBUG
dstrc_enter("f2_shearstress");
#endif

iel  = ele->numnp;
typ  = ele->distyp;
actmat  = ele->mat-1;
visc    = mat[actmat].m.fluid->viscosity;
fdyn    = alldyn[genprob.numff].fdyn;
/*-------------------------------------------------------------------------*/

for (node=0; node<iel; node++)
{
 xyze[0][node]=ele->node[node]->x[0];
 xyze[1][node]=ele->node[node]->x[1];
 actnode  = ele->node[node];
 for(i=0; i<2; i++)
 {
  eveln[i][node] = actnode->sol_increment.a.da[ipos->velnp][i];
 }
}

/*-------------------------------------------------------------------------*/
for (node=0; node<iel; node++)
 {
   actnode  = ele->node[node];
   actgnode = actnode->gnode;
   if (actgnode->dirich==NULL) continue;
   if(actgnode->dirich->dirich_onoff.a.iv[0]==1 && actgnode->dirich->dirich_onoff.a.iv[1]==1)
   {
    if(actnode->sol_increment.a.da[ipos->velnp][0] == 0.0 && actnode->sol_increment.a.da[ipos->velnp][1] == 0.0)
    {
/*---------------get local coordinates of nodes------------------------*/
     r = f2_rsn(node,0,iel);
     s = f2_rsn(node,1,iel);

/*--------------- get values of  shape functions and their derivatives */
      switch(typ)
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
       icode   = 3;
       ihoel   = 1;
       f2_rec(funct,deriv,deriv2,r,s,typ,icode);
      break;
      case tri3: case tri6:   /* --> tri - element */
        icode   = 3;
        ihoel   = 1;
	f2_tri(funct,deriv,deriv2,r,s,typ,icode);
      break;
      default:
         dserror("typ unknown!");
      } /* end switch(typ) */

/*-------------------------------------------- compute Jacobian matrix */
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);

/*------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);

/*---------------------------------- compute global derivates for vel. */
      f2_vder(vderxy,derxy,eveln,iel);

/*------------------------------------- compute dimensless shearstress
            (ATTENTION: you've to divide with the square of flow rate) */
      actnode->fluid_varia->c_f_shear += 2*visc*(vderxy[0][1]+vderxy[1][0])/actnode->numele;

/*------------------- compute shearvelocity for the scaned coordinates */
      if (FABS(actnode->x[0]-fdyn->coord_scale[0])<EPS7 && FABS(actnode->x[1]-fdyn->coord_scale[1])<EPS15)
       fdyn->washvel += sqrt(visc*FABS(vderxy[0][1]+vderxy[1][0]))/actnode->numele;
    }
   }
 }


#ifdef DEBUG
dstrc_exit();
#endif

return;
} /*end of f2_shearstress */

/*!---------------------------------------------------------------------
\brief control routine for curvature calculation

<pre>                                                         genk 02/03

</pre>
\param    *ele     ELEMENt        (i)    actual element
\param     imyrank INT            (i)    proc number
\return void

------------------------------------------------------------------------*/
void f2_curvature(
	           ELEMENT        *ele,
		   INT             imyrank
		 )
{
#ifdef D_FSI
INT      i;
INT      ngline,foundline;
DOUBLE **kappa;
GLINE    *gline[4];
FLUID_FREESURF_CONDITION *linefs[4];

#ifdef DEBUG
dstrc_enter("f2_curvature");
#endif

fdyn  = alldyn[genprob.numff].fdyn;

kappa=ele->e.f2->kappa_ND.a.da;
/*-------------------------------------------------- store old solution */
if (fdyn->fsstnif!=0)
for (i=0;i<ele->numnp;i++) kappa[i][0]=kappa[i][1];

/*------------------------------------- number of lines to this element */
ngline=ele->g.gsurf->ngline;
foundline=0;
/*---------- loop over lines, check for freesurface conditions on lines */
for (i=0; i<ngline; i++)
{
   gline[i] = ele->g.gsurf->gline[i];
   linefs[i] = gline[i]->freesurf;
   if(linefs[i]==NULL) continue;
   foundline++;
}

if (foundline==0) goto end;

switch(ele->distyp)
{
case quad4: case tri3:
   f2_calq4curv(gline,linefs,ele,foundline,ngline,xyze,deriv,
                kappa);
break;
case quad8: case quad9:
   f2_calq8curv(gline,linefs,ele,ngline,xyze,kappa);
break;
case tri6:
   dserror("calculation of curvature not implemented for tri elements!\n");
default:
   dserror("distyp not valid!\n");
break;
}

end:

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

#else
dserror("FSI-functions not compiled in!\n");
#endif

return;
} /* end of f2_curvature */
/*!---------------------------------------------------------------------
\brief evaluate height function for free surface

<pre>                                                        genk 04/04

in this routine the height function is solved for the position of the
free surface

</pre>
\param   *ele              ELEMENT              actual element
\param   *estif_global     ARRAY                ele stiffness matrix
\param   *eforce_global   ARRAY                ele rhs
\param   *ipos                           (i)   node array positions
\param   *container        CONTAINER            container
\return void

------------------------------------------------------------------------*/
void f2_heightfunc(
                   ELEMENT              *ele,
                   ARRAY                *estif_global,
		   ARRAY                *eforce_global,
                   ARRAY_POSITION       *ipos,
		   CONTAINER            *container
		   )
{
#ifdef D_FSI
INT       i;
INT       nir ,nil;
INT       ngline;
INT       foundline;
DIS_TYP   typ;
INT       iedgnod[MAXNOD_F2];
INT       line,ngnode;
DOUBLE    phimin,phimax;
DOUBLE    velint[2],vel2int[2];
GLINE    *gline[4];
FLUID_FREESURF_CONDITION *linefs[4];
FLUID_DYNAMIC *fdyn;
NODE      *actfnode;

#ifdef DEBUG
dstrc_enter("f2_heightfunc");
#endif

fdyn   = alldyn[genprob.numff].fdyn;

amzero(estif_global);
amzero(eforce_global);

/*--------------------------------------------- set element coordinates */
f2_alecoor(ele,xyze);

/*-------------------------------------------------- set element values */
for(i=0;i<ele->numnp;i++) /* loop nodes of element */
{
   actfnode=ele->node[i];
   if(actfnode->hfdof==NULL) continue;
   /*---------------------------------------- set element values at n+1 */
   evelng[0][i]   =actfnode->sol_increment.a.da[ipos->velnp][0];
   evelng[1][i]   =actfnode->sol_increment.a.da[ipos->velnp][1];
   ephing[i]      =actfnode->xfs[1];
   /*------------------------------------------ set element values at n */
   eveln[0][i]    =actfnode->sol_increment.a.da[ipos->veln][0];
   eveln[1][i]    =actfnode->sol_increment.a.da[ipos->veln][1];
   ephin[i]       =actfnode->sol_increment.a.da[ipos->veln][3];
} /* end of loop over nodes of element */

typ  = ele->distyp;
nir  = ele->e.f2->nGP[0];
foundline=0;

/*------------------------------------- number of lines to this element */
ngline=ele->g.gsurf->ngline;

/*---------- loop over lines, check for freesurface conditions on lines */
for (i=0; i<ngline; i++)
{
   gline[i] = ele->g.gsurf->gline[i];
   linefs[i] = gline[i]->freesurf;
   if(linefs[i]==NULL) continue;
   foundline++;
}

if (foundline!=1)
   dserror("no or too many element edges at free surface!\n");

/*------------------------------------------ set number of gauss points */
nil = IMAX(nir,2);

/*------------------------------------- loop over lines at free surface */
for (line=0; line<ngline; line++)
{
   if (linefs[line]==NULL) continue;
   /*------------------------------------ check number of nodes on line */
   ngnode = gline[line]->ngnode;
   /*--------------------------------------------------- get edge nodes */
   f2_iedg(iedgnod,ele,line,0);
   phimax = DMAX(ephing[iedgnod[0]],ephing[iedgnod[1]]);
   phimin = DMIN(ephing[iedgnod[0]],ephing[iedgnod[1]]);
   fdyn->dphimax = DMAX(phimax-phimin,EPS6);

   /*-------------------------------- calculate stabilisation parameter */
   if (fdyn->hf_stab==1)
   f2_stabpar_hfsep(ele,xyze,evelng,funct,velint,ZERO,ZERO,ZERO,ZERO,
   iedgnod,ngnode,typ);
   /*--------------------------------- integration loop on actual gline */
   f2_calint_hfsep(ele,funct,deriv,wa1,xyze,ngnode,nil,
                   iedgnod,velint,vel2int,evelng,eveln,ephing,ephin,derxy,typ,
                   estif,eforce);
} /* end of loop over glines */

/*------------------------------------------- copy iedgnod to container */
container->ngnode=ngnode;
container->iedgnod=iedgnod;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
#endif
return;
} /* end of f2_heightfunc */
/*!---------------------------------------------------------------------
\brief control routine for stabilisation calculation

<pre>                                                         genk 05/04

especially in the ALE-case it seems to be advantageous to calculate
the stabilisation parameters at the end of a time step in order to use
them as tau_? at time (n)

</pre>
\param      *ele     ELEMENT        the actual element
\param  *ipos                           (i)   node array positions
\return void

------------------------------------------------------------------------*/
void f2_calstab(ELEMENT *ele, ARRAY_POSITION *ipos)
{
INT             i;
NODE           *actfnode;
INT             test=0;
INT            *hasext;

#ifdef DEBUG
dstrc_enter("f2_calstab");
#endif

hasext=&test;

/*--------------------------------------------- get actual co-ordinates */
if (ele->e.f2->is_ale==0)
for(i=0;i<ele->numnp;i++)
{
   xyze[0][i]=ele->node[i]->x[0];
   xyze[1][i]=ele->node[i]->x[1];
}
else
   f2_alecoor(ele,xyze);

/*------------------------------------------------ get actual velocity */
for(i=0;i<ele->numnp;i++) /* loop nodes of element */
{
   actfnode=ele->node[i];
   evelng[0][i]   =actfnode->sol_increment.a.da[ipos->velnp][0];
   evelng[1][i]   =actfnode->sol_increment.a.da[ipos->velnp][1];
}

/*-------------------------- calculate element size and stab-parameter: */
f2_calelesize(ele,NULL,xyze,funct,deriv,deriv2,xjm,derxy,
              vderxy,evelng,NULL,NULL,NULL,&visc,ipos,1);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of f2_calstab */


/*!---------------------------------------------------------------------
\brief control routine for calculation of normal

<pre>                                                         genk 05/04

In this routine the mass consistent normal is calculated.
See Gresho & Sani pp.450

</pre>
/param      *ele     ELEMENT        the actual element

\return void

------------------------------------------------------------------------*/
void f2_calnormal(ELEMENT *ele)
{
INT      i,k;
INT      iel;
INT      lr,ls,nir,nis,intc;
DOUBLE    e1,e2,facr,facs,fac,det;
NODE    *actnode;
DIS_TYP  typ;
FLUID_DYNAMIC *fdyn;
FLUID_DATA    *data;

#ifdef DEBUG
dstrc_enter("f2_calnormal");
#endif

fdyn = alldyn[genprob.numff].fdyn;
data = fdyn->data;

/*--------------------------------------------- get actual co-ordinates */
if (ele->e.f2->is_ale==0)
for(i=0;i<ele->numnp;i++)
{
   xyze[0][i]=ele->node[i]->x[0];
   xyze[1][i]=ele->node[i]->x[1];
}
else
   f2_alecoor(ele,xyze);

typ  = ele->distyp;
iel=ele->numnp;

/*------- get integraton data and check if elements are "higher order" */
switch (typ)
{
case quad4: case quad8: case quad9:  /* --> quad - element */
   nir = ele->e.f2->nGP[0];
   nis = ele->e.f2->nGP[1];
break;
case tri3: case tri6: /* --> tri - element */
    /* initialise integration */
   nir  = ele->e.f2->nGP[0];
   nis  = 1;
   intc = ele->e.f2->nGP[1];
break;
default:
   dserror("typ unknown!");
} /* end switch(typ) */

/*----------------------------------------- now perform the integration */
for (lr=0;lr<nir;lr++)
{
   for (ls=0;ls<nis;ls++)
   {
/*--------------- get values of  shape functions and their derivatives */
      switch(typ)
      {
      case quad4: case quad8: case quad9:   /* --> quad - element */
         e1   = data->qxg[lr][nir-1];
         facr = data->qwgt[lr][nir-1];
         e2   = data->qxg[ls][nis-1];
         facs = data->qwgt[ls][nis-1];
         f2_rec(funct,deriv,deriv2,e1,e2,typ,2);
      break;
      case tri3: case tri6:   /* --> tri - element */
         e1   = data->txgr[lr][intc];
         facr = data->twgt[lr][intc];
         e2   = data->txgs[lr][intc];
         facs = ONE;
         f2_tri(funct,deriv,deriv2,e1,e2,typ,2);
      break;
      default:
         dserror("typ unknown!");
      } /* end switch(typ) */
/*--------------------------------------------- compute Jacobian matrix */
      f2_jaco(xyze,deriv,xjm,&det,iel,ele);
      fac = facr*facs*det;
/*-------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);
/*--------- compute normal according to Gresho & Sani chapter 3.13.1. e */
      for (k=0;k<ele->numnp;k++)
      {
         actnode=ele->node[k];
         if (actnode->actn==NULL) continue;
         actnode->actn[0]+=derxy[0][k]*fac;
         actnode->actn[1]+=derxy[1][k]*fac;
      }
   }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calnormal */


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
void f2_caleleres(
	           ELEMENT        *ele,
	           ARRAY          *eforce_global,
                   ARRAY_POSITION *ipos,
                   INT            *hasdirich,
                   INT            *hasext
	       )
{
INT             i;
DOUBLE          estress[3][MAXNOD_F2];
NODE           *actnode;

#ifdef DEBUG
dstrc_enter("f2_caleleres");
#endif

/*--------------------------------------------- initialise with ZERO ---*/
amzero(eforce_global);
*hasdirich=0;
*hasext=0;

switch(ele->e.f2->is_ale)
{
case 0:
/*---------------------------------------------------- set element data */
   f2_calset(ele,xyze,eveln,evelng,evhist,epren,edeadn,edeadng,ipos,hasext);

   /*---------------------------------------------- get viscosity ---*/
   visc = mat[ele->mat-1].m.fluid->viscosity;

   /*--------------------------------------------- stab-parameter ---*/
   f2_caltau(ele,xyze,funct,deriv,xjm,evelng,visc);

   /*-------------------------------- perform element integration ---*/
   f2_int_res(ele,hasext,eforce,xyze,funct,deriv,deriv2,xjm,derxy,
              derxy2,evelng,evhist,NULL,epren,edeadng,vderxy,
              vderxy2,visc,wa1,wa2,estress);
break;
case 1:
   /*---------------------------------------------- set element data ---*/
   f2_calseta(ele,xyze,eveln,evelng,evhist,ealecovn,
              ealecovng,egridv,epren,edeadn,edeadng,ekappan,ekappang,
	      ephin,ephing,evnng,evnn,ipos,hasext,0);

   /*------------------------------------------------- get viscosity ---*/
   visc = mat[ele->mat-1].m.fluid->viscosity;

   /*------------------------------------------------ stab-parameter ---*/
   f2_caltau(ele,xyze,funct,deriv,xjm,ealecovng,visc);
   /*---------- prepare second derivatives in stress projection case ---*/
   if (fdyn->stresspro)
   {
      for(i=0;i<ele->numnp;i++) /* loop nodes of element */
      {
         actnode=ele->node[i];
         /*--------------------------------- set nodal stress values ---*/
         estress[0][i]=actnode->sol_increment.a.da[ipos->stresspro][0];
         estress[1][i]=actnode->sol_increment.a.da[ipos->stresspro][1];
         estress[2][i]=actnode->sol_increment.a.da[ipos->stresspro][2];
      }
   }

   /*----------------------------------- perform element integration ---*/
   f2_int_res(ele,hasext,eforce,xyze,funct,deriv,deriv2,xjm,derxy,
              derxy2,evelng,evhist,ealecovng,epren,edeadng,
              vderxy,vderxy2,visc,wa1,wa2,estress);
break;
default:
   dserror("parameter is_ale not 0 or 1!\n");
} /*end switch */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_caleleres */

/*!---------------------------------------------------------------------
\brief control routine for integration of stress projection matrix and rhs

<pre>                                                        chfoe 11/04

This routine controls the integration of the elemental matrix and rhs
to serve an L2-stress projection. The projection matrix has the structure
of a mass matrix (i.e. int(N^T N)dx ) and the right hand side carries
weighted stress values.

For further comments see header of f2_int_stress_project

</pre>

\param  *ele	         ELEMENT	(i)   actual element
\param  *hasext          INT	        (o)   element flag
\param **estif_global    ARRAY	        (o)   ele projection matrix
\param  *eforce_global   ARRAY	        (o)   elemental rhs
\param  *ipos                           (i)   node array positions
\return void

------------------------------------------------------------------------*/
void f2_calstresspro(ELEMENT   *ele,
                     INT       *hasext,
                     ARRAY     *estif_global,
		     ARRAY     *eforce_global,
                     ARRAY_POSITION *ipos
  )
{
INT i;
NODE *actnode;

#ifdef DEBUG
dstrc_enter("f2_calstresspro");
#endif
/*------------------------------------------------ initialise with ZERO */
amzero(estif_global);
amzero(eforce_global);

switch(ele->e.f2->is_ale)
{
case 0:
/*-------------------------------------------- set element coordinates -*/
   for(i=0;i<ele->numnp;i++)
   {
      xyze[0][i]=ele->node[i]->x[0];
      xyze[1][i]=ele->node[i]->x[1];
   }
break;
case 1:
   f2_alecoor(ele,xyze);
break;
default:
   dserror("parameter is_ale not 0 or 1!\n");
} /*end switch */



for(i=0;i<ele->numnp;i++) /* loop nodes of element */
{
   actnode=ele->node[i];
/*----------------------------------- set element velocities (n+gamma) */
   evelng[0][i]=actnode->sol_increment.a.da[ipos->velnp][0];
   evelng[1][i]=actnode->sol_increment.a.da[ipos->velnp][1];
} /* end of loop over nodes of element */

if (ele->e.f2->stab_type == stab_usfem)
{
   /*-------------------------------- perform element integration ---*/
   f2_int_stress_project(ele,hasext,estif,eforce,xyze,funct,deriv,
                         xjm,derxy,evelng,vderxy);
}
else dserror("No stress projection with other stabilisation than USFEM");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
}





/*-----------------------------------------------------------------------*/
/*!
  \brief control function for error calculation of f2 elements


  \param ele        *ELEMENT        (i) the element
  \param container  *CONTAINER      (i) contains variables defined in container.h
  \param ipos       *ARRAY_POSITION (i)

  \return void

  \author mn
  \date   08/05
 */
/*-----------------------------------------------------------------------*/
void f2_calerr(
    ELEMENT          *ele,
    CONTAINER        *container,
    ARRAY_POSITION   *ipos)
{

  INT       i;
  DOUBLE    visc;
  NODE     *actnode;

#ifdef DEBUG
  dstrc_enter("f2_err");
#endif


  /* get viscosity */
  visc = mat[ele->mat-1].m.fluid->viscosity;


  switch(ele->e.f2->is_ale)
  {
    case 0:

      /* set element coordinates */
      for(i=0;i<ele->numnp;i++)
      {
        xyze[0][i]=ele->node[i]->x[0];
        xyze[1][i]=ele->node[i]->x[1];
      }
      break;


    case 1:
      f2_alecoor(ele,xyze);
      break;


    default:
      dserror("parameter is_ale not 0 or 1!\n");

  }  /* switch(ele->e.f2->is_ale) */


  /* loop nodes of element */
  for(i=0;i<ele->numnp;i++)
  {
    actnode=ele->node[i];

    /* set element velocities (n+gamma) */
    evelng[0][i]=actnode->sol_increment.a.da[ipos->velnp][0];
    evelng[1][i]=actnode->sol_increment.a.da[ipos->velnp][1];

    /* set pressures (n+1) */
    epren[i]   =actnode->sol_increment.a.da[ipos->velnp][2];

  }  /* for(i=0;i<ele->numnp;i++) */


  /* perform element integration */
  f2_int_kim_moin_err(ele,xyze,funct,deriv,xjm,evelng,visc,
      epren,container);

#ifdef DEBUG
  dstrc_exit();
#endif

} /* f2_calerr */


#endif
/*! @} (documentation module close)*/

