/*!----------------------------------------------------------------------
\file
\brief element control routine

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
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*/		    
static ARRAY     eveln_a;  /* element velocities at (n)                 */
static double  **eveln;
static ARRAY     evelng_a; /* element velocities at (n+gamma)           */
static double  **evelng;
static ARRAY     ealecovn_a;  /* element ale-convective velocities      */
static double  **ealecovn;    /* at (n)                                 */
static ARRAY     ealecovng_a; /* element ale-convective velocities      */
static double  **ealecovng;   /* at (n+gamma)                           */
static ARRAY     egridv_a; /* element grid velocity                     */
static double  **egridv;
static ARRAY     epren_a;  /* element pressures at (n)	                */
static double   *epren;
static ARRAY     edeadn_a; /* element dead load (selfweight)            */
static double   *edeadng;
static ARRAY     edeadng_a;/* element dead load (selfweight)            */
static double   *edeadn;
static ARRAY     funct_a;  /* shape functions                           */
static double   *funct;
static ARRAY     deriv_a;  /* first natural derivatives                 */
static double  **deriv;
static ARRAY     deriv2_a; /* second natural derivatives                */
static double  **deriv2;
static ARRAY     xyze_a;
static double  **xyze;   
static ARRAY     xyzen_a;
static double  **xyzen;   
static ARRAY     xjm_a;    /* jocobian matrix                           */
static double  **xjm;
static ARRAY     velint_a; /* velocities at integration point           */
static double   *velint;
static ARRAY     vel2int_a;/* velocities at integration point           */
static double   *vel2int;
static ARRAY     alecovint_a; /* ale-convective vel at integration point*/
static double   *alecovint;
static ARRAY     gridvint_a;  /* grid-vel at integration point          */
static double   *gridvint;
static ARRAY	 covint_a; /* convective velocities at integr. point    */
static double   *covint;
static ARRAY     vderxy_a; /* vel - derivatives                         */
static double  **vderxy;
static ARRAY     pderxy_a; /* pre -derivatives                          */
static double   *pderxy;
static ARRAY     vderxy2_a;/* vel - 2nd derivatives                     */
static double  **vderxy2;
static ARRAY     derxy_a;  /* coordinate - derivatives                  */
static double  **derxy;
static ARRAY     derxy2_a; /* 2nd coordinate - derivatives              */
static double  **derxy2;
static ARRAY     sigmaint_a; /* fluid stresses at integration point     */
static double  **sigmaint;
static ARRAY     ekappan_a; /* surface curvature at (n)                  */
static double   *ekappan;
static ARRAY     ekappang_a; /* surface curvature at (n+1)               */
static double   *ekappang;
static ARRAY     w1_a;     /* working array of arbitrary chosen size    */
static double  **wa1;      /* used in different element routines        */
static ARRAY     w2_a;     /* working array of arbitrary chosen size    */
static double  **wa2;      /* used in different element routines        */
static double  **estif;    /* pointer to global ele-stif                */
static double  **emass;    /* pointer to galerkin ele-stif              */
static double   *etforce;  /* pointer to Time RHS                       */
static double   *eiforce;  /* pointer to Iteration RHS                  */
static double   *edforce;  /* pointer to RHS due to dirichl. conditions */

static ARRAY     eddy_a;       /* element turbulent ken. energy  at (n+gamma)*/
static double   *eddy;

double           visc;
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
\param  *data	         FLUID_DATA     (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *ele	         ELEMENT	(i)   actual element
\param  *eleke	         ELEMENT	(i)   element for turbulence-model
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *etforce_global  ARRAY	        (o)   element time force
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param  *hasdirich       int	        (o)   element flag
\param  *hasext          int	        (o)   element flag         
\param   imyrank         int            (i)   proc number
\param   velgrad	         int	        (i)   flag
\param   init	         int	        (i)   init flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_calele(
                FLUID_DATA     *data, 
                FLUID_DYN_CALC *dynvar, 
	          ELEMENT        *ele,             
                ELEMENT        *eleke, 
                ARRAY          *estif_global,   
                ARRAY          *emass_global,   
	          ARRAY          *etforce_global,       
	          ARRAY          *eiforce_global, 
		    ARRAY          *edforce_global,		
		    int            *hasdirich,      
                int            *hasext,
                int             imyrank,
		    int             velgrad,            
		    int             init            
	       )
{

#ifdef DEBUG 
dstrc_enter("f2_calele");
#endif

if (init==1) /* allocate working arrays and set pointers */
{
   eveln     = amdef("eveln"    ,&eveln_a    ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   evelng    = amdef("evelng"   ,&evelng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   ealecovn  = amdef("ealecovn" ,&ealecovn_a ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   ealecovng = amdef("ealecovng",&ealecovng_a,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   egridv    = amdef("egridv"   ,&egridv_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   epren     = amdef("epren"    ,&epren_a    ,MAXNOD_F2,1,"DV");
   edeadn    = amdef("edeadn"   ,&edeadn_a   ,2,1,"DV");
   edeadng   = amdef("edeadng"  ,&edeadng_a  ,2,1,"DV");      
   funct     = amdef("funct"    ,&funct_a    ,MAXNOD_F2,1,"DV");
   deriv     = amdef("deriv"    ,&deriv_a    ,2,MAXNOD_F2,"DA");
   deriv2    = amdef("deriv2"   ,&deriv2_a   ,3,MAXNOD_F2,"DA");
   xjm       = amdef("xjm"      ,&xjm_a      ,2,2        ,"DA");
   xyze      = amdef("xyze"     ,&xyze_a     ,2,MAXNOD_F2,"DA");
   xyzen     = amdef("xyzen"    ,&xyzen_a    ,2,MAXNOD_F2,"DA");
   velint    = amdef("velint"   ,&velint_a   ,NUM_F2_VELDOF,1,"DV");
   vel2int   = amdef("vel2int"  ,&vel2int_a  ,NUM_F2_VELDOF,1,"DV");
   alecovint = amdef("alecovint",&alecovint_a,NUM_F2_VELDOF,1,"DV");
   gridvint  = amdef("gridvint" ,&gridvint_a ,NUM_F2_VELDOF,1,"DV");
   covint    = amdef("covint"   ,&covint_a   ,NUM_F2_VELDOF,1,"DV");
   vderxy    = amdef("vderxy"   ,&vderxy_a   ,2,2,"DA");
   pderxy    = amdef("pderxy"   ,&pderxy_a   ,2,1,"DV");
   vderxy2   = amdef("vderxy2"  ,&vderxy2_a  ,2,3,"DA");
   derxy     = amdef("derxy"    ,&derxy_a    ,2,MAXNOD_F2,"DA");
   derxy2    = amdef("derxy2"   ,&derxy2_a   ,3,MAXNOD_F2,"DA");
   sigmaint  = amdef("sigmaint" ,&sigmaint_a ,3,MAXGAUSS ,"DA");
   ekappan   = amdef("ekappan"  ,&ekappan_a  ,MAXNOD_F2,1 ,"DV");
   ekappang  = amdef("ekappang" ,&ekappang_a ,MAXNOD_F2,1 ,"DV");
   eddy      = amdef("eddy"     ,&eddy_a  ,MAXNOD_F2,1,"DV");
   wa1       = amdef("wa1"      ,&w1_a       ,50,50,"DA");
   wa2       = amdef("wa2"      ,&w2_a       ,50,50,"DA");  
/*                                               \- size is arbitrary chosen!  */
   estif   = estif_global->a.da;
   emass   = emass_global->a.da;
   eiforce = eiforce_global->a.dv;
   etforce = etforce_global->a.dv;
   edforce = edforce_global->a.dv;
   goto end;
} /* endif (init==1) */

/*------------------------------------------------ initialise with ZERO */
amzero(estif_global);
amzero(emass_global);
amzero(eiforce_global);
amzero(etforce_global);
amzero(edforce_global);
*hasdirich=0;
*hasext=0;

switch(ele->e.f2->is_ale)
{
case 0:
/*compute the shear stresses (only if k-omega or k-epsilon is activated)*/
if (velgrad==1) 
{
 f2_shearstress(ele,dynvar,evelng,vderxy,xjm,xyze,deriv,deriv2,derxy,funct);
 goto end;
}

/*---------------------------------------------------- set element data */
   f2_calset(dynvar,ele,xyze,eveln,evelng,epren,edeadn,edeadng,hasext);

/*-------------------------- calculate element size and stab-parameter: */
   f2_calelesize(ele,eleke,data,dynvar,xyze,funct,deriv,deriv2,xjm,derxy,
                 vderxy,evelng,velint,wa1,eddy,&visc);
/*-------------------------------- calculate element stiffness matrices */
/*                                            and element force vectors */
   f2_calint(data,ele,dynvar,hasext,
             estif,emass,etforce,eiforce,
	     xyze,funct,deriv,deriv2,xjm,derxy,derxy2,
	     eveln,evelng,epren,edeadn,edeadng,
	     velint,vel2int,covint,vderxy,pderxy,vderxy2,
	     wa1,wa2,visc);
break;
case 1:
/*---------------------------------------------------- set element data */
   f2_calseta(dynvar,ele,xyze,eveln,evelng,ealecovn,
               ealecovng,egridv,epren,edeadn,edeadng,ekappan,ekappang,hasext);
/*-------------------------- calculate element size and stab-parameter: */   
   f2_calelesize(ele,eleke,data,dynvar,xyze,funct,deriv,deriv2,xjm,derxy,
                 vderxy,ealecovng,velint,wa1,eddy,&visc);
/*-------------------------------- calculate element stiffness matrices */
/*                                            and element force vectors */
   f2_calinta(data,ele,dynvar,hasext,imyrank,
              estif,emass,etforce,eiforce,
	      xyze,funct,deriv,deriv2,xjm,derxy,derxy2,
	      eveln,evelng,ealecovn,ealecovng,egridv,epren,edeadn,edeadng,
	      velint,vel2int,alecovint,gridvint,covint,vderxy,pderxy,vderxy2,
	      ekappan,ekappang,wa1,wa2);  
break;
default:
   dserror("parameter is_ale not 0 or 1!\n");
} /*end switch */

switch(ele->e.f2->fs_on)
{
case 0: case 1: /* no or explict free surface */
   /*-------------- add emass and estif to estif and permute the matrix */
   f2_permestif(estif,emass,wa1,ele,dynvar);
   /*------------------------------ permute element load vector etforce */
   if (dynvar->nif!=0)
      f2_permeforce(etforce,wa1,ele->numnp);   
   /*------------------------------ permute element load vector eiforce */
   if (dynvar->nii+(*hasext)!=0)
      f2_permeforce(eiforce,wa1,ele->numnp);
break;
case 2: /* implict free surface */
   dsassert(ele->e.f2->is_ale!=0,"element at free surface has to be ALE!\n");
   f2_permestif_ifs(estif,emass,wa1,ele,dynvar);
   /*------------------------------ permute element load vector etforce */
   if (dynvar->nif!=0)
      f2_permeforce_ifs(etforce,wa1,ele);
   /*------------------------------ permute element load vector eiforce */
   if (dynvar->nii+(*hasext)!=0)
      f2_permeforce_ifs(eiforce,wa1,ele);
break;   
default:
   dserror("parameter fs_on out of range!\n");
}

/*------------------------------- calculate element load vector edforce */
fluid_caldirich(ele,edforce,estif,hasdirich);

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

\param     str     FLUID_STRESS   (i)    flag for stress calculation
\param     viscstr int            (i)    viscose stresses yes/no?
\param    *data    FLUID_DATA     (i)
\param    *ele     ELEMENt        (i)    actual element 
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_stress(FLUID_STRESS  str, 
               int           viscstr,
	       FLUID_DATA   *data, 
	       ELEMENT      *ele      )
{
int       i;        /* simply a counter                                 */
int       coupled;  /* flag for fsi interface element                   */
int       iel;      /* number of nodes per element                      */
int       actmat;   /* actual material number                           */
GNODE    *actgnode; /* actual gnode                                     */

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
      if (actgnode->mfcpnode[genprob.numsf]==NULL) continue;
      coupled=1;
      break;    
   }
   if (coupled==1) 
   f2_calfsistress(viscstr,data,ele,eveln,epren,funct,
                   deriv,derxy,vderxy,xjm,xyze,sigmaint);      
break;
#endif
case str_liftdrag:
   dserror("lift&drag computation not implemented yet\n");
break;
case str_all:
   dserror("stress computation for all elements not implemented yet\n");
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
\brief control routine for curvature calculation

<pre>                                                         genk 02/03  
			     			
</pre>
\param    *data    FLUID_DATA     (i)
\param    *dynvar  FLUID_DYN_CALC (i)
\param    *ele     ELEMENt        (i)    actual element 
\param     imyrank int            (i)    proc number
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_curvature(
	           FLUID_DATA     *data,
		   FLUID_DYN_CALC *dynvar,  
	           ELEMENT        *ele,
		   int             imyrank
		 )
{
#ifdef D_FSI
int      i;
int      ngline,foundline;
double **kappa;
GLINE    *gline[4];
FLUID_FREESURF_CONDITION *linefs[4];

#ifdef DEBUG 
dstrc_enter("f2_curvature");
#endif

kappa=ele->e.f2->kappa_ND.a.da;
/*-------------------------------------------------- store old solution */
if (dynvar->fsstnif!=0)
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
   f2_calq4curv(gline,linefs,dynvar,ele,foundline,ngline,xyze,deriv,
                kappa);
break;
case quad8: case quad9:
   f2_calq8curv(gline,linefs,dynvar,ele,ngline,xyze,kappa);
break;
case tri6:
   dserror("calculation of curvature not implemented for tri elements!\n");
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

#endif
/*! @} (documentation module close)*/
