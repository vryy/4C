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

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;

/*----------------------------------------------------------------------*/
static ARRAY     eveln_a;  /* element velocities at (n)                 */
static DOUBLE  **eveln;
/*NOTE: if there is no classic time rhs (as described in WAW) the array
         eveln is misused and does NOT contain the velocity at time (n)
	 but rather a linear combination of old velocities and 
	 accelerations depending upon the time integration scheme!!!!! */
static ARRAY     evelng_a; /* element velocities at (n+gamma)           */
static DOUBLE  **evelng;
static ARRAY     ealecovn_a;  /* element ale-convective velocities      */
static DOUBLE  **ealecovn;    /* at (n)                                 */
static ARRAY     ealecovng_a; /* element ale-convective velocities      */
static DOUBLE  **ealecovng;   /* at (n+gamma)                           */
static ARRAY     egridv_a; /* element grid velocity                     */
static DOUBLE  **egridv;
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
static ARRAY     xyze_a;
static DOUBLE  **xyze;   
static ARRAY     xyzen_a;
static DOUBLE  **xyzen;   
static ARRAY     xjm_a;    /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     velint_a; /* velocities at integration point           */
static DOUBLE   *velint;
static ARRAY     vel2int_a;/* velocities at integration point           */
static DOUBLE   *vel2int;
static ARRAY     alecovint_a; /* ale-convective vel at integration point*/
static DOUBLE   *alecovint;
static ARRAY     gridvint_a;  /* grid-vel at integration point          */
static DOUBLE   *gridvint;
static ARRAY	 covint_a; /* convective velocities at integr. point    */
static DOUBLE   *covint;
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
static DOUBLE  **sigmaint;
static ARRAY     ekappan_a; /* surface curvature at (n)                  */
static DOUBLE   *ekappan;
static ARRAY     ekappang_a; /* surface curvature at (n+1)               */
static DOUBLE   *ekappang;
static ARRAY     w1_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa1;      /* used in different element routines        */
static ARRAY     w2_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa2;      /* used in different element routines        */
static DOUBLE  **estif;    /* pointer to global ele-stif                */
static DOUBLE  **emass;    /* pointer to galerkin ele-stif              */
static DOUBLE   *etforce;  /* pointer to Time RHS                       */
static DOUBLE   *eiforce;  /* pointer to Iteration RHS                  */
static DOUBLE   *edforce;  /* pointer to RHS due to dirichl. conditions */

static ARRAY     eddy_a;       /* element turbulent ken. energy  at (n+gamma)*/
static DOUBLE   *eddy;

static DOUBLE    visc;
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
\param  *hasdirich       INT	        (o)   element flag
\param  *hasext          INT	        (o)   element flag         
\param   imyrank         INT            (i)   proc number
\param   velgrad         INT	        (i)   flag
\param   is_relax        INT            (i)   flag
\param   init	         INT	        (i)   init flag
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
		INT            *hasdirich,      
                INT            *hasext,
                INT             imyrank,
		INT             is_relax,
		INT             init            
	       )
{
INT		readfrom;	/* where to read dbc from 		*/

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
               ealecovng,egridv,epren,edeadn,edeadng,ekappan,ekappang,hasext,
	       is_relax);
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
if (is_relax)			/* calculation for relaxation parameter	*/
   readfrom = 7;
else				/* standard case			*/
   readfrom = 3;		

fluid_caldirich(ele,edforce,estif,hasdirich,readfrom);
   
/*------------------------- calculate emass * vel(n) for BDF2 method ---*/
if (dynvar->nim)
   f2_massrhs(ele,emass,eveln,eiforce); 
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
\param    *data      FLUID_DATA     (i)
\param    *ele       ELEMENt        (i)    actual element 
\param     is_relax  INT            (i)    flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_stress(FLUID_STRESS  str, 
               INT           viscstr,
	       FLUID_DATA   *data, 
	       ELEMENT      *ele,
	       INT           is_relax  )
{
INT       i;        /* simply a counter                                 */
INT       coupled;  /* flag for fsi interface element                   */
INT       iel;      /* number of nodes per element                      */
INT       actmat;   /* actual material number                           */
INT       ldflag;
GNODE    *actgnode; /* actual gnode                                     */
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
      if (actgnode->mfcpnode[genprob.numsf]==NULL) continue;
      coupled=1;
      break;    
   }
   if (coupled==1) 
   f2_calfsistress(viscstr,data,ele,eveln,epren,funct,
                   deriv,derxy,vderxy,xjm,xyze,sigmaint,is_relax);      
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
   f2_calelestress(viscstr,data,ele,eveln,epren,funct,
                   deriv,derxy,vderxy,xjm,xyze,sigmaint);      
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
\brief calculate shearstress for TURBULENCE MODEL

<pre>                                                         he  12/02

				      
</pre>
\param   *ele      ELEMENT	     (i)    actual element
\param  *dynvar    FLUID_DYN_CALC  (i)
\param  **evel     DOUBLE	     (i)    vel. at nodes
\param  **vderxy    DOUBLE	     (-)    vel. derivates
\param  **xjm       DOUBLE	     (-)    jacobi matrix
\param  **deriv     DOUBLE	     (-)    local deriv.
\param  **deriv2    DOUBLE	     (-)    local 2nd deriv.
\param  **derxy     DOUBLE	     (-)    global deriv.
\param  *funct      DOUBLE	     (-)    shape funct.
\return void                                                                       

------------------------------------------------------------------------*/
void f2_shearstress(
	           ELEMENT    *ele,
                 FLUID_DYN_CALC  *dynvar 
                 )
{
DOUBLE           det;    
INT              actmat;       /* material number of the element */
INT              iel,ntyp;
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
ntyp = ele->e.f2->ntyp; 
typ  = ele->distyp;
actmat  = ele->mat-1;
visc    = mat[actmat].m.fluid->viscosity;
/*-------------------------------------------------------------------------*/

for (node=0; node<iel; node++) 
{
 xyze[0][node]=ele->node[node]->x[0];
 xyze[1][node]=ele->node[node]->x[1];
 actnode  = ele->node[node];
 for(i=0; i<2; i++)
 {
  eveln[i][node] = actnode->sol_increment.a.da[3][i]; 
 }
}

/*-------------------------------------------------------------------------*/
for (node=0; node<iel; node++) 
 {
   actnode  = ele->node[node];
   actgnode = actnode->gnode;      
   if(actgnode->dirich->dirich_onoff.a.iv[0]==1 && actgnode->dirich->dirich_onoff.a.iv[1]==1)
   {
    if(actnode->sol_increment.a.da[3][0] == 0.0 && actnode->sol_increment.a.da[3][1] == 0.0)
    {
/*---------------get local coordinates of nodes------------------------*/
     r = f2_rsn(node,0,iel);
     s = f2_rsn(node,1,iel);

/*--------------- get values of  shape functions and their derivatives */
      switch(ntyp)  
      {
      case 1:   /* --> quad - element */
       icode   = 3;
       ihoel   = 1;
       f2_rec(funct,deriv,deriv2,r,s,typ,icode);
      break;
      case 2:   /* --> tri - element */              
       if (iel>3)
       {
        icode   = 3;
        ihoel   = 1;
       }
	 f2_tri(funct,deriv,deriv2,r,s,typ,icode);
      break; 
      default:
         dserror("ntyp unknown!");
      } /* end switch(ntyp) */
     
/*-------------------------------------------- compute Jacobian matrix */
      f2_jaco(xyze,funct,deriv,xjm,&det,iel,ele);

/*------------------------------------------- compute global derivates */
      f2_gder(derxy,deriv,xjm,det,iel);

/*---------------------------------- compute global derivates for vel. */
      f2_vder(vderxy,derxy,eveln,iel);	

/*------------------------------------- compute dimensless shearstress
            (ATTENTION: you've to divide with the square of flow rate) */
      actnode->fluid_varia->c_f_shear += 2*visc*(vderxy[0][1]+vderxy[1][0])/actnode->numele;       

/*------------------- compute shearvelocity for the scaned coordinates */
      if (FABS(actnode->x[0]-dynvar->coord_scale[0])<EPS7 && FABS(actnode->x[1]-dynvar->coord_scale[1])<EPS15)
       dynvar->washvel += sqrt(visc*FABS(vderxy[0][1]+vderxy[1][0]))/actnode->numele;       
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
\param    *data    FLUID_DATA     (i)
\param    *dynvar  FLUID_DYN_CALC (i)
\param    *ele     ELEMENt        (i)    actual element 
\param     imyrank INT            (i)    proc number
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_curvature(
	           FLUID_DATA     *data,
		   FLUID_DYN_CALC *dynvar,  
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
