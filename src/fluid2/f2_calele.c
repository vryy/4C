/*!----------------------------------------------------------------------
\file
\brief element control routine

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "../headers/solution.h"
#include "../fluid_full/fluid_prototypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
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
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *etforce_global  ARRAY	        (o)   element time force
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param*  hasdirich       int	        (o)   element flag
\param   init	         int	        (i)   init flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_calele(
                FLUID_DATA     *data, 
                FLUID_DYN_CALC *dynvar, 
	        ELEMENT        *ele,             
                ARRAY          *estif_global,   
                ARRAY          *emass_global,   
	        ARRAY          *etforce_global,       
	        ARRAY          *eiforce_global, 
		ARRAY          *edforce_global,	
		int            *hasdirich,      
		int             init            
	       )
{
static ARRAY     eveln_a;  /* element velocities at (n)                 */
static double  **eveln;
static ARRAY     evelng_a; /* element velocities at (n+gamma)           */
static double  **evelng;
static ARRAY     epren_a;  /* element pressures at (n)	                */
static double   *epren;
static ARRAY     funct_a;  /* shape functions                           */
static double   *funct;
static ARRAY     deriv_a;  /* first natural derivatives                 */
static double  **deriv;
static ARRAY     deriv2_a; /* second natural derivatives                */
static double  **deriv2;
static ARRAY     xjm_a;    /* jocobian matrix                           */
static double  **xjm;
static ARRAY     velint_a; /* velocities at integration point           */
static double   *velint;
static ARRAY     vel2int_a;/* velocities at integration point           */
static double   *vel2int;
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
static ARRAY     w1_a;     /* working array of arbitrary chosen size    */
static double  **wa1;      /* used in different element routines        */
static ARRAY     w2_a;     /* working array of arbitrary chosen size    */
static double  **wa2;      /* used in different element routines        */
static double  **estif;    /* pointer to global ele-stif                */
static double  **emass;    /* pointer to galerkin ele-stif              */
static double   *etforce;  /* pointer to Time RHS                       */
static double   *eiforce;  /* pointer to Iteration RHS                  */
static double   *edforce;  /* pointer to RHS due to dirichl. conditions */

#ifdef DEBUG 
dstrc_enter("f2_calele");
#endif

if (init==1) /* allocate working arrays and set pointers */
{
   eveln   = amdef("evln"   ,&eveln_a  ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   evelng  = amdef("evelng" ,&evelng_a ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   epren   = amdef("epren"  ,&epren_a  ,MAXNOD_F2,1,"DV");
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
   wa1     = amdef("wa1"    ,&w1_a     ,30,30,"DA");
   wa2     = amdef("wa2"    ,&w2_a     ,30,30,"DA");  
/*                                        \- size is arbitrary chosen!  */
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

/*---------------------------------------------------- set element data */
f2_calset(dynvar,ele,eveln,evelng,epren);

/*-------------------------- calculate element size and stab-parameter: */
f2_calelesize(ele,data,dynvar,funct,deriv,deriv2,xjm,evelng,velint,wa1);

/*-------------------------------- calculate element stiffness matrices */
/*                                            and element force vectors */
f2_calint(data,ele,dynvar,
          estif,emass,etforce,eiforce,
	  funct,deriv,deriv2,xjm,derxy,derxy2,
	  eveln,evelng,epren,
	  velint,vel2int,covint,vderxy,pderxy,vderxy2,
	  wa1,wa2);

/*----------------- add emass and estif to estif and permute the matrix */
f2_permestif(estif,emass,wa1,ele->numnp,dynvar);

/*--------------------------------- permute element load vector etforce */
if (dynvar->nif!=0)
   f2_permeforce(etforce,wa1,ele->numnp);

/*--------------------------------- permute element load vector eiforce */
if (dynvar->nii!=0)
   f2_permeforce(eiforce,wa1,ele->numnp);

/*------------------------------- calculate element load vector edforce */
fluid_caldirich(ele,edforce,estif,hasdirich);

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calele */


#endif
