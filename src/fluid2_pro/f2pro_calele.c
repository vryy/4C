/*!----------------------------------------------------------------------
\file
\brief element control routine

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
/*!---------------------------------------------------------------------                                         
\brief control routine for element integration of fluid2pro

<pre>                                                         basol 10/02

This routine controls the element evaluation:
-actual vel. and pres. variables are set
-element integration is performed --> element stiffness matrix and 
                                  --> element load vectors
-element load vector due to dirichlet conditions is calculated				      
			     
</pre>
\param  *data	         FLUID_DATA     (i)
\param  *dynvar	         FLUID_DYN_CALC (i)
\param  *elev	         ELEMENT	(i)   actual velocity element
\param  *elep	         ELEMENT	(i)   actual pressure element
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *lmass_global    ARRAY	        (o)   lumped mass matrix
\param  *gradopr_global  ARRAY	        (o)   gradient operator
\param  *etforce_global  ARRAY	        (o)   element time force 
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param*  hasdirich       int	        (o)   element flag
\param   init	         int	        (i)   init flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2pro_calele(
                FLUID_DATA     *data, 
                FLUID_DYN_CALC *dynvar, 
	        ELEMENT        *elev, 
	        ELEMENT        *elep,            
                ARRAY          *estif_global,   
                ARRAY          *emass_global,
		ARRAY          *lmass_global,
		ARRAY          *gradopr_global,   
	        ARRAY          *etforce_global,
		ARRAY          *eiforce_global, 
		ARRAY          *edforce_global,
		INT             *hasdirich,      
		INT             init            
	       )
{
static ARRAY     eveln_a;    /* element velocities at (n)               */
static DOUBLE  **eveln;
static ARRAY     epren_a;    /* element pressures at (n)	        */
static DOUBLE  *epren;
static ARRAY     funct_a;    /* shape functions for velocity            */
static DOUBLE  *funct;
static ARRAY     functpr_a;  /* shape functions for pressure            */
static DOUBLE  *functpr;
static ARRAY     deriv_a;    /* first natural derivatives for velocity shape fnct.  */
static DOUBLE  **deriv;
static ARRAY     derivpr_a;  /* first natural derivatives for pressure shape fnct.  */
static DOUBLE  **derivpr;
static ARRAY     deriv2_a;   /* second natural derivatives for velocity shape fnct. */
static DOUBLE   **deriv2;
static ARRAY     xjm_a;      /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     velint_a;   /* velocities at integration point           */
static DOUBLE  *velint;
static ARRAY	 covint_a;   /* convective velocities at integr. point    */
static DOUBLE  *covint;
static ARRAY     vderxy_a;   /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     pderxy_a;   /* pre -derivatives                          */
static DOUBLE  *pderxy;
static ARRAY     derxy_a;    /* coordinate - derivatives for velocity shape functions */
static DOUBLE  **derxy;
static ARRAY     derxypr_a;  /* coordinate - derivatives for pressure shape functions */
static DOUBLE  **derxypr;
static ARRAY     xyze_a;     /* coordinates of the element */
static DOUBLE  **xyze; 
static ARRAY     w1_a;         /* working array of arbitrary chosen size     */
static DOUBLE  **wa1;          /* used in different element routines         */
static DOUBLE  **estif;        /* pointer to global ele-stif                 */
static DOUBLE  **emass;        /* pointer to galerkin ele-stif               */
static DOUBLE  **lmass;        /* pointer to galerkin lumped ele mass matrix */ 
static DOUBLE  **gradopr;      /* pointer to gradient operator               */
static DOUBLE  *etforce ;      /* pointer to Time RHS                        */
static DOUBLE  *eiforce;       /* pointer to Iteration RHS                   */
static DOUBLE  *edforce;       /* pointer to RHS due to dirichl. conditions  */
DOUBLE  dirich[MAXDOFPERELE];  /* dirichlet values of act. ele               */
INT dirich_onoff[MAXDOFPERELE]; /* dirichlet flags of act. ele               */ 
INT i,j,iel;
iel=elev->numnp;
#ifdef DEBUG 
dstrc_enter("f2pro_calele");
#endif
if (init==1) /* allocate working arrays and set pointers */
{
   eveln   = amdef("eveln"   ,&eveln_a   ,2,9,"DA");
   epren   = amdef("epren"   ,&epren_a   ,4,1,"DV");
   funct   = amdef("funct"   ,&funct_a   ,9,1,"DV");
   functpr = amdef("functpr" ,&functpr_a ,4,1,"DV");
   deriv   = amdef("deriv"   ,&deriv_a   ,2,9,"DA");
   deriv2  = amdef("deriv2"  ,&deriv2_a  ,3,9,"DA");
   derivpr = amdef("derivpr" ,&derivpr_a ,2,4,"DA");
   xjm     = amdef("xjm"     ,&xjm_a     ,2,2,"DA");
   xyze    = amdef("xyze"    ,&xyze_a    ,2,9,"DA");
   velint  = amdef("velint"  ,&velint_a  ,2,1,"DV");
   covint  = amdef("covint"  ,&covint_a  ,2,1,"DV");
   vderxy  = amdef("vderxy"  ,&vderxy_a  ,2,2,"DA");
   pderxy  = amdef("pderxy"  ,&pderxy_a  ,2,1,"DV");
   derxy   = amdef("derxy"   ,&derxy_a   ,2,9,"DA");
   derxypr = amdef("derxypr" ,&derxypr_a ,2,4,"DA");
   wa1     = amdef("wa1"    ,&w1_a     ,30,30,"DA");
   estif   = estif_global->a.da;
   emass   = emass_global->a.da;
   lmass   = lmass_global->a.da;
   gradopr = gradopr_global->a.da;
   eiforce = eiforce_global->a.dv;
   etforce = etforce_global->a.dv;
   edforce = edforce_global->a.dv;
   goto end;
} /* endif (init==1) */
/*------------------------------------------------ initialise with ZERO */
amzero(estif_global);
amzero(emass_global);
amzero(lmass_global);
amzero(gradopr_global);
amzero(eiforce_global);
amzero(etforce_global);
amzero(edforce_global);
*hasdirich=0;
/*---------------------------------------------------- set element data */
if (dynvar->pro_calveln==1) 
   f2pro_calset(elev,elep,xyze,eveln,epren);
/*-------------------------------- calculate element stiffness matrices */
/*                                            and element force vectors */
f2pro_calint(data,elev,elep,dynvar,
             estif,emass,gradopr,etforce,eiforce,
             xyze,funct,functpr,deriv,derivpr,xjm,derxy,derxypr,
             eveln,epren,
             velint,covint,vderxy,pderxy,wa1,
             dirich,deriv2,dirich_onoff);     
/*-------------- calculate the element lumped mass matrix just one time */
if (dynvar->pro_lum==1)	  
   f2pro_lmass(lmass,emass,elev->numnp);
/*------------------------------- calculate element load vector edforce */
if (dynvar->pro_caldirich==1)
{
   switch (dynvar->pro_profile)
   { 
   case 1:
      fluid_pm_caldirich_parabolic(elev,edforce,estif,emass,dynvar->dta,dynvar->theta,hasdirich);
   break;
   case 2:
      fluid_pm_caldirich(elev,edforce,estif,emass,dynvar->dta,dynvar->theta,hasdirich);
   break;
   case 3:
      fluid_pm_caldirich_cyl(elev,edforce,estif,emass,dynvar->dta,dynvar->theta,hasdirich);
   break;
   default:
      dserror("unknown velocity profile!\n");
   } /*end of switch (dynvar->pro_profile)*/
}/*end of if (dynvar->pro_kvv==1)*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2pro_calele */

#endif
/*! @} (documentation module close)*/
