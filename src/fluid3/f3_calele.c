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
#include "../fluid_full/fluid_prototypes.h"
#include "fluid3_prototypes.h"
#include "fluid3.h"
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
static ARRAY     eveln_a;  /* element velocities at (n) 		*/
static DOUBLE  **eveln;
/*NOTE: if there is no classic time rhs (as described in WAW) the array
         eveln is misused and does NOT contain the velocity at time (n)
	 but rather a linear combination of old velocities and 
	 accelerations depending upon the time integration scheme!!!!! */
static ARRAY     evelng_a; /* element velocities at (n+gamma)		*/
static DOUBLE  **evelng;
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
static ARRAY     velint_a; /* velocities at integration point		*/
static DOUBLE   *velint;
static ARRAY     vel2int_a;/* velocities at integration point		*/
static DOUBLE   *vel2int;
static ARRAY	 covint_a; /* convective velocities at integr. point	*/
static DOUBLE   *covint;
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
static ARRAY     w1_a;     /* working array of arbitrary chosen size	*/
static DOUBLE  **wa1;      /* used in different element routines	*/
static ARRAY     w2_a;     /* working array of arbitrary chosen size	*/
static DOUBLE  **wa2;      /* used in different element routines	*/
static DOUBLE  **estif;    /* pointer to global ele-stif		*/
static DOUBLE  **emass;    /* pointer to galerkin ele-stif		*/
static DOUBLE   *etforce;  /* pointer to Time RHS			*/
static DOUBLE   *eiforce;  /* pointer to Iteration RHS  		*/
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
\param  *data	         FLUID_DATA     (i)
\param  *ele	         ELEMENT	(i)   actual element
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *etforce_global  ARRAY	        (o)   element time force
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param  *hasdirich       INT	        (o)   element flag
\param  *hasext          INT	        (o)   element flag
\param   init	         INT	        (i)   init flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f3_calele(
                FLUID_DATA     *data, 
	        ELEMENT        *ele,
                ARRAY          *estif_global,
                ARRAY          *emass_global, 
	        ARRAY          *etforce_global,
	        ARRAY          *eiforce_global,
		ARRAY          *edforce_global,	
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
   eveln   = amdef("evln"    ,&eveln_a   ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
   evelng  = amdef("evelng"  ,&evelng_a  ,NUM_F3_VELDOF,MAXNOD_F3,"DA");
   epren   = amdef("epren"   ,&epren_a   ,MAXNOD_F3,1,"DV");
   edeadn  = amdef("edeadn"  ,&edeadn_a  ,3,1,"DV");
   edeadng = amdef("edeadng" ,&edeadng_a ,3,1,"DV"); 
   funct   = amdef("funct"   ,&funct_a   ,MAXNOD_F3,1,"DV");
   deriv   = amdef("deriv"   ,&deriv_a   ,3,MAXNOD_F3,"DA");
   deriv2  = amdef("deriv2"  ,&deriv2_a  ,6,MAXNOD_F3,"DA");
   xyze    = amdef("xyze"    ,&xyze_a    ,3,MAXNOD_F3,"DA");
   xjm     = amdef("xjm"     ,&xjm_a     ,3,3        ,"DA");
   velint  = amdef("velint"  ,&velint_a  ,NUM_F3_VELDOF,1,"DV");
   vel2int = amdef("vel2int" ,&vel2int_a ,NUM_F3_VELDOF,1,"DV");
   covint  = amdef("covint"  ,&covint_a  ,NUM_F3_VELDOF,1,"DV");
   vderxy  = amdef("vderxy"  ,&vderxy_a  ,3,3,"DA");
   pderxy  = amdef("pderxy"  ,&pderxy_a  ,3,1,"DV");
   vderxy2 = amdef("vderxy2" ,&vderxy2_a ,3,6,"DA");
   derxy   = amdef("derxy"   ,&derxy_a   ,3,MAXNOD_F3,"DA");
   derxy2  = amdef("derxy2"  ,&derxy2_a  ,6,MAXNOD_F3,"DA");
   sigmaint= amdef("sigmaint",&sigmaint_a,MAXGAUSS ,6,"DA");
   wa1     = amdef("wa1"     ,&w1_a      ,300,300        ,"DA");
   wa2     = amdef("wa2"     ,&w2_a      ,300,300        ,"DA");  
/*                                        \- size is chosen arbitrarily! */
   estif   = estif_global->a.da;
   emass   = emass_global->a.da;
   eiforce = eiforce_global->a.dv;
   etforce = etforce_global->a.dv;
   edforce = edforce_global->a.dv;
   
   fdyn    = alldyn[genprob.numff].fdyn;
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

/*---------------------------------------------------- set element data */
f3_calset(ele,eveln,evelng,epren,edeadn,edeadng,hasext);

/*------------------------- calculate element size and stab-parameter: */
f3_calelesize(ele,data,funct,deriv,deriv2,derxy,xjm,evelng,velint,wa1);

/*------------------------------- calculate element stiffness matrices */
/*                                           and element force vectors */
f3_calint(data,ele,hasext,
          estif,emass,etforce,eiforce,
	  funct,deriv,deriv2,xjm,derxy,derxy2,
	  eveln,evelng,epren,edeadn,edeadng,
	  velint,vel2int,covint,vderxy,pderxy,vderxy2,
	  wa1,wa2);

/*----------------- add emass and estif to estif and permute the matrix */
f3_permestif(estif,emass,wa1,ele->numnp);

/*--------------------------------- permute element load vector etforce */
if (fdyn->nif!=0)
   f3_permeforce(etforce,wa1,ele->numnp);

/*--------------------------------- permute element load vector eiforce */
if (fdyn->nii+(*hasext)!=0)
   f3_permeforce(eiforce,wa1,ele->numnp);

/*------------------------------- calculate element load vector edforce */
fluid_caldirich(ele,edforce,estif,hasdirich,3);

/*----------------------------------------- calculate emass * vec(n) ---*/
if (fdyn->nim)
   f3_massrhs(ele,emass,eveln,eiforce); 
   
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f3_calele */





/*!--------------------------------------------------------------------- 
  \brief control routine for stress calculation in f3

  <pre>                                                         mn 03/04

  </pre>

  \param     str       FLUID_STRESS   (i)    flag for stress calculation
  \param     viscstr   INT            (i)    viscose stresses yes/no?
  \param    *data      FLUID_DATA     (i)
  \param    *ele       ELEMENt        (i)    actual element 
  \param     is_relax  INT            (i)    flag

  \return void                                               

  ------------------------------------------------------------------------*/
void f3_stress(
    FLUID_STRESS  str, 
    INT           viscstr,
    FLUID_DATA   *data, 
    ELEMENT      *ele,
    INT           is_relax
    )
{
  INT       i;
  INT       coupled;      /* flag for fsi interface element */
  INT       iel;          /* number of nodes per element */
  GNODE    *actgnode;     /* actual gnode */


#ifdef DEBUG 
  dstrc_enter("f3_stress");
#endif

  switch(str)
  {
    case str_none: /* do nothing */
      break;


#if 0
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
        f3_calfsistress(viscstr,data,ele,eveln,epren,funct,
            deriv,derxy,vderxy,xjm,xyze,sigmaint,is_relax);      
      break;
#endif
#endif


    case str_liftdrag:
    case str_all:
      f3_calelestress(viscstr,data,ele,eveln,epren,funct,
          deriv,derxy,vderxy,xjm,wa1,xyze,sigmaint);      
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


#endif
