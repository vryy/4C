/*!----------------------------------------------------------------------
\file
\brief element control routine

<pre>
Maintainer: Thomas Hettich
            hettich@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hettich/
            0771 - 685-6575
</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID2 
#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"
#include "fluid2_prototypes.h"
#include "fluid2.h"
/*!---------------------------------------------------------------------                                         
\brief control routine for element integration of fluid2

<pre>                                                        he  12/02

This routine controls the element evaluation:
-actual kapeps variables are set
-stabilisation parameters are calculated
-element integration is performed --> element stiffness matrix and 
                                  --> element load vectors
-element load vector due to dirichlet conditions is calculated				      
			     
</pre>
\param  *data	       FLUID_DATA         (i)
\param  *dynvar	       FLUID_DYN_CALC     (i)
\param  *eleke	       ELEMENT	        (i)   actual element
\param  *elev	       ELEMENT	        (i)   actual element for velocity
\param  *estif_global    ARRAY	        (o)   ele stiffnes matrix
\param  *emass_global    ARRAY	        (o)   ele mass matrix
\param  *etforce_global  ARRAY	        (o)   element time force
\param  *eiforce_global  ARRAY	        (o)   ele iteration force
\param  *edforce_global  ARRAY	        (o)   ele dirichlet force
\param  *eproforce_global  ARRAY	        (o)   ele production force
\param  *hasdirich       INT	              (o)   element flag
\param  *hasext          INT	              (o)   element flag
\param   init	        INT	              (i)   init flag
\return void                                               
                                 
------------------------------------------------------------------------*/
void f2_calele_tu(
                FLUID_DATA     *data, 
                FLUID_DYN_CALC *dynvar, 
	          ELEMENT        *eleke,             
                ELEMENT        *elev,          
                ARRAY          *estif_global,   
                ARRAY          *emass_global,   
	          ARRAY          *etforce_global,       
	          ARRAY          *eiforce_global, 
                ARRAY          *edforce_global,		
                ARRAY          *eproforce_global,		
                INT            *hasdirich,      
                INT            *hasext,
		    INT             init            
	         )
{
INT              hasdead;
static ARRAY     kapepsn_a;     /* element turbulent ken. energy  at (n)    */
static DOUBLE   *kapepsn;
static ARRAY     kapepsg_a;     /* element turbulent ken. energy  at (n+gamma)*/
static DOUBLE   *kapepsg; 
static ARRAY     kapepspro_a;   /* kappa for Production-term              */
static DOUBLE   *kapepspro;
static ARRAY     kappa_a;       /* kappa for epsilon equation                */
static DOUBLE   *kappa;
static ARRAY     kappan_a;      /* kappan for epsilon equation               */
static DOUBLE   *kappan;
static ARRAY     epsilon_a;     /* epsilon for LOW-REYNOLD's MODEL          */
static DOUBLE   *epsilon;
static ARRAY     eddyg_a;       /* element turbulent ken. energy  at (n+gamma)*/
static DOUBLE   *eddyg;
static ARRAY     eddypro_a;     /* element turbulent ken. energy  for prod. term */
static DOUBLE   *eddypro;
static ARRAY     kapepsderxy_a; /* kapeps - derivatives                         */
static DOUBLE   *kapepsderxy;
static ARRAY     kapepsderxy2_a;/* kapeps - 2nd derivatives                     */
static DOUBLE   *kapepsderxy2;
static ARRAY     velint_dc_a;   /* element element velocity for DISC. CAPT. at integ.point */
static DOUBLE   *velint_dc;

static ARRAY     evel_a;        /* element velocities                        */
static DOUBLE  **evel;
static ARRAY     velint_a;      /* element element velocity  at integ.point */
static DOUBLE   *velint;
static ARRAY     vderxy_a;      /* vel  derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     vderxy2_a;     /* 2nd vel  derivatives                         */
static DOUBLE  **vderxy2;

static ARRAY     funct_a;       /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;       /* first natural derivatives                 */
static DOUBLE  **deriv;
static ARRAY     deriv2_a;      /* second natural derivatives                */
static DOUBLE  **deriv2;
static ARRAY     xyze_a;
static DOUBLE  **xyze;   
static ARRAY     xjm_a;         /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     derxy_a;       /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     derxy2_a;      /* 2nd coordinate - derivatives              */
static DOUBLE  **derxy2;
static ARRAY     w1_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa1;      /* used in different element routines        */
static ARRAY     w2_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa2;      /* used in different element routines        */
static DOUBLE  **estif;    /* pointer to global ele-stif                */
static DOUBLE  **emass;    /* pointer to galerkin ele-stif              */
static DOUBLE   *etforce;  /* pointer to Time RHS                       */
static DOUBLE   *eiforce;  /* pointer to Iteration RHS                  */
static DOUBLE   *edforce;  /* pointer to RHS due to dirichl. conditions */
static DOUBLE   *eproforce;  /* pointer to RHS due to dirichl. conditions */


#ifdef DEBUG 
dstrc_enter("f2_calele_tu");
#endif

if (init==1) /* allocate working arrays and set pointers */
{
   evel         = amdef("evel"   ,&evel_a  ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
   kapepsn      = amdef("kapepsn"   ,&kapepsn_a  ,MAXNOD_F2,1,"DV");
   kapepsg      = amdef("kapepsg"   ,&kapepsg_a  ,MAXNOD_F2,1,"DV");
   kapepspro    = amdef("kapepspro"   ,&kapepspro_a  ,MAXNOD_F2,1,"DV");
   kappa        = amdef("kappa"   ,&kappa_a  ,MAXNOD_F2,1,"DV");
   kappan       = amdef("kappan"   ,&kappan_a  ,MAXNOD_F2,1,"DV");
   epsilon      = amdef("epsilon"   ,&epsilon_a  ,MAXNOD_F2,1,"DV");
   eddyg        = amdef("eddyg"   ,&eddyg_a  ,MAXNOD_F2,1,"DV");
   eddypro      = amdef("eddypro"   ,&eddypro_a  ,MAXNOD_F2,1,"DV");
   funct        = amdef("funct"  ,&funct_a  ,MAXNOD_F2,1,"DV");
   deriv        = amdef("deriv"  ,&deriv_a  ,2,MAXNOD_F2,"DA");
   deriv2       = amdef("deriv2" ,&deriv2_a ,3,MAXNOD_F2,"DA");
   xjm          = amdef("xjm"    ,&xjm_a    ,2,2        ,"DA");
   xyze         = amdef("xyze"     ,&xyze_a     ,2,MAXNOD_F2,"DA");
   velint       = amdef("velint" ,&velint_a ,1,1,"DV");
   velint_dc    = amdef("velint_dc" ,&velint_dc_a ,1,1,"DV");
   kapepsderxy  = amdef("kapepsderxy" ,&kapepsderxy_a ,2,1,"DV");
   kapepsderxy2 = amdef("kapepsderxy2",&kapepsderxy2_a,3,1,"DV");
   derxy   = amdef("derxy"  ,&derxy_a  ,2,MAXNOD_F2,"DA");
   derxy2  = amdef("derxy2" ,&derxy2_a ,3,MAXNOD_F2,"DA");
   vderxy  = amdef("vderxy" ,&vderxy_a ,2,MAXNOD_F2,"DA");
   vderxy2 = amdef("vderxy2" ,&vderxy2_a ,3,MAXNOD_F2,"DA");
   wa1     = amdef("wa1"    ,&w1_a     ,30,30,"DA");
   wa2     = amdef("wa2"    ,&w2_a     ,30,30,"DA");  
/*                                        \- size is arbitrary chosen!  */
   estif   = estif_global->a.da;
   emass   = emass_global->a.da;
   eiforce = eiforce_global->a.dv;
   etforce = etforce_global->a.dv;
   edforce = edforce_global->a.dv;
   eproforce = eproforce_global->a.dv;

   goto end;
} /* endif (init==1) */

/*------------------------------------------------ initialise with ZERO */
amzero(estif_global);
amzero(emass_global);
amzero(eiforce_global);
amzero(etforce_global);
amzero(edforce_global);
amzero(eproforce_global);
*hasdirich=0;
*hasext=0;

/*---------------------------------------------------- set element data */
f2_calset_tu(dynvar,data,eleke,elev,kapepsn,kapepsg,kapepspro,eddyg,eddypro,
             kappa,kappan,epsilon,evel,xyze);

/*-------------- calculate stab-parameter and parameter for DISC. CAPT. */
f2_calelesize_tu(eleke,elev,dynvar,data,funct,deriv,deriv2,evel,eddyg,velint,
                 velint_dc,kapepsn,xjm,xyze,derxy,kapepsderxy,wa1);

/*-------------------------------- calculate element stiffness matrices */
/*                                            and element force vectors */
f2_calint_tu(data,eleke,elev,dynvar,estif,emass,etforce,eiforce,eproforce,
	       funct,deriv,deriv2,xjm,xyze,derxy,derxy2,kapepsn,kapepsg,eddyg,eddypro,
             kappa,kappan,epsilon,kapepspro,kapepsderxy,kapepsderxy2,velint,
             velint_dc,evel,vderxy,vderxy2,wa1,wa2);

/*--------------------------------------- add emass and estif to estif  */

f2_estifadd_tu(estif,emass,wa1,eleke->numnp,dynvar);

/*------------------------------- calculate element load vector edforce */
fluid_caldirich_tu(dynvar,eleke,edforce,estif,hasdirich);

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of f2_calele */


#endif
