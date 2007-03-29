/*!----------------------------------------------------------------------
\file
\brief service routines for fluid2_TDS element

<pre>
Maintainer: 
            
            
            
</pre>

------------------------------------------------------------------------*/


/*!
\addtogroup FLUID2
*//*! @{ (documentation module open)*/
#ifdef D_FLUID2
#ifdef D_FLUID2_TDS
#include "../headers/standardtypes.h"
#include "fluid2_TDS_prototypes.h" 

#include "../fluid2/fluid2.h"
#include "../fluid2/fluid2_prototypes.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of material laws                                              |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _MATERIAL  *mat;
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

/*!---------------------------------------------------------------------
\brief set data arrays for element calculation

<pre>                                                        gammi 12/06

get the element velocities, accelerations, coordinates at different times
</pre>
\param   *ele      ELEMENT          (i)  actual element
\param  **xyze     DOUBLE           (o)  nodal coordinates
\param  **eaccng   DOUBLE           (o)  ele accs at time n+alpha_M
\param  **evelng   DOUBLE           (o)  ele vels at time n+alpha_F
\param   *epreng   DOUBLE           (o)  ele pres at time n+1
\param   *edeadng  DOUBLE           (o)  ele dead load at n+g (selfweight)
\param   *ipos     ARRAY_POSITION   (i)  node array positions
\param    visc     DOUBLE           (o)  kinematic viscosity
\return void

------------------------------------------------------------------------*/
void f2_inc_gen_alpha_calset(
	        ELEMENT         *ele,
                DOUBLE         **xyze,
                DOUBLE         **eaccng,
	        DOUBLE         **evelng,
		DOUBLE          *epreng,
		DOUBLE          *edeadng,
                ARRAY_POSITION *ipos,
		DOUBLE         *visc
	      )
{
INT    i;           /* simply some counters                             */
INT    actcurve;    /* actual time curve                                */

int    velnp;
int    velnm;
int    accnm;

DOUBLE acttimefac;  /* time factor from actual curve                    */
DOUBLE acttimefacn; /* time factor at time (n)                          */
DOUBLE acttime;
NODE  *actnode;     /* actual node                                      */
GSURF *actgsurf;


FLUID_DYNAMIC *fdyn;

#ifdef DEBUG
dstrc_enter("f2_inc_gen_alpha_calset");
#endif

/*------------------------------------------------------- initialise ---*/
fdyn  = alldyn[genprob.numff].fdyn;

velnp = ipos->velnp;
velnm = ipos->velnm;
accnm = ipos->accnm;

/*---------------------------------------------------- get viscosity ---*/
visc[0] = mat[ele->mat-1].m.fluid->viscosity;



for(i=0;i<ele->numnp;i++) /* loop nodes of element */
{
   actnode=ele->node[i];

   /*------------------------------------------ set element coordinates */
   xyze[0][i]  =actnode->x[0];
   xyze[1][i]  =actnode->x[1];
   
   /*---------------------------- set element accelerations (n+alpha_M) */
   eaccng[0][i]=actnode->sol_increment.a.da[accnm][0];
   eaccng[1][i]=actnode->sol_increment.a.da[accnm][1];
   
   /*------------------------------- set element velocities (n+alpha_F) */
   evelng[0][i]=actnode->sol_increment.a.da[velnm][0];
   evelng[1][i]=actnode->sol_increment.a.da[velnm][1];

   /*------------------------------- set element pressure (n+1) -------*/
   epreng   [i]=actnode->sol_increment.a.da[velnp][2];
} /* end of loop over nodes of element */


/*------------------------------------------------ check for dead load */
actgsurf = ele->g.gsurf;
if (actgsurf->neum!=NULL)
{
   actcurve = actgsurf->neum->curve-1;
   if (actcurve<0)
   {
      acttimefac =ONE;
      acttimefacn=ONE;
   }
   else
   {
      acttime = fdyn->acttime;
      dyn_facfromcurve(actcurve,acttime,&acttimefac) ;
      acttime = fdyn->acttime-fdyn->dta;
      dyn_facfromcurve(actcurve,acttime,&acttimefacn) ;
   }
   for (i=0;i<2;i++)
   {

      if (actgsurf->neum->neum_onoff.a.iv[i]==0)
      {
         edeadng[i] = ZERO;
      }

      

      if (actgsurf->neum->neum_type==neum_dead  &&
          actgsurf->neum->neum_onoff.a.iv[i]!=0)
      {
         edeadng[i] = actgsurf->neum->neum_val.a.dv[i]*acttimefac;
      }
   }


}
else
{
   for (i=0;i<2;i++)
   {
       edeadng[i] = ZERO;
   } 
}

/*---------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of f2_calset */






/*!---------------------------------------------------------------------
\brief storing results in solution history

<pre>                                                        gammi 12/06

This is the nonlinear iteration update

                        i -> i+1

  o Add incremental accelerations in DIST_VECTOR to the node
    accleration (time n+1)
  o Update acceleration  at time n+alpha_M
  o Update velocity at time time n+1 and n+alpha_F
  o Update pressure at time time n+1

in this routine the results in the DIST_VECTOR are put to the nodes in
a certain place in ARRAY sol_increment.
Result has to be allreduced and is put to the whole field on each proc.
If necassary the norms for the iteration convergence check of the
nonlinear iteration scheme are calculated.

</pre>
\param **actfield      FIELD	      (i)    actual field
\param  *actintra      INTRA	      (i)    actual intra comm.
\param	*sol 	       DIST_VECTOR    (i)    solution vector
\param  *ipos          ARRAY_POSITION (i)    where to store acc, velo etc
\param	*sysarray      SPARSE_ARRAY   (i)
\param	*sysarray_typ  SPARSE_TYP     (i)
\param	*vrat          DOUBLE	      (o)    vel.  conv. ratio
\param	*prat          DOUBLE	      (o)    pre.  conv. ratio
\param  *grat          DOUBLE         (o)    grid  conv. ratio
\return void

------------------------------------------------------------------------*/
void fluid_result_incre_for_genalpha(
                          FIELD             *actfield,
                          INT                disnum,
                          INTRA             *actintra,
			  DIST_VECTOR       *sol,
			  DIST_VECTOR       *rhs,
                          ARRAY_POSITION    *ipos,
			  SPARSE_ARRAY      *sysarray,
			  SPARSE_TYP        *sysarray_typ,
			  DOUBLE            *vrat,
			  DOUBLE            *prat,
                          DOUBLE            *grat
		       )
{
INT      i,j;          /* simply some counters                         */
INT      dof;          /* actual dof number                            */
INT      predof;       /* number of pressure dof (2=2D; 3=3D)          */
INT      numeq_total;  /* total number of equations                    */
int      veln ,accn ;
int      velnm,accnm;
int      velnp,accnp;
DOUBLE   dvnorm=0;
DOUBLE   dpnorm=0;
DOUBLE   dgnorm=0;
DOUBLE    vnorm=0;
DOUBLE    pnorm=0;  /* values for norm calculation                  */
DOUBLE    gnorm=0;
NODE    *actnode;      /* actual node                                  */
static ARRAY    result_a;
static DOUBLE  *result;       /* redundandent result vector                   */

FLUID_DYNAMIC *fdyn;

double  alpha_M,alpha_F,theta;
double  dt;

#ifdef DEBUG
dstrc_enter("fluid_result_incre_or_genalpha");
#endif

fdyn = alldyn[genprob.numff].fdyn;
/*----------------------------------------------------------------------*/
numeq_total = sol->numeq_total;
predof      = fdyn->numdf-1;

veln        = ipos->veln;
accn        = ipos->accn;

velnm       = ipos->velnm;
accnm       = ipos->accnm;

velnp       = ipos->velnp;
accnp       = ipos->accnp;

theta       = fdyn->theta;
alpha_M     = fdyn->alpha_m;
alpha_F     = fdyn->alpha_f;

dt          = fdyn->dt;

/*------------------------- allocate space to allreduce the DIST_VECTOR */
if(result_a.Typ==cca_XX)
   result = amdef("result",&result_a,numeq_total,1,"DV");

if (result_a.fdim < numeq_total) {
  amdel(&result_a);
  result = amdef("result",&result_a,numeq_total,1,"DV");
}
amzero(&result_a);

/*------------------ copy distributed result to redundant result vector */
solserv_reddistvec(
                      sol,
                      sysarray,
                      sysarray_typ,
                      result,
                      sol->numeq_total,
                      actintra
                     );


switch (fdyn->itnorm) /* switch to norm */
{
/*-------------------------------------------------------------------------*/
case fncc_L2: /* L_2 norm */
   /*-----------  loop nodes and put the result back to the node structure */
   for (i=0; i<actfield->dis[disnum].numnp; i++)
   {
      actnode = &(actfield->dis[disnum].node[i]);
      /*------------------------------ enlarge sol_increment, if necessary */
      for (j=0; j<actnode->numdf; j++) /* loop dofs and calculate the norms */
      {
         dof = actnode->dof[j];
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
         if (actnode->gnode->dirich!=NULL &&
             actnode->gnode->dirich->dirich_onoff.a.iv[j]!=0)
           continue;
#else
         if (dof>=numeq_total)
	 {
	     continue;
	 }
#endif
         if (j==predof) /* pressure dof */
         {
            pnorm  += DSQR(result[dof]+actnode->sol_increment.a.da[velnp][j]);
            dpnorm += DSQR(result[dof]);

            /*
	      o Update pressure at time time n+1
	    */

	    actnode->sol_increment.a.da[velnp][j] += result[dof];
	 } /* endif pressure dof */
         else /* acceleration and vel - dof */
         {

	    vnorm  += DSQR(theta*dt*result[dof]+actnode->sol_increment.a.da[velnp][j]);
            dvnorm += DSQR(theta*dt*result[dof]);
	    	    
            /*
	      o Add incremental accelerations in DIST_VECTOR to the node
	        accleration (time n+1)
	    */

            actnode->sol_increment.a.da[accnp][j] += result[dof];
	    
            /*
	      o Update acceleration  at time n+alpha_M
	    */

            actnode->sol_increment.a.da[accnm][j] =
		actnode->sol_increment.a.da[accn ][j]
		+
		(alpha_M  )*
		(actnode->sol_increment.a.da[accnp][j]
		 -
		 actnode->sol_increment.a.da[accn ][j]);
	    	    
            /*
	      o Update velocity at time time n+1
	    */
	    
	    actnode->sol_increment.a.da[velnp][j]+=
		theta*dt* result[dof];
		/*actnode->sol_increment.a.da[veln][j]
		+
		dt*actnode->sol_increment.a.da[accn][j]*/
		/*(actnode->sol_increment.a.da[accnp][j]
		 -
		 actnode->sol_increment.a.da[accn][j])*/

            /*
	      o Update velocity at time time n+alpha_F
	    */
	    
	    actnode->sol_increment.a.da[velnm][j]=
		actnode->sol_increment.a.da[veln ][j]
		+
		(alpha_F  )*
		(actnode->sol_increment.a.da[velnp][j]
		 -
		 actnode->sol_increment.a.da[veln ][j]);

	    
         } /* endif acceleration and vel - dof */
      } /* end of loop over dofs */
   } /* end of loop over nodes */
   dvnorm = sqrt(dvnorm);
    vnorm = sqrt( vnorm);
   dpnorm = sqrt(dpnorm);
    pnorm = sqrt( pnorm);
   dgnorm = sqrt(dgnorm);
    gnorm = sqrt( gnorm);
break;
/*-------------------------------------------------------------------*/
default:
   dserror("unknown norm for convergence check!\n");
} /* end of switch(fdyn->itnorm) */

/*------------------------------------------- check for "ZERO-field" */
if (vnorm<EPS5)
{
   vnorm = ONE;
#ifdef DEBUG
    printf("ATTENTION: zero vel field - norm <= 1.0e-5 set to 1.0!! \n");
#endif
}
/* If prat==NULL we are not interessted in the pressure norm. This is
 * the case in the projection method where we use discontinous
 * pressures. */
if (prat!=NULL)
{
  if (pnorm<EPS5)
  {
    pnorm = ONE;
#ifdef DEBUG
    printf("ATTENTION: zero pre field - norm <= 1.0e-5 set to 1.0!! \n");
#endif
  }
}
if (gnorm<EPS5)
{
   gnorm = ONE;
}
/*------------------------------------- set final convergence ratios */
*vrat = dvnorm/vnorm;
if (prat!=NULL)
  *prat = dpnorm/pnorm;


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif

return;
} /* end of  fluid_result_incre_for_genalpha */

/*!---------------------------------------------------------------------
\brief estimate new trial values for new timestep

<pre>                                                        gammi 11/06

Assuming constant velocity and pressure, we have the new trial values

              +-              -+   +-      -+
              | ux^{n+1}_{(0)} |   | ux^{n} |
              |                | = |        |
              | uy^{n+1}_{(0)} |   | uy^{n} |
              +-              -+   +-      -+
              
                 p^{n+1}_{(0)}=p^{n}


and according to this we have for the accelerations


              +-              -+            +-      -+
              | ax^{n+1}_{(0)} |   gamma-1  | ax^{n} |
              |                | = -------  |        |
              | ay^{n+1}_{(0)} |    gamma   | ay^{n} |
              +-              -+            +-      -+
                 

</pre>
\param  *actpart       PARTITION        (i)
\param  *actintra      INTRA            (i)
\param  *actfield      FIELD            (i)
\param  *ipos          ARRAY_POSITION   (i)
\param   disnum_calc   INT              (i)                
\return void

------------------------------------------------------------------------*/

void f2_estimate_new_trial_values_for_inc_gen_alpha(
     PARTITION      *actpart,
     INTRA          *actintra,
     FIELD          *actfield,
     ARRAY_POSITION *ipos,
     INT             disnum_calc)
{
INT      nn,ne,j,i;       /* simply some counters                    */

/* element related data */
ELEMENT *actele;       /* actual element                               */

int      nir=0,nis=0;       /* gausspoint numbers */
int      icode=0,intc=0;    /* required for tri elements            */
int      ihoel=0;           /* higher order element flag            */
double   det;

DOUBLE    e1,e2;      /* natural coordinates of integr. point           */

int      lr,ls;       /* Gausspoint counters */

double   visc;

FLUID_DATA      *data;

double  fact;

/* node related data */

NODE    *actnode;      /* actual node                                  */

/* positions of values in sol_increment */
int      veln ,accn ; 
int      velnm,accnm;
int      velnp,accnp;

int      predof;      /* pressure degree of freedom number */

/* gauss point related data */

/* the residual of the momentum equation without body force */
double   res_mod[2];

/* the time increment of the subscale acceleration */
double   tin[2]; 


/* higher order terms */
DOUBLE  hot   [2];

/* pressure gradient */ 
DOUBLE  gradp [2];

/* intermediate acceleration (n+alpha_M) and velocities (n+alpha_F)    */
double   accint[2],velint[2];

/* divergence of velocity */
double   divu;

/* time integration */
FLUID_DYNAMIC   *fdyn;

double  alpha_M,alpha_F,theta;

double  dt;

double  aftdt;

/* stabilisation parameters */
double   tau_M;
double   tau_C;


/* abbreviations containing the stabilisation parameters */
double  C1,C2,C3;
double  facM;
double  amtau_M;

/* the old subscale velocities and accelerations */
double   sv_old    [2];
double   sv_acc_old[2];

/* the new subscale velocities and accelerations */
double   sv_new        [2];
double   sv_acc_new    [2];
double   sv_acc_mod_new[2];

DOUBLE  viscs2[2][2*MAXNOD]; /* viscous term incluiding 2nd derivatives */


static ARRAY     xyze_a;   /* actual element coordinates                */
static DOUBLE  **xyze;
static ARRAY     xjm_a;    /* jocobian matrix                           */
static DOUBLE  **xjm;
static ARRAY     funct_a;  /* shape functions                           */
static DOUBLE   *funct;
static ARRAY     deriv_a;  /* first natural derivatives                 */
static DOUBLE  **deriv;
static ARRAY     eveln_a;  /* element velocities at (n)                 */
static DOUBLE  **eveln;
static ARRAY     evelng_a; /* element velocities at (n+gamma)           */
static DOUBLE  **evelng;
static ARRAY     eaccng_a;  /* element accelerations at (n+alpha_M)     */
static DOUBLE  **eaccng;
static ARRAY     derxy_a;  /* coordinate - derivatives                  */
static DOUBLE  **derxy;
static ARRAY     vderxy_a; /* vel - derivatives                         */
static DOUBLE  **vderxy;
static ARRAY     deriv2_a; /* second natural derivatives                */
static DOUBLE  **deriv2;
static ARRAY     w1_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa1;      /* used in different element routines        */
static ARRAY     w2_a;     /* working array of arbitrary chosen size    */
static DOUBLE  **wa2;      /* used in different element routines        */
static ARRAY     derxy2_a; /* 2nd coordinate - derivatives              */
static DOUBLE  **derxy2;
static ARRAY     vderxy2_a;/* vel - 2nd derivatives                     */
static DOUBLE  **vderxy2;
static ARRAY     epreng_a; /* element pressures at (n)	                */
static DOUBLE   *epreng;
static DOUBLE   *edeadng;
static ARRAY     edeadng_a;/* element dead load (selfweight)            */

#ifdef DEBUG
dstrc_enter("f2_estimate_new_trial_values_for_inc_gen_alpha");
#endif


fdyn = alldyn[genprob.numff].fdyn;

data   = fdyn->data;

/*----------------------------------------------------------------------*/
predof      = 2;


veln        = ipos->veln;
accn        = ipos->accn;

velnm       = ipos->velnm;
accnm       = ipos->accnm;

velnp       = ipos->velnp;
accnp       = ipos->accnp;

theta       = fdyn->theta;
alpha_M     = fdyn->alpha_m;
alpha_F     = fdyn->alpha_f;

fact        = (theta-1)/theta;

dt          = fdyn->dt;
aftdt       = alpha_F*dt*theta;


xyze      = amdef("xyze"     ,&xyze_a     ,2,MAXNOD_F2,"DA");
xjm       = amdef("xjm"      ,&xjm_a      ,2,2        ,"DA");
funct     = amdef("funct"    ,&funct_a    ,MAXNOD_F2,1,"DV");
deriv     = amdef("deriv"    ,&deriv_a    ,2,MAXNOD_F2,"DA");
eveln     = amdef("eveln"    ,&eveln_a    ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
evelng    = amdef("evelng"   ,&evelng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
eaccng    = amdef("eaccng"   ,&eaccng_a   ,NUM_F2_VELDOF,MAXNOD_F2,"DA");
vderxy    = amdef("vderxy"   ,&vderxy_a   ,2,2,"DA");
derxy     = amdef("derxy"    ,&derxy_a    ,2,MAXNOD_F2,"DA");
deriv2    = amdef("deriv2"   ,&deriv2_a   ,3,MAXNOD_F2,"DA");
wa1       = amdef("wa1"      ,&w1_a       ,MAXDOFPERELE,MAXDOFPERELE,"DA");
wa2       = amdef("wa2"      ,&w2_a       ,MAXDOFPERELE,MAXDOFPERELE,"DA");
derxy2    = amdef("derxy2"   ,&derxy2_a   ,3,MAXNOD_F2,"DA");
vderxy2   = amdef("vderxy2"  ,&vderxy2_a  ,2,3,"DA");
epreng    = amdef("epreng"   ,&epreng_a   ,MAXNOD_F2,1,"DV");
edeadng   = amdef("edeadng"  ,&edeadng_a  ,2,1,"DV");

/*----------- set initial trial values for the large scale quantities */
/*                                      assuming constant velocities| */

for (nn=0;nn<actfield->dis[disnum_calc].numnp;nn++)
{
    actnode=&actfield->dis[disnum_calc].node[nn];

    for (j=0;j<2;j++)
    {
	/*----------------------------------------- estimate new velocities */
	/* the intermediate solution is the same since we estimate
	 * a constant velocity --- 'cept for dirichlet boundaries!!!!       */
	actnode->sol_increment.a.da[velnm][j]=
	    actnode->sol_increment.a.da[veln ][j]
	    +
	    alpha_F*(actnode->sol_increment.a.da[velnp][j]
		     -
		     actnode->sol_increment.a.da[veln ][j]);


	
	/* estimate new accelerations  */
	actnode->sol_increment.a.da[accnp][j]
	    =(actnode->sol_increment.a.da[velnp][j]-actnode->sol_increment.a.da[veln][j])/(theta*dt)
	    +fact*actnode->sol_increment.a.da[accn][j];

	
	/* calculate the estimated intermediate accelerations */
	actnode->sol_increment.a.da[accnm][j]=
	    actnode->sol_increment.a.da[accn ][j]
	    +
	    alpha_M*(actnode->sol_increment.a.da[accnp][j]
		     -
		     actnode->sol_increment.a.da[accn ][j]);
	
    }

    /*----------------------------------------- estimate new pressure */
    actnode->sol_increment.a.da[velnp][predof]=
	actnode->sol_increment.a.da[veln ][predof];
}


/*-------------- set initial trial values for the subscale quantities */

for (ne=0;ne<actpart->pdis[disnum_calc].numele;ne++)
{
    actele=actpart->pdis[disnum_calc].element[ne];

        /*---------------------------------------------- get viscosity ---*/
    visc = mat[actele->mat-1].m.fluid->viscosity;
    
    /*------- get integraton data and check if elements are "higher order" */
    switch (actele->distyp)
    {
	case quad4: case quad8: case quad9:  /* --> quad - element */
	    icode   = 3; /* flag for higher order elements                 */
	    ihoel   = 1; /* flag for eveluation of shape functions         */
	    /* initialise integration */
	    nir = actele->e.f2->nGP[0];
	    nis = actele->e.f2->nGP[1];
	    break;
	case tri6: /* --> tri - element */
	    icode   = 3; /* flag for higher order elements                 */
	    ihoel   = 1; /* flag for eveluation of shape functions         */
	    /* initialise integration */
	    nir  = actele->e.f2->nGP[0];
	    nis  = 1;
	    break;
	case tri3:
	    ihoel  =0;  /* flag for higher order elements                 */
	    icode  =2;  /* flag for eveluation of shape functions         */
	    /* initialise integration */
	    nir  = actele->e.f2->nGP[0];
	    nis  = 1;
	    break;
	default:
	    dserror("typ unknown!");
    } /* end switch(typ) */


    /*------------------------------------------ set element coordinates -*/
    f2_inc_gen_alpha_calset(
	actele,
	xyze,
	eaccng,
	evelng,
	epreng,
	edeadng,
	ipos,
	&visc
	);


    /*--------------------------------------------- stab-parameter ---*/
    f2_get_time_dependent_sub_tau(actele,xyze,funct,deriv,evelng,NULL,visc);
    
    tau_M  = fdyn->tau[0];

    tau_C  = fdyn->tau[2];

    C1 = (alpha_M-theta)*dt*(tau_C/(alpha_M*tau_C+aftdt));
    C2 = ((alpha_F-1.0)*dt*theta+alpha_M*tau_C)/(alpha_M*tau_C+aftdt);   
    C3 = (tau_C/(alpha_M*tau_C+aftdt))*theta*dt;

    amtau_M = alpha_M*tau_M;
    facM    = 1./(amtau_M+aftdt);


    
    
    /* loop Gausspoints */
    for (lr=0;lr<nir;lr++)
    {
	for (ls=0;ls<nis;ls++)
	{
	    /* set subscale data at gausspoints*/

	    
            /*------- get integraton data and check if elements are
	                                                 "higher order" */
	    switch (actele->distyp)
	    {
		case quad4: case quad8: case quad9:  /* --> quad - element */
		    e1   = data->qxg[lr][nir-1];
		    e2   = data->qxg[ls][nis-1];
		    f2_rec(funct,deriv,deriv2,e1,e2,actele->distyp,icode);
		    break;
		case tri6: /* --> tri - element */
		    e1   = data->txgr[lr][intc];
		    e2   = data->txgs[lr][intc];
		    f2_tri(funct,deriv,deriv2,e1,e2,actele->distyp,icode);
		    break;
		case tri3:
		    e1   = data->txgr[lr][intc];
		    e2   = data->txgs[lr][intc];
		    f2_tri(funct,deriv,deriv2,e1,e2,actele->distyp,icode);
		    break;
		default:
		    dserror("typ unknown!");
	    } /* end switch(typ) */
   
	    for(i=0;i<2;i++) 
	    {
		sv_old    [i]=actele->e.f2->sub_vel.a.da    [i][lr*nis+ls];

		sv_acc_old[i]=actele->e.f2->sub_vel_acc.a.da[i][lr*nis+ls]
		    +1./alpha_M*edeadng[i];
	    }
   
	    /*------------------ compute Jacobian matrix at time n+1 ---*/
	    f2_jaco(xyze,deriv,xjm,&det,actele->numnp,actele);

	    /*----------------------------- compute global derivates ---*/
	    f2_gder(derxy,deriv,xjm,det,actele->numnp);

	    /*------ get velocities (n+alpha_F) at integration point ---*/
	    f2_veci(velint,funct,evelng,actele->numnp);

	    /*--- get accelerations (n+alpha_M) at integration point ---*/
	    f2_veci(accint,funct,eaccng,actele->numnp);
	    
	    /*--- get velocity (n+alpha_F,i) derivatives at integration
	                                                          point */
	    f2_vder(vderxy,derxy,evelng,actele->numnp);

	    if (ihoel!=0)
	    {
		f2_gder2(xyze,xjm,wa1,wa2,derxy,derxy2,deriv2,actele->numnp);
		f2_vder2(vderxy2,derxy2,evelng,actele->numnp);
	    }

	    /* compute residual of continuity equation */
	    divu     = vderxy[0][0] + vderxy[1][1];

	    /*------------------------------ get higher order terms ---*/
	    if(ihoel!=0)
	    {
		hot    [0]=0.5 * (2.0*vderxy2[0][0]
				  +
				  (vderxy2[0][1] + vderxy2[1][2]));
		hot    [1]=0.5 * (2.0*vderxy2[1][1]
				  +
				  (vderxy2[1][0] + vderxy2[0][2]));
	    }	
	    else
	    {
		hot    [0]=0;
		hot    [1]=0;
	    }

	    
	    /*------------------------------- get pressure gradients ---*/
	    gradp[0] = gradp[1] = 0.0;
	    
	    for (i=0; i<actele->numnp; i++)
	    {
		gradp[0] += derxy[0][i] * epreng[i];
		gradp[1] += derxy[1][i] * epreng[i];
	    }

	    /*---------------------------------- calculate the residual */
	    for(i=0;i<2;i++) 
	    {
		/* residual without body force */
		res_mod[i] = accint[i] + 0 - 2*visc*hot[i]
		    + gradp[i];

		
		sv_acc_mod_new[i] = - facM * sv_old[i]
		    - facM * (tau_M+alpha_F*dt-amtau_M-aftdt) * sv_acc_old[i]
		    - facM * (tau_M*res_mod[i]+aftdt/alpha_M*edeadng[i]);

		actele->e.f2->sub_vel_acc_trial.a.da[i][lr*nis+ls]
		    =sv_acc_mod_new[i];

		sv_acc_new[i] = sv_acc_mod_new[i]+ 1./alpha_M*edeadng[i];


		sv_new[i] = sv_old[i] + dt*sv_acc_old[i]
		    + dt*theta *(sv_acc_new[i]-sv_acc_old[i]);
		
		actele->e.f2->sub_vel_trial.a.da[i][lr*nis+ls]=sv_new [i];
	    }

	    /* estimate subscale pressure (accelerations) according to
	     * the prescribed values of the large scale quantities
	     * (constant velocity/pressure estimate or values from
	     * dirichlet boundary conditions) */

	    actele->e.f2->sub_pres_trial.a.dv[lr*nis+ls]
		=
		C1*actele->e.f2->sub_pres_acc.a.dv[lr*nis+ls]
		+
		C2*actele->e.f2->sub_pres.a.dv[lr*nis+ls]
		-
		C3*divu;

	    actele->e.f2->sub_pres_acc_trial.a.dv[lr*nis+ls]=
		(theta-1.)/theta*actele->e.f2->sub_pres_acc.a.dv[lr*nis+ls]
		+
		(actele->e.f2->sub_pres_trial.a.dv[lr*nis+ls]
		 -
		 actele->e.f2->sub_pres      .a.dv[lr*nis+ls])/(theta*dt);


	    /* estimate subscale velocities (accelerations) according
	     * to the prescribed values of the large scale quantities
	     * (constant velocity/pressure estimate or values from
	     * dirichlet boundary conditions) */

	    for(i=0;i<2;i++)
	    {
		actele->e.f2->sub_vel_acc_trial.a.da[i][lr*nis+ls]
		    =
		    sv_acc_mod_new[i];
		
		actele->e.f2->sub_vel_trial.a.da[i][lr*nis+ls]
		    =
		    sv_new[i];

	    		
	    }
	}
    }
    
}


/* clean up ------------------------------------------------------------*/

amdel(&xyze_a);
amdel(&xjm_a);
amdel(&funct_a);
amdel(&deriv_a);
amdel(&eveln_a);
amdel(&evelng_a);
amdel(&eaccng_a);
amdel(&derxy_a);
amdel(&vderxy_a);
amdel(&deriv2_a);
amdel(&w1_a);
amdel(&w2_a);
amdel(&derxy2_a);
amdel(&vderxy2_a);
amdel(&epreng_a);
amdel(&edeadng_a);



#ifdef DEBUG
dstrc_exit();
#endif

return;

} /* end of f2_estimate_new_trial_values_for_inc_gen_alpha */

#endif /* D_FLUID2_TDS */
#endif /* D_FLUID2     */
/*! @} (documentation module close)*/
