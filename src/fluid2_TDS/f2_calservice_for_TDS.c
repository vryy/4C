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
DOUBLE   dvnorm=ZERO;
DOUBLE   dpnorm=ZERO;
DOUBLE   dgnorm=ZERO;
DOUBLE    vnorm=ZERO;
DOUBLE    pnorm=ZERO;  /* values for norm calculation                  */
DOUBLE    gnorm=ZERO;
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
         if (dof>=numeq_total) continue;
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
	    
	    actnode->sol_increment.a.da[velnp][j]=
		actnode->sol_increment.a.da[veln][j]
		+
		dt*actnode->sol_increment.a.da[accn][j]
		+
		theta*dt* 
		(actnode->sol_increment.a.da[accnp][j]
		 -
		 actnode->sol_increment.a.da[accn][j]);

            /*
	      o Update velocity at time time n+alpha_F
	    */
	    
	    actnode->sol_increment.a.da[velnm][j]=
		actnode->sol_increment.a.da[veln ][j]
		+
		(alpha_F  )*
		(actnode->sol_increment.a.da[velnp][j]
		 -
		 actnode->sol_increment.a.da[veln][j]);

	    
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

#endif /* D_FLUID2_TDS */
#endif /* D_FLUID2     */
/*! @} (documentation module close)*/
