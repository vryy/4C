/*!----------------------------------------------------------------------
\file
\brief implicit and semi-implicit time integration algorithm for fsi

------------------------------------------------------------------------*/
/*! 
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "../fluid_full/fluid_prototypes.h"
#include "fsi_prototypes.h"

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h                                                  
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;                      
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;   
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | int                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern int            numcurve;
extern struct _CURVE *curve;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];
 
/*!---------------------------------------------------------------------                                         
\brief implicit and semi-implicit algorithms for 
       multifield fluid problems

<pre>                                                         genk 09/02

this functions solves the fluid within a multifield problem in an
ALE-framework. The mesh velocity is determined based on the displacements
of the mesh (fsi_ale()).
			     
</pre>   

\param *fsidyn   FSI_DYNAMIC                                 (i)
\param *fdyn	 FLUID_DYNAMIC                               (i)
\param *actfield FIELD            actual field               (i)
\param  mctrl    int              evaluation flag            (i)
\param  numff    int              number of fluid field      (i)
\return void 
\warning up to now only the One-Step-Theta scheme combined with a
fixed-point-like iteration scheme is tested! 

------------------------------------------------------------------------*/
void fsi_fluid(
                       FSI_DYNAMIC    *fsidyn,
		       FLUID_DYNAMIC  *fdyn, 
		       FIELD          *actfield, 
		       int             mctrl,
		       int             numff
	      )
{
static int             itnum;              /* counter for nonlinear iteration  */
int                    i;                  /* simply a counter                 */
static int             numeq;              /* number of equations on this proc */
static int             numeq_total;        /* total number of equations        */
int                    init;               /* flag for solver_control call     */
static int             nsysarray=1;        /* one system matrix                */
static int             actsysarray=0;      /* number of actual sysarray        */
static int             outstep;            /* counter for time integration     */
static int             pssstep;
static int             nfrastep;           /* number of steps for fractional-
                                              step-theta procedure             */
int                    calstress=1;        /* flag for stress calculation      */
int                    converged=0;        /* convergence flag                 */
int                    steady=0;           /* flag for steady state            */
static int             actpos;             /* actual position in sol. history  */
double                 vrat,prat;
static double          grat;               /* convergence ratios               */
double                 t1,ts,te;           /*				       */
static double          tes=ZERO;           /*				       */
static double          tss=ZERO;           /*				       */
static SOLVAR         *actsolv;            /* pointer to active sol. structure */
static PARTITION      *actpart;            /* pointer to active partition      */
static INTRA          *actintra;           /* pointer to active intra-communic.*/
static CALC_ACTION    *action;             /* pointer to the cal_action enum   */
static FLUID_DYN_CALC *dynvar;             /* pointer to fluid_dyn_calc        */

static ARRAY           ftimerhs_a;
static double         *ftimerhs;	   /* time - RHS		       */
static ARRAY           fiterhs_a;
static double         *fiterhs;	           /* iteration - RHS  		       */
static ARRAY           time_a;             /* stored time                      */
static ARRAY           totarea_a;
static double        **totarea;
static CONTAINER       container;          /* variables for calelm             */
static FLUID_STRESS    str;  

#ifdef DEBUG 
dstrc_enter("fsi_fluid");
#endif

switch (mctrl)
{
/*======================================================================*
 |                      I N I T I A L I S A T I O N                     |
 *======================================================================*/
case 1: 

fdyn->dt=fsidyn->dt;
fdyn->maxtime=fsidyn->maxtime;
fdyn->nstep=fsidyn->nstep;
grat=ZERO;
dsassert(fdyn->iop==4,"TIMEINTEGR for fluid: only ONE-STEP-THETA implemented!\n"); 
/*------------------------------------------ initialiase some counters */
outstep=0;
pssstep=0;  
/*--------------------------------------------------- set some pointers */
/*---------------------------- only valid for single field problem !!!! */
actsolv     = &(solv[numff]);
actpart     = &(partition[numff]);
action      = &(calc_action[numff]);
dynvar      = &(fdyn->dynvar);
container.fieldtyp = actfield->fieldtyp;
container.actndis  = 0;
container.turbu    = fdyn->turbu;
str         = str_fsicoupling;
dynvar->acttime=ZERO;
/*---------------- if we are not parallel, we have to allocate an alibi * 
  ---------------------------------------- intra-communicator structure */
#ifdef PARALLEL 
actintra    = &(par.intra[numff]);
#else
actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = fluid;
actintra->intra_rank     = 0;
actintra->intra_nprocs   = 1;
#endif

/*- there are only procs allowed in here, that belong to the fluid -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != fluid) break;

/*------------------------------- init the dist sparse matrices to zero */
solserv_zero_mat(
                 actintra,
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray])
                );

/*---------------------------- get global and local number of equations */
solserv_getmatdims(actsolv->sysarray[actsysarray],
                   actsolv->sysarray_typ[actsysarray],
                   &numeq,
                   &numeq_total);

/*--------------------------------------- allocate 1 dist. vector 'rhs' */
actsolv->nrhs = 1;
solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
solserv_zero_vec(&(actsolv->rhs[0]));

/*------------------------------------------- allocate 1 dist. solution */		       
actsolv->nsol= 1;
solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
solserv_zero_vec(&(actsolv->sol[0]));
                                     
/*--------------- allocate one redundant vector ftimerhs of full lenght */
/*        this is used by the element routines to assemble the  Time RHS*/
ftimerhs = amdef("ftimerhs",&ftimerhs_a,numeq_total,1,"DV");

/*---------------  allocate one redundant vector fiterhs of full lenght */
/*   this is used by the element routines to assemble the  Iteration RHS*/
fiterhs = amdef("fiterhs",&fiterhs_a,numeq_total,1,"DV");

/*--------------------------- allocate one vector for storing the time */
if (par.myrank==0) amdef("time",&time_a,1000,1,"DV");
/*--------------------------- allocate one vector for storing the area */
if (fdyn->checkarea>0)
{ 
   totarea = amdef("area",&totarea_a,fsidyn->nstep,fdyn->itemax,"DA");
   amzero(&totarea_a);
}
/*--------------------------------------------- initialise fluid field */
fluid_init(actfield,fdyn,7,str);		     
actpos=0;   

/*-------------------------------------- init the dirichlet-conditions */
fluid_initdirich(actfield, fdyn);

/*---------------------------------- initialize solver on all matrices */
/*
NOTE: solver init phase has to be called with each matrix one wants to 
      solve with. Solver init phase has to be called with all matrices
      one wants to do matrix-vector products and matrix scalar products.
      This is not needed by all solver libraries, but the solver-init 
      phase is cheap in computation (can be costly in memory)
*/
/*--------------------------------------------------- initialize solver */
init=1;
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(actsolv->sol[0]),
               &(actsolv->rhs[0]),
               init);
	       
/*------------------------------------- init the assembly for stiffness */
init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);
	       	       
/*------------------------------- init the element calculating routines */
*action = calc_fluid_init;
calinit(actfield,actpart,action,&container);

/*-------------------------------------------- initialise energey check */
if (fsidyn->ichecke>0) fsi_dyneint(NULL,NULL,NULL,fdyn,2);

/*------------------------------------------------ output to the screen */
if (par.myrank==0) printf("\n\n");
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
printf("PROC  %3d | FIELD FLUID     | number of equations      : %10d \n", 
        par.myrank,numeq);
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
if (par.myrank==0)
printf("          | FIELD FLUID     | total number of equations: %10d \n",numeq_total);
if (par.myrank==0) printf("\n\n");

/*-------------------------------- calculate curvature at the beginning */
if (fdyn->surftens!=0)
{
   fluid_tcons(fdyn,dynvar);
   *action = calc_fluid_curvature;
   fluid_curvature(actfield,actpart,actintra,action);
}

/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
monitoring(actfield,numff,actpos,0,fdyn->time);
/*-------------------------------------- print out initial data to .out */
out_sol(actfield,actpart,actintra,fdyn->step,actpos);
actpos++;
break;

/*======================================================================*
 |                     S O L U T I O N    P H A S E                     |
 *======================================================================*/
/* nodal solution history fluid field:                                  *
 * sol[0][j]           ... initial data 				*
 * sol[1...actpos][j]  ... solution for visualisation (real pressure)	*
 * sol_increment[0][j] ... solution at time (n-1)			*
 * sol_increment[1][j] ... solution at time (n) 			*
 * sol_increment[2][j] ... solution at time (n+g)			*
 * sol_increment[3][j] ... solution at time (n+1)			*
 * sol_increment[4][i] ... grid velocity time (n) -> (n+1)		*
 * sol_increment[5][i] ... convective velocity at time (n)		*
 * sol_increment[6][i] ... convective velocity at time (n+1)	        *
 * sol_mf[0][j]        ... solution at time (n+1)			*
 * sol_mf[1][j]        ... nodal stresses at FS-interface at time (n+1) *
 *======================================================================*/
case 2:

/*------- copy solution from sol_increment[1][j] to sol_increment[3][j] */
if (fsidyn->ifsi>=4) solserv_sol_copy(actfield,0,1,1,1,3);
/*- there are only procs allowed in here, that belong to the fluid -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != fluid) break;
/*------------------------------------------ check (starting) algorithm */
if (fdyn->step<=(fdyn->nums+1))
   fluid_startproc(fdyn,&nfrastep);

/*------------------------------ calculate constants for time algorithm */
fluid_tcons(fdyn,dynvar);
dynvar->acttime=fdyn->time;

/*------------------------------------------------ output to the screen */
if (par.myrank==0) printf("Solving FLUID by One-Step-Theta ...\n"); 
/*--------------------------------------------------------- ALE-PHASE I */
if (fsidyn->iale!=0)
{
   /*--------------------------------------------- get the gridvelocity */
   if (fsidyn->iale>0) fsi_alecp(actfield,dynvar,fdyn->numdf,1);
   else  dserror("ALE field by function not implemented yet!\n");
  /*----------------------------------------------- change element flag */
  dynvar->ishape=1;
  /*------------------- calculate ALE-convective velocities at time (n) */
  fsi_aleconv(actfield,fdyn->numdf,5,1);     
}
 
/*--------------------- set dirichlet boundary conditions for  timestep */
fluid_setdirich(actfield,fdyn);

/*-------------------------------------------------- initialise timerhs */
amzero(&ftimerhs_a);

if (fdyn->itchk!=0 && par.myrank==0)
{
   printf("---------------------------------------------------------------- \n");
   printf("|- step/max -|-  tol     [norm] -|- vel. error -|- pre. error -|\n");
}
itnum=1;
/*======================================================================* 
 |               N O N L I N E A R   I T E R A T I O N                  |
 *======================================================================*/
nonlniter:

/*------------------------- calculate constants for nonlinear iteration */
fluid_icons(fdyn,dynvar,itnum);

/*-------------------------------------------------------- ALE-PHASE II */
if (fsidyn->iale!=0)
{
   /*---- for implicit free surface we have to update the grid velocity */
   if (fdyn->freesurf==2)
   {
      if (fsidyn->iale>0) fsi_alecp(actfield,dynvar,fdyn->numdf,2);   
      else  dserror("ALE field by function not implemented yet!\n");
     /*-------------------------------------------- change element flag */
      dynvar->ishape=1;   
   }
   /*--------------- calculate ale-convective velocities at  time (n+1) */
   fsi_aleconv(actfield,fdyn->numdf,6,3);
}

/*--------------------------------- calculate curvature at free surface */
if (dynvar->surftens!=0) 
{
   *action = calc_fluid_curvature;
   fluid_curvature(actfield,actpart,actintra,action);
}
/*---------------------------- intitialise global matrix and global rhs */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_mat(actintra,&(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray]));

/*------------------------------------------- initialise iterations-rhs */
amzero(&fiterhs_a);

/*-------------- form incremental matrices, residual and element forces */
*action = calc_fluid;
t1=ds_cputime();
container.dvec         = NULL;
container.ftimerhs     = ftimerhs;
container.fiterhs      = fiterhs;
container.global_numeq = numeq_total;
container.nii          = dynvar->nii;
container.nif          = dynvar->nif;
container.kstep        = 0;
container.fieldtyp     = actfield->fieldtyp;
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
       &container,action);
te=ds_cputime()-t1;
tes+=te;	     

/*--------------------------------------------------------------------- *
 | build the actual rhs-vector:                                         |
 |        rhs = ftimerhs + fiterhs                                      |
 *----------------------------------------------------------------------*/
/* add time-rhs: */
assemble_vec(actintra,
             &(actsolv->sysarray_typ[actsysarray]),
             &(actsolv->sysarray[actsysarray]),
             &(actsolv->rhs[0]),
             ftimerhs,
             1.0
             );
/* add iteration-rhs: */
assemble_vec(actintra,
             &(actsolv->sysarray_typ[actsysarray]),
             &(actsolv->sysarray[actsysarray]),
             &(actsolv->rhs[0]),
             fiterhs,
             1.0
             );

/*-------------------------------------------------------- solve system */
init=0;
t1=ds_cputime();
solver_control(actsolv, actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(actsolv->sol[0]),
               &(actsolv->rhs[0]),
               init);
ts=ds_cputime()-t1;
tss+=ts;

/*-- set flags for stability parameter evaluation and convergence check */
dynvar->ishape=0;

/*--- return solution to the nodes and calculate the convergence ratios */
fluid_result_incre(actfield,actintra,&(actsolv->sol[0]),3,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray]),
		     &vrat,&prat,&grat,fdyn);	     

/*---------------------------------------------------- store total area */
if (fdyn->checkarea>0)
{
   if (totarea_a.fdim>=fsidyn->step && totarea_a.sdim>=itnum)
   totarea[fsidyn->step-1][itnum-1] = dynvar->totarea;
}

/*----------------------------------------- iteration convergence check */
converged = fluid_convcheck(fdyn,vrat,prat,grat,itnum,te,ts);

/*--------------------- check if nonlinear iteration has to be finished */
if (converged==0)
{
   itnum++;
   goto nonlniter;
}
/*----------------------------------------------------------------------*
 | -->  end of nonlinear iteration                                      |
 *----------------------------------------------------------------------*/
 
/*-------------------------------------------------- steady state check */
 /* no steady state check for fsi-problems!!!                           */


/*------------------------- calculate stresses transferred to structure */
if (fsidyn->ifsi>0)
{
   *action = calc_fluid_stress;
   container.nii= 0;
   container.nif= 0;
   container.str=str;
   calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
          &container,action);

   /* since stresses are stored locally at the element it's necassary to 
      reduce them to all procs! ----------------------------------------*/
   dsassert(actsolv->parttyp==cut_elements,"Stress reduction for 'cut_nodes' not possible\n");
   fluid_reducestress(actintra,actfield,fdyn->numdf,str);
   /*----------------------------------------- store stresses in sol_mf */
#if 0
   solserv_zerosol(actfield,0,3,1);
#endif
   solserv_sol_zero(actfield,0,3,1);
   fsi_fluidstress_result(actfield,fdyn->numdf);
}

if (fsidyn->ifsi>=4)
break;

/*======================================================================* 
 |                       F I N A L I S I N G                            |
 *======================================================================*/
case 3:

/*---------------------- for multifield fluid problems with freesurface */
/*-------------- copy solution from sol_increment[3][j] to sol_mf[0][j] */
/* check this for FSI with free surface!!! */
if (fdyn->freesurf>0) solserv_sol_copy(actfield,0,1,3,3,0);

/*-------- copy solution from sol_increment[3][j] to sol_increment[1[j] */
solserv_sol_copy(actfield,0,1,1,3,1);

/*---------------------------------------------- finalise this timestep */
outstep++;
pssstep++;

/*-------- copy solution from sol_increment[3][j] to sol_[actpos][j]   
           and transform kinematic to real pressure --------------------*/
fluid_sol_copy(actfield,0,1,0,3,actpos,fdyn->numdf);


if (outstep==fdyn->upout && ioflags.fluid_sol_file==1)
{
   outstep=0;
   out_sol(actfield,actpart,actintra,fdyn->step,actpos);
}
/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
monitoring(actfield,numff,actpos,fdyn->step,fdyn->time);

if (pssstep==fsidyn->uppss && ioflags.fluid_vis_file==1 && par.myrank==0)
{
   pssstep=0;
   /*--------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim)
      amredef(&(time_a),time_a.fdim+100,1,"DV");
   time_a.a.dv[actpos] = fdyn->time;   
   actpos++;
   fsidyn->actpos = actpos;
} 
/*--------------------------------------------------------------------- */
break;   

/*======================================================================* 
 |                C L E A N I N G   U P   P H A S E                     |
 *======================================================================*/ 
case 99:
/*- there are only procs allowed in here, that belong to the fluid -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != fluid) break;
if (pssstep==0) actpos--;
/*---------------------------------------------- print out to .mon file */
if (ioflags.monitor==1 && par.myrank==0)
out_monitor(actfield,numff);

/*------------------------------------- print out solution to .out file */
if (outstep!=0 && ioflags.fluid_sol_file==1)
out_sol(actfield,actpart,actintra,fdyn->step,actpos);

/*------------------------------------ print out solution to 0.pss file */
if (ioflags.fluid_vis_file==1 && par.myrank==0)
{
   if (pssstep!=0)
   {
      /*------------------------------------------ store time in time_a */
      if (actpos >= time_a.fdim)
      amredef(&(time_a),time_a.fdim+1000,1,"DV");
      time_a.a.dv[actpos] = fdyn->time;   
   }   
   visual_writepss(actfield,actpos+1,&time_a);
}

/*-------------------------------------- output of area to monitor file */
if (fdyn->checkarea>0) out_area(totarea_a);

/*---------------------------------- print total CPU-time to the screen */
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
for (i=0;i<par.nprocs;i++)
{
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
if (par.myrank==i)
{
printf("\n");
printf("PROC  %3d | FIELD FLUID     | total time element for calculations: %10.3E \n", 
        par.myrank,tes);
printf("PROC  %3d | FIELD FLUID     | total time for solver              : %10.3E \n", 
        par.myrank,tss);
}
}
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif

/*------------------------------------------------------------- tidy up */   
amdel(&ftimerhs_a);
amdel(&fiterhs_a);
if (par.myrank==0) amdel(&time_a);
if (fdyn->checkarea>0) amdel(&totarea_a);
solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);

#ifndef PARALLEL 
CCAFREE(actintra);
#endif

break;
default:
   dserror("Parameter out of range: mctrl \n");
} /* end switch (mctrl) */

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_isi */ 

#endif

/*! @} (documentation module close)*/
