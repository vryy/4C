/*!----------------------------------------------------------------------
\file
\brief implicit and semi-implicit time integration algorithm for fluid

------------------------------------------------------------------------*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "../headers/solution_mlpcg.h"
#include "../headers/solution.h"
#include "fluid_prototypes.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
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
 
/*----------------------------------------------------------------------*
 | routine to control implicit and semi-implicit algorithms for fluid   |
 | problems combined with Newton and fixed point iteration schemes.     | 
 | IN PROGRESS: ONE-STEP-THETA                                          |
 |              fixed point like iteration                              |
 |              only Euler no ALE                                       |
 |                                                          genk  03/02 |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------                                         
\brief implicit and semi-implicit algorithms for fluid problems

<pre>                                                         genk 03/02

this routine conrols the implicit and semi-implicit algorithms for fluid
problems combined with different nonlinear iteration schemes.

time-discretisation:
fdyn->iop=2: Semi-Implicit-One-Step Method
fdyn->iop=3: Semi-Implicit-Two-Step Method
fdyn->iop=4: One-Step-Theta Scheme
fdyn->iop=5: Fractional-Step-Theta Scheme 

see dissertation of W.A. Wall chapter 4.2 'Zeitdiskretisierung'

nonlinear iteration scheme:
fdyn->ite=0: no nonlinear iteration
fdyn->ite=1: fixed-point-like iteration
fdyn->ite=2: Newton iteration
fdyn->ite=3: fixed-point iteration

see dissertation chapter 4.3 'Linearisierung und Iteratonsverfahren'.
			     
</pre>   
\param *fdyn	 FLUID_DYNAMIC (i)    

\return void 
\warning up to now only the One-Step-Theta scheme combined with a
fixed-point-like iteration scheme is tested! 

------------------------------------------------------------------------*/
void fluid_isi(FLUID_DYNAMIC *fdyn)
{
int             itnum;              /* counter for nonlinear iteration  */
int             i;                  /* simply a counter                 */
int             numeq;              /* number of equations on this proc */
int             numeq_total;        /* total number of equations        */
int             init;               /* flag for solver_control call     */
int             nsysarray=1;        /* one system matrix                */
int             actsysarray=0;      /* number of actual sysarray        */
int             istep=0;            /* counter for time integration     */
int             iststep=0;          /* counter for time integration     */
int             nfrastep;           /* number of steps for fractional-
                                       step-theta procedure             */
int             actcurve;           /* actual timecurve                 */
int             converged=0;        /* convergence flag                 */
int             steady=0;           /* flag for steady state            */
int             actpos;             /* actual position in sol. history  */
double          vrat,prat;          /* convergence ratios               */
double          t1,ts,te;	    /*					*/
double          tfs=0.0;            /*					*/
double          tes=0.0;            /*					*/
double          tas=0.0;            /*					*/
double          tss=0.0;            /*					*/
double          tds=0.0;            /*					*/
double          tcs=0.0;            /*					*/
double          tsts=0.0;           /* variables for time tracing	*/

SOLVAR         *actsolv;            /* pointer to active sol. structure */
PARTITION      *actpart;            /* pointer to active partition      */
FIELD          *actfield;           /* pointer to active field          */
INTRA          *actintra;           /* pointer to active intra-communic.*/
CALC_ACTION    *action;             /* pointer to the cal_action enum   */
FLUID_DYN_CALC *dynvar;             /* pointer to fluid_dyn_calc        */

ARRAY           ftimerhs_a;
double         *ftimerhs;	    /* time - RHS			*/
ARRAY           fiterhs_a;
double         *fiterhs;	    /* iteration - RHS  		*/
ARRAY           time_a;             /* stored time                      */

#ifdef DEBUG 
dstrc_enter("fluid_isi");
#endif

/*--------------------------------------------------- set some pointers */
/*---------------------------- only valid for single field problem !!!! */
actfield    = &(field[0]);
actsolv     = &(solv[0]);
actpart     = &(partition[0]);
action      = &(calc_action[0]);
dynvar      = &(fdyn->dynvar);

/*---------------- if we are not parallel, we have to allocate an alibi * 
  ---------------------------------------- intra-communicator structure */
#ifdef PARALLEL 
actintra    = &(par.intra[0]);
#else
actintra    = (INTRA*)CALLOC(1,sizeof(INTRA));
if (!actintra) dserror("Allocation of INTRA failed");
actintra->intra_fieldtyp = fluid;
actintra->intra_rank     = 0;
actintra->intra_nprocs   = 1;
#endif

/*- there are only procs allowed in here, that belong to the fluid -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != fluid) goto end;

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

/*------------------------------------------------ output to the screen */
#ifdef PARALLEL
printf("number of eqations on PROC %3d : %10d \n", 
        par.myrank,numeq);
if (par.myrank==0)
   printf("total number of equations: %10d \n",numeq_total);	
#else
printf("total number of equations: %10d \n",numeq_total);
#endif

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
amdef("time",&time_a,1000,1,"DV");

/*--------------------------------------------- initialise fluid field */
fluid_init(actfield, fdyn);		     
actpos=0;

/*--------------------------------------- init all applied time curves */
for (actcurve=0; actcurve<numcurve; actcurve++)
   dyn_init_curve(actcurve,fdyn->nstep,fdyn->dt,fdyn->maxtime);   

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
init_assembly(actpart,actsolv,actintra,actfield,actsysarray);
	       	       
/*------------------------------- init the element calculating routines */
*action = calc_fluid_init;
calinit_fluid(actfield,actpart,action);

/*-------------------------------------- print out initial data to .out */
out_sol(actfield,actpart,actintra,fdyn->step,actpos);

/*----------------------------------------------------------------------* 
 | do the time loop                                                     |
 *----------------------------------------------------------------------*/
timeloop:
fdyn->step++;
iststep++;

/*------------------------------------------ check (starting) algorithm */
if (fdyn->step<=(fdyn->nums+1))
   fluid_startproc(fdyn,&nfrastep);

/*------------------------------ calculate constants for time algorithm */
fluid_tcons(fdyn,dynvar);
fdyn->time += dynvar->dta; 

/*------------------------------------------------ output to the screen */
if (par.myrank==0) fluid_algoout(fdyn,dynvar);

/*--------------------- set dirichlet boundary conditions for  timestep */
t1=ds_cputime();
fluid_setdirich(actfield,fdyn);
tds+=ds_cputime()-t1;

/*-------------------------------------------------- initialise timerhs */
amzero(&ftimerhs_a);

if (fdyn->itchk!=0 && par.myrank==0)
{
   printf(" _____________________________________________________________ \n");
   printf("|- step/max -|-  tol     [norm] -|- vel. error -|- pre. error-| \n");
}
itnum=1;
/*----------------------------------------------------------------------*
 | do the nonlinear iteration                                           |
 *----------------------------------------------------------------------*/
nonlniter:

/*------------------------- calculate constants for nonlinear iteration */
fluid_icons(fdyn,dynvar,itnum);

/*---------------------------- intitialise global matrix and global rhs */
solserv_zero_vec(&(actsolv->rhs[0]));
solserv_zero_mat(actintra,&(actsolv->sysarray[actsysarray]),
                 &(actsolv->sysarray_typ[actsysarray]));

/*------------------------------------------- initialise iterations-rhs */
amzero(&fiterhs_a);

/*-------------- form incremental matrices, residual and element forces */
*action = calc_fluid;
t1=ds_cputime();
calelm_fluid(actfield,actsolv,actpart,actintra,actsysarray,-1,
             ftimerhs,fiterhs,numeq_total,dynvar->nii,
	     dynvar->nif,0,action);
te=ds_cputime()-t1;
tes+=te;	     

/*--------------------------------------------------------------------- *
 | build the actual rhs-vector:                                         |
 |        rhs = ftimerhs + fiterhs                                      |
 *----------------------------------------------------------------------*/
/* add time-rhs: */
t1=ds_cputime();
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
tas+=ds_cputime()-t1;

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
t1=ds_cputime();
fluid_result_incre(actfield,actintra,&(actsolv->sol[0]),3,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray]),
		     &vrat,&prat,fdyn);
tcs+=ds_cputime()-t1;		     

/*----------------------------------------- iteration convergence check */
converged = fluid_convcheck(fdyn,vrat,prat,itnum,te,ts);

/*--------------------- check if nonlinear iteration has to be finished */
if (converged==0)
{
   itnum++;
   goto nonlniter;
}

/*----------------------------------------------------------------------*
 | -->  end of nonlinear iteration                                      |
 *----------------------------------------------------------------------*/
if (fdyn->itchk!=0 && par.myrank==0)
{
   printf("|____________|___________________|______________|_____________| \n");
   printf("\n"); 
}   

/*-------------------------------------------------- steady state check */
t1=ds_cputime();
if (fdyn->stchk==iststep)
{
   iststep=0;
   steady = fluid_steadycheck(fdyn,actfield,numeq_total);
}
tsts+=ds_cputime()-t1;

/*--------------------------------- copy solution at (n+1) to place (n) *
                                      in solution history sol_increment */
fluid_copysol(fdyn,actfield,3,1,0);
/*---------------------------------------------- finalise this timestep */
istep++;
t1=ds_cputime();

/*------------------------------------------ store results in node->sol *
 * on position actpos; actpos=0: initial data;                          *
 *                     actpos=1... stored timesteps --------------------*/
/*#######################################################################
  at the moment kinematic pressure is stored and written to output file
  this has to be changed!!!!!
  maybe the whole element formulation will be
  transformed to real pressure!!!
  #####################################################################*/

if (istep==fdyn->idisp)
{
   istep=0;
   actpos++;
   fluid_copysol(fdyn,actfield,3,actpos,1);
   /*--------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim)
   {
      amredef(&(time_a),time_a.fdim+1000,1,"DV");
   }
   time_a.a.dv[actpos] = fdyn->time;
   /*------------------------------------------------ print out to .out */
   out_sol(actfield,actpart,actintra,fdyn->step,actpos);
} 
tfs+=ds_cputime()-t1;  

/*--------------------- check time and number of steps and steady state */
if (fdyn->step < fdyn->nstep && fdyn->time <= fdyn->maxtime && steady==0)
   goto timeloop;

/*------------------------------------ print out solution to 0.pss file */
if (ioflags.fluid_vis_file==1 && par.myrank==0) 
   fluid_writevispss(actfield,actpos+1,&time_a);

/*---------------------------------- print total CPU-time to the screen */
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
for (i=0;i<par.nprocs;i++)
{
if (par.myrank==i)
{
printf("\n");
printf("PROC %3d: TOTAL TIME for setting dirichlet values: %10.3#E \n", par.myrank,tds);
printf("PROC %3d: TOTAL TIME for element calculations: %10.3#E \n", par.myrank,tes);
printf("PROC %3d: TOTAL TIME for assembly of RHS: %10.3#E \n", par.myrank,tas);
printf("PROC %3d: TOTAL TIME for SOLVER: %10.3#E \n", par.myrank,tss);
printf("PROC %3d: TOTAL TIME for calculation of convergence ratios: %10.3#E \n", par.myrank,tcs);
printf("PROC %3d: TOTAL TIME for steady state check: %10.3#E \n", par.myrank,tsts);
printf("PROC %3d: TOTAL TIME finalising of time step: %10.3#E \n", par.myrank,tfs);
}
}
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
end:

/*--------------------------------------------------- cleaning up phase */
amdel(&ftimerhs_a);
amdel(&fiterhs_a);
solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
solserv_del_vec(&(actsolv->sol),actsolv->nsol);

/*----------------------------------------------------------------------*/
#ifndef PARALLEL 
FREE(actintra);
#endif

#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_isi */ 

#endif
