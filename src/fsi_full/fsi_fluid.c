/*!----------------------------------------------------------------------
\file
\brief implicit and semi-implicit time integration algorithm for fsi

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

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
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
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
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT            numcurve;
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
\param *actfield FIELD            actual field               (i)
\param  mctrl    INT              evaluation flag            (i)
\param  numff    INT              number of fluid field      (i)
\return void 
\warning up to now only the One-Step-Theta scheme combined with a
fixed-point-like iteration scheme is tested! 

------------------------------------------------------------------------*/
void fsi_fluid(
		       FIELD          *actfield, 
		       INT             mctrl,
		       INT             numff
	      )
{
static INT             itnum;              /* counter for nonlinear iteration  */
INT                    i,j;		   /* counters				*/
static INT             numeq;              /* number of equations on this proc */
static INT             numeq_total;        /* total number of equations        */
INT                    init;               /* flag for solver_control call     */
static INT             nsysarray=1;        /* one system matrix                */
static INT             actsysarray=0;      /* number of actual sysarray        */
static INT             outstep;            /* counter for output control       */
static INT             pssstep;            /* counter for output control       */
static INT             restartstep;        /* counter for restart control      */
static INT             nfrastep;           /* number of steps for fractional-
                                              step-theta procedure             */
static INT             restart;
INT                    calstress=1;        /* flag for stress calculation      */
INT                    converged=0;        /* convergence flag                 */
INT                    steady=0;           /* flag for steady state            */

DOUBLE			liftdrag[3];	/* array with lift & drag coeff.*/
DOUBLE			recv[3];	

DLINE		      *actdline;

static INT             actpos;             /* actual position in sol. history  */
DOUBLE                 vrat,prat;
static DOUBLE          grat;               /* convergence ratios               */
DOUBLE                 t1,ts,te;           /*				       */
static DOUBLE          tes=ZERO;           /*				       */
static DOUBLE          tss=ZERO;           /*				       */
static SOLVAR         *actsolv;            /* pointer to active sol. structure */
static PARTITION      *actpart;            /* pointer to active partition      */
static INTRA          *actintra;           /* pointer to active intra-communic.*/
static CALC_ACTION    *action;             /* pointer to the cal_action enum   */

static ARRAY           ftimerhs_a;
static DOUBLE         *ftimerhs;	   /* time - RHS		       */
static ARRAY           fiterhs_a;
static DOUBLE         *fiterhs;	           /* iteration - RHS  		       */
static ARRAY           time_a;             /* stored time                      */
static ARRAY           totarea_a;
static DOUBLE        **totarea;
static CONTAINER       container;          /* variables for calelm             */
static FLUID_STRESS    str;  
static FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
static FSI_DYNAMIC    *fsidyn;             /* fsi dynamic variables     */

#ifdef DEBUG 
dstrc_enter("fsi_fluid");
#endif

switch (mctrl)
{
/*======================================================================*
 |                      I N I T I A L I S A T I O N                     |
 *======================================================================*/
case 1: 
fdyn   = alldyn[genprob.numff].fdyn;
fsidyn = alldyn[genprob.numaf+1].fsidyn;

fdyn->dt=fsidyn->dt;
fdyn->maxtime=fsidyn->maxtime;
fdyn->nstep=fsidyn->nstep;
grat=ZERO;
dsassert(fdyn->iop==4,"TIMEINTEGR for fluid: only ONE-STEP-THETA implemented!\n"); 
/*------------------------------------------ initialiase some counters */
outstep=0;
pssstep=0;  
restartstep=0;
/*--------------------------------------------------- set some pointers */
/*---------------------------- only valid for single field problem !!!! */
actsolv     = &(solv[numff]);
actpart     = &(partition[numff]);
action      = &(calc_action[numff]);
restart     = genprob.restart;
container.fieldtyp = actfield->fieldtyp;
container.actndis  = 0;
container.turbu    = fdyn->turbu;
str         = str_fsicoupling;
fdyn->acttime=ZERO;

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
amdef("time",&time_a,1000,1,"DV");
/*--------------------------- allocate one vector for storing the area */
if (fdyn->checkarea>0)
{ 
   totarea = amdef("area",&totarea_a,fsidyn->nstep,fdyn->itemax,"DA");
   amzero(&totarea_a);
}
/*--------------------------------------------- initialise fluid field */
if (restart>0)
{
   if (fdyn->init>0)
      dserror("Initial field either by restart, or by function or from file ...\n");
   else
   {
      fdyn->resstep=genprob.restart;
      fdyn->init=2;
   }
}
fluid_init(actpart,actintra,actfield,action,&container,7,str);		     
actpos=0;   

/*-------------------------------------- init the dirichlet-conditions */
fluid_initdirich(actfield);

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

/*---------------------------------- allocate fluid integration data ---*/
alldyn[genprob.numff].fdyn->data = (FLUID_DATA*)CCACALLOC(1,sizeof(FLUID_DATA));

/*------------------------------- init the element calculating routines */
*action = calc_fluid_init;
calinit(actfield,actpart,action,&container);

/*-------------------------------------------- initialise energey check */
if (fsidyn->ichecke>0) fsi_dyneint(NULL,2);

/*------------------------------------------------ output to the screen */
if (par.myrank==0) printf("\n\n");
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
for (i=0;i<par.nprocs;i++)
if (par.myrank==i)
printf("PROC  %3d | FIELD FLUID     | number of equations      : %10d \n", 
        par.myrank,numeq);
#ifdef PARALLEL
MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
if (par.myrank==0)
printf("          | FIELD FLUID     | total number of equations: %10d \n",numeq_total);
if (par.myrank==0) printf("\n\n");

/*----------------------------------------- calculate initial curvature */
if (fdyn->surftens!=0 && restart==0)
{
   fluid_tcons();
   *action = calc_fluid_curvature;
   fluid_curvature(actfield,actpart,actintra,action);
}

/*------------------------- init lift&drag calculation real FSI-problem */
if (fdyn->liftdrag==3)

/*----------------------------------------init lift&drag calculation ---*/
fluid_liftdrag(0,NULL,NULL,NULL,NULL,NULL,NULL);

/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
monitoring(actfield,numff,actpos,0,fdyn->acttime);
/*-------------------------------------- print out initial data to .out */
out_sol(actfield,actpart,actintra,fdyn->step,actpos);
/*---------- calculate time independent constants for time algorithm ---*/
fluid_cons();
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
   fluid_startproc(&nfrastep,0);

/*------------------------------ calculate constants for time algorithm */
fluid_tcons();

/*------------------------------------------------ output to the screen */
if (par.myrank==0) printf("Solving FLUID by One-Step-Theta ...\n"); 
/*--------------------------------------------------------- ALE-PHASE I */
if (fsidyn->iale!=0)
{
   /*--------------------------------------------- get the gridvelocity */
   if (fsidyn->iale>0) fsi_alecp(actfield,fdyn->numdf,1);
   else  dserror("ALE field by function not implemented yet!\n");
  /*----------------------------------------------- change element flag */
  fdyn->ishape=1;
  /*------------------- calculate ALE-convective velocities at time (n) */
  fsi_aleconv(actfield,fdyn->numdf,5,1);     
}
 
/*--------------------- set dirichlet boundary conditions for  timestep */
fluid_setdirich(actfield,3);

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
fluid_icons(itnum);

/*-------------------------------------------------------- ALE-PHASE II */
if (fsidyn->iale!=0)
{
   /*---- for implicit free surface we have to update the grid velocity */
   if (fdyn->freesurf==2)
   {
      if (fsidyn->iale>0) fsi_alecp(actfield,fdyn->numdf,2);   
      else  dserror("ALE field by function not implemented yet!\n");
     /*-------------------------------------------- change element flag */
      fdyn->ishape=1;   
   }
   /*--------------- calculate ale-convective velocities at  time (n+1) */
   fsi_aleconv(actfield,fdyn->numdf,6,3);
}

/*--------------------------------- calculate curvature at free surface */
if (fdyn->surftens!=0) 
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
container.nii          = fdyn->nii;
container.nif          = fdyn->nif;
container.kstep        = 0;
container.fieldtyp     = actfield->fieldtyp;
container.is_relax     = 0;
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
fdyn->ishape=0;

/*--- return solution to the nodes and calculate the convergence ratios */
fluid_result_incre(actfield,actintra,&(actsolv->sol[0]),3,
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sysarray_typ[actsysarray]),
		     &vrat,&prat,&grat);	     

/*---------------------------------------------------- store total area */
if (fdyn->checkarea>0)
{
   if (totarea_a.fdim>=fsidyn->step && totarea_a.sdim>=itnum)
   totarea[fsidyn->step-1][itnum-1] = fdyn->totarea;
}

/*----------------------------------------- iteration convergence check */
converged = fluid_convcheck(vrat,prat,grat,itnum,te,ts);

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
   container.is_relax = 0;
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

/*----------------------------------------------- lift&drag computation */
if (fdyn->liftdrag>0)
{
   *action = calc_fluid_liftdrag;
   container.str=str_liftdrag;
   fluid_liftdrag(fdyn->liftdrag,action,&container,actfield,
                  actsolv,actpart,actintra);
}

/*---------------------- for multifield fluid problems with freesurface */
/*-------------- copy solution from sol_increment[3][j] to sol_mf[0][j] */
/* check this for FSI with free surface!!! */
if (fdyn->freesurf>0) solserv_sol_copy(actfield,0,1,3,3,0);

/*-------- copy solution from sol_increment[3][j] to sol_increment[1[j] */
solserv_sol_copy(actfield,0,1,1,3,1);

/*---------------------------------------------- finalise this timestep */
outstep++;
pssstep++;
restartstep++;

if (pssstep==fsidyn->uppss && ioflags.fluid_vis_file==1)
{
   pssstep=0;
   /*--------------------------------------------- store time in time_a */
   if (actpos >= time_a.fdim)
      amredef(&(time_a),time_a.fdim+100,1,"DV");
   time_a.a.dv[actpos] = fdyn->acttime;   
   actpos++;
} 

/*-------- copy solution from sol_increment[3][j] to sol_[actpos][j]   
           and transform kinematic to real pressure --------------------*/
solserv_sol_copy(actfield,0,1,0,3,actpos);
fluid_transpres(actfield,0,0,actpos,fdyn->numdf-1,0);

if (outstep==fdyn->upout && ioflags.fluid_sol_file==1)
{
   outstep=0;
   out_sol(actfield,actpart,actintra,fdyn->step,actpos);
}
/*---------------------------------------------------------- monitoring */
if (ioflags.monitor==1)
monitoring(actfield,numff,actpos,fdyn->step,fdyn->acttime);

fsidyn->actpos = actpos;

/*------------------------------------------- write restart to pss file */
if (restartstep==fsidyn->res_write_evry)
{
   restartstep=0;
   restart_write_fluiddyn(fdyn,actfield,actpart,actintra,action,&container);   
}

/*--------------------------------------------------------------------- */
break;   

/*======================================================================*
 |     S O L U T I O N    F O R    S T E E P E S T    D E S C E N T     |
 |                        E V A L U A T I O N                           |
 *======================================================================*/
/* nodal solution history fluid field:                                  *
 * sol[0][j]           ... initial data 				*
 * sol[1...actpos][j]  ... solution for visualisation (real pressure)	*
 * sol_increment[0][j] ... solution at time (n-1)			*
 * sol_increment[1][j] ... solution at time (n) 			*
 * sol_increment[2][j] ... solution at time (n+g)			*
 * sol_increment[3][j] ... solution at time (n+1)			*
 * sol_increment[4][i] ... grid velocity time (n) -> (n+1) #		*
 * sol_increment[5][i] ... convective velocity at time (n)		*
 * sol_increment[6][i] ... convective velocity at time (n+1) #	        *
 * sol_increment[7][i] ... fluid solution for relax.-param. steep. desc.*
 * #: these vectors also used for steepest descent calculation		*
 * sol_mf[0][j]        ... solution at time (n+1)			*
 * sol_mf[1][j]        ... nodal stresses at FS-interface at time (n+1) *
 *======================================================================*/
case 6:

if (fsidyn->ifsi != 6) 
dserror("No auxiliary fluid solution within this coupling scheme");

/*- there are only procs allowed in here, that belong to the fluid -----*/
/* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
if (actintra->intra_fieldtyp != fluid) break;

/*------------------------------ calculate constants for time algorithm */
fluid_tcons();          

/*------------------------------------------------ output to the screen */
if (par.myrank==0) printf("          - Solving FLUID ...\n"); 
/*--------------------------------------------------------- ALE-PHASE I */
if (fsidyn->iale!=0)
{
  /*----------------------------------------------- change element flag */
  fdyn->ishape=1;
  /*------------------- calculate ALE-convective velocities at time (n) */
  fsi_aleconv(actfield,fdyn->numdf,6,3);     
}
 
/*----------------------------------- set dirichlet boundary conditions */
fluid_setdirich_sd(actfield);

/*------------------------- calculate constants for nonlinear iteration */
/*
nik <->  EVALUATION OF LHS-MATRICES (w/o NONLINEAR TERM)    
nic <->  EVALUATION OF NONLINEAR LHS N-CONVECTIVE	    
nir <->  EVALUATION OF NONLINEAR LHS N-REACTION 	    
nie <->  EVALUATE ONLY LHS-TERMS FOR EXPLICIT VELOCITY      
nil <->  EVALUATION OF LUMPED MASS MATRIX (Mvv-lumped)      
nif <->  EVALUATION OF "TIME - RHS" (F-hat)		    
nii <->  EVALUATION OF "ITERATION - RHS"		    
nis <->  STATIONARY CASE (NO TIMEDEPENDENT TERMS) */
fdyn->nik = 1;
fdyn->nic = 2;
fdyn->nir = 0;
fdyn->nie = 0;
fdyn->nil = 0;
fdyn->nif = 0;
fdyn->nii = 0;
fdyn->nis = 0;
/*--------------------------------- calculate curvature at free surface */
if (fdyn->surftens!=0) 
{
   dserror("steepest descent method with free surface not yet implemented.");
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
container.ftimerhs     = ftimerhs;  /* not used here */
container.fiterhs      = fiterhs;
container.global_numeq = numeq_total;
container.nii          = fdyn->nii;
container.nif          = fdyn->nif;
container.kstep        = 0;
container.fieldtyp     = actfield->fieldtyp;
container.is_relax     = 1;
calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
       &container,action);
te=ds_cputime()-t1;
tes+=te;	     

/*--------------------------------------------------------------------- *
 | build the actual rhs-vector:                                         |
 |        rhs = fiterhs                                      |
 *----------------------------------------------------------------------*/
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
fdyn->ishape=0;

/*--------- return solution to the nodes to increment vector place 7 ---*/
solserv_result_incre(
		     actfield,
		     actintra,
		     &(actsolv->sol[actsysarray]),
		     7,
		     &(actsolv->sysarray[actsysarray]),
		     &(actsolv->sysarray_typ[actsysarray]),
		     0);     

/*------------------------- calculate stresses transferred to structure */
if (fsidyn->ifsi>0)
{
   *action = calc_fluid_stress;
   container.nii= 0;
   container.nif= 0;
   container.str=str;
   container.is_relax = 1;
   calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
          &container,action);

   /* since stresses are stored locally at the element it's necassary to 
      reduce them to all procs! ----------------------------------------*/
   dsassert(actsolv->parttyp==cut_elements,"Stress reduction for 'cut_nodes' not possible\n");
   fluid_reducestress(actintra,actfield,fdyn->numdf,str);

   /*----------------------------------------- store stresses in sol_mf */
   solserv_sol_zero(actfield,0,3,1);
   fsi_fluidstress_result(actfield,fdyn->numdf);
}

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
if (ioflags.fluid_vis_file==1)
{
   if (pssstep!=0)
   {
      /*------------------------------------------ store time in time_a */
      if (actpos >= time_a.fdim)
      amredef(&(time_a),time_a.fdim+1000,1,"DV");
      time_a.a.dv[actpos] = fdyn->acttime;   
   }   
   if (par.myrank==0) visual_writepss(actfield,actpos+1,&time_a);
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
amdel(&time_a);
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
