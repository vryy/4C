/*!----------------------------------------------------------------------
\file
\brief ls_fluid.c

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../fluid_full/fluid_prototypes.h"
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct  _FIELD          *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct  _GENPROB         genprob;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct  _SOLVAR         *solv;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct  _PARTITION      *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct  _IO_FLAGS        ioflags;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct  _PAR             par;
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
extern INT                      numcurve;
extern struct _CURVE           *curve;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum   _CALC_ACTION      calc_action[MAXFIELD];
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES            allfiles;



static INT             numff;              /* field number for fluid field     */
static INT             numls;              /* field number for level set field */
static INT             numeq;              /* number of equations on this proc */
static INT             numeq_total;        /* total number of equations        */
static INT             init;               /* flag for solver_control call     */
static INT             nsysarray = 1;      /* one system matrix                */
static INT             actsysarray = 0;    /* number of actual sysarray        */
static INT             outstep;            /* counter for output control       */
static INT             resstep = 0;        /* counter for output control       */
static INT             pssstep;            /* counter for output control       */
static INT             restartstep = 0;    /* counter for restart control      */
static INT             nfrastep;           /* number of steps for fractional-
                                              step-theta procedure             */
static INT             restart;

static INT             actpos;             /* actual position in sol. history  */
static DOUBLE          grat;               /* convergence ratios               */
static DOUBLE          tes = ZERO;         /*				       */
static DOUBLE          tss = ZERO;         /*				       */
static DOUBLE          tts = ZERO;         /*				       */
static FIELD          *actfield;           /* pointer to active field          */
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
static FILE           *out;
static FLUID_STRESS    str;
static FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables          */
static LS_DYNAMIC     *lsdyn;              /* ls dynamic variables             */



/*!----------------------------------------------------------------------
\brief control subroutine for fluid sub-problem in coupled level set /
extended finite element problem

<pre>                                                            irhan 05/04
control subroutine for fluid sub-problem in coupled level set /
extended finite element problem
</pre>

*----------------------------------------------------------------------*/
void ls_fluid(
  INT     eflag
  )
{
#ifdef DEBUG
  dstrc_enter("ls_fluid");
#endif
/*----------------------------------------------------------------------*/

  switch (eflag)
  {
/*-------------------------------------------- initialize fluid problem */
      case 1:
        ls_fluid_init();
        break;
/*------------------------------------------------- solve fluid problem */
      case 2:
        ls_fluid_solv();
        break;
/*---------------------------------------------- finalize fluid problem */
      case 3:
        ls_fluid_fina();
        break;
/*--------------------------------------------------------------- clean */
      case 99:
        ls_fluid_clea();
        break;
/*------------------------------------------------------------- default */
      default:
        dserror("action unknown\n");
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_fluid */



/*!----------------------------------------------------------------------
\brief initialization of sub-problem fluid

<pre>                                                            irhan 05/04
initialization of sub-problem fluid
</pre>

*----------------------------------------------------------------------*/
void ls_fluid_init()
{
  INT     i;
  INT     actcurve;

#ifdef DEBUG
  dstrc_enter("ls_fluid_init");
#endif
/*----------------------------------------------------------------------*/

  /* I N I T I A L I S A T I O N */
  /* set */
  numff = genprob.numff;
  numls = genprob.numls;
  fdyn   = alldyn[numff].fdyn;
  lsdyn = alldyn[numls].lsdyn;

  fdyn->dt = lsdyn->dt;
  fdyn->maxtime = lsdyn->maxtime;
  fdyn->nstep = lsdyn->nstep;
  grat = ZERO;
  dsassert(fdyn->iop==4,"TIMEINTEGR for fluid: only ONE-STEP-THETA implemented!\n");

  /* initialiase some counters */
  outstep=0;
  pssstep=0;
  restartstep=0;

  /* set file pointer out */
  out = allfiles.out_out;

  /* set pointers */
  actfield = &(field[numff]);
  actsolv = &(solv[numff]);
  actpart = &(partition[numff]);
  action = &(calc_action[numff]);
  restart = genprob.restart;
  container.actndis = 0;
  container.turbu = fdyn->turbu;
  container.fieldtyp = actfield->fieldtyp;

  /* set flag for stress evaluation */
  str = str_none;

  if (fdyn->liftdrag!=0)
    str = str_liftdrag;

/**************************BE CAREFUL************************************/
/**************************BE CAREFUL************************************/
/**************************BE CAREFUL************************************/
  if (ioflags.fluid_stress_gid==1)
  {
    dserror("ioflags.fluid_stress_gid=1 but not implemented!");
    /* str = str_all;*/
  }
/**************************BE CAREFUL************************************/
/**************************BE CAREFUL************************************/
/**************************BE CAREFUL************************************/

  fdyn->acttime = ZERO;

  /*
   * if we are not parallel, we have to allocate an alibi
   * intra-communicator structure
   */
#ifdef PARALLEL
  actintra    = &(par.intra[numff]);
#else
  actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  actintra->intra_fieldtyp = fluid;
  actintra->intra_rank     = 0;
  actintra->intra_nprocs   = 1;
#endif
  /*
   * there are only procs allowed in here, that belong to the fluid
   * intracommunicator (in case of nonlinear fluid. dyn., this should be all)
   */
  if (actintra->intra_fieldtyp != fluid) goto end;

  /* init the dist sparse matrices to zero */
  solserv_zero_mat(
    actintra,
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->sysarray_typ[actsysarray])
    );

  /* get global and local number of equations */
  solserv_getmatdims(
    &(actsolv->sysarray[actsysarray]),
    actsolv->sysarray_typ[actsysarray],
    &numeq,
    &numeq_total
    );

  /* output to the screen */
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

  /* write general data to .out */
  if (par.myrank==0)
  {
    fprintf(out,"max. values for:\n");
    fprintf(out," step |    time    | ite|   ite tol  | steady tol |\n");
    fprintf(out,"---------------------------------------------------\n");
    fprintf(out,"%5d | %10.3E | %2d | %10.3E | %10.3E |\n",
            fdyn->nstep,fdyn->maxtime,fdyn->itemax,fdyn->ittol,fdyn->sttol);
    fprintf(out,"---------------------------------------------------\n");
    fprintf(out,"\n");
    fprintf(out," step |    time    | ite| vel. error | pre. error |\n");
    fprintf(out,"---------------------------------------------------\n");
  }

  /* allocate dist. vectors 'rhs' */
  if (fdyn->iop == 7) actsolv->nrhs = 2; /* two dist. rhs vecs for BDF2	*/
  else actsolv->nrhs = 1;
  solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
  solserv_zero_vec(&(actsolv->rhs[0]));

  /* allocate dist. solution vectors */
  actsolv->nsol = 1;
  solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
  solserv_zero_vec(&(actsolv->sol[0]));

  /*
   * allocate one redundant vector ftimerhs of full lenght
   * this is used by the element routines to assemble the  Time RHS
   */
  ftimerhs = amdef("ftimerhs",&ftimerhs_a,numeq_total,1,"DV");

  /*
   * allocate one redundant vector fiterhs of full lenght
   * this is used by the element routines to assemble the  Iteration RHS
   */
  fiterhs = amdef("fiterhs",&fiterhs_a,numeq_total,1,"DV");

  /* allocate one vector for storing the time */
  if (ioflags.fluid_vis_file==1 )
    amdef("time",&time_a,1000,1,"DV");

  /* initialise fluid field */
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

  fluid_init(actpart,actintra,actfield,action,&container,4,str);
  actpos=0;

  /* initialize multilevel algorithm */
#if defined(FLUID2_ML) || defined(FLUID3_ML)
  if (fdyn->mlfem==1) fluid_ml_init(actfield);
#endif

  /* init all applied time curves -*/
  for (actcurve=0; actcurve<numcurve; actcurve++)
    dyn_init_curve(actcurve,fdyn->nstep,fdyn->dt,fdyn->maxtime);

  /* init the dirichlet-conditions -*/
  fluid_initdirich(actfield);

  /* initialize solver on all matrices */

  /*
    NOTE => solver init phase has to be called with each matrix one wants
    to solve with. Solver init phase has to be called with all matrices
    one wants to do matrix-vector products and matrix scalar products.
    This is not needed by all solver libraries, but the solver-init
    phase is cheap in computation (can be costly in memory)
  */

  /* initialize solver */
  init=1;
  solver_control(
    actsolv, actintra,
    &(actsolv->sysarray_typ[actsysarray]),
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->sol[0]),
    &(actsolv->rhs[0]),
    init
    );

  /* init the assembly for stiffness */
  init_assembly(actpart,actsolv,actintra,actfield,actsysarray,0);

  /* allocate fluid integration data */
  alldyn[genprob.numff].fdyn->data = (FLUID_DATA*)CCACALLOC(1,sizeof(FLUID_DATA));

  /* init the element calculating routines */
  *action = calc_fluid_init;
  calinit(actfield,actpart,action,&container);

  /* print out initial data to .out */
  out_sol(actfield,actpart,actintra,fdyn->step,actpos);

  /* print out initial data to .flavia.res */
  if (ioflags.fluid_sol_gid==1 && par.myrank==0)
  {
    out_gid_sol("velocity",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
    out_gid_sol("pressure",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
  }
  if (ioflags.fluid_stress_gid==1 && par.myrank==0)
  {
    out_gid_sol("stress",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
  }

  /* calculate time independent constants for time algorithm */
  fluid_cons();

 end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_fluid_init */



/*!----------------------------------------------------------------------
\brief solution of sub-problem fluid

<pre>                                                            irhan 05/04
solution of sub-problem fluid
</pre>

*----------------------------------------------------------------------*/
void ls_fluid_solv()
{
  INT        converged = 0;        /* convergence flag */
  INT        itnum = 0;            /* iteration number */
  DOUBLE     vrat,prat;
  DOUBLE     t1,ts,te;

#ifdef DEBUG
  dstrc_enter("ls_fluid_solv");
#endif
/*----------------------------------------------------------------------*/

  /* check (starting) algorithm */
  if (fdyn->step<=(fdyn->nums+1)) fluid_startproc(&nfrastep,0);

  /* calculate constants for time algorithm */
  fluid_tcons();
  /* set new absolute time */
  if (fdyn->iop == 1) /* generalised alpha is solved for time n+alpha_f */
    fdyn->acttime += fdyn->dta * fdyn->alpha_f;
  else
    fdyn->acttime += fdyn->dta;

  /* output to the screen */
  if (par.myrank==0) fluid_algoout();

  /* set dirichlet boundary conditions to sol_increment[3] */
  fluid_setdirich(actfield,3);

  /* initialise timerhs */
  amzero(&ftimerhs_a);

  /* prepare time rhs in mass form */
  if(fdyn->time_rhs==0)
    fluid_prep_rhs(actfield);

  /* start time step on the screen */
  if (fdyn->itnorm!=fncc_no && par.myrank==0)
  {
    printf("----------------------------------------------------------------\n");
    printf("|- step/max -|-  tol     [norm] -|- vel. error -|- pre. error -| \n");
  }
  itnum=1;

/*======================================================================*
 |           N O N L I N E A R   I T E R A T I O N                      |
 *======================================================================*/
 nonlniter:

  /* calculate constants for nonlinear iteration */
  fluid_icons(itnum);

  /* intitialise global matrix and global rhs */
  solserv_zero_vec(&(actsolv->rhs[0]));
  solserv_zero_mat(
    actintra,&(actsolv->sysarray[actsysarray]),
    &(actsolv->sysarray_typ[actsysarray])
    );

  /* initialise iterations-rhs */
  amzero(&fiterhs_a);

  /* form incremental matrices, residual and element forces */
  *action = calc_fluid;
  t1=ds_cputime();
  container.dvec         = NULL;
  container.ftimerhs     = ftimerhs;
  container.fiterhs      = fiterhs;
  container.global_numeq = numeq_total;
  container.nii          = fdyn->nii;
  container.nif          = fdyn->nif;
  container.nim          = fdyn->nim;
  container.kstep        = 0;
  container.is_relax     = 0;
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
         &container,action);
  te=ds_cputime()-t1;
  tes+=te;

/*--------------------------------------------------------------------- *
 | build the actual rhs-vector:                                         |
 |        rhs = ftimerhs + fiterhs 				        |
 *----------------------------------------------------------------------*/
  /* add time-rhs: */
  if (fdyn->time_rhs)
    assemble_vec(
      actintra,
      &(actsolv->sysarray_typ[actsysarray]),
      &(actsolv->sysarray[actsysarray]),
      &(actsolv->rhs[0]),
      ftimerhs,
      1.0
      );
  /* add iteration-rhs: */
  assemble_vec(
    actintra,
    &(actsolv->sysarray_typ[actsysarray]),
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->rhs[0]),
    fiterhs,
    1.0
    );

  /* solve system */
  init=0;
  t1=ds_cputime();
  solver_control(
    actsolv, actintra,
    &(actsolv->sysarray_typ[actsysarray]),
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->sol[0]),
    &(actsolv->rhs[0]),
    init
    );
  ts=ds_cputime()-t1;
  tss+=ts;

  /* set flags for stability parameter evaluation and convergence check */
  fdyn->ishape=0;

  /* return solution to the nodes and calculate the convergence ratios */
  fluid_result_incre(
    actfield,actintra,&(actsolv->sol[0]),3,
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->sysarray_typ[actsysarray]),
    &vrat,&prat,NULL
    );

  /* iteration convergence check */
  converged = fluid_convcheck(vrat,prat,ZERO,itnum,te,ts);

  /* check if nonlinear iteration has to be finished */
  if (converged==0)
  {
    itnum++;
    goto nonlniter;
  }
/*----------------------------------------------------------------------*
 | -->  end of nonlinear iteration                                      |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_fluid_solv */



/*!----------------------------------------------------------------------
\brief finalization of sub-problem fluid

<pre>                                                            irhan 05/04
finalization of sub-problem fluid
</pre>

*----------------------------------------------------------------------*/
void ls_fluid_fina()
{
  DOUBLE     t2,tt;

#ifdef DEBUG
  dstrc_enter("ls_fluid_fina");
#endif
/*----------------------------------------------------------------------*/

  /* extrapolate from n+alpha_f to n+1 for generalised alpha */
  if (fdyn->iop == 1)
  {
    solserv_sol_zero(actfield,0,1,2);
    solserv_sol_add(actfield,0,1,1,3,2,1.0/fdyn->alpha_f);
    solserv_sol_add(actfield,0,1,1,1,2,1.0-1.0/fdyn->alpha_f);
    solserv_sol_copy(actfield,0,1,1,2,3);

    fdyn->acttime += fdyn->dta * (1.0 - fdyn->alpha_f);
  }

  /* steady state check */
/*if (fdyn->stchk==iststep)
  {
    iststep=0;
    steady = fluid_steadycheck(actfield,numeq_total);
    }*/

  /* lift&drag computation */
  if (fdyn->liftdrag>0)
  {
    container.str = str;
    *action = calc_fluid_liftdrag;
    fluid_liftdrag(
      1,action,&container,actfield,
      actsolv,actpart,actintra
      );
  }

  /* stress computation */
  if (ioflags.fluid_stress_gid==1)
  {
    container.str = str;
    *action = calc_fluid_stress;
    calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
           &container,action);
  }

  /*
    copy solution from sol_increment[1][j] to sol_increment[0][j]
    -> prev. solution becomes (n-1)-solution of next time step
  */
  solserv_sol_copy(actfield,0,1,1,1,0);

  /*
    copy solution from sol_increment[3][j] to sol_increment[1][j]
    -> actual solution becomes previous solution of next time step
  */
  solserv_sol_copy(actfield,0,1,1,3,1);

  /* finalise this timestep */
  outstep++;
  pssstep++;
  resstep++;
  restartstep++;

  /* write solution to .pss */
  if (pssstep==fdyn->uppss && ioflags.fluid_vis_file==1)
  {
    pssstep=0;
    /* store time in time_a */
    if (actpos >= time_a.fdim)
      amredef(&(time_a),time_a.fdim+1000,1,"DV");
    time_a.a.dv[actpos] = fdyn->acttime;
    actpos++;
  }

  /*
   * copy solution from sol_increment[3][j] to sol_[actpos][j]
   * and transform kinematic to real pressure
   */
  solserv_sol_copy(actfield,0,1,0,3,actpos);
/**************************BE CAREFUL************************************/
/**************************BE CAREFUL************************************/
/**************************BE CAREFUL************************************/

/*  fluid_transpres(actfield,0,0,actpos,fdyn->numdf-1,0);*/

/**************************BE CAREFUL************************************/
/**************************BE CAREFUL************************************/
/**************************BE CAREFUL************************************/

  /* copy solution on level 2 at (n+1) to place (n) for multi-level FEM */
#if defined(FLUID2_ML) || defined(FLUID3_ML)
  if (fdyn->mlfem==1) fluid_smcopy(actpart);
#endif

  /* write solution to .flavia.res */
  if (resstep==fdyn->upres && par.myrank==0)
  {
    resstep=0;
    /*out_checkfilesize(1);*/

    if(ioflags.fluid_sol_gid==1)
    {
      out_gid_sol("velocity",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
      out_gid_sol("pressure",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
    }
    if(ioflags.fluid_stress_gid==1)
    {
      out_gid_sol("stress",actfield,actintra,fdyn->step,actpos,fdyn->acttime);
    }
  }

  /* write solution to .out */
  if (outstep==fdyn->upout && ioflags.fluid_sol_file==1)
  {
    outstep=0;
    out_sol(actfield,actpart,actintra,fdyn->step,actpos);
  }

  /* write restart to pss file */
  if (restartstep==fdyn->uprestart)
  {
    restartstep=0;
    restart_write_fluiddyn(fdyn,actfield,actpart,actintra,action,&container);
  }

  tt=ds_cputime()-t2;
  tts+=tt;
  printf("PROC  %3d | total time for this time step: %10.3e \n",par.myrank,tt);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_fluid_solv */



/*!----------------------------------------------------------------------
\brief end of sub-problem fluid

<pre>                                                            irhan 05/04
end of sub-problem fluid
</pre>

*----------------------------------------------------------------------*/
void ls_fluid_clea()
{
  INT     i;

#ifdef DEBUG
  dstrc_enter("ls_fluid_clea");
#endif
/*----------------------------------------------------------------------*/

  if (pssstep==0) actpos--;
  /* print out solution to .out file */
  if (outstep!=0 && ioflags.fluid_sol_file==1)
    out_sol(actfield,actpart,actintra,fdyn->step,actpos);

  /* print out solution to 0.pss file */
  if (ioflags.fluid_vis_file==1)
  {
    if (pssstep!=0)
    {
      /* store time in time_a */
      if (actpos >= time_a.fdim)
        amredef(&(time_a),time_a.fdim+1000,1,"DV");
      time_a.a.dv[actpos] = fdyn->acttime;
    }
    if (par.myrank==0) visual_writepss(actfield,actpos+1,&time_a);
  }

  /* print total CPU-time to the screen */
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
  for (i=0;i<par.nprocs;i++)
  {
    if (par.myrank==i)
    {
      printf("\n");
      printf("PROC  %3d | FIELD FLUID     | total time element for calculations: %10.3E \n",
             par.myrank,tes);
      printf("PROC  %3d | FIELD FLUID     | total time for solver              : %10.3E \n",
             par.myrank,tss);
      printf("PROC  %3d | FIELD FLUID     | total time for time loop           : %10.3E \n",
             par.myrank,tts);
    }

    /* write total cpu-time to .out */
    if (par.myrank==0)
    {
      fprintf(out,"\n");
      fprintf(out," total time element for calculations: %10.3E \n", tes);
      fprintf(out," total time for solver              : %10.3E \n", tss);
      fprintf(out," total time for time loop           : %10.3E \n", tts);
    }
  }
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif

  /* cleaning up phase */
  amdel(&ftimerhs_a);
  amdel(&fiterhs_a);
  if (ioflags.fluid_vis_file==1 )
    amdel(&time_a);
  solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
  solserv_del_vec(&(actsolv->sol),actsolv->nsol);

/*----------------------------------------------------------------------*/
#ifndef PARALLEL
  CCAFREE(actintra);
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_fluid_clea */
/*! @} (documentation module close)*/
#endif
