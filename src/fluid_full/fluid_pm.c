/*!
\file
\brief projection method algorithm for fluid

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

This is a new implementation of the projection method as described in
GRESHO, Philip M. "On the theory of Semi-Implicit Projection Methods
for viscous incompressible flow and its implementation via a finite
element method that also introduces a nearly consistent mass matrix
part 1 & part 2" International Journal for Numerical Methods in
Fluids, VOL 11, pp. 587-659, (1990).

We use discontinous pressure here with pressure dofs connected to the
elements. The elements themselves are inf-sub stable, no pressure
stabilization needed. Thus there is only one discretization. The
pressure equation therefore cannot be attached to a matching
discretization but must be expressed with stand alone sparse matrices.
Additional these matrices must support the weird C^T*ML^-1*C
product. This motivated a new sparse matrix type. The existing ones
are not suitable because those are directly connected to
discretizations, elements, nodes and such.

One issue with the algorithm is parallelization. It does work. And it
does make things complicated. In the end, because the pressure values
belong to the elements here, we have to communicate element
values. This differs from the normal behavior of ccarat where we
communicate nodal data.

The current implementation is a working version, however it is not
optimized for speed. The fluid momentum equation is treated fully
implicit. Work has to be done in this area to obtain a very fast
algorithm.

*/
/*!
\addtogroup FLUID_PM
*//*! @{ (documentation module open)*/

#ifdef D_FLUID_PM

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../solver/solver_sparse.h"
#include "fluid_prototypes.h"
#include "fluid_pm_prototypes.h"
#include "../io/io.h"

#include "../fluid3_pro/fluid3pro.h"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
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

/*----------------------------------------------------------------------*
 | routine to control implicit and semi-implicit algorithms for fluid   |
 | problems combined with Newton and fixed point iteration schemes.     |
 | IN PROGRESS: ONE-STEP-THETA                                          |
 |              fixed point like iteration                              |
 |              only Euler no ALE                                       |
 |                                                          genk  03/02 |
 *----------------------------------------------------------------------*/
/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;




/*----------------------------------------------------------------------*/
/*!
  \brief the projection method

  In essence a copy of the global fluid algorithm.

  We use discontinous pressure and just one discretization. The
  pressure dofs are stored by the elements.

  The old version done by basol has been merged with the new code.

  \author u.kue
  \date 12/05
 */
/*----------------------------------------------------------------------*/
void fluid_pm()
{
  INT             itnum;        /* counter for nonlinear iteration  */
  INT             i;            /* simply a counter                 */
  INT             numeq;        /* number of equations on this proc */
  INT             numeq_total;  /* total number of equations        */
  INT             init;         /* flag for solver_control call     */
  INT             actsysarray=0; /* number of actual sysarray        */
  INT             outstep=0;    /* counter for output control       */
  INT             resstep=0;    /* counter for output control       */
  INT             pssstep=0;    /* counter for output control	*/
  INT             restartstep=0;
  INT             iststep=0;    /* counter for time integration     */
  INT             nfrastep;     /* number of steps for fractional-
                                   step-theta procedure             */
  INT             actcurve;     /* actual timecurve                 */
  INT             converged=0;  /* convergence flag                 */
  INT             steady=0;     /* flag for steady state            */
  INT             actpos;       /* actual position in sol. history  */
  INT             restart;
  INT		repeat=0;       /* flag, repeat time step?		*/
  INT		repeated=0;     /* has time step been repeated?	*/
  DOUBLE          vrat;         /* convergence ratio                */
  DOUBLE          t1,t2,ts,te,tt;
  DOUBLE          tes=0.0;
  DOUBLE          tss=0.0;
  DOUBLE          tts=0.0;
  DOUBLE          fact1,fact2;
  FLUID_STRESS    str;

  SOLVAR         *actsolv;      /* pointer to active sol. structure */
  PARTITION      *actpart;      /* pointer to active partition      */
  FIELD          *actfield;     /* pointer to active field          */
  INTRA          *actintra;     /* pointer to active intra-communic.*/
  CALC_ACTION    *action;       /* pointer to the cal_action enum   */

  ARRAY           frhs_a;
  DOUBLE         *frhs;         /* iteration - RHS                  */
  ARRAY           fgradprhs_a;
  DOUBLE         *fgradprhs;

  CONTAINER       container;    /* contains variables defined in container.h */
  FILE           *out = allfiles.out_out;

  ARRAY_POSITION *ipos;
  FLUID_DYNAMIC  *fdyn;

  INT             disnum_calc;
  INT             disnum_io;

  INT numpdof;

  PARALLEL_SPARSE grad;
  PARALLEL_SPARSE pmat;

  SOLVAR presolv;
  INT pnumeq_total;
  INT pnumeq;

  DOUBLE* lmass;

  INT stiff_array;              /* indice of the active system sparse matrix */
  INT mass_array;               /* indice of the active system sparse matrix */

#ifdef BINIO
  BIN_OUT_FIELD   out_context;
  BIN_OUT_FIELD   restart_context;
#endif

#ifdef DEBUG
  dstrc_enter("fluid_pm");
#endif

  /*                    I N I T I A L I S A T I O N                       */

  fdyn = alldyn[genprob.numff].fdyn;

  actfield    = &(field[0]);
  actsolv     = &(solv[0]);
  actpart     = &(partition[0]);
  action      = &(calc_action[0]);
  restart     = genprob.restart;

#ifdef SUBDIV
  if (actfield->subdivide > 0)
  {
    disnum_calc = 1;
    disnum_io   = ioflags.output_dis;
  }
  else
#endif
    disnum_calc = disnum_io = 0;

  /*
   * This is again highly dangerous! We need an additional mass matrix
   * and assume we can enlarge the solver's matrix array by one to put
   * it there. We also assume that the current stiffness matrix is the
   * last matrix in that array.
   *
   * Okay. The stiffness matrix used to be the only one in fluid
   * calculations. No need to worry. And the subdiscretization stuff
   * seems to put the io discretization (that does not need an mass
   * matrix, no matrix actually) before the calculation
   * discretization. Seems fine. But don't mess around with these
   * discretizations any further.
   *
   * These data structures need a clean-up badly! */
  actsysarray = disnum_calc;
  stiff_array = disnum_calc;
  mass_array  = disnum_calc+1;

  container.disnum    = disnum_calc;
  container.turbu     = fdyn->turbu;
  container.fieldtyp  = actfield->fieldtyp;

  ipos = &(actfield->dis[disnum_calc].ipos);

  /* set flag for stress evaluation */
  str         = str_none;

  fdyn->acttime=0.;

  /*---------------- if we are not parallel, we have to allocate an alibi *
    ---------------------------------------- intra-communicator structure */
#ifdef PARALLEL
  actintra    = &(par.intra[0]);
#else
  actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  actintra->intra_fieldtyp = fluid;
  actintra->intra_rank     = 0;
  actintra->intra_nprocs   = 1;
#endif

  /*- there are only procs allowed in here, that belong to the fluid -----*/
  /* intracommunicator (in case of nonlinear fluid. dyn., this should be all)*/
  if (actintra->intra_fieldtyp != fluid) goto end;

  /* ------------------------------------------------ */
  /* Enumerate discontinous pressure dofs */

  numpdof = pm_assign_press_dof(actfield, actpart, disnum_calc, actintra);

  /* ------------------------------------------------ */
  /* Allocate space for the additional mass matrix
   * This follows the structure code. The same ugliness here. */

  /* stiff_array already exists, so copy the mask of it to mass_array
   * reallocate the vector of sparse matrices and the vector of there
   * types formerly lenght 1, now lenght 2
   * (Not quite right because of the subdiscretizations.) */
  actsolv->nsysarray += 1;
  actsolv->sysarray_typ = (SPARSE_TYP*)CCAREALLOC(actsolv->sysarray_typ,
                                                  actsolv->nsysarray*sizeof(SPARSE_TYP));
  actsolv->sysarray = (SPARSE_ARRAY*)CCAREALLOC(actsolv->sysarray,
                                                actsolv->nsysarray*sizeof(SPARSE_ARRAY));

  /* copy the matrices sparsity mask from stiff_array to mass_array */
  solserv_alloc_cp_sparsemask(actintra,
                              &(actsolv->sysarray_typ[stiff_array]),
                              &(actsolv->sysarray[stiff_array]),
                              &(actsolv->sysarray_typ[mass_array]),
                              &(actsolv->sysarray[mass_array]));

  /* inverted lumed masses in vector form */
  lmass = (DOUBLE*)CCAMALLOC(actfield->dis[disnum_calc].numeq*sizeof(DOUBLE));

  /* ------------------------------------------------ */
  /* The gradient matrix G */
  /* This is a parallel matrix, each processor owns a separate slice. */

  /* the discontinuous pressure version */
  parallel_sparse_init(&grad,
                       actfield->dis[disnum_calc].numeq,
                       numpdof*actfield->dis[disnum_calc].numele,
                       sd_slice_vertical,
                       numpdof*actpart->pdis[disnum_calc].numlele);

#ifdef PARALLEL

  /* We need to provide the update array ourselves. */
  pm_fill_gradient_update(actpart, disnum_calc, actintra, numpdof, &grad);

#endif

  /* Now build the mask. */
  pm_gradient_mask_mat(actfield, actpart, disnum_calc, actintra, numpdof, &grad);

  /* ------------------------------------------------ */
  /* The pressure matrix D*ML^-1*G */
  /* With D^T = G */

  /* This one is horizontally sliced. In the huge matrix there are as
   * many rows and columns as there are elements times the number of
   * pressure dofs per element. No dirichlet conditions on the
   * pressure supported right now. */
  parallel_sparse_init(&pmat,
                       numpdof*actpart->pdis[disnum_calc].numlele,
                       numpdof*actfield->dis[disnum_calc].numele,
                       sd_slice_horizontal,
                       numpdof*actfield->dis[disnum_calc].numele);

#ifdef PARALLEL

  /*
   * Here we slice horizontally and thus we have to remember the
   * global row numbers */

  /*
   * We have the same slicing here as in the gradient matrix, only the
   * direction changed. No need to repeat the above loop. */
  memcpy(pmat.update, grad.update, pmat.slice.rows*sizeof(INT));

#endif

  /* create a matrix suitable for a solver */
  /* We need two copies anyway, because our solvers are allowed to
   * destroy their matrix on solving. So in order to keep the
   * matrix-matrix multiplication simple and independent of the solver
   * we first create an independent sparse and convert it afterwards. */

  pm_build_pmat_sparse_mask(actfield, actpart, disnum_calc, actintra, numpdof, &pmat);

  /* Convert the sparse matrix to something solvable. :) */
  /* Concerns the sparse mask only. */
  parallel_sparse_convert_aztec(&pmat, &presolv);

  /* ------------------------------------------------ */

  /*---------------------------- get global and local number of equations */
  solserv_getmatdims(&(presolv.sysarray[0]),
                     presolv.sysarray_typ[0],
                     &pnumeq,
                     &pnumeq_total);

  /*---------------------------------------- allocate dist. vectors 'rhs' */
  presolv.nrhs = 1;
  solserv_create_vec(&(presolv.rhs),presolv.nrhs,pnumeq_total,pnumeq,"DV");
  solserv_zero_vec(&(presolv.rhs[0]));

  /*---------------------------------- allocate dist. solution vectors ---*/
  presolv.nsol = 1;
  solserv_create_vec(&(presolv.sol),presolv.nsol,pnumeq_total,pnumeq,"DV");
  solserv_zero_vec(&(presolv.sol[0]));

  /* Setup the solver for the pressure matrix. That is again special
   * because this matrix is not (directly) connected to a
   * discretization. */
  solver_control(actfield,disnum_calc,
                 &presolv, actintra,
                 &(presolv.sysarray_typ[0]),
                 &(presolv.sysarray[0]),
                 &(presolv.sol[0]),
                 &(presolv.rhs[0]),
                 1);

  /*------------------------------------ prepare lift&drag calculation ---*/
  if (fdyn->liftdrag==ld_stress)
    str         = str_liftdrag;
  if (fdyn->liftdrag==ld_nodeforce)
    fluid_liftdrag(-1,action,&container,actfield,disnum_calc,
                   actsolv,actpart,actintra,ipos);

  if (ioflags.fluid_stress==1)
    str         = str_all;

  /*------------------------------- init the dist sparse matrices to zero */
  solserv_zero_mat(
    actintra,
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->sysarray_typ[actsysarray])
    );

  /*---------------------------- get global and local number of equations */
  solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                     actsolv->sysarray_typ[actsysarray],
                     &numeq,
                     &numeq_total);

  /*---------------------------------------- allocate dist. vectors 'rhs' */
  actsolv->nrhs = 2;
  solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
  solserv_zero_vec(&(actsolv->rhs[0]));

  /*---------------------------------- allocate dist. solution vectors ---*/
  actsolv->nsol = 1;
  solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
  solserv_zero_vec(&(actsolv->sol[0]));

  /*---------------  allocate one redundant vector frhs of full lenght ---*/
  /* this is used by the element routines to assemble the Iteration RHS */
  frhs = amdef("frhs",&frhs_a,numeq_total,1,"DV");
  fgradprhs = amdef("fgradprhs",&fgradprhs_a,numeq_total,1,"DV");

  /* screen output */
  pm_out_screen_header(numeq, numeq_total, actintra, out, fdyn);

  /*--------------------------------------------- initialise fluid field */
  if (restart != 0)
  {
    if (fdyn->init>0)
      /*dserror("Initial field either by restart, or by function or from file ...\n");*/
      printf("Restart: Initial field not used!!!\n\n\n");

    fdyn->resstep=genprob.restart;
    fdyn->init=2;
  }

  fluid_init_pos_euler(ipos);
  fluid_init(actpart, actintra, actfield, disnum_calc, disnum_io, action, &container, 8,ipos,str);
  actpos=0;

  /*------------------------------------ initialize multilevel algorithm */
#if defined(FLUID2_ML) || defined(FLUID3_ML)
  if (fdyn->mlfem==1)
  {
    dserror("no multi level during projection method");
  }
#endif

  /*--------------------------------------- init all applied time curves -*/
  for (actcurve=0; actcurve<numcurve; actcurve++)
    dyn_init_curve(actcurve,fdyn->nstep,fdyn->dt,fdyn->maxtime);

  /*-------------------------------------- init the dirichlet-conditions -*/
  fluid_initdirich(actfield, disnum_calc, ipos);

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
  solver_control(actfield,disnum_calc,actsolv, actintra,
                 &(actsolv->sysarray_typ[actsysarray]),
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sol[0]),
                 &(actsolv->rhs[0]),
                 init);
  /* initialize the mass matrix, too, so we can use the global
   * matrix-vector product */
  solver_control(actfield,disnum_calc,actsolv, actintra,
                 &(actsolv->sysarray_typ[mass_array]),
                 &(actsolv->sysarray[mass_array]),
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

#ifdef BINIO

  /* initialize binary output
   * It's important to do this only after all the node arrays are set
   * up because their sizes are used to allocate internal memory. */
  init_bin_out_field(&out_context,
                     &(actsolv->sysarray_typ[actsysarray]),
                     &(actsolv->sysarray[actsysarray]),
                     actfield, actpart, actintra, disnum_io);

  if (disnum_io != disnum_calc)
    init_bin_out_field(&restart_context,
                       &(actsolv->sysarray_typ[actsysarray]),
                       &(actsolv->sysarray[actsysarray]),
                       actfield, actpart, actintra, disnum_calc);

#endif

  /*--------------------------------------------- calculate nodal normals */
  /* this is fluid2 only... ??? */
  /*fluid_cal_normal(actfield, disnum_calc, 1, action);*/

  /*------------------------------------------------- define local co-sys */
  fluid_locsys(actfield, disnum_calc, fdyn);

#ifdef SUBDIV
  /* transfer the solution to the nodes of the master-dis */
  if (actfield->subdivide > 0)
  {
    solserv_sol_trans(actfield, disnum_calc, node_array_sol, actpos);
  }
#endif

  /*-------------------------------------- print out initial data to .out */
  if (ioflags.output_out==1 && ioflags.fluid_sol==1 && par.myrank==0)
    out_sol(actfield,actpart,disnum_io,actintra,fdyn->step,actpos);

  /*------------------------------- print out initial data to .flavia.res */
  if (ioflags.output_gid==1 && par.myrank==0)
  {
    if (ioflags.fluid_sol==1)
    {
      out_gid_sol("velocity",actfield,disnum_io,actintra,fdyn->step,actpos,fdyn->acttime);
      out_gid_sol("average_pressure",actfield,disnum_io,actintra,fdyn->step,actpos,fdyn->acttime);
    }
    if (ioflags.fluid_stress==1)
    {
      out_gid_sol("stress",actfield,disnum_io,actintra,fdyn->step,actpos,fdyn->acttime);
    }
  }

  /*-------------------------------------------- write solution to binary */
#ifdef BINIO
  if (ioflags.output_bin==1)
  {
    if (ioflags.fluid_sol==1)
    {
      out_results(&out_context, fdyn->acttime, fdyn->step, actpos, OUTPUT_VELOCITY | OUTPUT_AV_PRESSURE);
    }
    if (ioflags.fluid_stress==1)
    {
      out_results(&out_context, fdyn->acttime, fdyn->step, actpos, OUTPUT_STRESS);
    }
  }
#endif

  /*---------------------------------------------------------- monitoring */
  if (ioflags.monitor==1)
  {
    out_monitor(actfield,genprob.numff,ZERO,1);
    monitoring(actfield,disnum_io,genprob.numff,actpos,fdyn->acttime);
  }

  /*---------- calculate time independent constants for time algorithm ---*/
  fluid_cons();

  dsmemreport();

  /* -------------------------------------------------------------------- */
  /* Now setup is done and we actually calculate the values of G and
   * M. */

  /* ------------------------------------------------ */
  /* Calculate gradient and mass matrices, including the inverted
   * lumped mass matrix. Velocity dirichlet conditions are observed. */
  solserv_zero_mat(actintra,&(actsolv->sysarray[mass_array]),
                   &(actsolv->sysarray_typ[mass_array]));
  sparse_zero(&grad.slice);
  pm_calelm(actfield, actpart, disnum_calc,
            actsolv, mass_array,
            actintra, ipos, numpdof, &grad, lmass);

  /* Calculate global mass matrix */

  /* matrix-matrix-matrix multiplication. */
  /* Ok. We have sparse matrices and the inner one is diagonal. */

  parallel_sparse_pm_matmat(&pmat, &grad, lmass, actintra);

  /*======================================================================*
   |                         T I M E L O O P                              |
   *======================================================================*/
  /* nodal solution history fluid field:                                  *
   * sol[0][j]           ... initial data                                 *
   * sol[0...actpos][j]  ... solution for visualisation (real pressure)   *
   * sol_increment[position][j] ... solution at some time level time      *
   * position flags:                                                      *
   *  ipos->velnp ... velocity at time (n+1)                              *
   *  ipos->veln  ... velocity at time (n)                                *
   *  ipos->velnm ... velocity at time (n-1)                              *
   *  ipos->accn  ... acceleration at time (n)                            *
   *  ipos->accnm ... acceleration at time (n-1)                          *
   *  ipos->hist  ... old solution data depending on time integration     *
   *  ipos->pred  ... predicted solution for new time level               *
   *  ipos->terr  ... local truncation error                              *
   *======================================================================*/

  /* timeloop */
  for (;;)
  {
    t2=ds_cputime();

    fdyn->step++;
    iststep++;

    if (par.myrank==0)
      fprintf(out,"%5d | %10.3E |",fdyn->step,fdyn->acttime);

    /*------------------------------------------ check (starting) algorithm */
    if (fdyn->step<=(fdyn->nums+1)) fluid_startproc(&nfrastep,0);

    /*------------------------------ calculate constants for time algorithm */
    fluid_tcons();

    /*-------------------------------------------- set new absolute time ---*/
    if (fdyn->iop == 1)/* generalised alpha is solved for time n+alpha_f 	*/
      fdyn->acttime += fdyn->dta * fdyn->alpha_f;
    else
      fdyn->acttime += fdyn->dta;

    /*------------------------------------------------ output to the screen */
    if (par.myrank==0) fluid_algoout();

    /*------------------------ predictor step for adaptive time stepping ---*/
    if (fdyn->adaptive && fdyn->step > 1)
    {
      fluid_predictor(actfield, disnum_calc, ipos,fdyn->iop);
    }
    else
    {
      /* do explicit predictor step to start iteration from better value */
      fact1 = fdyn->dta*(1.0+fdyn->dta/fdyn->dtp);
      fact2 = DSQR(fdyn->dta/fdyn->dtp);
      solserv_sol_add(actfield,0,node_array_sol_increment,
                      node_array_sol_increment,
                      ipos->accn,ipos->velnp, fact1);
      solserv_sol_add(actfield,0,node_array_sol_increment,
                      node_array_sol_increment,
                      ipos->veln,ipos->velnp,-fact2);
      solserv_sol_add(actfield,0,node_array_sol_increment,
                      node_array_sol_increment,
                      ipos->velnm,ipos->velnp,fact2);
    }
    /*-------- set dirichlet boundary conditions to sol_increment[velnp] ---*/
    fluid_setdirich(actfield, disnum_calc, ipos,ipos->velnp);

    /*------------------------------------ prepare time rhs in mass form ---*/
    fluid_prep_rhs(actfield, disnum_calc, ipos);

    /*------------------------------------- start time step on the screen---*/
    if (fdyn->itnorm!=fncc_no && par.myrank==0)
    {
      printf("-------------------------------------------------\n");
      printf("|- step/max -|-  tol     [norm] -|- vel. error -|\n");
    }
    itnum=1;

    /*======================================================================*
     |           N O N L I N E A R   I T E R A T I O N                      |
     *======================================================================*/

    /* nonlniter */
    for (;;)
    {
      /*------------------------- calculate constants for nonlinear iteration */
      /*fluid_icons(itnum); this is not needed any more */

      /*---------------------------- intitialise global matrix and global rhs */
      solserv_zero_vec(&(actsolv->rhs[0]));
      solserv_zero_vec(&(actsolv->rhs[1]));
      solserv_zero_mat(actintra,&(actsolv->sysarray[actsysarray]),
                       &(actsolv->sysarray_typ[actsysarray]));

      /*------------------------------------------- re-initialise neumann bcs */
      inherit_design_dis_neum(&(actfield->dis[disnum_calc]));

      /*------------------------------------------- initialise iterations-rhs */
      amzero(&frhs_a);
      amzero(&fgradprhs_a);

      /*-------------- form incremental matrices, residual and element forces */
#ifdef PERF
      perf_begin(81);
#endif

      *action = calc_fluid;
      t1=ds_cputime();
      container.dvec         = NULL;
      container.frhs         = frhs;
      container.fgradprhs    = fgradprhs;
      container.global_numeq = numeq_total;
      container.nii          = fdyn->nii;
      container.kstep        = 0;
      container.is_relax     = 0;
      calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
             &container,action);
      te=ds_cputime()-t1;
      tes+=te;

#ifdef PERF
      perf_end(81);
#endif

      /* Add the pressure term M*ML^-1*G*p(n) to the rhs. The vector
       * G*p(n) has already been calculated by the element call. */
      /* (This term could be done just once for each nonlinear
       * iteration, but maybe we'll do the convection explicit and
       * circumvent the iteration altogether.) */

      /*
       * We do ML^-1*G*p(n) at first because that results in a vector
       * and we are left with an ordinary matrix-vector
       * multiplication.
       *
       * We do the multiplication on the global level for
       * simplicity. But still we only handle those dofs that belong
       * to the local processor. */
      for (i=0; i<actpart->pdis[disnum_calc].numnp; ++i)
      {
        INT dof;
        NODE* actnode;
        actnode = actpart->pdis[disnum_calc].node[i];
        if (actnode->proc == actintra->intra_rank)
        {
          for (dof=0; dof<actnode->numdf; ++dof)
          {
            INT gdof;
            gdof = actnode->dof[dof];
            if (gdof < actfield->dis[disnum_calc].numeq)
            {
              fgradprhs[gdof] *= lmass[gdof];
            }
          }
        }
      }

      /*
       * Ok. We waste time and space here. The total vector fgradprhs
       * (only the values that belong to this processor are filled!)
       * is copied to a distributed one. The matrix-vector
       * multiplication allreduces it, does its work and distributes
       * it back.
       *
       * If this takes too much time we can optimize it, but not for
       * all matrix types supported here. */
      assemble_vec(actintra,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   &(actsolv->rhs[1]),
                   fgradprhs,
                   -fdyn->thsr
        );

      solserv_sparsematvec(actintra,
                           &(actsolv->rhs[0]),
                           &(actsolv->sysarray[mass_array]),
                           &(actsolv->sysarray_typ[mass_array]),
                           &(actsolv->rhs[1]));

      /*--------------------------------------------------------- add rhs: ---*/
      /* The pressure part is already there. Add the rhs terms from
       * the element calls. */
      assemble_vec(actintra,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   &(actsolv->rhs[0]),
                   frhs,
                   1.0
        );

      /*-------------------------------------------------------- solve system */
#ifdef PERF
      perf_begin(80);
#endif

      init=0;
      t1=ds_cputime();
      solver_control(actfield,disnum_calc,actsolv, actintra,
                     &(actsolv->sysarray_typ[actsysarray]),
                     &(actsolv->sysarray[actsysarray]),
                     &(actsolv->sol[0]),
                     &(actsolv->rhs[0]),
                     init);
      ts=ds_cputime()-t1;
      tss+=ts;

#ifdef PERF
      perf_end(80);
#endif

      /*-- set flags for stability parameter evaluation and convergence check */
      fdyn->ishape=0;

      /*--- return solution to the nodes and calculate the convergence ratios */
      /* note: fdyn->numdf is always genprob.ndim+1. This is
       * ridiculous if we use discontinuous pressure because there are
       * only genprob.ndim dofs per node. But it allows us to use the
       * ordinary fluid service functions. */
      fluid_result_incre(actfield, 0,actintra,&(actsolv->sol[0]),
                         ipos->velnp,
                         &(actsolv->sysarray[actsysarray]),
                         &(actsolv->sysarray_typ[actsysarray]),
                         &vrat,NULL,NULL);

      /*----------------------------------------- iteration convergence check */
      converged = fluid_convcheck(vrat,0.,0.,itnum,te,ts);

      /*--------------------- check if nonlinear iteration has to be finished */
      if (converged==0)
      {
        itnum++;
      }
      else
      {
        break;
      }
    }

    /*----------------------------------------------------------------------*
     | -->  end of nonlinear iteration                                      |
     *----------------------------------------------------------------------*/

    /* Now we need to solve the pressure equation. */
    solserv_zero_vec(&(presolv.rhs[0]));
    parallel_sparse_copy_aztec(&pmat, &presolv);

    /* build up the rhs */
    pm_calprhs(actfield, actpart, disnum_calc, actintra, ipos, numpdof, &(presolv.rhs[0]));

    /* solve for the pressure increment */
    solver_control(actfield,disnum_calc,&presolv, actintra,
                   &(presolv.sysarray_typ[0]),
                   &(presolv.sysarray[0]),
                   &(presolv.sol[0]),
                   &(presolv.rhs[0]),
                   0);

    /* update pressure */
    pm_press_update(actfield, actpart, disnum_calc, actintra, ipos, numpdof, &(presolv.sol[0]), fdyn->dta);

    /* update velocity */
    pm_vel_update(actfield, actpart, disnum_calc, actintra, ipos, lmass, frhs, fgradprhs);

    /*---------- extrapolate from n+alpha_f to n+1 for generalised alpha ---*/
    if (fdyn->iop == 1)
    {
      solserv_sol_zero(actfield,0,node_array_sol_increment,ipos->accnm);
      solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,ipos->velnp,ipos->accnm,1.0/fdyn->alpha_f);
      solserv_sol_add(actfield,0,node_array_sol_increment,node_array_sol_increment,ipos->veln,ipos->accnm,1.0-1.0/fdyn->alpha_f);
      solserv_sol_copy(actfield,0,node_array_sol_increment,node_array_sol_increment,ipos->accnm,ipos->velnp);

      fdyn->acttime += fdyn->dta * (1.0 - fdyn->alpha_f);
    }
    /*-------------------------- check time step size in adaptive regime ---*/
    if(fdyn->adaptive)
    {
      if (fdyn->step > 1)
      {
        /*------------------------ evaluate local truncation error (LTE) ---*/
        fluid_lte(actfield, disnum_calc, ipos,fdyn->iop);

        /*------------------ evaluate norm of LTE and new time step size ---*/
        fluid_lte_norm(actpart, disnum_calc, actintra,ipos,
                       &iststep,&repeat,&repeated,itnum);
        if (repeat)
        {
          repeated++;
          continue;
        }
      }
      else
      {
        repeated = 0;
        fdyn->dt_prop = fdyn->dta;
      }
    }

    /*---------------------------------------------- update acceleration ---*/
    /*----- copy solution from sol_increment[5][j] to sol_increment[4][j] */
    /*--- -> prev. acceleration becomes (n-1)-accel. of next time step ---*/
    if (fdyn->step > 1)
      solserv_sol_copy(actfield,0,node_array_sol_increment,
                       node_array_sol_increment,
                       ipos->accn,ipos->accnm);
    /*------------------------ evaluate acceleration in this time step ---*
     *-------------------------------- depending on integration method ---*/
    if (fdyn->step == 1)
    { /* do just a linear interpolation within the first timestep */
      solserv_sol_zero(actfield,0,node_array_sol_increment,ipos->accnm);
      solserv_sol_add(actfield,0,node_array_sol_increment,
                      node_array_sol_increment,
                      ipos->velnp,ipos->accn, 1.0/fdyn->dta);
      solserv_sol_add(actfield,0,node_array_sol_increment,
                      node_array_sol_increment,
                      ipos->veln,ipos->accn,-1.0/fdyn->dta);
      solserv_sol_copy(actfield,0,node_array_sol_increment,
                       node_array_sol_increment,
                       ipos->accn,ipos->accnm);
    }
    else
      fluid_acceleration(actfield, disnum_calc, ipos,fdyn->iop);

    /*------------------------------------------- update time step sizes ---*/
    fdyn->dtp = fdyn->dta;
    if (fdyn->adaptive)
    {
      fdyn->dta = fdyn->dt_prop;
    }

    /*-------------------------------------------------- steady state check */
    if (fdyn->stchk==iststep)
    {
      iststep=0;
      steady = fluid_steadycheck(actfield, disnum_calc, ipos,numeq_total);
    }
    else
    {
      if (par.myrank==0)
        fprintf(out,"            |            |");
    }

    /*----------------------------------------------- lift&drag computation */
    if (fdyn->liftdrag==ld_stress || fdyn->liftdrag==ld_nodeforce)
    {
      container.str = str;
      *action = calc_fluid_liftdrag;
      fluid_liftdrag(1,action,&container,actfield, disnum_calc,
                     actsolv,actpart,actintra,ipos);
    }

    /* stress computation */
    if (ioflags.fluid_stress==1)
    {
      container.str = str;
      *action = calc_fluid_stress;
      calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
             &container,action);
    }

    /* error calculation for beltrami and kim-moin */
    if (fdyn->init==8 || fdyn->init==9)
    {
      container.vel_error  = 0.0;
      container.pre_error  = 0.0;
      container.error_norm = 2;

      *action = calc_fluid_error;
      calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
             &container,action);


      container.vel_error = sqrt(container.vel_error);
      container.pre_error = sqrt(container.pre_error);

      printf("\n  L2_err  :  velocity %f  pressure %f \n\n",
             container.vel_error,container.pre_error);
    }

    /*------- copy solution from sol_increment[1][j] to sol_increment[0][j] */
    /*------- -> prev. solution becomes (n-1)-solution of next time step ---*/
    solserv_sol_copy(actfield,0,node_array_sol_increment,
                     node_array_sol_increment,ipos->veln,ipos->velnm);

    /* copy solution from sol_increment[velnp][j] to sol_increment[veln][j] */
    /*--- -> actual solution becomes previous solution of next time step ---*/
    solserv_sol_copy(actfield,0,node_array_sol_increment,
                     node_array_sol_increment,ipos->velnp,ipos->veln);

    /*-------- copy solution from sol_increment[ipos->velnp][j] to sol_[actpos][j] */
    solserv_sol_copy(actfield,0,node_array_sol_increment,
                     node_array_sol,ipos->veln,actpos);

    /*---------------------------------------------- finalise this timestep */
    outstep++;
    pssstep++;
    resstep++;
    restartstep++;

    /*------------------------------------------- write restart to pss file */
    if (restartstep==fdyn->uprestart)
    {
      restartstep=0;
#ifdef BINIO
      if (disnum_io != disnum_calc)
        restart_write_bin_fluiddyn(&restart_context,fdyn);
      else
        restart_write_bin_fluiddyn(&out_context,fdyn);
#else
      restart_write_fluiddyn(fdyn,actfield,actpart,actintra,action,&container);
#endif
    }

#ifdef SUBDIV
    /* transfer the solution to the nodes of the master-dis */
    if (actfield->subdivide > 0)
    {
      solserv_sol_trans(actfield, disnum_calc, node_array_sol, actpos);
    }
#endif

    /*--------------------------------------- write solution to .flavia.res */
    if (resstep==fdyn->upres && par.myrank==0 && ioflags.output_gid==1)
    {
      if(ioflags.fluid_sol==1)
      {
        out_gid_sol("velocity",actfield,disnum_io,actintra,fdyn->step,actpos,fdyn->acttime);
        out_gid_sol("average_pressure",actfield,disnum_io,actintra,fdyn->step,actpos,fdyn->acttime);
      }
      if(ioflags.fluid_stress==1)
      {
        out_gid_sol("stress",actfield,disnum_io,actintra,fdyn->step,actpos,fdyn->acttime);
      }
    }

    /*--------------------------------------- write solution to binary */
#ifdef BINIO
    if (resstep==fdyn->upres && ioflags.output_bin==1)
    {
      if (ioflags.fluid_sol==1)
      {
        out_results(&out_context, fdyn->acttime, fdyn->step, actpos, OUTPUT_VELOCITY | OUTPUT_AV_PRESSURE);
      }
      if (ioflags.fluid_stress==1)
      {
        out_results(&out_context, fdyn->acttime, fdyn->step, actpos, OUTPUT_STRESS);
      }
    }
#endif

    if (resstep==fdyn->upres)
    {
      resstep=0;
    }

    /*---------------------------------------------- write solution to .out */
    if (outstep==fdyn->upout && ioflags.output_out==1 && ioflags.fluid_sol==1)
    {
      outstep=0;
      out_sol(actfield,actpart, disnum_io, actintra,fdyn->step,actpos);
    }

    /*---------------------------------------------------------- monitoring */
    if (ioflags.monitor==1)
      monitoring(actfield, disnum_io, genprob.numff,actpos,fdyn->acttime);

    tt=ds_cputime()-t2;
    tts+=tt;
    printf("PROC  %3d | total time for this time step: %10.3e \n",par.myrank,tt);


    if (par.myrank==0)
    {
      fprintf(out,"            |            |");
      fprintf(out," %10.3E |\n",tt);
      fflush(out);
    }

    dsmemreport();

    /*--------------------- check time and number of steps and steady state */
    if (!(fdyn->step < fdyn->nstep && fdyn->acttime <= fdyn->maxtime && steady==0))
    {
      break;
    }
  }

  /*----------------------------------------------------------------------*
   | -->  end of timeloop                                                 |
   *----------------------------------------------------------------------*/

  /*======================================================================*
   |                      F I N A L I S I N G                             |
   *======================================================================*/
  if (pssstep==0) actpos--;
  /*------------------------------------- print out solution to .out file */
  if (outstep!=0 && ioflags.output_out==1 && ioflags.fluid_sol==1)
    out_sol(actfield,actpart, disnum_io, actintra,fdyn->step,actpos);

  /*---------------------------------- print total CPU-time to the screen */
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
      fprintf(out,"\n\n\n");
      fprintf(out," total time element for calculations: %10.3E \n", tes);
      fprintf(out," total time for solver              : %10.3E \n", tss);
      fprintf(out," total time for time loop           : %10.3E \n", tts);
    }
  }
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
end:

#ifdef BINIO
  destroy_bin_out_field(&out_context);
  if (disnum_io != disnum_calc)
    destroy_bin_out_field(&restart_context);
#endif

  /*--------------------------------------------------- cleaning up phase */
  parallel_sparse_destroy(&grad);
  parallel_sparse_destroy(&pmat);

  CCAFREE(lmass);

  amdel(&frhs_a);
  amdel(&fgradprhs_a);
  if (ioflags.fluid_vis==1 )
    solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
  solserv_del_vec(&(actsolv->sol),actsolv->nsol);

#ifndef PARALLEL
  CCAFREE(actintra);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


#endif
/*! @} (documentation module close)*/
