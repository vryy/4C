/*!
\file
\brief projection method algorithm for fluid

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

This is the projection method as described in GRESHO, Philip M. "On
the theory of Semi-Implicit Projection Methods for viscous
incompressible flow and its implementation via a finite element method
that also introduces a nearly consistent mass matrix part 1 & part 2"
International Journal for Numerical Methods in Fluids, VOL 11,
pp. 587-659, (1990).

Here we have a continous pressure discretication.

*/
/*!
\addtogroup FLUID_PM
*//*! @{ (documentation module open)*/

#ifdef D_FLUID_PM

#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../solver/solver_sparse.h"
#include "../solver/solver_trilinos_service.H"
#include "fluid_prototypes.h"
#include "fluid_pm_prototypes.h"
#include "../io/io.h"

#ifdef D_FLUID2_PRO
#include "../fluid2_pro/fluid2pro.h"
#include "../fluid2_pro/fluid2pro_prototypes.h"
#endif

#ifdef D_FLUID3_PRO
#include "../fluid3_pro/fluid3pro.h"
#include "../fluid3_pro/fluid3pro_prototypes.h"
#endif

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


  \author u.kue
  \date 10/06
 */
/*----------------------------------------------------------------------*/
void fluid_pm_cont()
{
#if defined(TRILINOS_PACKAGE) && defined(PM_TRILINOS)
  INT             itnum;        /* counter for nonlinear iteration  */
  INT             i;            /* simply a counter                 */
  INT             numeq;        /* number of equations on this proc */
  INT             numeq_total;  /* total number of equations        */

  INT             pnumeq;
  INT             pnumeq_total;

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
  SOLVAR         *pressolv;
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
  INT             press_dis;

  TRILINOSMATRIX grad;
  TRILINOSMATRIX lmass;

  DIST_VECTOR   *press_rhs;
  DIST_VECTOR   *press_sol;

  INT stiff_array;              /* indice of the active system sparse matrix */
  INT mass_array;               /* indice of the active system sparse matrix */
  INT press_array;

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
  pressolv    = &(solv[1]);
  actpart     = &(partition[0]);
  action      = &(calc_action[0]);
  restart     = genprob.restart;

#ifdef SUBDIV
  if (actfield->subdivide > 0)
  {
    dserror("not supported");
    disnum_calc = 1;
    disnum_io   = ioflags.output_dis;
  }
  else
#endif
    disnum_calc = disnum_io = 0;
  press_dis = 1;

  /*
   * This is again highly dangerous! We need an additional mass matrix
   * and assume we can enlarge the solver's matrix array by one to put
   * it there.
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
  press_array = stiff_array+1;
  mass_array  = press_array+1;

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

  /* change discretization type in main discretization */
  for (i=0; i<actfield->dis[disnum_calc].numele; ++i)
  {
    ELEMENT* actele;
    actele = &actfield->dis[disnum_calc].element[i];
    switch (actele->eltyp)
    {
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:
    {
      FLUID2_PRO* f2pro;
      f2pro = actele->e.f2pro;
      switch (f2pro->dm)
      {
      case dm_q2pm1:
	/* switch to taylor-hood */
	f2pro->dm = dm_q2q1;
	f2pro->other->e.f2pro->dm = dm_q2q1;
	break;
      case dm_q1p0:
	f2pro->dm = dm_q1q1;
	f2pro->other->e.f2pro->dm = dm_q1q1;
	break;
      default:
        dserror("discretization mode %d currently unsupported", f2pro->dm);
      }
      break;
    }
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:
    {
      FLUID3_PRO* f3pro;
      f3pro = actele->e.f3pro;
      switch (f3pro->dm)
      {
      case dm_q2pm1:
	/* switch to "taylor-hood" */
	f3pro->dm = dm_q2q1;
	f3pro->other->e.f3pro->dm = dm_q2q1;
	break;
      case dm_q1p0:
	f3pro->dm = dm_q1q1;
	f3pro->other->e.f3pro->dm = dm_q1q1;
	break;
      default:
        dserror("discretization mode %d currently unsupported", f3pro->dm);
      }
      break;
    }
#endif
    default:
      dserror("element type %d undefined",actele->eltyp);
    }
  }

  /* ------------------------------------------------ */
  /* The gradient matrix G */
  /* This is a parallel matrix, each processor owns a separate slice. */

  /* the continuous pressure version */

  if (!genprob.usetrilinosalgebra)
  {
    dserror("trilinos algebra required");
  }

  {
    TRILINOSMATRIX* global_stiff;
    TRILINOSMATRIX* global_press;

    /* we take the update map from our global solvers */
    if (actsolv->sysarray_typ[stiff_array]!=trilinos)
      dserror("fatal: trilinos solver expected");
    if (actsolv->sysarray_typ[press_array]!=trilinos)
      dserror("fatal: trilinos solver expected");

    global_stiff = actsolv->sysarray[stiff_array].trilinos;
    global_press = actsolv->sysarray[press_array].trilinos;

    /* Let the distinguished pressure solver point to the pressure's
     * matrix. This way we can use the pressure solver's flags for
     * solving the pressure equation. Yeah! What a hack! */

    pressolv->sysarray_typ = (SPARSE_TYP*)  CCACALLOC(1,sizeof(SPARSE_TYP));
    pressolv->sysarray     = (SPARSE_ARRAY*)CCACALLOC(1,sizeof(SPARSE_ARRAY));
    pressolv->sysarray_typ[0] = trilinos;
    pressolv->sysarray[0].trilinos = global_press;

    /* setup diagonal mass matrix */
    memset(&lmass,0,sizeof(TRILINOSMATRIX));
    lmass.numeq_total = global_stiff->numeq_total;
    lmass.numeq       = global_stiff->numeq;

    am_alloc_copy(&global_stiff->update,&lmass.update);
    construct_trilinos_diagonal_matrix(actintra,&lmass);

    /* setup rectangular gradient matrix (FillComplete with arguments
     * needed.) */
    memset(&grad,0,sizeof(TRILINOSMATRIX));
    grad.numeq_total = global_press->numeq_total;
    grad.numeq       = global_press->numeq;

    am_alloc_copy(&global_press->update,&grad.update);
    construct_trilinos_matrix(actintra,&grad);
  }

  /*---------------------------- get global and local number of equations */
  solserv_getmatdims(&(actsolv->sysarray[press_array]),
                     actsolv->sysarray_typ[press_array],
                     &pnumeq,
                     &pnumeq_total);

  /*---------------------------------------- allocate dist. vectors 'rhs' */
  solserv_create_vec(&press_rhs,1,pnumeq_total,pnumeq,"DV");
  solserv_zero_vec(press_rhs);

  /*---------------------------------- allocate dist. solution vectors ---*/
  solserv_create_vec(&press_sol,1,pnumeq_total,pnumeq,"DV");
  solserv_zero_vec(press_sol);

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

  /* we need two entries on the pressure discretization */
  /*
   * sol_increment[0]  ... phi
   * sol_increment[1]  ... press
   */
  solserv_sol_zero(actfield,press_dis,node_array_sol_increment,1);

  solserv_sol_zero(actfield,press_dis,node_array_sol,0);

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
  /* initialize pressure matrix */
  /* We initialize the matrix in actsolv using the parameters from
   * pressolv. */
  solver_control(actfield,disnum_calc,pressolv, actintra,
                 &(actsolv->sysarray_typ[press_array]),
                 &(actsolv->sysarray[press_array]),
                 NULL,
                 NULL,
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
      out_gid_sol("projected_pressure",actfield,press_dis,actintra,fdyn->step,actpos,fdyn->acttime);
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
      out_results(&out_context, fdyn->acttime, fdyn->step, actpos, OUTPUT_VELOCITY);
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

  pm_calelm_cont(actfield, actpart, disnum_calc, press_dis,
		 actsolv, mass_array,
		 actintra, ipos, &grad, &lmass);

  {
    TRILINOSMATRIX* global_press;

    if (actsolv->sysarray_typ[press_array]!=trilinos)
      dserror("fatal: trilinos solver expected");

    global_press = actsolv->sysarray[press_array].trilinos;
    mult_trilinos_mmm_cont(global_press,&grad,0,&lmass,0,&grad,1);
  }

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

/* #define PM_ITERATION */
#ifdef PM_ITERATION
    {
      INT loop;
    for (loop=0; loop<10; ++loop)
    {
#else
#endif
    
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

#ifdef QUASI_NEWTON
      fdyn->itnum = itnum;
#endif

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

      assemble_vec(actintra,
                   &(actsolv->sysarray_typ[actsysarray]),
                   &(actsolv->sysarray[actsysarray]),
                   &(actsolv->rhs[0]),
                   fgradprhs,
                   -fdyn->thsl
        );

      matvec_trilinos(&(actsolv->rhs[1]),
		      &(actsolv->rhs[0]),
		      &lmass);

      solserv_zero_vec(&(actsolv->rhs[0]));

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
      fluid_result_incre(actfield, 0,actintra,&(actsolv->sol[0]),&(actsolv->rhs[0]),
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

#if 1
    if (actintra->intra_rank==0)
    {
      fprintf(allfiles.gidres,"RESULT \"vel_star\" \"ccarat\" %d VECTOR ONNODES\n"
              "RESULTRANGESTABLE \"standard_fluid    \"\n"
              "COMPONENTNAMES \"e1\" \"e2\"\n"
              "VALUES\n",fdyn->step);

      for (i=0; i<actfield->dis[0].numnp; ++i)
      {
        NODE* n = &actfield->dis[0].node[i];
        fprintf(allfiles.gidres,"%d %e %e\n",n->Id+1,
                n->sol_increment.a.da[ipos->velnp][0],
                n->sol_increment.a.da[ipos->velnp][1]);
      }

      fprintf(allfiles.gidres,"END VALUES\n");
    }
#endif
    
    solserv_zero_vec(press_rhs);
    amzero(&frhs_a);
    amzero(&fgradprhs_a);

    /* build up the rhs */
    pm_calprhs_cont(actfield, actpart, disnum_calc, actintra, ipos, press_rhs, fgradprhs, frhs);

    assemble_vec(actintra,
		 &(actsolv->sysarray_typ[press_array]),
		 &(actsolv->sysarray[press_array]),
		 press_rhs,
		 frhs,
		 1.0
      );

    /* solve for the pressure increment */
    solver_control(actfield,press_dis,pressolv,actintra,
                   &(actsolv->sysarray_typ[press_array]),
                   &(actsolv->sysarray[press_array]),
                   press_sol,
                   press_rhs,
                   0);

    /* update pressure */
    /* solserv_result_incre does not work due to misleading dirichlet
       conditions. The conditions on the velocity interfere. */
    {
      INT i;
      amzero(&frhs_a);
      solserv_reddistvec(press_sol,
			 &(actsolv->sysarray[press_array]),
			 &(actsolv->sysarray_typ[press_array]),
			 frhs,
			 press_sol->numeq_total,
			 actintra);
      for (i=0; i<actfield->dis[press_dis].numnp; i++)
      {
	NODE* actnode;
	actnode = &(actfield->dis[press_dis].node[i]);
	actnode->sol_increment.a.da[0][0]  = frhs[actnode->dof[0]];
	/* actnode->sol_increment.a.da[1][0] += 2./fdyn->dta*frhs[actnode->dof[0]]; */
	/*actnode->sol_increment.a.da[1][0] += frhs[actnode->dof[0]];*/
        actnode->sol_increment.a.da[1][0] += frhs[actnode->dof[0]] / fdyn->thsl;
      }
    }

#ifdef PM_ITERATION
    }}
#else
    
    /* update velocity */
    pm_vel_update(actfield, actpart, disnum_calc, actintra, ipos,
		  &lmass, actsolv, actsysarray,
		  frhs, fgradprhs);

#endif
    
#if 1
    if (actintra->intra_rank==0)
    {
      fprintf(allfiles.gidres,"RESULT \"press_rhs\" \"ccarat\" %d SCALAR ONNODES\n"
              "RESULTRANGESTABLE \"standard_fluid    \"\n"
              "COMPONENTNAMES \"pressure\"\n"
              "VALUES\n",fdyn->step);

      for (i=0; i<actfield->dis[press_dis].numnp; ++i)
      {
        NODE* n = &actfield->dis[press_dis].node[i];
        fprintf(allfiles.gidres,"%d %f\n",
                n->Id+1-genprob.nodeshift,
                /* frhs[n->dof[0]]); */
                press_rhs->vec.a.dv[n->dof[0]]);
      }

      fprintf(allfiles.gidres,"END VALUES\n");

      fprintf(allfiles.gidres,"RESULT \"press_sol\" \"ccarat\" %d SCALAR ONNODES\n"
              "RESULTRANGESTABLE \"standard_fluid    \"\n"
              "COMPONENTNAMES \"pressure\"\n"
              "VALUES\n",fdyn->step);

      for (i=0; i<actfield->dis[press_dis].numnp; ++i)
      {
        NODE* n = &actfield->dis[press_dis].node[i];
        fprintf(allfiles.gidres,"%d %f\n",
                n->Id+1-genprob.nodeshift,
                /* frhs[n->dof[0]]); */
                press_sol->vec.a.dv[n->dof[0]]);
      }

      fprintf(allfiles.gidres,"END VALUES\n");
    }
#endif

#if 0
    if (actintra->intra_rank==0)
    {
      fprintf(allfiles.gidres,"RESULT \"elementcall\" \"ccarat\" %d VECTOR ONNODES\n"
              "RESULTRANGESTABLE \"standard_fluid    \"\n"
              "COMPONENTNAMES \"e1\" \"e2\"\n"
              "VALUES\n",fdyn->step);

      for (i=0; i<actfield->dis[0].numnp; ++i)
      {
        DOUBLE d1,d2;
        NODE* n = &actfield->dis[0].node[i];

        d1 = frhs[n->dof[0]];
        d2 = frhs[n->dof[1]];

        if (n->gnode->dirich && n->gnode->dirich->dirich_onoff.a.iv[0])
          d1 = 0;
        if (n->gnode->dirich && n->gnode->dirich->dirich_onoff.a.iv[1])
          d2 = 0;

        fprintf(allfiles.gidres,"%d %e %e\n",n->Id+1,d1,d2);
      }

      fprintf(allfiles.gidres,"END VALUES\n");


      fprintf(allfiles.gidres,"RESULT \"elementcall2\" \"ccarat\" %d VECTOR ONNODES\n"
              "RESULTRANGESTABLE \"standard_fluid    \"\n"
              "COMPONENTNAMES \"e1\" \"e2\"\n"
              "VALUES\n",fdyn->step);

      for (i=0; i<actfield->dis[0].numnp; ++i)
      {
        DOUBLE d1,d2;
        NODE* n = &actfield->dis[0].node[i];

        d1 = fgradprhs[n->dof[0]];
        d2 = fgradprhs[n->dof[1]];

        if (n->gnode->dirich && n->gnode->dirich->dirich_onoff.a.iv[0])
          d1 = 0;
        if (n->gnode->dirich && n->gnode->dirich->dirich_onoff.a.iv[1])
          d2 = 0;

        fprintf(allfiles.gidres,"%d %e %e\n",n->Id+1,d1,d2);
      }

      fprintf(allfiles.gidres,"END VALUES\n");
    }
#endif

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

    solserv_sol_copy(actfield,press_dis,
		     node_array_sol_increment,
                     node_array_sol,
		     1,
		     actpos);

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
	out_gid_sol("projected_pressure",actfield,press_dis,actintra,fdyn->step,actpos,fdyn->acttime);
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
        out_results(&out_context, fdyn->acttime, fdyn->step, actpos, OUTPUT_VELOCITY);
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
#else
  dserror("TRILINOS_PACKAGE and PM_TRILINOS required");
#endif
}


#endif
/*! @} (documentation module close)*/
