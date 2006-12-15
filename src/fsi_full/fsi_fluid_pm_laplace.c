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
#ifdef D_FLUID_PM


#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../fluid_full/fluid_prototypes.h"
#include "../fluid_full/fluid_pm_prototypes.h"
#include "fsi_prototypes.h"
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


/*----------------------------------------------------------------------*/
/*!
  \brief setup fsi fluid algorithm
 */
/*----------------------------------------------------------------------*/
void fsi_fluid_pm_laplace_setup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  )
{
#if defined(TRILINOS_PACKAGE) && defined(PM_TRILINOS)
  INT             numff;              /* actual number of fluid field     */
  INT             i;                  /* counters                         */
  INT             numeq;              /* number of equations on this proc */
  INT             numeq_total;        /* total number of equations        */
  INT             pnumeq;
  INT             pnumeq_total;
  INT             init;               /* flag for solver_control call     */
  INT             actsysarray=0;      /* number of actual sysarray        */
  INT             restart;

  DOUBLE          grat;               /* convergence ratios               */
  SOLVAR         *actsolv;            /* pointer to active sol. structure */
  SOLVAR         *pressolv;
  PARTITION      *actpart;            /* pointer to active partition      */
  INTRA          *actintra;           /* pointer to active intra-communic.*/
  CALC_ACTION    *action;             /* pointer to the cal_action enum   */

  DOUBLE         *frhs;               /* iteration - RHS                  */
  DOUBLE         *fgradprhs;

#ifdef PARALLEL
  INT             numddof;            /* number of Dirichlet dofs         */
  DOUBLE         *fcouple;            /* to store fsi coupling forces     */
  DOUBLE         *recvfcouple;        /* to communicate fsi coupling forces*/
#endif

  DOUBLE         *totarea;
  CONTAINER       container;          /* variables for calelm             */
  FLUID_STRESS    str=str_none;
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  FSI_DYNAMIC    *fsidyn;             /* fsi dynamic variables     */
  ARRAY_POSITION *ipos;

  INT stiff_array;              /* indice of the active system sparse matrix */
  INT mass_array;               /* indice of the active system sparse matrix */
  INT press_array;
  INT             press_dis;

  perf_begin(46);

  /* initialiase some counters */
  work->outstep=0;
  work->restartstep=0;

#ifndef PARALLEL
  work->dummy_intra.intra_fieldtyp = fluid;
  work->dummy_intra.intra_rank     = 0;
  work->dummy_intra.intra_nprocs   = 1;
#endif

  /****************************************/

  press_dis = 1;
  
  numff         = genprob.numff;
  fdyn          = alldyn[numff].fdyn;
  fsidyn        = alldyn[genprob.numaf+1].fsidyn;

  actsysarray   = disnum_calc;
  stiff_array = disnum_calc;
  press_array = stiff_array+1;
  mass_array  = press_array+1;

  fdyn->dt      = fsidyn->dt;
  fdyn->maxtime = fsidyn->maxtime;
  fdyn->nstep   = fsidyn->nstep;
  grat          = ZERO;

  /* set some pointers */
  /* -> only valid for single field problem !!!! */
  actsolv     = &(solv[numff]);
  pressolv    = &(solv[numff+3]);
  actpart     = &(partition[numff]);
  action      = &(calc_action[numff]);
  restart     = genprob.restart;

  memset(&container,0,sizeof(CONTAINER));
  container.fieldtyp = actfield->fieldtyp;
  container.disnum   = disnum_calc;
  container.turbu    = fdyn->turbu;

  if (genprob.probtyp == prb_fsi)
    str = str_fsicoupling;

  fdyn->acttime=ZERO;

  if (fdyn->freesurf == 5)
    fdyn->hf_stab=0;

  ipos = &(actfield->dis[disnum_calc].ipos);

  /* if we are not parallel, we have to allocate an alibi *
   * intra-communicator structure */
#ifdef PARALLEL
  actintra    = &(par.intra[numff]);
#else
  actintra    = &work->dummy_intra;
#endif

  if (actintra->intra_fieldtyp != fluid)
    dserror("only fluid allowed");

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
  }
  
  /*---------------------------- get global and local number of equations */
  solserv_getmatdims(&(actsolv->sysarray[press_array]),
                     actsolv->sysarray_typ[press_array],
                     &pnumeq,
                     &pnumeq_total);

  /*---------------------------------------- allocate dist. vectors 'rhs' */
  solserv_create_vec(&work->press_rhs,1,pnumeq_total,pnumeq,"DV");
  solserv_zero_vec(work->press_rhs);

  /*---------------------------------- allocate dist. solution vectors ---*/
  solserv_create_vec(&work->press_sol,1,pnumeq_total,pnumeq,"DV");
  solserv_zero_vec(work->press_sol);
  
  /****************************************/


  /* init the dist sparse matrices to zero */
  solserv_zero_mat(
    actintra,
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->sysarray_typ[actsysarray])
    );


  /* get global and local number of equations */
  solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
		     actsolv->sysarray_typ[actsysarray],
		     &numeq,
		     &numeq_total);


  /* allocate 1 dist. vector 'rhs' */
  actsolv->nrhs = 2;
  solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
  solserv_zero_vec(&(actsolv->rhs[0]));


  /* allocate dist. solution vectors */
  actsolv->nsol = 1;

  solserv_create_vec(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");
  for (i=0; i<actsolv->nsol; i++)
    solserv_zero_vec(&(actsolv->sol[i]));


  /* allocate one redundant vector frhs of full lenght */
  /* -> this is used by the element routines to assemble the Iteration RHS*/
  frhs = amdef("frhs",&work->frhs_a,numeq_total,1,"DV");
  fgradprhs = amdef("fgradprhs",&work->fgradprhs_a,numeq_total,1,"DV");


  /* allocate one vector for storing the area */
  if (fdyn->checkarea > 0)
  {
    totarea = amdef("area",&work->totarea_a,fdyn->itemax,1,"DV");
    amzero(&work->totarea_a);
  }


  /* init lift&drag calculation real FSI-problem */
  if (fdyn->liftdrag == ld_stress)
  {
    str = str_liftdrag;
    if(genprob.numfld == 3) fluid_liftdrag(0,NULL,NULL,NULL,0,NULL,NULL,NULL,NULL);
  }
  if (fdyn->liftdrag == ld_nodeforce)
    fluid_liftdrag(-1,action,&container,actfield,disnum_calc,
		   actsolv,actpart,actintra,ipos);



  /* init fsi coupling via consistent nodal forces */
  if (fsidyn->coupforce == cf_nodeforce)
  {
    /* Scheme for parallelisation of fsi coupling forces:
       a) Solution for all dofs (including Dirichlet boundary conditions):
       - allocate vector of full length on every proc
       - collect nodal force values on the respective proc
       - allreduce the full vector
       - distribute the values to the nodal field sol_mf[1]
       b) Solution without Dbcs
       - allocate vector of length of number of Dirichlet dofs
       - collect nodal force vectors and allreduce this (small) vector
       - distribute the values to the nodal field sol_mf[1]
    */
#ifdef PARALLEL

#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
    /*-----------  allocate vectors for coupling force of full lenght ---*/
    fcouple = amdef("fcouple",&work->fcouple_a,numeq_total,1,"DV");
    recvfcouple = amdef("recvfcouple",&work->recvfcouple_a,numeq_total,1,"DV");
    numddof = 0;
#else
    numddof = actfield->dis[disnum_calc].numdf-numeq_total;
    fcouple = amdef("fcouple",&work->fcouple_a,numddof,1,"DV");
    recvfcouple = amdef("recvfcouple",&work->recvfcouple_a,numddof,1,"DV");
#endif /* SOLVE_DIRICH */

    fsi_cbf(&(actpart->pdis[disnum_calc]),fcouple,ipos,numeq_total,1);
#else
    fsi_cbf(&(actpart->pdis[disnum_calc]),NULL,ipos,0,1);
    solserv_sol_zero(actfield,disnum_calc,node_array_sol_mf,1);
#endif /* PARALLEL */
  }



  /* initialise fluid field */
  if (restart > 0)
  {
    if (fdyn->init>0)
      dserror("Initial field either by restart, or by function or from file ...\n");
    else
    {
      fdyn->resstep=genprob.restart;
      fdyn->init=2;
    }
  }

  fluid_init_pos_ale(actfield,disnum_calc);

  if(fdyn->iop == 4) ipos->numincr = 9;
  else               ipos->numincr = 7;
  if ((fsidyn->ifsi == fsi_iter_stagg_steep_desc) ||
      (fsidyn->ifsi == fsi_iter_stagg_steep_desc_force))
    /* ... steepest descent relaxation ... */
    ipos->numincr++;            /* ... each need one more solution field. */

  ipos->nummf = 2;
  ipos->mf_velnp = 0;
  ipos->mf_forcenp = 1;

  if (fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_force)
  {
    ipos->nummf += 1;
    ipos->mf_forcen = 2;
    solserv_sol_zero(actfield, disnum_calc, node_array_sol_mf, ipos->mf_forcen);
  }
  else if (fsidyn->ifsi==fsi_iter_stagg_steep_desc_force)
  {
    ipos->nummf += 4;
    ipos->mf_forcen = 2;
    ipos->mf_forcecpy = 3;
    ipos->mf_sd_g = 4;
    ipos->mf_velcpy = 5;
    solserv_sol_zero(actfield, disnum_calc, node_array_sol_mf, ipos->mf_velcpy);
  }

  fluid_init(actpart,actintra,actfield,disnum_calc,disnum_io,action,
	     &container,ipos->numincr,ipos,str);

  /* we need two entries on the pressure discretization */
  /*
   * sol_increment[0]  ... phi
   * sol_increment[1]  ... press
   */
  solserv_sol_zero(actfield,press_dis,node_array_sol_increment,1);

  solserv_sol_zero(actfield,press_dis,node_array_sol,0);
  
  /* we need two entries on the pressure discretization */
  /*
   * sol_increment[0]  ... phi
   * sol_increment[1]  ... press
   */
  solserv_sol_zero(actfield,press_dis,node_array_sol_increment,1);

  solserv_sol_zero(actfield,press_dis,node_array_sol,0);

  /* init the dirichlet-conditions */
  fluid_initdirich(actfield, disnum_calc, ipos);



  /* initialize solver on all matrices:
   * ----------------------------------
   * NOTE: solver init phase has to be called with each matrix one wants to
   * solve with. Solver init phase has to be called with all matrices
   * one wants to do matrix-vector products and matrix scalar products.
   * This is not needed by all solver libraries, but the solver-init
   * phase is cheap in computation (can be costly in memory)
   */


  /* initialize solver */
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



  /* init the assembly for stiffness */
  init_assembly(actpart,actsolv,actintra,actfield,actsysarray,disnum_calc);


  /* allocate fluid integration data */
  alldyn[numff].fdyn->data = (FLUID_DATA*)CCACALLOC(1,sizeof(FLUID_DATA));


  /* init the element calculating routines */
  *action = calc_fluid_init;
  calinit(actfield,actpart,action,&container);


  /* initialise energey check */
  if (fsidyn->ichecke>0) fsi_dyneint(NULL,disnum_calc,2);


  /* output to the screen */
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



  /* initialise height function solution */
  if (fdyn->freesurf == 3)
    fluid_heightfunc(1,&grat,actfield,actpart,actintra,action,
		     &container,ipos);


  /* calculate curvature at the beginning */
  if (fdyn->surftens != 0)
  {
    fluid_tcons();
    *action = calc_fluid_curvature;
    fluid_curvature(actfield,actpart,actintra,action);
  }


  /* calculate nodal normals */
  fluid_cal_normal(actfield,disnum_calc,1,action);


  /* define local co-sys */
  fluid_locsys(actfield,disnum_calc,fdyn);


  /* predictor for free surface at the beginning */
  if (fdyn->freesurf > 0)
    fluid_updfscoor(actfield, fdyn, fdyn->dt, ipos, -1);


  /* monitoring */
  if (ioflags.monitor == 1)
  {
    out_monitor(actfield,numff,ZERO,1);
    monitoring(actfield,disnum_calc,numff,0,fdyn->acttime);
  }


  /* init area monitoring */
  if (fdyn->checkarea > 0) out_area(work->totarea_a,fdyn->acttime,0,1);


  /* print out initial data to .out */
  if (ioflags.output_out == 1  &&  ioflags.fluid_sol == 1)
    out_sol(actfield,actpart,disnum_io,actintra,fdyn->step,0);


  /* calculate time independent constants for time algorithm */
  fluid_cons();


#ifdef D_MORTAR
  /* redefine the size of sol_mf from 2 to 3, the third field is necessary*/
  /* to store the nodal forces due to fsi */
  solserv_sol_zero(actfield, disnum_calc, node_array_sol_mf, 3);
#endif


#ifdef BINIO
  /* initialize binary output
   * It's important to do this only after all the node arrays are set
   * up because their sizes are used to allocate internal memory. */
  init_bin_out_field(&work->out_context,
		     &(actsolv->sysarray_typ[actsysarray]),
		     &(actsolv->sysarray[actsysarray]),
		     actfield, actpart, actintra, disnum_io);

  if (disnum_io != disnum_calc)
    init_bin_out_field(&work->restart_context,
                       &(actsolv->sysarray_typ[actsysarray]),
                       &(actsolv->sysarray[actsysarray]),
                       actfield, actpart, actintra, disnum_calc);
#endif


#ifdef PERF
  perf_end(46);
#endif

#else
  dserror("TRILINOS_PACKAGE and PM_TRILINOS required");
#endif
}



      /*======================================================================*
        |                     S O L U T I O N    P H A S E                   |
       *======================================================================*/
      /* nodal solution history fluid field:                                  *
       * sol[0][j]           ... initial data                                 *
       * sol[1...0][j]  ... solution for visualisation (real pressure)   *
       * sol_increment[flag][j] ... solution value needed further             *
       * sol_mf[0][j]        ... solution at time (n+1)                       *
       * sol_mf[1][j]        ... nodal stresses at FS-interface at time (n+1) *
       * in mortar cases only:                                                *
       * sol_mf[2][j]        ... nodal forces at FS-interface at time (n+1)   *
       * sol_increment flags:                                                 *
       *  velnm  ...  nodal solution at time (n-1)                            *
       *  veln   ...  nodal solution at time (n)                              *
       *  velnp  ...  nodal solution at time (n+1)                            *
       *  accnm  ...  nodal acceleration at time (n-1)                        *
       *  accn   ...  nodal acceleration at time (n)                          *
       *  hist   ...  linear combination of history values needed for rhs     *
       *  gridv  ...  nodal grid velocity within actual time step             *
       *  convn  ...  nodal convective velocity at time (n)                   *
       *  convnp ...  nodal convective velocity at time (n+1)                 *
       *======================================================================*/

/*----------------------------------------------------------------------*/
/*!
  \brief fsi fluid nonlinear loop

  Calculate on time step but do not update. We might want to do it
  again.

 */
/*----------------------------------------------------------------------*/
void fsi_fluid_pm_laplace_calc(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io,
  FIELD          *alefield,
  INT             adisnum_calc
  )
{
#if defined(TRILINOS_PACKAGE) && defined(PM_TRILINOS)
  INT             itnum;              /* counter for NR-Iterations        */
  INT             numff;              /* actual number of fluid field     */
  INT             numeq;              /* number of equations on this proc */
  INT             numeq_total;        /* total number of equations        */
  INT             init;               /* flag for solver_control call     */
  INT             actsysarray=0;      /* number of actual sysarray        */
  INT             nfrastep;           /* number of steps for fractional-
					 step-theta procedure             */
  INT             restart;
  INT             converged;          /* convergence flag                 */

  DOUBLE          vrat,prat;
  DOUBLE          grat;               /* convergence ratios               */
  DOUBLE          t1,ts,te;
  SOLVAR         *actsolv;            /* pointer to active sol. structure */
  SOLVAR         *pressolv;
  PARTITION      *actpart;            /* pointer to active partition      */
  INTRA          *actintra;           /* pointer to active intra-communic.*/
  CALC_ACTION    *action;             /* pointer to the cal_action enum   */

  DOUBLE         *frhs;               /* iteration - RHS                  */
  DOUBLE         *fgradprhs;

#ifdef PARALLEL
  INT             numddof;            /* number of Dirichlet dofs         */
  DOUBLE         *fcouple;            /* to store fsi coupling forces     */
  DOUBLE         *recvfcouple;        /* to communicate fsi coupling forces*/
#endif

  DOUBLE         *totarea;
  CONTAINER       container;          /* variables for calelm             */
  FLUID_STRESS    str=str_none;
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  FSI_DYNAMIC    *fsidyn;             /* fsi dynamic variables     */
  ARRAY_POSITION *ipos;

  INT stiff_array;              /* indice of the active system sparse matrix */
  INT mass_array;               /* indice of the active system sparse matrix */
  INT press_array;
  INT             press_dis;
  
#ifdef PERF
  perf_begin(47);
#endif

  /****************************************/

  frhs = work->frhs_a.a.dv;
  fgradprhs = work->fgradprhs_a.a.dv;
  
#ifdef PARALLEL
  fcouple = work->fcouple_a.a.dv;
  recvfcouple = work->recvfcouple_a.a.dv;
#endif

  totarea = work->totarea_a.a.dv;

  numff         = genprob.numff;
  fdyn          = alldyn[numff].fdyn;
  fsidyn        = alldyn[genprob.numaf+1].fsidyn;

  actsysarray   = disnum_calc;
  stiff_array = disnum_calc;
  press_array = stiff_array+1;
  mass_array  = press_array+1;


  fdyn->dt      = fsidyn->dt;
  fdyn->maxtime = fsidyn->maxtime;
  fdyn->nstep   = fsidyn->nstep;
  grat          = ZERO;

  press_dis = 1;

  /* set some pointers */
  /* -> only valid for single field problem !!!! */
  actsolv     = &(solv[numff]);
  pressolv    = &(solv[numff+3]);
  actpart     = &(partition[numff]);
  action      = &(calc_action[numff]);
  restart     = genprob.restart;

  memset(&container,0,sizeof(CONTAINER));
  container.fieldtyp = actfield->fieldtyp;
  container.disnum   = disnum_calc;
  container.turbu    = fdyn->turbu;

  if (genprob.probtyp == prb_fsi)
    str = str_fsicoupling;

  ipos = &(actfield->dis[disnum_calc].ipos);

  /* if we are not parallel, we have to allocate an alibi *
   * intra-communicator structure */
#ifdef PARALLEL
  actintra    = &(par.intra[numff]);
#else
  actintra    = &work->dummy_intra;
#endif

  if (actintra->intra_fieldtyp != fluid)
    dserror("only fluid allowed");

  /* -------------------------------------------------------------------- */
  /* Create the pressure Laplace matrix. */

  trilinos_zero_matrix(actsolv->sysarray[mass_array].trilinos);
  trilinos_zero_matrix(actsolv->sysarray[press_array].trilinos);
  
  pm_calelm_laplace(actfield,actpart,disnum_calc,press_dis,
		    actsolv,press_array,mass_array,
		    actintra,ipos);
  
  /****************************************/

  /* get global and local number of equations */
  solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
		     actsolv->sysarray_typ[actsysarray],
		     &numeq,
		     &numeq_total);

#ifdef PARALLEL
#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
  numddof = 0;
#else
  numddof = actfield->dis[disnum_calc].numdf-numeq_total;
#endif
#endif

  grat = ZERO;

  /* check (starting) algorithm */
  if (fdyn->step<=(fdyn->nums+1))
    fluid_startproc(&nfrastep,0);


  /* calculate constants for time algorithm */
  fluid_tcons();


  /* output to the screen */
  if (par.myrank==0)
  {
    if (fdyn->iop == 4) printf("Solving FLUID by One-Step-Theta ...\n");
    else if (fdyn->iop == 7) printf("Solving FLUID by BDF2 ...\n");
    else dserror("wrong time integration scheme");
  }


  /* ALE-PHASE I */
  if (fsidyn->iale == 1)
  {
    /* get the gridvelocity */
    fsi_alecp(actfield,disnum_calc,alefield,adisnum_calc,fdyn->dta,fdyn->numdf,1);

    /* change element flag */
    fdyn->ishape=1;

    /* calculate ALE-convective velocities at time (n) */
    fsi_aleconv(actfield,disnum_calc,fdyn->numdf,ipos->convn,ipos->veln);
  }

  else  dserror("ALE field by function not implemented yet!\n");


  /* set dirichlet boundary conditions for  timestep */
  fluid_setdirich(actfield,disnum_calc,ipos,ipos->velnp);


  /* prepare time rhs in mass form */
  fluid_prep_rhs(actfield, disnum_calc, ipos);


  /* start time step for fluid on the screen */
  if (fdyn->itnorm!=fncc_no && par.myrank==0)
  {
    if (fdyn->freesurf>1)
    {
      printf("------------------------------------------------------------------------------- \n");
      printf("|- step/max -|-  tol     [norm] -|- vel. error -|- pre. error -|-  fs error  -|\n");
    }
    else
    {
      printf("---------------------------------------------------------------- \n");
      printf("|- step/max -|-  tol     [norm] -|- vel. error -|- pre. error -|\n");
    }
  }
  itnum=1;

#ifdef PERF
  perf_end(47);
#endif

  /*======================================================================*
   |               N O N L I N E A R   I T E R A T I O N                  |
   *======================================================================*/
nonlniter:

#ifdef PERF
  perf_begin(48);
#endif

  fdyn->itnum=itnum;

  /* calculate constants for nonlinear iteration */
  if(fdyn->freesurf) fluid_icons(itnum);

  /* ALE-PHASE II */
  if (fsidyn->iale == 1)
  {
    /* for implicit free surface we have to update the
       grid velocity during the iteration */
    if (fdyn->freesurf>1 && itnum>1)
    {
      fsi_alecp(actfield,disnum_calc,alefield,adisnum_calc,fdyn->dta,fdyn->numdf,fdyn->freesurf);

      /* change element flag */
      fdyn->ishape=1;
    }

    /* calculate ale-convective velocities at  time (n+1) */
    fsi_aleconv(actfield,disnum_calc,fdyn->numdf,ipos->convnp,ipos->velnp);
  }
  else
    dserror("ALE field by function not implemented yet!\n");


  /* calculate curvature at free surface */
  if (fdyn->surftens!=0)
  {
    *action = calc_fluid_curvature;
    fluid_curvature(actfield,actpart,actintra,action);
  }

  /* intitialise global matrix and global rhs */
  solserv_zero_vec(&(actsolv->rhs[0]));
  solserv_zero_mat(actintra,&(actsolv->sysarray[actsysarray]),
		   &(actsolv->sysarray_typ[actsysarray]));

  /* re-initialise neumann bcs */
  inherit_design_dis_neum(&(actfield->dis[disnum_calc]));

  /* initialise rhs */
  amzero(&work->frhs_a);
  amzero(&work->fgradprhs_a);

  /* form incremental matrices, residual and element forces */
  *action = calc_fluid;
  t1=ds_cputime();
  container.dvec         = NULL;
  container.frhs         = frhs;
  container.global_numeq = numeq_total;
  container.nii          = fdyn->nii;
  container.kstep        = 0;
  container.fieldtyp     = actfield->fieldtyp;
  container.is_relax     = 0;
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
	 &container,action);
  te=ds_cputime()-t1;
  work->tes+=te;

  /* Add the pressure term M*ML^-1*G*p(n) to the rhs. The vector
   * G*p(n) has already been calculated by the element call. */
  /* (This term could be done just once for each nonlinear iteration,
   * but maybe we'll do the convection explicit and circumvent the
   * iteration altogether.) */
  
  assemble_vec(actintra,
               &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(actsolv->rhs[0]),
               frhs,
               -fdyn->thsl
    );

  /*--------------------------------------------------------- add rhs: ---*/
  /* The pressure part is already there. Add the rhs terms from the
   * element calls. */
  assemble_vec(actintra,
	       &(actsolv->sysarray_typ[actsysarray]),
	       &(actsolv->sysarray[actsysarray]),
	       &(actsolv->rhs[0]),
	       frhs,
	       1.0
    );


#ifdef PERF
  perf_end(48);
#endif


  /* solve system */

#ifdef PERF
  perf_begin(49);
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
  work->tss+=ts;

#ifdef PERF
  perf_end(49);
#endif

  /* set flags for stability parameter evaluation and convergence check */
  fdyn->ishape=0;

#ifdef PERF
  perf_begin(50);
#endif

  /* return solution to the nodes and calculate the convergence ratios */
  fluid_result_incre(actfield, disnum_calc,actintra,&(actsolv->sol[0]),&(actsolv->rhs[0]),ipos->velnp,
		     &(actsolv->sysarray[actsysarray]),
		     &(actsolv->sysarray_typ[actsysarray]),
		     &vrat,&prat,&grat);

  /* store total area */
  if (fdyn->checkarea > 0)
  {
    if (work->totarea_a.fdim < itnum)
      dserror("cannot store totarea!\n");
    totarea[itnum-1] = fdyn->totarea;
  }

  /* solve heightfunction seperately */
  if (fdyn->freesurf == 3)
    fluid_heightfunc(2,&grat,actfield,actpart,actintra,action,
		     &container,ipos);

  /* update coordinates at free surface */
  if (fdyn->freesurf > 1)
    fluid_updfscoor(actfield, fdyn, fdyn->dta, ipos, 1);

  /* based on the new position calculate normal at free surface */
  if (itnum == 1) fluid_cal_normal(actfield,disnum_calc,0,action);

  /* iteration convergence check */
  converged = fluid_convcheck(vrat,prat,grat,itnum,te,ts);

#ifdef PERF
  perf_end(50);
#endif


  /* check if nonlinear iteration has to be finished */
  if (converged == 0)
  {
    itnum++;
    goto nonlniter;
  }
  /*----------------------------------------------------------------------*
   * -->  end of nonlinear iteration                                      *
   *----------------------------------------------------------------------*/


    solserv_zero_vec(work->press_rhs);
    amzero(&work->frhs_a);
    amzero(&work->fgradprhs_a);

    /* build up the rhs */
    pm_calprhs_cont(actfield, actpart, disnum_calc, actintra, ipos, work->press_rhs, fgradprhs, frhs);

    assemble_vec(actintra,
		 &(actsolv->sysarray_typ[press_array]),
		 &(actsolv->sysarray[press_array]),
		 work->press_rhs,
		 frhs,
		 1.0
      );

#if 0
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
		frhs[n->dof[0]]);
      }

      fprintf(allfiles.gidres,"END VALUES\n");
    }
#endif

#ifdef PERF
    perf_begin(80);
#endif
    /* solve for the pressure increment */
    solver_control(actfield,press_dis,pressolv,actintra,
                   &(actsolv->sysarray_typ[press_array]),
                   &(actsolv->sysarray[press_array]),
                   work->press_sol,
                   work->press_rhs,
                   0);
#ifdef PERF
    perf_end(80);
#endif

    /* update pressure */
    /* solserv_result_incre does not work due to misleading dirichlet
       conditions. The conditions on the velocity interfere. */
    {
      INT i;
      amzero(&work->frhs_a);
      solserv_reddistvec(work->press_sol,
			 &(actsolv->sysarray[press_array]),
			 &(actsolv->sysarray_typ[press_array]),
			 frhs,
			 work->press_sol->numeq_total,
			 actintra);
      for (i=0; i<actfield->dis[press_dis].numnp; i++)
      {
	NODE* actnode;
	actnode = &(actfield->dis[press_dis].node[i]);
	actnode->sol_increment.a.da[0][0]  = frhs[actnode->dof[0]];
	actnode->sol_increment.a.da[1][0] += frhs[actnode->dof[0]]/fdyn->thsl;
      }
    }

#if 0
    if (actintra->intra_rank==0)
    {
      fprintf(allfiles.gidres,"RESULT \"press_sol\" \"ccarat\" %d SCALAR ONNODES\n"
              "RESULTRANGESTABLE \"standard_fluid    \"\n"
              "COMPONENTNAMES \"pressure\"\n"
              "VALUES\n",fdyn->step);

      for (i=0; i<actfield->dis[press_dis].numnp; ++i)
      {
        NODE* n = &actfield->dis[press_dis].node[i];
        fprintf(allfiles.gidres,"%d %f\n",
                n->Id+1-genprob.nodeshift,
		frhs[n->dof[0]]);
      }

      fprintf(allfiles.gidres,"END VALUES\n");
    }
#endif

#if 0
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

    /* update velocity */

    solserv_zero_vec(&actsolv->sol[0]);
    solserv_zero_vec(&actsolv->rhs[0]);
    amzero(&work->frhs_a);

    /* build up rhs off velocity corretion */
    pm_calvrhs(actfield, actpart, disnum_calc, actintra, ipos, &actsolv->rhs[0], frhs);

    assemble_vec(actintra,
		 &(actsolv->sysarray_typ[mass_array]),
		 &(actsolv->sysarray[mass_array]),
		 &actsolv->rhs[0],
		 frhs,
		 1
      );

    /*-------------------------------------------------------- solve system */
#ifdef PERF
    perf_begin(80);
#endif

    init=0;
    t1=ds_cputime();
    solver_control(actfield,disnum_calc,actsolv, actintra,
		   &(actsolv->sysarray_typ[mass_array]),
		   &(actsolv->sysarray[mass_array]),
		   &(actsolv->sol[0]),
		   &(actsolv->rhs[0]),
                     init);
    ts=ds_cputime()-t1;
    work->tss+=ts;

#ifdef PERF
    perf_end(80);
#endif

    solserv_sol_zero(actfield,disnum_calc,node_array_sol_increment,ipos->velnm);

    fluid_result_incre(actfield,disnum_calc,actintra,&(actsolv->sol[0]),&(actsolv->rhs[0]),
		       ipos->velnm,
		       &(actsolv->sysarray[actsysarray]),
		       &(actsolv->sysarray_typ[actsysarray]),
		       &vrat,NULL,NULL);

    solserv_sol_add(actfield,disnum_calc,
		    node_array_sol_increment,
		    node_array_sol_increment,
		    ipos->velnm,
		    ipos->velnp,1);


  /* output of area to monitor file */
  if (fdyn->checkarea>0) out_area(work->totarea_a,fdyn->acttime,itnum,0);


#ifdef PERF
  perf_begin(51);
#endif

  /* calculate stresses transferred to structure */
  if (fsidyn->ifsi==fsi_basic_sequ_stagg ||
      fsidyn->ifsi==fsi_sequ_stagg_pred ||
      fsidyn->ifsi==fsi_sequ_stagg_shift ||
      fsidyn->ifsi==fsi_iter_stagg_fixed_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_steep_desc ||
      fsidyn->ifsi==fsi_iter_stagg_CHEB_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_Newton_FD ||
      fsidyn->ifsi==fsi_iter_stagg_Newton_I ||
      fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_force ||
      fsidyn->ifsi==fsi_iter_stagg_steep_desc_force)
  {
    if (fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_force ||
	fsidyn->ifsi==fsi_iter_stagg_steep_desc_force)
    {
      /* Store the old nodal forces */
      solserv_sol_copy(actfield,disnum_calc,
		       node_array_sol_mf,
		       node_array_sol_mf,
		       ipos->mf_forcenp,
		       ipos->mf_forcen);
    }

    if(fsidyn->coupforce == cf_nodeforce)
    {
      solserv_sol_zero(actfield,disnum_calc,node_array_sol_mf,ipos->mf_forcenp);
#ifdef PARALLEL

      amzero(&work->fcouple_a);
      fsi_cbf(&(actpart->pdis[disnum_calc]),fcouple,ipos,numeq_total,0);

      fsi_allreduce_coupforce(fcouple,recvfcouple,numeq_total,numddof,
			      actintra,actfield,disnum_calc);

#else

      fsi_cbf(&(actpart->pdis[disnum_calc]),NULL,ipos,0,0);

#endif  /* PARALLEL */
    }
    else
    {
      *action = calc_fluid_stress;
      container.nii= 0;
      container.str=str;
      container.is_relax = 0;
      calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
	     &container,action);

      /* since stresses are stored locally at the element it's necassary to
	 reduce them to all procs! */
      dsassert(actsolv->parttyp==cut_elements,
	       "Stress reduction for 'cut_nodes' not possible\n");

      fluid_reducestress(actintra,actpart,actfield,disnum_calc,fdyn->numdf,str);


      /* store stresses in sol_mf */
      solserv_sol_zero(actfield,disnum_calc,node_array_sol_mf,ipos->mf_forcenp);
      fsi_fluidstress_result(actfield,disnum_calc,fdyn->numdf);
    }
  }

#ifdef PERF
  perf_end(51);
#endif


#ifdef D_MORTAR
  if(fsidyn->coupmethod == 0) /* mortar method */
  {
    /* redefine the size of sol_mf from 2 to 3, the third field is */
    /* necessary to store the nodal forces due to fsi */
    solserv_sol_zero(actfield, disnum_calc, node_array_sol_mf, 3);
  }
#endif

#else
  dserror("TRILINOS_PACKAGE and PM_TRILINOS required");
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief fsi fluid finalize time step.
 */
/*----------------------------------------------------------------------*/
void fsi_fluid_pm_laplace_final(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  )
{
#if defined(TRILINOS_PACKAGE) && defined(PM_TRILINOS)
  INT             numff;              /* actual number of fluid field     */
  INT             actsysarray=0;      /* number of actual sysarray        */
  INT             restart;

  DOUBLE          grat;               /* convergence ratios               */
  SOLVAR         *actsolv;            /* pointer to active sol. structure */
  PARTITION      *actpart;            /* pointer to active partition      */
  INTRA          *actintra;           /* pointer to active intra-communic.*/
  CALC_ACTION    *action;             /* pointer to the cal_action enum   */

  DOUBLE         *frhs;               /* iteration - RHS                  */
  DOUBLE         *fgradprhs;

#ifdef PARALLEL
  DOUBLE         *fcouple;            /* to store fsi coupling forces     */
  DOUBLE         *recvfcouple;        /* to communicate fsi coupling forces*/
#endif

  DOUBLE         *totarea;
  CONTAINER       container;          /* variables for calelm             */
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  FSI_DYNAMIC    *fsidyn;             /* fsi dynamic variables     */
  ARRAY_POSITION *ipos;
  INT             press_dis;

#ifdef PERF
  perf_begin(52);
#endif

  /****************************************/

  frhs = work->frhs_a.a.dv;
  fgradprhs = work->fgradprhs_a.a.dv;
  
  press_dis = 1;
  
#ifdef PARALLEL
  fcouple = work->fcouple_a.a.dv;
  recvfcouple = work->recvfcouple_a.a.dv;
#endif

  totarea = work->totarea_a.a.dv;

  actsysarray   = disnum_calc;

  numff         = genprob.numff;
  fdyn          = alldyn[numff].fdyn;
  fsidyn        = alldyn[genprob.numaf+1].fsidyn;

  fdyn->dt      = fsidyn->dt;
  fdyn->maxtime = fsidyn->maxtime;
  fdyn->nstep   = fsidyn->nstep;
  grat          = ZERO;

  /* set some pointers */
  /* -> only valid for single field problem !!!! */
  actsolv     = &(solv[numff]);
  actpart     = &(partition[numff]);
  action      = &(calc_action[numff]);
  restart     = genprob.restart;

  memset(&container,0,sizeof(CONTAINER));
  container.fieldtyp = actfield->fieldtyp;
  container.disnum   = disnum_calc;
  container.turbu    = fdyn->turbu;

  ipos = &(actfield->dis[disnum_calc].ipos);

  /* if we are not parallel, we have to allocate an alibi *
   * intra-communicator structure */
#ifdef PARALLEL
  actintra    = &(par.intra[numff]);
#else
  actintra    = &work->dummy_intra;
#endif

  if (actintra->intra_fieldtyp != fluid)
    dserror("only fluid allowed");

  /****************************************/

  /* lift&drag computation */
  if (fdyn->liftdrag>0)
  {
    *action = calc_fluid_liftdrag;
    container.str=str_liftdrag;
    fluid_liftdrag(genprob.numfld,action,&container,actfield,disnum_calc,
		   actsolv,actpart,actintra,ipos);
  }


  /* update acceleration */
  if (fdyn->iop==4) /* for step theta */
  {

    /* evaluate acceleration in this time step *
     * depending on integration method */
    if (fdyn->step == 1)
    {
      /* do just a linear interpolation within the first timestep */
      solserv_sol_zero(actfield,disnum_calc,node_array_sol_increment,ipos->accn);
      solserv_sol_add(actfield,disnum_calc,
                      node_array_sol_increment,
		      node_array_sol_increment,
		      ipos->velnp,
                      ipos->accn,
                      1.0/fdyn->dta);
      solserv_sol_add(actfield,disnum_calc,
                      node_array_sol_increment,
		      node_array_sol_increment,
		      ipos->veln,
                      ipos->accn,
                      -1.0/fdyn->dta);
      solserv_sol_copy(actfield,disnum_calc,
                       node_array_sol_increment,
		       node_array_sol_increment,
                       ipos->accn,
                       ipos->accnm);
    }
    else
    {
      /* previous acceleration becomes (n-1)-acceleration of next step    */
      solserv_sol_copy(actfield,disnum_calc,
                       node_array_sol_increment,
		       node_array_sol_increment,
		       ipos->accn,
                       ipos->accnm);
      /* the following alternative to the previous copy command interferes
	 with the restart! */
      /*   leftspace = ipos->accnm;
	   ipos->accnm = ipos->accn;
	   ipos->accn  = leftspace;    */
      fluid_acceleration(actfield,disnum_calc,ipos,fdyn->iop);
    }
  }



  /* make predictor at free surface */
  if (fdyn->freesurf>0)
    fluid_updfscoor(actfield, fdyn, fdyn->dta, ipos, 0);


  /* based on the predictor calculate new normal at free surface */
  fluid_cal_normal(actfield,disnum_calc,2,action);


  /* copy solution from sol_increment[veln][j] to sol_increment[velnm][j] */
  /* -> prev. solution becomes (n-1)-solution of next time step */
  solserv_sol_copy(actfield,disnum_calc,
		   node_array_sol_increment,
		   node_array_sol_increment,
		   ipos->veln,
		   ipos->velnm);

  /* copy solution from sol_increment[velnp][j] to sol_increment[veln][j] */
  /* -> actual solution becomes previous solution of next time step */
  solserv_sol_copy(actfield,disnum_calc,
		   node_array_sol_increment,
		   node_array_sol_increment,
		   ipos->velnp,
		   ipos->veln);

  solserv_sol_copy(actfield,press_dis,
                   node_array_sol_increment,
                   node_array_sol,
                   1,
                   0);
  
  /* alternative to the above copys which does not work with restart: */
  /*-------------------------- shift position of old velocity solution ---*/
  /* leftspace = ipos->velnm;
     ipos->velnm = ipos->veln; */

  /*--------------------- shift position of previous velocity solution ---*/
  /* ipos->veln = ipos->velnp; */

  /*--------- set place for new solution to be solved in the next step ---*/
  /* ipos->velnp = leftspace; */

  /*---------- it is however necessary to have the newest solution ...
    ... still on ipos->velnp ---*/
  /*
    solserv_sol_copy(actfield,0,node_array_sol_increment,
    node_array_sol_increment,ipos->veln,ipos->velnp);
  */

  /* for multifield fluid problems with freesurface */
  /* copy solution from sol_increment[ipos->veln][j] to sol_mf[0][j] */
  /* check this for FSI with free surface!!! */
  if (fdyn->freesurf>0)
    solserv_sol_copy(actfield,disnum_calc,
		     node_array_sol_increment,
		     node_array_sol_mf,
		     ipos->velnp,
		     ipos->mf_velnp);

  /* finalise this timestep */
  work->outstep++;
  work->restartstep++;

  /*copy solution from sol_increment[ipos->veln][j] to sol[0][j]
    and transform kinematic to real pressure */
  solserv_sol_copy(actfield,disnum_calc,node_array_sol_increment,
		   node_array_sol,ipos->velnp,0);
  fluid_transpres(actfield,disnum_calc,0,0,fdyn->numdf-1,0);


#ifdef SUBDIV
  /* transfer the solution to the nodes of the master-dis */
  if (actfield->subdivide > 0)
  {
    solserv_sol_trans(actfield, disnum_calc, node_array_sol, 0);
  }
#endif



  if (work->outstep==fdyn->upout && ioflags.output_out==1 && ioflags.fluid_sol==1)
  {
    work->outstep=0;
    out_sol(actfield,actpart,disnum_io,actintra,fdyn->step,0);
  }


  /* write restart to pss file */
  if (work->restartstep==fsidyn->uprestart)
  {
    work->restartstep=0;
#ifdef BINIO
  if(disnum_io != disnum_calc)
    restart_write_bin_fluiddyn(&work->restart_context, fdyn);
  else
    restart_write_bin_fluiddyn(&work->out_context, fdyn);
#else
    restart_write_fluiddyn(fdyn,actfield,actpart,actintra,action,&container);
#endif
  }


  /* monitoring */
  if (ioflags.monitor==1)
    monitoring(actfield,disnum_calc,numff,0,fdyn->acttime);

#ifdef PERF
  perf_end(52);
#endif
#else
  dserror("TRILINOS_PACKAGE and PM_TRILINOS required");
#endif
}



    /*======================================================================*
      |     S O L U T I O N    F O R    S T E E P E S T    D E S C E N T     |
      |                        E V A L U A T I O N                           |
     *======================================================================*/
    /* nodal solution history fluid field:                                  *
     * sol[0][j]           ... initial data                                 *
     * sol[1...0][j]  ... solution for visualisation (real pressure)   *
     * sol_increment[0][j] ... solution at time (n-1)                       *
     * sol_increment[1][j] ... solution at time (n)                         *
     * sol_increment[2][j] ... solution at time (n+g)                       *
     * sol_increment[3][j] ... solution at time (n+1)                       *
     * sol_increment[4][i] ... grid velocity time (n) -> (n+1) #            *
     * sol_increment[5][i] ... convective velocity at time (n)              *
     * sol_increment[6][i] ... convective velocity at time (n+1) #          *
     * sol_increment[7][i] ... fluid solution for relax.-param. steep. desc.*
     * #: these vectors also used for steepest descent calculation          *
     * sol_mf[0][j]        ... solution at time (n+1)                       *
     * sol_mf[1][j]        ... nodal stresses at FS-interface at time (n+1) *
     *======================================================================*/


/*----------------------------------------------------------------------*/
/*!
  \brief fsi fluid steepest decent relaxation parameter calculation.

  Do the required fluid sensitivity calculation.
 */
/*----------------------------------------------------------------------*/
void fsi_fluid_pm_laplace_sd(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  )
{
#if defined(TRILINOS_PACKAGE) && defined(PM_TRILINOS)
  INT             numff;              /* actual number of fluid field     */
  INT             numeq;              /* number of equations on this proc */
  INT             numeq_total;        /* total number of equations        */
  INT             init;               /* flag for solver_control call     */
  INT             actsysarray=0;      /* number of actual sysarray        */
  INT             restart;

  DOUBLE          grat;               /* convergence ratios               */
  DOUBLE          t1,ts,te;
  SOLVAR         *actsolv;            /* pointer to active sol. structure */
  PARTITION      *actpart;            /* pointer to active partition      */
  INTRA          *actintra;           /* pointer to active intra-communic.*/
  CALC_ACTION    *action;             /* pointer to the cal_action enum   */

  DOUBLE         *frhs;               /* iteration - RHS                  */
  DOUBLE         *fgradprhs;

#ifdef PARALLEL
  INT             numddof;            /* number of Dirichlet dofs         */
  DOUBLE         *fcouple;            /* to store fsi coupling forces     */
  DOUBLE         *recvfcouple;        /* to communicate fsi coupling forces*/
#endif

  DOUBLE         *totarea;
  CONTAINER       container;          /* variables for calelm             */
  FLUID_STRESS    str=str_none;
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  FSI_DYNAMIC    *fsidyn;             /* fsi dynamic variables     */
  ARRAY_POSITION *ipos;

  /****************************************/

  frhs = work->frhs_a.a.dv;
  fgradprhs = work->fgradprhs_a.a.dv;
  
  totarea = work->totarea_a.a.dv;

  actsysarray   = disnum_calc;

  numff         = genprob.numff;
  fdyn          = alldyn[numff].fdyn;
  fsidyn        = alldyn[genprob.numaf+1].fsidyn;

  fdyn->dt      = fsidyn->dt;
  fdyn->maxtime = fsidyn->maxtime;
  fdyn->nstep   = fsidyn->nstep;
  grat          = ZERO;

  /* set some pointers */
  /* -> only valid for single field problem !!!! */
  actsolv     = &(solv[numff]);
  actpart     = &(partition[numff]);
  action      = &(calc_action[numff]);
  restart     = genprob.restart;

  memset(&container,0,sizeof(CONTAINER));
  container.fieldtyp = actfield->fieldtyp;
  container.disnum   = disnum_calc;
  container.turbu    = fdyn->turbu;

  if (genprob.probtyp == prb_fsi)
    str = str_fsicoupling;

  ipos = &(actfield->dis[disnum_calc].ipos);

  /* if we are not parallel, we have to allocate an alibi *
   * intra-communicator structure */
#ifdef PARALLEL
  actintra    = &(par.intra[numff]);
#else
  actintra    = &work->dummy_intra;
#endif

  if (actintra->intra_fieldtyp != fluid)
    dserror("only fluid allowed");

  /* get global and local number of equations */
  solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
		     actsolv->sysarray_typ[actsysarray],
		     &numeq,
		     &numeq_total);

#ifdef PARALLEL
  fcouple = work->fcouple_a.a.dv;
  recvfcouple = work->recvfcouple_a.a.dv;

#if defined(SOLVE_DIRICH) || defined(SOLVE_DIRICH2)
    numddof = 0;
#else
    numddof = actfield->dis[disnum_calc].numdf-numeq_total;
#endif /* SOLVE_DIRICH */
#endif

  dsassert(fsidyn->ifsi == fsi_iter_stagg_steep_desc ||
           fsidyn->ifsi == fsi_iter_stagg_steep_desc_force,
	   "No auxiliary fluid solution within this coupling scheme");

  /* calculate constants for time algorithm */
  fluid_tcons();

  /* output to the screen */
  if (par.myrank==0) printf("          - Solving FLUID ...\n");

  /* ALE-PHASE I */
  if (fsidyn->iale!=0)
  {
    /* change element flag */
    fdyn->ishape=1;

    /* calculate ALE-convective velocities at time (n) */
    fsi_aleconv(actfield,disnum_calc,fdyn->numdf,ipos->convnp,ipos->velnp);
  }

  /* set dirichlet boundary conditions */
  fluid_setdirich_sd(actfield,disnum_calc, ipos);

  /* calculate constants for nonlinear iteration */
  /*
    nir <->  EVALUATION OF NONLINEAR LHS N-REACTION
    nil <->  EVALUATION OF LUMPED MASS MATRIX (Mvv-lumped)
    nii <->  EVALUATION OF "ITERATION - RHS"
    nis <->  STATIONARY CASE (NO TIMEDEPENDENT TERMS) */
  fdyn->nir = 0;
  fdyn->nil = 0;
  fdyn->nii = 0;
  fdyn->nis = 0;

  /* calculate curvature at free surface */
  if (fdyn->surftens!=0)
  {
    dserror("steepest descent method with free surface not yet implemented.");
    *action = calc_fluid_curvature;
    fluid_curvature(actfield,actpart,actintra,action);
  }

  /* intitialise global matrix and global rhs */
  solserv_zero_vec(&(actsolv->rhs[0]));
  solserv_zero_mat(actintra,&(actsolv->sysarray[actsysarray]),
		   &(actsolv->sysarray_typ[actsysarray]));

  /* initialise iterations-rhs */
  amzero(&work->frhs_a);

  /* form incremental matrices, residual and element forces */
  *action = calc_fluid;
  t1=ds_cputime();
  container.dvec         = NULL;
  container.frhs         = frhs;
  container.global_numeq = numeq_total;
  container.nii          = fdyn->nii;
  container.kstep        = 0;
  container.fieldtyp     = actfield->fieldtyp;
  container.is_relax     = 1;
  calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
	 &container,action);
  te=ds_cputime()-t1;
  work->tes+=te;

  /* build rhs vector */
  assemble_vec(actintra,
	       &(actsolv->sysarray_typ[actsysarray]),
	       &(actsolv->sysarray[actsysarray]),
	       &(actsolv->rhs[0]),
	       frhs,
	       1.0
    );

  /* solve system */
  init=0;
  t1=ds_cputime();
  solver_control(actfield,disnum_calc,actsolv, actintra,
		 &(actsolv->sysarray_typ[actsysarray]),
		 &(actsolv->sysarray[actsysarray]),
		 &(actsolv->sol[0]),
		 &(actsolv->rhs[0]),
		 init);
  ts=ds_cputime()-t1;
  work->tss+=ts;

  /* set flags for stability parameter evaluation and convergence check */
  fdyn->ishape=0;

  /* return solution to the nodes to increment vector */
  solserv_result_incre(
    actfield,
    disnum_calc,
    actintra,
    &(actsolv->sol[actsysarray]),
    ipos->relax,
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->sysarray_typ[actsysarray]));

  /* calculate stresses transferred to structure */
  if (fsidyn->ifsi==fsi_basic_sequ_stagg ||
      fsidyn->ifsi==fsi_sequ_stagg_pred ||
      fsidyn->ifsi==fsi_sequ_stagg_shift ||
      fsidyn->ifsi==fsi_iter_stagg_fixed_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_steep_desc ||
      fsidyn->ifsi==fsi_iter_stagg_CHEB_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_force ||
      fsidyn->ifsi==fsi_iter_stagg_Newton_FD ||
      fsidyn->ifsi==fsi_iter_stagg_Newton_I ||
      fsidyn->ifsi==fsi_iter_stagg_steep_desc_force)
  {
    if (fsidyn->ifsi==fsi_iter_stagg_steep_desc_force)
    {
      /*
       * Hack: the force calculation expects to use the normal velocity. */
      /* backup */
      solserv_sol_copy(actfield, disnum_calc,
		       node_array_sol_increment,
		       node_array_sol_mf,
		       ipos->velnp,
		       ipos->mf_velcpy);
      solserv_sol_copy(actfield, disnum_calc,
		       node_array_sol_increment,
		       node_array_sol_increment,
		       ipos->relax,
		       ipos->velnp);
    }

    if (fsidyn->coupforce == cf_nodeforce)
    {
      dserror("node force coupling with steepest decent not yet working");
      solserv_sol_zero(actfield,disnum_calc,node_array_sol_mf,1);

#ifdef PARALLEL

      amzero(&work->fcouple_a);
      fsi_cbf(&(actpart->pdis[disnum_calc]),fcouple,ipos,numeq_total,0);

      fsi_allreduce_coupforce(fcouple,recvfcouple,numeq_total,numddof,
			      actintra,actfield,disnum_calc);

#else

      fsi_cbf(&(actpart->pdis[disnum_calc]),NULL,ipos,0,0);

#endif  /* PARALLEL */
    }
    else
    {
      *action = calc_fluid_stress;
      container.nii= 0;
      container.str=str;
      container.is_relax = 1;
      calelm(actfield,actsolv,actpart,actintra,actsysarray,-1,
	     &container,action);

      /* since stresses are stored locally at the element it's necassary to
	 reduce them to all procs! */
      dsassert(actsolv->parttyp==cut_elements,
	       "Stress reduction for 'cut_nodes' not possible\n");
      fluid_reducestress(actintra,actpart,actfield,disnum_calc,fdyn->numdf,str);

      /* store stresses in sol_mf */
      solserv_sol_zero(actfield,disnum_calc,node_array_sol_mf,ipos->mf_forcenp);
      fsi_fluidstress_result(actfield,disnum_calc,fdyn->numdf);
    }

    if (fsidyn->ifsi==fsi_iter_stagg_steep_desc_force)
    {
      /* restore */
      solserv_sol_copy(actfield, disnum_calc,
		       node_array_sol_mf,
		       node_array_sol_increment,
		       ipos->mf_velcpy,
		       ipos->velnp);
    }
  }
#else
  dserror("TRILINOS_PACKAGE and PM_TRILINOS required");
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief setup fsi fluid binary output
 */
/*----------------------------------------------------------------------*/
void fsi_fluid_pm_laplace_output(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  )
{
#if defined(TRILINOS_PACKAGE) && defined(PM_TRILINOS)
#ifdef BINIO
  if (ioflags.output_bin==1 && ioflags.fluid_sol==1)
  {
    INT             numff;              /* actual number of fluid field     */
    FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */

    numff         = genprob.numff;
    fdyn          = alldyn[numff].fdyn;

    out_results(&work->out_context, fdyn->acttime, fdyn->step, 0, OUTPUT_VELOCITY);
    out_results(&work->out_context, fdyn->acttime, fdyn->step, 0, OUTPUT_PRESSURE);
  }
#endif
#else
  dserror("TRILINOS_PACKAGE and PM_TRILINOS required");
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief cleanup fsi fluid algorithm
 */
/*----------------------------------------------------------------------*/
void fsi_fluid_pm_laplace_cleanup(
  FSI_FLUID_WORK *work,
  FIELD          *actfield,
  INT             disnum_calc,
  INT             disnum_io
  )
{
#if defined(TRILINOS_PACKAGE) && defined(PM_TRILINOS)
  INT             numff;              /* actual number of fluid field     */
  INT             i;                  /* counters                         */
  INT             actsysarray=0;      /* number of actual sysarray        */

  SOLVAR         *actsolv;            /* pointer to active sol. structure */
  PARTITION      *actpart;            /* pointer to active partition      */
  INTRA          *actintra;           /* pointer to active intra-communic.*/

#ifdef PARALLEL
  DOUBLE         *fcouple;            /* to store fsi coupling forces     */
  DOUBLE         *recvfcouple;        /* to communicate fsi coupling forces*/
#endif

  DOUBLE         *totarea;
  FLUID_DYNAMIC  *fdyn;               /* fluid dynamic variables   */
  FSI_DYNAMIC    *fsidyn;             /* fsi dynamic variables     */
  ARRAY_POSITION *ipos;

  /****************************************/

#ifdef PARALLEL
  fcouple = work->fcouple_a.a.dv;
  recvfcouple = work->recvfcouple_a.a.dv;
#endif

  totarea = work->totarea_a.a.dv;

  actsysarray   = disnum_calc;

  numff         = genprob.numff;
  fdyn          = alldyn[numff].fdyn;
  fsidyn        = alldyn[genprob.numaf+1].fsidyn;

  /* set some pointers */
  /* -> only valid for single field problem !!!! */
  actsolv     = &(solv[numff]);
  actpart     = &(partition[numff]);

  ipos = &(actfield->dis[disnum_calc].ipos);

  /* if we are not parallel, we have to allocate an alibi *
   * intra-communicator structure */
#ifdef PARALLEL
  actintra    = &(par.intra[numff]);
#else
  actintra    = &work->dummy_intra;
#endif

  if (actintra->intra_fieldtyp != fluid)
    dserror("only fluid allowed");

  /****************************************/

  /* print out solution to .out file */
  if (work->outstep!=0 && ioflags.output_out==1 && ioflags.fluid_sol==1)
    out_sol(actfield,actpart,disnum_io,actintra,fdyn->step,0);

  /* print total CPU-time to the screen */
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
	     par.myrank,work->tes);
      printf("PROC  %3d | FIELD FLUID     | total time for solver              : %10.3E \n",
	     par.myrank,work->tss);
    }
  }
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif

  /* tidy up */
  amdel(&work->frhs_a);
#ifdef PARALLEL
  if (fsidyn->coupforce == cf_nodeforce)
  {
    amdel(&work->fcouple_a);
    amdel(&work->recvfcouple_a);
  }
#endif
  if (fdyn->checkarea>0) amdel(&work->totarea_a);
  solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
  solserv_del_vec(&(actsolv->sol),actsolv->nsol);

#ifdef BINIO
  destroy_bin_out_field(&work->out_context);
  if(disnum_io != disnum_calc)
    destroy_bin_out_field(&work->restart_context);
#endif
#else
  dserror("TRILINOS_PACKAGE and PM_TRILINOS required");
#endif
}


#endif
#endif  /* ifdef D_FSI */


/*! @} (documentation module close)*/
