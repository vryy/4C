/*!----------------------------------------------------------------------
\file
\brief ale part of fsi-problems

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
/*!
\addtogroup FSI
*//*! @{ (documentation module open)*/
#ifdef D_FSI
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../ale3/ale3.h"
#include "../ale2/ale2.h"
#include "fsi_prototypes.h"
#include "../io/io.h"
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | struct _GENPROB       genprob; defined in global_control.c           |
 *----------------------------------------------------------------------*/
extern struct _GENPROB       genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;
/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate static variables if needed                       |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _STATIC_VAR  *statvar;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum _CALC_ACTION calc_action[MAXFIELD];
/*----------------------------------------------------------------------*
 |                                                          mn 06/02    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
 extern struct _IO_FLAGS     ioflags;
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
 |                                                       m.gee 06/01    |
 | structure allfiles, which holds all file pointers                    |
 | is defined in input_control_global.c
 *----------------------------------------------------------------------*/
 extern struct _FILES  allfiles;


void fsi_ale_nln_setup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  INT        	i;		  /* a counter  		                        */
  INT      numaf;            /* actual number of ale field                         */
  INT 	numeq;  	  /* number of equations on this proc                   */
  INT 	numeq_total;	  /* total number of equations over all procs           */
  INT        	init;		  /* init flag for solver                               */
  INT 	actsysarray;	  /* active sparse system matrix in actsolv->sysarray[] */
  INT      constsysarray;    /* second system matrix for steepest descent 		*/
  INT      numsys;           /* number of system matrices				*/
  SOLVAR       *actsolv;     /* pointer to the fields SOLVAR structure             */
  PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
  INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
  CALC_ACTION  *action;      /* pointer to the structures cal_action enum          */
  DOUBLE       *dirich;
  CONTAINER     container;   /* contains variables defined in container.h          */
  SPARSE_TYP    array_typ;   /* type of psarse system matrix                       */

  FSI_DYNAMIC      *fsidyn;
  ALE_DYNAMIC      *adyn;
  ARRAY_POSITION *ipos;

#ifndef PARALLEL
  work->dummy_intra.intra_fieldtyp = ale;
  work->dummy_intra.intra_rank     = 0;
  work->dummy_intra.intra_nprocs   = 1;
#endif

  numaf  = genprob.numaf;
  adyn   = alldyn[numaf].adyn;
  fsidyn = alldyn[numaf+1].fsidyn;
  ipos   = &(actfield->dis[disnum_calc].ipos);

  /*====================================================================*
   * nodal solution history ale field:                                  *
   * sol[1...0][j]  ... solution for visualisation (real pressure)	*
   * sol_mf[0][i]        ... displacements at (n)			*
   * sol_mf[1][i]        ... displacements at (n+1) 		        *
   * sol_mf[2][i]        ... used in fsi_gradient.c                     *
   * sol_mf[3][i]        ... displacements at (n+1) stored in fsi_calc_disp4ale*
   *====================================================================*/
  ipos->nummf = 4;
  ipos->mf_dispn = 0;
  ipos->mf_dispnp = 1;
  ipos->mf_dispnm = 2;
  ipos->mf_posnp = 2;

  ipos->numincr = 2;
  ipos->dispnp = 0;
  ipos->dispn = 1;

  adyn->dt=fsidyn->dt;
  adyn->maxtime=fsidyn->maxtime;
  adyn->nstep=fsidyn->nstep;

  container.isdyn   = 0;
  container.disnum  = disnum_calc;
  container.pos     = 0;

  work->outstep=0;
  work->restartstep=0;

  /*------------ the distributed system matrix, which is used for solving */
  actsysarray= disnum_calc;
  numsys = 1;  /* the number of system matrices */

  /*--------------------------------------------------- set some pointers */
  actsolv     = &(solv[numaf]);
  actpart     = &(partition[numaf]);
  action      = &(calc_action[numaf]);
  container.fieldtyp  = actfield->fieldtyp;

#ifdef PARALLEL
  actintra    = &(par.intra[numaf]);
#else
  actintra    = &work->dummy_intra;
#endif

  /*- there are only procs allowed in here, that belong to the structural */
  /*    intracommunicator (in case of linear statics, this should be all) */
  if (actintra->intra_fieldtyp != ale)
    dserror("only ale allowed");

  /*------------------------------------------------ typ of global matrix */
  array_typ   = actsolv->sysarray_typ[actsysarray];

  /*---------------------------- get global and local number of equations */
  /*   numeq equations are on this proc, the total number of equations is */
  /*                                                          numeq_total */
  solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                     actsolv->sysarray_typ[actsysarray],
                     &numeq,
                     &numeq_total);

  /*--- in the case of relaxation parameter calculations we would like to
    have a second system matrix which remains constant over time.     */
  if (fsidyn->ifsi == fsi_iter_stagg_steep_desc)
  {
    /*----------- one sysarray already exists, so copy the mask of it to */
    /* reallocate the vector of sparse matrices and the vector of there types
       formerly lenght 1, now lenght 2 */
    numsys++;
    actsolv->nsysarray=2;
    constsysarray = actsysarray+1;
    actsolv->sysarray_typ =
      (SPARSE_TYP*)CCAREALLOC(actsolv->sysarray_typ,actsolv->nsysarray*sizeof(SPARSE_TYP));
    actsolv->sysarray =
      (SPARSE_ARRAY*)CCAREALLOC(actsolv->sysarray,actsolv->nsysarray*sizeof(SPARSE_ARRAY));
    /*-copy the matrices sparsity mask */
    solserv_alloc_cp_sparsemask(  actintra,
                                  &(actsolv->sysarray_typ[actsysarray]),
                                  &(actsolv->sysarray[actsysarray]),
                                  &(actsolv->sysarray_typ[constsysarray]),
                                  &(actsolv->sysarray[constsysarray]));
  }

/*------------------------------------------------ output to the screen */
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
  for (i=0;i<par.nprocs;i++)
    if (par.myrank==i)
      printf("PROC  %3d | FIELD ALE       | number of equations      : %10d \n",
             par.myrank,numeq);
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
  if (par.myrank==0)
    printf("          | FIELD ALE       | total number of equations: %10d \n",numeq_total);
  if (par.myrank==0) printf("\n\n");

/*---------------------------------- number of rhs and solution vectors */
  actsolv->nrhs=1;
  actsolv->nsol=1;
  solserv_create_vec(&(actsolv->rhs),1,numeq_total,numeq,"DV");
  solserv_create_vec(&(actsolv->sol),1,numeq_total,numeq,"DV");
/*------------------------------ init the created dist. vectors to zero */
  solserv_zero_vec(&(actsolv->rhs[0]));
  solserv_zero_vec(&(actsolv->sol[0]));
/*----------- create a vector of full length for dirichlet part of rhs */
  dirich = amdef("intforce",&work->dirich_a,numeq_total,1,"DV");
/*--------------------------------------------------- initialize solver */
  init=1;
  for (i = actsysarray; i<actsysarray+numsys; i++)
  {
    solver_control(
      actfield,
      disnum_calc,
      actsolv,
      actintra,
      &(actsolv->sysarray_typ[i]),
      &(actsolv->sysarray[i]),
      &(actsolv->sol[0]),
      &(actsolv->rhs[0]),
      init
      );
    /*------------------------------ init the dist sparse matrix to zero */
    /*            NOTE: Has to be called after solver_control(init=1) */
    solserv_zero_mat(
      actintra,
      &(actsolv->sysarray[i]),
      &(actsolv->sysarray_typ[i])
      );
    /*---------------------------- init the assembly for sparse matrices */
    init_assembly(actpart,actsolv,actintra,actfield,i,disnum_calc);
  }

/*------------------------------- init the element calculating routines */
  *action = calc_ale_init_nln;
  calinit(actfield,actpart,action,&container);

/*--------------------------------------------------- init ale field ---*/
  fsi_init_ale(actfield,ipos->numincr);

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

/*--------------------------------------------------- check for restart */
  if (genprob.restart!=0)
  {
#if defined(BINIO)
    restart_read_bin_aledyn(adyn,
                            &(actsolv->sysarray_typ[i]),
                            &(actsolv->sysarray[i]),
                            actfield,
                            actpart,
                            disnum_io,
                            actintra,
                            genprob.restart);
#else
    restart_read_aledyn(genprob.restart,adyn,actfield,actpart,actintra);
#endif
  }

/*---------------------------------------------------------- monitoring */
  if (ioflags.monitor==1)
  {
    out_monitor(actfield,numaf,ZERO,1);
    monitoring(actfield,disnum_calc,numaf,0,adyn->time);
  }

/*------------------------------------------- print out results to .out */
#ifdef PARALLEL
  if (ioflags.ale_disp==1 && par.myrank==0)
    out_gid_domains(actfield,disnum_io);
#endif

  if (fsidyn->ifsi == fsi_iter_stagg_steep_desc) /* if we do steepest descent method! */
  {
    /* init linear element routines also */
    *action = calc_ale_init;
    calinit(actfield,actpart,action,&container);
    /* evaluate linear element stiffnesses (remain constant) */
    *action = calc_ale_stiff;
    container.dvec         = NULL;
    container.dirich       = NULL;
    container.global_numeq = 0;
    container.kstep        = 0;
    calelm(actfield,actsolv,actpart,actintra,constsysarray,-1,&container,action);
  }
}


void fsi_ale_nln_calc(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  INT      numaf;            /* actual number of ale field                         */
  INT 	numeq;  	  /* number of equations on this proc                   */
  INT 	numeq_total;	  /* total number of equations over all procs           */
  INT        	init;		  /* init flag for solver                               */
  INT 	actsysarray;	  /* active sparse system matrix in actsolv->sysarray[] */
  INT      numsys;           /* number of system matrices				*/
  SOLVAR       *actsolv;     /* pointer to the fields SOLVAR structure             */
  PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
  INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
  CALC_ACTION  *action;      /* pointer to the structures cal_action enum          */
  DOUBLE       *dirich;
  CONTAINER     container;   /* contains variables defined in container.h          */
  SPARSE_TYP    array_typ;   /* type of psarse system matrix                       */

  FSI_DYNAMIC      *fsidyn;
  ALE_DYNAMIC      *adyn;
  ARRAY_POSITION *ipos;

  dirich = work->dirich_a.a.dv;

  numaf  = genprob.numaf;
  adyn   = alldyn[numaf].adyn;
  fsidyn = alldyn[numaf+1].fsidyn;
  ipos   = &(actfield->dis[disnum_calc].ipos);

  container.isdyn   = 0;
  container.disnum  = disnum_calc;
  container.pos     = 0;

  /*------------ the distributed system matrix, which is used for solving */
  actsysarray= disnum_calc;
  numsys = 1;  /* the number of system matrices */

  /*--------------------------------------------------- set some pointers */
  actsolv     = &(solv[numaf]);
  actpart     = &(partition[numaf]);
  action      = &(calc_action[numaf]);
  container.fieldtyp  = actfield->fieldtyp;

#ifdef PARALLEL
  actintra    = &(par.intra[numaf]);
#else
  actintra    = &work->dummy_intra;
#endif

  /*- there are only procs allowed in here, that belong to the structural */
  /*    intracommunicator (in case of linear statics, this should be all) */
  if (actintra->intra_fieldtyp != ale)
    dserror("only ale allowed");

  /*------------------------------------------------ typ of global matrix */
  array_typ   = actsolv->sysarray_typ[actsysarray];

  /*---------------------------- get global and local number of equations */
  /*   numeq equations are on this proc, the total number of equations is */
  /*                                                          numeq_total */
  solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                     actsolv->sysarray_typ[actsysarray],
                     &numeq,
                     &numeq_total);

  if (par.myrank==0)
  {
    printf("Solving ALE (Je_stiffened)...\n");
    printf("\n");
  }

/*--------------------------------------- sequential staggered schemes: */
/*------------------------ copy from nodal sol_mf[1][j] to sol_mf[0][j] */
  if (fsidyn->ifsi==fsi_basic_sequ_stagg ||
      fsidyn->ifsi==fsi_sequ_stagg_pred)
    solserv_sol_copy(actfield,disnum_calc,
		     node_array_sol_mf,
		     node_array_sol_mf,
		     ipos->mf_dispnp,
		     ipos->mf_dispn);

  dsassert(fsidyn->ifsi!=fsi_sequ_stagg_shift,
           "ale-solution handling not implemented for algo with DT/2-shift!\n");

/*------------------------------ init the created dist. vectors to zero */
  solserv_zero_vec(&(actsolv->rhs[0]));
  solserv_zero_vec(&(actsolv->sol[0]));
/*--------------------------------------------------------------------- */
  amzero(&work->dirich_a);
/*-------------------------set dirichlet boundary conditions on at time */
  ale_setdirich_increment_fsi(actfield,disnum_calc,adyn,0);
/*----------------------------------------------------------------------*/
  solserv_zero_mat(actintra,
                   &(actsolv->sysarray[actsysarray]),
                   &(actsolv->sysarray_typ[actsysarray])
    );
/*--- call element routines to calculate & assemble stiffness matrix ---*/
  if (adyn->step <= adyn->num_initstep)
    *action = calc_ale_stiff_stress;
  else
    *action = calc_ale_stiff_nln;
  container.dvec         = NULL;
  container.dirich       = dirich;
  container.global_numeq = numeq_total;
  container.isdyn        = 1;
  switch (adyn->measure_quality)
  {
    /*--------------------------------*/
  case no_quality:
    container.quality = 0;
    break;
    /*--------------------------------*/
  case aspect_ratio:
    container.quality = 1;
    break;
    /*--------------------------------*/
  case corner_angle:
    container.quality = 2;
    break;
    /*--------------------------------*/
  case min_detF:
    container.quality = 3;
    break;
  default:
    dserror("element quality unknown");
    break;
  }
  calelm(actfield,actsolv,actpart,actintra,
         actsysarray,-1,&container,action);
/*------------------------ add rhs from prescribed displacements to rhs */
  assemble_vec(actintra,&(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),&(actsolv->rhs[0]),
               dirich,1.0);

/*--------------------------------------------------------- call solver */
  init=0;
  solver_control(
    actfield,
    disnum_calc,
    actsolv,
    actintra,
    &(actsolv->sysarray_typ[actsysarray]),
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->sol[0]),
    &(actsolv->rhs[0]),
    init
    );
/*---------------------allreduce the result and put it to sol_increment */
  solserv_result_incre(
    actfield,
    disnum_calc,
    actintra,
    &(actsolv->sol[0]),
    ipos->dispnp,
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->sysarray_typ[actsysarray])
    );

/*---------- add actual solution increment to sol (to serve output): ---*/
/* step 1: */
  /* copy prev. solution form sol_increment[1][j] to sol[0][j] ---*/
  solserv_sol_copy(actfield,disnum_calc,
		   node_array_sol_increment,
		   node_array_sol,
		   ipos->dispn,
		   0);
/* step 2: */
  /*----------- add actual solution increment to sol[0][j] ---*/
  solserv_sol_add(actfield,disnum_calc,
		  node_array_sol_increment,
		  node_array_sol,
		  ipos->dispnp,
		  0,
		  1.0);

/* save actual solution to be the previous one in the next time step ---*/
/*--------- copy actual solution from sol[0][i] to sol_mf[1][i] ---*/
  solserv_sol_copy(actfield,disnum_calc,
		   node_array_sol,
		   node_array_sol_mf,
		   0,
		   ipos->mf_dispnp);
}


void fsi_ale_nln_final(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  INT      numaf;            /* actual number of ale field                         */
  INT 	numeq;  	  /* number of equations on this proc                   */
  INT 	numeq_total;	  /* total number of equations over all procs           */
  INT 	actsysarray;	  /* active sparse system matrix in actsolv->sysarray[] */
  INT      numsys;           /* number of system matrices				*/
  SOLVAR       *actsolv;     /* pointer to the fields SOLVAR structure             */
  PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
  INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
  CALC_ACTION  *action;      /* pointer to the structures cal_action enum          */
  DOUBLE       *dirich;
  CONTAINER     container;   /* contains variables defined in container.h          */
  SPARSE_TYP    array_typ;   /* type of psarse system matrix                       */

  FSI_DYNAMIC      *fsidyn;
  ALE_DYNAMIC      *adyn;
  ARRAY_POSITION *ipos;

  dirich = work->dirich_a.a.dv;

  numaf  = genprob.numaf;
  adyn   = alldyn[numaf].adyn;
  fsidyn = alldyn[numaf+1].fsidyn;
  ipos   = &(actfield->dis[disnum_calc].ipos);

  container.isdyn   = 0;
  container.disnum  = disnum_calc;
  container.pos     = 0;

  /*------------ the distributed system matrix, which is used for solving */
  actsysarray= disnum_calc;
  numsys = 1;  /* the number of system matrices */

  /*--------------------------------------------------- set some pointers */
  actsolv     = &(solv[numaf]);
  actpart     = &(partition[numaf]);
  action      = &(calc_action[numaf]);
  container.fieldtyp  = actfield->fieldtyp;

#ifdef PARALLEL
  actintra    = &(par.intra[numaf]);
#else
  actintra    = &work->dummy_intra;
#endif

  /*- there are only procs allowed in here, that belong to the structural */
  /*    intracommunicator (in case of linear statics, this should be all) */
  if (actintra->intra_fieldtyp != ale)
    dserror("only ale allowed");

  /*------------------------------------------------ typ of global matrix */
  array_typ   = actsolv->sysarray_typ[actsysarray];

  /*---------------------------- get global and local number of equations */
  /*   numeq equations are on this proc, the total number of equations is */
  /*                                                          numeq_total */
  solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                     actsolv->sysarray_typ[actsysarray],
                     &numeq,
                     &numeq_total);

  /*------------------------------------ for iterative staggared schemes: */
  /*------------------------ copy from nodal sol_mf[1][j] to sol_mf[0][j] */
  if (fsidyn->ifsi==fsi_iter_stagg_fixed_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_steep_desc ||
      fsidyn->ifsi==fsi_iter_stagg_CHEB_rel_param ||
      fsidyn->ifsi==fsi_iter_stagg_AITKEN_rel_force ||
      fsidyn->ifsi==fsi_iter_stagg_steep_desc_force ||
      fsidyn->ifsi==fsi_iter_stagg_Newton_FD ||
      fsidyn->ifsi==fsi_iter_stagg_Newton_I ||
      fsidyn->ifsi==fsi_coupling_freesurface)
  {
    solserv_sol_copy(actfield,disnum_calc,
		     node_array_sol_mf,
		     node_array_sol_mf,
		     ipos->mf_dispn,
		     ipos->mf_dispnm);
    solserv_sol_copy(actfield,disnum_calc,
		     node_array_sol_mf,
		     node_array_sol_mf,
		     ipos->mf_dispnp,
		     ipos->mf_dispn);
  }

/*--------------------- to get the corrected free surface position copy
  --------------------------------- from sol_mf[1][j] to sol[0][j] */
  solserv_sol_copy(actfield,disnum_calc,
		   node_array_sol_mf,
		   node_array_sol,
		   ipos->mf_dispn,
		   0);

/*---------------------------- from sol_mf[1][j] to sol_increment[1][j] */
  solserv_sol_copy(actfield,disnum_calc,
		   node_array_sol_mf,
		   node_array_sol_increment,
		   ipos->mf_dispnp,
		   ipos->dispn);

#ifdef SUBDIV
  /* transfer the solution to the nodes of the master-dis */
  if (actfield->subdivide > 0)
  {
    solserv_sol_trans(actfield, disnum_calc, node_array_sol, 0);
  }
#endif


/*------------------------------------------- print out results to .out */
  work->outstep++;
  work->restartstep++;

  if (work->outstep==adyn->updevry_disp && ioflags.ale_disp==1 && ioflags.output_out==1)
  {
    work->outstep=0;
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,0);
  }


/*---------------------------------------------------------- monitoring */
  if (ioflags.monitor==1)
    monitoring(actfield,disnum_calc,numaf,0,adyn->time);

/*------------------------------------------------- write restart data */
  if (work->restartstep==fsidyn->uprestart)
  {
    work->restartstep=0;
#ifdef BINIO
    if(disnum_io != disnum_calc)
      restart_write_bin_aledyn(&work->restart_context, adyn);
    else
      restart_write_bin_aledyn(&work->out_context, adyn);
#else
    restart_write_aledyn(adyn,actfield,actpart,actintra);
#endif
  }

/*--------------------------------------- do mesh quality statistics ---*/
  if (adyn->measure_quality!=no_quality)
    plot_ale_quality(actfield,disnum_calc,fsidyn->step,(fsidyn->time-fsidyn->dt),actintra,actpart);
/*--------------------------------------------------------------------- */
}


void fsi_ale_nln_sd(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  FIELD             *structfield,
  INT                sdisnum
  )
{
  INT      numaf;            /* actual number of ale field                         */
  INT 	numeq;  	  /* number of equations on this proc                   */
  INT 	numeq_total;	  /* total number of equations over all procs           */
  INT        	init;		  /* init flag for solver                               */
  INT 	actsysarray;	  /* active sparse system matrix in actsolv->sysarray[] */
  INT      constsysarray;    /* second system matrix for steepest descent 		*/
  INT      numsys;           /* number of system matrices				*/
  SOLVAR       *actsolv;     /* pointer to the fields SOLVAR structure             */
  PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
  INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
  CALC_ACTION  *action;      /* pointer to the structures cal_action enum          */
  DOUBLE       *dirich;
  CONTAINER     container;   /* contains variables defined in container.h          */
  SPARSE_TYP    array_typ;   /* type of psarse system matrix                       */

  FSI_DYNAMIC      *fsidyn;
  ALE_DYNAMIC      *adyn;
  ARRAY_POSITION *ipos;

  dirich = work->dirich_a.a.dv;

  numaf  = genprob.numaf;
  adyn   = alldyn[numaf].adyn;
  fsidyn = alldyn[numaf+1].fsidyn;
  ipos   = &(actfield->dis[disnum_calc].ipos);

  container.isdyn   = 0;
  container.disnum  = disnum_calc;
  container.pos     = 0;

  /*------------ the distributed system matrix, which is used for solving */
  actsysarray= disnum_calc;
  numsys = 1;  /* the number of system matrices */
  constsysarray = actsysarray+1;

  /*--------------------------------------------------- set some pointers */
  actsolv     = &(solv[numaf]);
  actpart     = &(partition[numaf]);
  action      = &(calc_action[numaf]);
  container.fieldtyp  = actfield->fieldtyp;

#ifdef PARALLEL
  actintra    = &(par.intra[numaf]);
#else
  actintra    = &work->dummy_intra;
#endif

  /*- there are only procs allowed in here, that belong to the structural */
  /*    intracommunicator (in case of linear statics, this should be all) */
  if (actintra->intra_fieldtyp != ale)
    dserror("only ale allowed");

  /*------------------------------------------------ typ of global matrix */
  array_typ   = actsolv->sysarray_typ[actsysarray];

  /*---------------------------- get global and local number of equations */
  /*   numeq equations are on this proc, the total number of equations is */
  /*                                                          numeq_total */
  solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                     actsolv->sysarray_typ[actsysarray],
                     &numeq,
                     &numeq_total);

  if (par.myrank==0)
  {
    printf("          - Solving ALE ... \n");
  }

  dsassert(fsidyn->ifsi!=fsi_sequ_stagg_shift,
           "ale-solution handling not implemented for algo with DT/2-shift!\n");

/*------------------------------ init the created dist. vectors to zero */
  solserv_zero_vec(&(actsolv->rhs[0]));
  solserv_zero_vec(&(actsolv->sol[0]));
/*--------------------------------------------------------------------- */
  amzero(&work->dirich_a);

/*-------------------------set dirichlet boundary conditions on at time */
/* note: ale Dirichlet boundary conditions different from ZERO are not
   helpful for fsi coupling problems and would cause trouble here.
   But there's no test to set all ordinary dbc = 0.0 !!!
   The required Dirichlet boundary conditions from fsi coupling
   are calculated here.						*/
  ale_setdirich(actfield,disnum_calc,adyn,structfield->dis[sdisnum].ipos.mf_sd_g);

/*------------------------------- call element-routines to assemble rhs */
  *action = calc_ale_rhs;
  ale_rhs(actsolv,actpart,disnum_calc,actintra,constsysarray,-1,dirich,
          numeq_total,&container,action);

/*------ add rhs from fsi coupling (-> prescribed displacements) to rhs */
  assemble_vec(actintra,&(actsolv->sysarray_typ[constsysarray]),
               &(actsolv->sysarray[constsysarray]),&(actsolv->rhs[0]),
               dirich,1.0);

/*--------------------------------------------------------- call solver */
/*---- system matrix has been calculated within initialisation phase ---*/
  init=0;
  solver_control(
    actfield,
    disnum_calc,
    actsolv,
    actintra,
    &(actsolv->sysarray_typ[constsysarray]),
    &(actsolv->sysarray[constsysarray]),
    &(actsolv->sol[0]),
    &(actsolv->rhs[0]),
    init
    );

/*-------------- allreduce the result and put it to sol_increment[0][i] */
  solserv_result_incre(
    actfield,
    disnum_calc,
    actintra,
    &(actsolv->sol[0]),
    ipos->dispnp,
    &(actsolv->sysarray[constsysarray]),
    &(actsolv->sysarray_typ[constsysarray])
    );
}


void fsi_ale_nln_output(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
#ifdef BINIO
  if (ioflags.output_bin)
    if (ioflags.ale_disp==1)
    {
      ALE_DYNAMIC      *adyn;
      adyn   = alldyn[genprob.numaf].adyn;
      out_results(&work->out_context, adyn->time, adyn->step, 0, OUTPUT_DISPLACEMENT);
    }
#endif
}


void fsi_ale_nln_cleanup(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io
  )
{
  INT      numaf;            /* actual number of ale field                         */
  INT 	numeq;  	  /* number of equations on this proc                   */
  INT 	numeq_total;	  /* total number of equations over all procs           */
  INT 	actsysarray;	  /* active sparse system matrix in actsolv->sysarray[] */
  INT      numsys;           /* number of system matrices				*/
  SOLVAR       *actsolv;     /* pointer to the fields SOLVAR structure             */
  PARTITION    *actpart;     /* pointer to the fields PARTITION structure          */
  INTRA        *actintra;    /* pointer to the fields intra-communicator structure */
  CALC_ACTION  *action;      /* pointer to the structures cal_action enum          */
  DOUBLE       *dirich;
  CONTAINER     container;   /* contains variables defined in container.h          */
  SPARSE_TYP    array_typ;   /* type of psarse system matrix                       */

  FSI_DYNAMIC      *fsidyn;
  ALE_DYNAMIC      *adyn;

  dirich = work->dirich_a.a.dv;

  numaf  = genprob.numaf;
  adyn   = alldyn[numaf].adyn;
  fsidyn = alldyn[numaf+1].fsidyn;

  container.isdyn   = 0;
  container.disnum  = disnum_calc;
  container.pos     = 0;

  /*------------ the distributed system matrix, which is used for solving */
  actsysarray= disnum_calc;
  numsys = 1;  /* the number of system matrices */

  /*--------------------------------------------------- set some pointers */
  actsolv     = &(solv[numaf]);
  actpart     = &(partition[numaf]);
  action      = &(calc_action[numaf]);
  container.fieldtyp  = actfield->fieldtyp;

#ifdef PARALLEL
  actintra    = &(par.intra[numaf]);
#else
  actintra    = &work->dummy_intra;
#endif

  /*- there are only procs allowed in here, that belong to the structural */
  /*    intracommunicator (in case of linear statics, this should be all) */
  if (actintra->intra_fieldtyp != ale)
    dserror("only ale allowed");

  /*------------------------------------------------ typ of global matrix */
  array_typ   = actsolv->sysarray_typ[actsysarray];

  /*---------------------------- get global and local number of equations */
  /*   numeq equations are on this proc, the total number of equations is */
  /*                                                          numeq_total */
  solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                     actsolv->sysarray_typ[actsysarray],
                     &numeq,
                     &numeq_total);

/*------------------------------------------- print out results to .out */
  if (work->outstep!=0 && ioflags.ale_disp==1 && ioflags.output_out==1)
    out_sol(actfield,actpart,disnum_io,actintra,adyn->step,0);

/*------------------------------------------------------------- tidy up */
  solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
  solserv_del_vec(&(actsolv->sol),actsolv->nsol);

#ifdef BINIO
  destroy_bin_out_field(&work->out_context);
  if(disnum_io != disnum_calc)
    destroy_bin_out_field(&work->restart_context);
#endif
}


#if 0
/*!----------------------------------------------------------------------
\brief  solving for mesh displacements considering previous displacement

<pre>                                                            ck 06/03

in this function the mesh of a multifield problem is solved as a
pseudo structure. The stiffness is calculated newly in every timestep and
takes the actual element shape into account. Nonuniform additional
stiffening factors depending on the minimal Jacobian in the element are
used. Displacements are prescribed at the fluid structure
interface and at the free surface as Dirichlet boundary conditions

</pre>

\param  *actfield   FIELD          (i)     ale field
\param   mctrl      INT            (i)     control flag
\warning
\return void

\sa   calling: calelm(), monitoring(), ale_setdirich_increment_fsi(),
               plot_ale_quality()
      called by: fsi_ale()

*----------------------------------------------------------------------*/
void fsi_ale_nln(
  FSI_ALE_WORK      *work,
  FIELD             *actfield,
  INT                disnum_calc,
  INT                disnum_io,
  INT                mctrl
  )
{
#ifdef DEBUG
dstrc_enter("fsi_ale_nln");
#endif

switch (mctrl)
{
/*======================================================================*
 |                      I N I T I A L I S A T I O N                     |
 *======================================================================*/
case 1:
  fsi_ale_nln_setup(work,actfield,disnum_calc,disnum_io);
  break;

/*======================================================================*
 |                     S O L U T I O N    P H A S E                     |
 *======================================================================*
 * nodal solution history ale field:                                    *
 * sol[1...0][j]  ... solution for visualisation       		*
 * sol_mf[0][i]        ... displacements at (n)			        *
 * sol_mf[1][i]        ... displacements at (n+1) 		        *
 * sol_increment[0][i] ... working array (actual sol. increment)	*
 * sol_increment[1][i] ... solution of previous time step 		*
 *======================================================================*/
case 2:
  fsi_ale_nln_calc(work,actfield,disnum_calc,disnum_io);
  if (alldyn[genprob.numaf+1].fsidyn->ifsi>=fsi_iter_stagg_fixed_rel_param)
    break;

/*======================================================================*
 |                       F I N A L I S I N G                            |
 *======================================================================*/
case 3:
  fsi_ale_nln_final(work,actfield,disnum_calc,disnum_io);
  break;

/*======================================================================*
 |                A D D I T I O N A L   S O L U T I O N                 |
 |    for determination of relaxation parameter via steepest descent    |
 |                            l i n e a r !                             |
 *======================================================================*
 * nodal solution history ale field:                                    *
 * sol[1...0][j]  ... solution for visualisation (real pressure)	*
 * sol_mf[0][i]        ... displacements at (n)			        *
 * sol_mf[1][i]        ... displacements at (n+1) 		        *
 * sol_mf[2][i]        ... grid position in relaxation parameter calc.  *
 * sol_increment[0][i] ... displacement used to determine omega_RELAX   *
 *======================================================================*/
case 6:
  fsi_ale_nln_sd(work,actfield,disnum_calc,disnum_io);
  break;

/*======================================================================*
                            Binary Output
 *======================================================================*/
case 98:
  fsi_ale_nln_output(work,actfield,disnum_calc,disnum_io);
  break;

/*======================================================================*
 |                C L E A N I N G   U P   P H A S E                     |
 *======================================================================*/
case 99:
  fsi_ale_nln_cleanup(work,actfield,disnum_calc,disnum_io);
  break;
default:
   dserror("Parameter out of range: mctrl \n");
} /* end switch (mctrl) */

#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of fsi_ale_nln */
#endif

#endif
/*! @} (documentation module close)*/
