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
#include "ls_prototypes.h"
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct  _FIELD         *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern         ALLDYNA        *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB         genprob;
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
extern struct _PARTITION      *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS        ioflags;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR             par;
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern struct _CURVE          *curve;
extern INT                     numcurve;
/*----------------------------------------------------------------------*
 | enum _CALC_ACTION                                      m.gee 1/02    |
 | command passed from control routine to the element level             |
 | to tell element routines what to do                                  |
 | defined globally in global_calelm.c                                  |
 *----------------------------------------------------------------------*/
extern enum   _CALC_ACTION     calc_action[MAXFIELD];


extern struct _LS_DYNAMIC     *lsdyn;           /* ls dynamic variables */
extern struct _FIELD          *lsetfield;



static INT             numls;            /* field number for level set field */
static INT             init;             /* flag for solver_control call     */
static INT             numeq;            /* number of equations on this proc */
static INT             numeq_total;      /* total number of equations        */
static INT             actsysarray = 0;  /* number of actual sysarray        */

static DOUBLE          lrat;
static DOUBLE          t1,ts,te;
static DOUBLE          tes = ZERO;
static DOUBLE          tss = ZERO;
static SOLVAR         *actsolv;          /* pointer to active sol. structure */
static PARTITION      *actpart;          /* pointer to active partition      */
static INTRA          *actintra;         /* pointer to active intra-communic.*/
static CALC_ACTION    *action;           /* pointer to the cal_action enum   */

static INT             restart;

static ARRAY           ftimerhs_a;
static DOUBLE         *ftimerhs;	 /* time - RHS		             */
static ARRAY           fiterhs_a;
static DOUBLE         *fiterhs;	         /* iteration - RHS  		     */
static CONTAINER       container;        /* variables for calelm             */



/*!----------------------------------------------------------------------
\brief control subroutine for level set sub-problem in coupled level set /
extended finite element problem

<pre>                                                            irhan 05/04
control subroutine for level set sub-problem in coupled level set /
extended finite element problem
</pre>

*----------------------------------------------------------------------*/
void ls_levelset(
  INT     eflag
  )
{
#ifdef DEBUG
  dstrc_enter("ls_levelset");
#endif
/*----------------------------------------------------------------------*/

  switch (eflag)
  {
/*----------------------------------------- initialize levelset problem */
      case 1:
        ls_levelset_init_data();
        ls_levelset_init_levelset();
        break;
/*---------------------------------------------- solve levelset problem */
      case 2:
        ls_levelset_solv();
        break;
/*------------------------------------------- finalize levelset problem */
      case 3:
        ls_levelset_fina();
        break;
/*--------------------------------------------------- initialize solver */
      case 4:
        ls_levelset_init_solver();
        break;
/*-------------------------------------- write restart data to pss file */
      case 5:
        ls_levelset_write_restart();
        break;
/*------------------------------------- read restart data from pss file */
      case 6:
        ls_levelset_read_restart();
        break;
/*-------------------- reset level set solution before reinitialization */
      case 7:
        ls_levelset_reset_sol_increment();
        break;
/*--------------------------------------------------------------- clean */
      case 99:
        ls_levelset_clea();
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
\brief initialization of sub-problem level set

<pre>                                                            irhan 05/04
initialization of sub-problem level set (data)
</pre>

*----------------------------------------------------------------------*/
void ls_levelset_init_data()
{
#ifdef DEBUG
  dstrc_enter("ls_levelset_init_data");
#endif
/*----------------------------------------------------------------------*/

  /* I N I T I A L I Z A T I O N */
  numls = genprob.numls;
  lrat = ZERO;
  /* set some pointers */
  actsolv = &(solv[numls]);
  actpart = &(partition[numls]);
  action = &(calc_action[numls]);
  restart = genprob.restart;
  container.fieldtyp = lsetfield->fieldtyp;
  container.actndis = 0;
  /*
   * if we are not parallel, we have to allocate an alibi
   * intra-communicator structure
   */
#ifdef PARALLEL
  actintra    = &(par.intra[numls]);
#else
  actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  if (!actintra) dserror("Allocation of INTRA failed");
  actintra->intra_fieldtyp = levelset;
  actintra->intra_rank     = 0;
  actintra->intra_nprocs   = 1;
#endif
  /*
   * there are only procs allowed in here, that belong to the levelset
   * intracommunicator (in case of nonlinear levelset dyn., this should be all)
   */
  if (actintra->intra_fieldtyp!=levelset) dserror("fieldtyp != levelset");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_levelset_init_data */



/*!----------------------------------------------------------------------
\brief initialization of sub-problem level set

<pre>                                                            irhan 05/04
initialization of sub-problem level set (levelset)
</pre>

*----------------------------------------------------------------------*/
void ls_levelset_init_levelset()
{
#ifdef DEBUG
  dstrc_enter("ls_levelset_init_levelset");
#endif
/*----------------------------------------------------------------------*/

  /* initialise level set field */
  if (genprob.restart>0)
  {
    /* turn on restart flag */
    lsdyn->resstep = genprob.restart;
    lsdyn->init = 2;
  }
  else if (genprob.restart==0)
  {
    lsdyn->init = 3;
  }
  ls_init(lsetfield,lsdyn,5);

  /* init the element calculating routines */
  *action = calc_ls_init;
  calinit(lsetfield,actpart,action,&container);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_levelset_init_levelset */



/*!----------------------------------------------------------------------
\brief initialization of sub-problem level set

<pre>                                                            irhan 05/04
initialization of sub-problem level set (solver)
</pre>

*----------------------------------------------------------------------*/
void ls_levelset_init_solver()
{
  INT     i;

#ifdef DEBUG
  dstrc_enter("ls_levelset_init_solver");
#endif
/*----------------------------------------------------------------------*/

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
  /* allocate 1 dist. vector 'rhs' */
  actsolv->nrhs = 1;
  solserv_create_vec(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");
  solserv_zero_vec(&(actsolv->rhs[0]));
  /* allocate 1 dist. solution */
  actsolv->nsol= 1;
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
  /* init the dirichlet conditions */
  ls_initdirich(lsetfield);
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
  init_assembly(actpart,actsolv,actintra,lsetfield,actsysarray,0);
  /* output to the screen */
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
  for (i=0;i<par.nprocs;i++)
    if (par.myrank==i)
      printf("PROC  %3d | FIELD LEVELSET     | number of equations      : %10d \n",
             par.myrank,numeq);
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
  if (par.myrank==0)
    printf("          | FIELD LEVELSET     | total number of equations: %10d \n",numeq_total);
  if (par.myrank==0) printf("\n\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_levelset_init_solver */



/*!----------------------------------------------------------------------
\brief solution of sub-problem level set

<pre>                                                            irhan 05/04
solution of sub-problem level set
</pre>

*----------------------------------------------------------------------*/
void ls_levelset_solv()
{
  INT        itnum;
  INT        converged;
  INT        PLACE;

#ifdef DEBUG
  dstrc_enter("ls_levelset_solv");
#endif
/*----------------------------------------------------------------------*/

  /* S O L U T I O N    P H A S E */

  /*
   * nodal solution history levelset field =>
   * sol_increment[0][j] ... solution at time (n)
   * sol_increment[1][j] ... solution at time (n+1)
   */

  /*
    there are only procs allowed in here, that belong to the levelset
    intracommunicator (in case of nonlinear levelset dyn., this should be all)
  */
  if (actintra->intra_fieldtyp!=levelset) dserror("fieldtyp!=levelset");
  /* output to the screen */
  if (par.myrank==0) ls_algout();
  /* set dirichlet boundary conditions for timestep */
  ls_setdirich(lsetfield,1);
  /* initialise timerhs */
  amzero(&ftimerhs_a);

  if (lsdyn->itchk!=0 && par.myrank==0)
  {
    printf("|------------------------------------------------|\n");
    printf("|- step/max -|-  tol     [norm]  -|- lset error -|\n");
  }
  itnum=1;
  lsdyn->nif = 1;
  /* N O N L I N E A R   I T E R A T I O N */
 nonlniter:
  /* intitialise global matrix and global rhs */
  solserv_zero_vec(&(actsolv->rhs[0]));
  solserv_zero_vec(&(actsolv->sol[0]));
  solserv_zero_mat(
    actintra,&(actsolv->sysarray[actsysarray]),
    &(actsolv->sysarray_typ[actsysarray])
    );
  /* initialise iterations-rhs */
  amzero(&fiterhs_a);
  /* form incremental matrices, residual and element forces */
  *action = calc_ls;
  t1=ds_cputime();
  container.dvec = NULL;
  container.ftimerhs = ftimerhs;
  container.fiterhs = fiterhs;
  container.global_numeq = numeq_total;
  container.nii = lsdyn->nii;
  container.nif = lsdyn->nif;
  container.kstep = 0;
  container.is_relax = 0;
  calelm(
    lsetfield,actsolv,actpart,actintra,actsysarray,actsysarray,
    &container,action
    );
  lsdyn->nif = 0;
  te=ds_cputime()-t1;
  tes+=te;
  /*
   * build the actual rhs-vector =>
   * rhs = ftimerhs + fiterhs
   */
  /* add time-rhs => */
  assemble_vec(
    actintra,
    &(actsolv->sysarray_typ[actsysarray]),
    &(actsolv->sysarray[actsysarray]),
         &(actsolv->rhs[0]),
    ftimerhs,
    1.0
    );
  /* add iteration-rhs => */
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
  /* set place */
  if (lsdyn->lsdata->ls2_calc_flag==ls2_advection ||
      lsdyn->lsdata->ls2_calc_flag==ls2_reinitialization)
  {
    PLACE = 1;
  }
  else if (lsdyn->lsdata->ls2_calc_flag==ls2_gradient)
  {
    PLACE = 2 + lsdyn->lsdata->gradient_direction;
  }
  /* return solution to the nodes and calculate the convergence ratios */
  ls_result_incre(
    lsetfield,actintra,&(actsolv->sol[0]),PLACE,
    &(actsolv->sysarray[actsysarray]),
    &(actsolv->sysarray_typ[actsysarray]),
    &lrat,lsdyn
    );
  /* iteration convergence check */
  converged = ls_convcheck(lsdyn,lrat,itnum,te,ts);
  /* check if nonlinear iteration has to be finished */
  if (converged==0)
  {
    itnum++;
    goto nonlniter;
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_levelset_solv */



/*!----------------------------------------------------------------------
\brief finalization of sub-problem level set

<pre>                                                            irhan 05/04
finalization of sub-problem level set
</pre>

*----------------------------------------------------------------------*/
void ls_levelset_fina()
{
#ifdef DEBUG
  dstrc_enter("ls_levelset_fina");
#endif
/*----------------------------------------------------------------------*/

  /* F I N A L I Z E */
  if (lsdyn->lsdata->ls2_calc_flag == ls2_advection)
  {
    /* copy solution from sol_increment[1][j] to sol_increment[0][j] */
    solserv_sol_copy(lsetfield,0,node_array_sol_increment,node_array_sol_increment,1,0);
  }
  else if (lsdyn->lsdata->ls2_calc_flag == ls2_reinitialization)
  {
    /* copy solution from sol_increment[1][j] to sol_increment[0][j] */
    solserv_sol_copy_reinit(lsetfield,0,node_array_sol_increment,node_array_sol_increment,1,0);
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_levelset_fina */



/*!----------------------------------------------------------------------
\brief end of sub-problem level set

<pre>                                                            irhan 05/04
end of sub-problem level set
</pre>

*----------------------------------------------------------------------*/
void ls_levelset_clea()
{
  INT     i;
#ifdef DEBUG
  dstrc_enter("ls_levelset_clea");
#endif
/*----------------------------------------------------------------------*/

  /* C L E A N I N G   U P   P H A S E */
  /*
   * there are only procs allowed in here, that belong to the levelset
   * intracommunicator (in case of nonlinear levelset. dyn., this should be all)
   */
  if (actintra->intra_fieldtyp!=levelset) dserror("fieldtyp!=levelset");
  /* print total CPU-time to the screen */
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
  for (i=0;i<par.nprocs;i++)
  {
    if (par.myrank==i)
    {
      printf("\n");
      printf("PROC  %3d | FIELD LEVELSET     | total time element for calculations: %10.3E \n",
             par.myrank,tes);
      printf("PROC  %3d | FIELD LEVELSET     | total time for solver              : %10.3E \n",
             par.myrank,tss);
    }
  }
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
  /* tidy up */
  amdel(&ftimerhs_a);
  amdel(&fiterhs_a);
  solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
  solserv_del_vec(&(actsolv->sol),actsolv->nsol);
#ifndef PARALLEL
  CCAFREE(actintra);
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_levelset_clea */



/*!----------------------------------------------------------------------
\brief write restart data for levelset

<pre>                                                            irhan 05/04
write restart data for levelset
</pre>

*----------------------------------------------------------------------*/
void ls_levelset_write_restart()
{
#ifdef DEBUG
  dstrc_enter("ls_levelset_write_restart");
#endif
/*----------------------------------------------------------------------*/

  /* write restart data into pss file */
  restart_write_levelset(
    lsdyn,lsetfield,actpart,
    actintra,action,&container
    );

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_levelset_write_restart */



/*!----------------------------------------------------------------------
\brief read restart data for levelset

<pre>                                                            irhan 05/04
read restart data for levelset
</pre>

*----------------------------------------------------------------------*/
void ls_levelset_read_restart()
{
#ifdef DEBUG
  dstrc_enter("ls_levelset_read_restart");
#endif
/*----------------------------------------------------------------------*/

  /* read restart data from pss file */
  restart_read_levelset(
    lsdyn->resstep,lsdyn,lsetfield,actpart,
    actintra,action,&container
    );

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_levelset_read_restart */



/*!----------------------------------------------------------------------
\brief reset level set solution before reinitialization

<pre>                                                            irhan 08/04
write reset level set solution before reinitialization
</pre>

*----------------------------------------------------------------------*/
void ls_levelset_reset_sol_increment()
{
  INT       i;
  NODE     *actnode;

#ifdef DEBUG
  dstrc_enter("ls_levelset_reset_sol_increment");
#endif
/*----------------------------------------------------------------------*/

  for (i=0; i<lsetfield->dis[0].numnp; i++)
  {
    /* access to the node */
    actnode = &(lsetfield->dis[0].node[i]);
    if (actnode->gnode->is_node_active==0)
    {
      /* reset solution */
      amzero(&(actnode->sol_increment));
    }
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of ls_levelset_reset_sol_increment */
/*! @} (documentation module close)*/
#endif
