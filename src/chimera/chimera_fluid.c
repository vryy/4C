#ifdef D_CHIMERA
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../fluid_full/fluid_prototypes.h"
#include "chimera.h"



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA           *alldyn;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB    genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD     *field;
/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR        par;
extern struct _PARTITION        *partition;
/*----------------------------------------------------------------------*
 |                                                       m.gee 02/02    |
 | number of load curves numcurve                                       |
 | vector of structures of curves                                       |
 | defined in input_curves.c                                            |
 | INT                   numcurve;                                      |
 | struct _CURVE      *curve;                                           |
 *----------------------------------------------------------------------*/
extern INT                numcurve;
extern struct _CURVE     *curve;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS          ioflags;
/*----------------------------------------------------------------------*
  | global variable *solv, vector of lenght numfld of structures SOLVAR  |
  | defined in solver_control.c                                          |
  |                                                                      |
  |                                                       m.gee 11/00    |
  *----------------------------------------------------------------------*/
extern struct _SOLVAR           *solv;
/*-----------------------------------------------------------------------*
  | enum _CALC_ACTION                                      m.gee 1/02    |
  | command passed from control routine to the element level             |
  | to tell element routines what to do                                  |
  | defined globally in global_calelm.c                                  |
  *----------------------------------------------------------------------*/
extern enum   _CALC_ACTION       calc_action[MAXFIELD];
extern struct _CHIMERA_DATA     *chm_data;





static FLUID_DYNAMIC  *fdyn;               /* dynamic control for fluid field  */
static DOUBLE          grat;               /* convergence ratios               */
static FIELD          *actfield;           /* pointer to fluid field struct    */
static SOLVAR         *actsolv;            /* pointer to active sol. structure */
static PARTITION      *actpart;            /* pointer to active partition      */
static INTRA          *actintra;           /* pointer to active intra-communic.*/
static DISCRET        *actdis;             /* pointer to active dis. structure */

static CONTAINER       container;          /* variables for calelm             */
static FLUID_STRESS    str;                /* variable for stress calculation  */
static CALC_ACTION    *action;             /* pointer to the cal_action enum   */


static ARRAY           ftimerhs_a;
static DOUBLE         *ftimerhs;	   /* time - RHS		       */
static ARRAY           fiterhs_a;
static DOUBLE         *fiterhs;	           /* iteration - RHS  	       	       */

static INT            *numeq;              /* number of equations on this proc */
static INT            *numeq_total;        /* total number of equations        */
static INT             neq_max;            /* maximum number of equations
                                            * over the whole discretizations   */

static DOUBLE          tes=ZERO;           /*				       */
static DOUBLE          tss=ZERO;           /*				       */


static FILE           *velx;
static FILE           *vely;
static FILE           *pres;



void chimera_initialize()
{
  INT     i;
  INT     init;
  INT     ctrl;
  INT     neq;
  INT     neq_total;
  INT     option=0;

#ifdef DEBUG
  dstrc_enter("chimera_initialize");
#endif
/*----------------------------------------------------------------------*/

  /*plaausibility check */
  dsassert(genprob.numfld==1,"**ERROR** Number of fields is to be ONE!\n");
  dsassert(genprob.numff==0,"**ERROR** Field number is to be ZERO!\n");

  /* access to the fluid field */
  actfield = &(field[genprob.numff]);
  dsassert(actfield->fieldtyp==fluid,"FIELD 0 has to be fluid\n");
  dsassert(actfield->ndis==2,"**ERROR** Number of discretization is to be TWO!\n");

  /* initialise the corresponding fluid dynamic struct */
  fdyn = alldyn[genprob.numff].fdyn;
  grat = ZERO;

  /* set some pointers */
  actpart            = &(partition[genprob.numff]);
  actsolv            = &(solv[genprob.numff]);
  action             = &(calc_action[genprob.numff]);
  str                = str_none;

  /*---------------- if we are not parallel, we have to allocate an alibi
    -------------------------------------- intra-communicator structure */
#ifdef PARALLEL
  actintra    = &(par.intra[genprob.numff]);
#else
  actintra    = (INTRA*)CCACALLOC(1,sizeof(INTRA));
  if (!actintra) dserror("Allocation of INTRA failed");
  actintra->intra_fieldtyp = fluid;
  actintra->intra_rank     = 0;
  actintra->intra_nprocs   = 1;
#endif

  /*- there are only procs allowed in here, that belong to the fluid -- */
  /* intracomm. (in case of nonlin. fluid. dyn., this should be all)    */
  if (actintra->intra_fieldtyp!=fluid) dserror("**ERROR** field != fluid");

  dsassert(actsolv->nsysarray==2, "two discretizations expected");

  /* initialize equation numbers */
  numeq       = (INT*)CCACALLOC(actsolv->nsysarray,sizeof(INT));
  numeq_total = (INT*)CCACALLOC(actsolv->nsysarray,sizeof(INT));

  for (i=0; i<actsolv->nsysarray; i++)
  {
    /* init the dist sparse matrices to zero */
    solserv_zero_mat(
      actintra,
      &(actsolv->sysarray[i]),
      &(actsolv->sysarray_typ[i])
      );
    /* get global and local number of equations */
    solserv_getmatdims(
      &(actsolv->sysarray[i]),
      actsolv->sysarray_typ[i],
      &neq,
      &neq_total
      );
    numeq[i] = neq;
    numeq_total[i] = neq_total;
  }

  /* allocate 2 dist. vector 'rhs' */
  actsolv->nrhs = actsolv->nsysarray;

  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */

  solserv_create_vec2(&(actsolv->rhs),actsolv->nrhs,numeq_total,numeq,"DV");

  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */

  /* initialize 'rhs' */
  for (i=0; i<actsolv->nsysarray; i++)
    solserv_zero_vec(&(actsolv->rhs[i]));

  /* allocate 2 dist. vector 'lhs' */
  actsolv->nsol = actsolv->nsysarray;
  solserv_create_vec2(&(actsolv->sol),actsolv->nsol,numeq_total,numeq,"DV");

  /* initialize 'lhs' */
  for (i=0; i<actsolv->nsysarray; i++)
    solserv_zero_vec(&(actsolv->sol[i]));

  /* find maximum number of equations */
  neq_max = 0;
  ctrl    = 0;
  for (i=0; i<actfield->ndis; i++)
  {
    ctrl = actfield->dis[i].numeq;
    if (ctrl<neq_max) continue;
    neq_max = ctrl;
  }

  /*------------- allocate one redundant vector ftimerhs of full lenght */
  /*        it is used by the element routines to assemble the Time RHS */
  ftimerhs = amdef("ftimerhs",&ftimerhs_a,neq_max,1,"DV");

  /*-------------  allocate one redundant vector fiterhs of full lenght */
  /*   it is used by the element routines to assemble the Iteration RHS */
  fiterhs = amdef("fiterhs",&fiterhs_a,neq_max,1,"DV");

  /* loop over the discretizations */
  for (i=0; i<actfield->ndis; i++)
  {
    actdis = &(actfield->dis[i]);
  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */

    /* initialize the fluid discretization */
    /*fluid_init_discretization(actpart,actintra,actdis,action,&container,SIX,str);*/
    fluid_init(actpart,actintra,actfield, i,action,&container,SIX,str);

    /* init the dirichlet-conditions for discretization */
    fluid_initdirich_discretization(actdis);

  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */

  }

  /*--------------------------------- initialize solver on all matrices */
  /*
    NOTE: solver init phase has to be called with each matrix one wants to
    solve with. Solver init phase has to be called with all matrices
    one wants to do matrix-vector products and matrix scalar products.
    This is not needed by all solver libraries, but the solver-init
    phase is cheap in computation (can be costly in memory)
  */
  /* initialize solver */
  init=1;
  for (i=0; i<actsolv->nsysarray; i++)
    solver_control(
      actsolv, actintra,
      &(actsolv->sysarray_typ[i]),
      &(actsolv->sysarray[i]),
      &(actsolv->sol[i]),
      &(actsolv->rhs[i]),
      init
      );

  /* init the assembly for stiffness */
  for (i=0; i<actsolv->nsysarray; i++)
  {
    init_assembly(actpart,actsolv,actintra,actfield,i,i);
  }

  /* allocate fluid integration data ---*/
  alldyn[genprob.numff].fdyn->data = (FLUID_DATA*)CCACALLOC(1,sizeof(FLUID_DATA));

  /* init the element calculating routines */
  *action = calc_fluid_init;
  calinit(actfield,actpart,action,&container);

  /* output to the screen */
  for (i=0; i<actsolv->nsysarray; i++)
  {
    if (par.myrank==0) printf("\n\n");
#ifdef PARALLEL
    MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
    for (i=0; i<par.nprocs; i++)
      if (par.myrank==i)
        printf("PROC  %3d | FIELD DISCRETIZATION[%1d] | number of equations      : %10d\n",
               par.myrank,i,numeq[i]);
#ifdef PARALLEL
    MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
    if (par.myrank==0)
      printf("PROC  %3d | FIELD DISCRETIZATION[%1d] | number of equations      : %10d\n",
             par.myrank,i,numeq_total[i]);
    if (par.myrank==0) printf("\n\n");
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of chimera_initialize */



void chimera_solve(
  INT      numdis
  )
{
  INT        i,j;
  INT        itnum;                /* counter for nonlinear iteration  */
  INT        iststep=0;            /* counter for time integration     */
  INT        init;

  INT        steady=0;             /* flag for steady state            */
  INT        converged=0;          /* convergence flag                 */

  DOUBLE     t1,ts,te;

  DOUBLE     vrat=0.0;
  DOUBLE     prat=0.0;             /* convergence ratios               */

  DISCRET    *actdis;
  NODE       *actnode;

#ifdef DEBUG
  dstrc_enter("chimera_solve");
#endif
/*----------------------------------------------------------------------*/

/*======================================================================*
 |                     S O L U T I O N    P H A S E                     |
 *======================================================================*/
/* nodal solution history fluid field:                                  *
 * sol[0][j]           ... initial data 				*
 * sol[1...actpos][j]  ... solution for visualisation (real pressure)	*
 * sol_increment[0][j] ... solution at time 0        			*
 * sol_increment[1][j] ... temporaer                                    *
 * sol_increment[2][j] ... solution at time (n+g)			*
 * sol_increment[3][j] ... solution at subdomain iteration (m+1)        *
 * sol_increment[4][j] ... solution at subdomain iteration (m)          *
 * sol_increment[5][j] ... solution at time (n)                         *
 *======================================================================*/

  /*- there are only procs allowed in here, that belong to the fluid ---*/
  /* intracomm. (in case of nonlin. fluid. dyn., this should be all)    */
  if (actintra->intra_fieldtyp!=fluid) dserror("**ERROR** field != fluid");

  /* access to active discretization */
  actdis = &(actfield->dis[numdis]);

  /* calculate constants for time algorithm */
  fluid_tcons();

  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */

  /* set dirichlet boundary conditions for  timestep */
  fluid_setdirich_discretization(actfield, actdis, 3);

  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */

/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/
  if (chm_data->interact_dis_on_off==1)
  {
    /* enforce mass conservation */
    chimera_continuity_interpolation(numdis);
  }
/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/

  /* initialise timerhs */
  amzero(&ftimerhs_a);

  if (fdyn->itnorm!=fncc_no && par.myrank==0)
  {
    printf("----------------------------------------------------------------------------- \n");
    printf("|- step/max -|-  tol     [norm] -|- vel. error -|- pre. error -|\n");
  }

  itnum=1;
/*=======================================================================*
 |                N O N L I N E A R   I T E R A T I O N                  |
 *=======================================================================*/
 nonlniter:

  /* calculate constants for nonlinear iteration */
  fluid_icons(itnum);

  /* intitialise global matrix and global rhs */
  solserv_zero_vec(&(actsolv->rhs[numdis]));
  solserv_zero_mat(
    actintra,&(actsolv->sysarray[numdis]),
    &(actsolv->sysarray_typ[numdis])
    );

  /* initialize iterations rhs */
  amzero(&fiterhs_a);

  /* form incremental matrices, residual and element forces */
  *action = calc_fluid;
  t1      = ds_cputime();

/********************** BE CAREFUL ! ********************************   */
/********************** BE CAREFUL ! ********************************   */
/********************** BE CAREFUL ! ********************************   */
/********************** BE CAREFUL ! ********************************   */
/********************** BE CAREFUL ! ********************************   */
  container.fieldtyp     = actfield->fieldtyp;
  container.actndis      = numdis; /*       => processing disc. numdis! */
  container.dvec         = NULL;
  container.ftimerhs     = ftimerhs;
  container.fiterhs      = fiterhs;
  container.global_numeq = numeq_total[numdis];
  container.nii          = fdyn->nii;
  container.nif          = fdyn->nif;
  container.nim          = fdyn->nim;
  container.kstep        = 0;
  container.is_relax     = 0;
  container.turbu    = fdyn->turbu;
/********************** BE CAREFUL ! ********************************   */
/********************** BE CAREFUL ! ********************************   */
/********************** BE CAREFUL ! ********************************   */
/********************** BE CAREFUL ! ********************************   */
/********************** BE CAREFUL ! ********************************   */

  /* compute element contributions */
  calelm(
    actfield,
    actsolv,
    actpart,
    actintra,
    numdis,        /*        => processing discretization numdis! */
    -1,
    &container,
    action
    );

  te   = ds_cputime()-t1;
  tes += te;

/*----------------------------------------------------------------------*
  | build the actual rhs-vector:                                        |
  |        rhs = ftimerhs + fiterhs                                     |
  *---------------------------------------------------------------------*/
  /* add time right hand side */
  assemble_vec(
    actintra,
    &(actsolv->sysarray_typ[numdis]),
    &(actsolv->sysarray[numdis]),
    &(actsolv->rhs[numdis]),
    ftimerhs,
    1.0
    );

  /* add iteration right hand side */
  assemble_vec(
    actintra,
    &(actsolv->sysarray_typ[numdis]),
    &(actsolv->sysarray[numdis]),
    &(actsolv->rhs[numdis]),
    fiterhs,
    1.0
    );

  init = 0;
  t1   = ds_cputime();
  /* solve system */
  solver_control(
    actsolv, actintra,
    &(actsolv->sysarray_typ[numdis]),
    &(actsolv->sysarray[numdis]),
    &(actsolv->sol[numdis]),
    &(actsolv->rhs[numdis]),
    init
    );
  ts   = ds_cputime()-t1;
  tss += ts;

  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */

  /* return solution to the nodes and calculate the convergence ratios */
  fluid_result_incre(
    actfield,
    numdis,
    actintra,
    &(actsolv->sol[numdis]),
    3,
    &(actsolv->sysarray[numdis]),
    &(actsolv->sysarray_typ[numdis]),
    &vrat,&prat,&grat
    );

  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */
  /* THE FOLLOWING SUBROUTINE IS NEW */

  /* iteration convergence check */
  converged = fluid_convcheck(vrat,prat,grat,itnum,te,ts);

  /* check if nonlinear iteration has to be finished */
  if (converged==0)
  {
    itnum++;
    goto nonlniter;
  }
/*-----------------------------------------------------------------------*
  | -->  end of nonlinear iteration                                      |
  *----------------------------------------------------------------------*/
  /* steady state check */
  if (fdyn->stchk==iststep)
  {
    iststep = 0;
    steady  = fluid_steadycheck(actfield,numeq_total[numdis]);
  }

  if (chm_data->interact_dis_on_off==1)
  {
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/
/****************************CHIMERA REGION START************************/

    /* copy solution from sol_increment[3][j] to sol_increment[4][j] */
    solserv_sol_copy(actfield,numdis,1,1,3,4);

/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/
/****************************CHIMERA REGION END**************************/
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of chimera_solve */



void chimera_finalize(
  INT      disnum
  )
{

#ifdef DEBUG
  dstrc_enter("chimera_finalize");
#endif
/*----------------------------------------------------------------------*/

  /* copy solution from sol_increment[3][j] to sol_increment[1[j] */
  solserv_sol_copy(actfield,disnum,1,1,3,1);

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of chimera_finalize */



void chimera_clean(
  )
{
  INT     i;

#ifdef DEBUG
  dstrc_enter("chimera_clean");
#endif
/*----------------------------------------------------------------------*/

  /*- there are only procs allowed in here, that belong to the fluid ---*/
  /* intracomm. (in case of nonlin. fluid. dyn., this should be all)    */
  if (actintra->intra_fieldtyp!=fluid) dserror("**ERROR** field != fluid");
  /*-------------------------------- print total CPU-time to the screen */
#ifdef PARALLEL
  MPI_Barrier(actintra->MPI_INTRA_COMM);
#endif
  for (i=0; i<par.nprocs; i++)
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

  /* tidy up */
  amdel(&ftimerhs_a);
  amdel(&fiterhs_a);
  solserv_del_vec(&(actsolv->rhs),actsolv->nrhs);
  solserv_del_vec(&(actsolv->sol),actsolv->nsol);

  if (ioflags.chimera_to_matlab==1)
  {
    fclose(velx);
    fclose(vely);
    fclose(pres);
  }

#ifndef PARALLEL
  CCAFREE(actintra);
#endif

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of chimera_clean */



void chimera_to_matlab()
{
  INT          i,j,k;
  ELEMENT     *actele;
  NODE        *actnode;
  FILE        *f1;

#ifdef DEBUG
  dstrc_enter("chimera_to_matlab");
#endif
/*----------------------------------------------------------------------*/

  /* I need the connectivity */
  f1 = fopen("src/chimera/to_matlab/chimera_to_matlab/chimera_to_matlab_connectivity","w");

  /* write connectivity to file */
  for (i=0; i<actfield->ndis; i++)
  {
    for (j=0; j<actfield->dis[i].numele; j++)
    {
      actele = &(actfield->dis[i].element[j]);
      for (k=0; k<actele->numnp; k++)
      {
        fprintf(f1,"%5d ",actele->lm[k]+1);
      }
      fprintf(f1,"\n");
    }
    fprintf(f1,"\n\n\n\n");
  }
  fclose(f1);

  /* I need the node coordinates */
  f1 = fopen("src/chimera/to_matlab/chimera_to_matlab/chimera_to_matlab_nodecoordinates","w");

  /* write node coordinates to file */
  for (i=0; i<actfield->ndis; i++)
  {
    for (j=0; j<actfield->dis[i].numnp; j++)
    {
      actnode = &(actfield->dis[i].node[j]);
      for (k=0; k<genprob.ndim; k++)
      {
        fprintf(f1,"%10.5E ",actnode->x[k]);
      }
      fprintf(f1,"\n");
    }
    fprintf(f1,"\n\n\n\n");
  }
  fclose(f1);

  f1 = fopen("src/chimera/to_matlab/chimera_to_matlab/chimera_to_matlab_nodeIds","w");

  /* write node Ids to file */
  for (i=0; i<actfield->ndis; i++)
  {
    for (j=0; j<actfield->dis[i].numnp; j++)
    {
      actnode = &(actfield->dis[i].node[j]);
      fprintf(f1,"%5d\n",actnode->Id+1);
    }
    fprintf(f1,"\n\n\n\n");
  }
  fclose(f1);

  /* I need some general information */
  f1 = fopen("src/chimera/to_matlab/chimera_to_matlab/chimera_to_matlab_general","w");

  fprintf(f1,"%5d\n",genprob.ndim);
  fprintf(f1,"%5d\n",actfield->ndis);
  fprintf(f1,"%5d\n",fdyn->nstep);
  /* write some general information to file */
  for (i=0; i<actfield->ndis; i++)
  {
    fprintf(f1,"%5d\n",actfield->dis[i].numnp);
    fprintf(f1,"%5d\n",actfield->dis[i].numele);
    fprintf(f1,"%5d\n",actfield->dis[i].element[0].numnp);
    fprintf(f1,"\n\n\n\n");
  }
  fclose(f1);

  /* open files to write solution history */
  velx = fopen("src/chimera/to_matlab/chimera_to_matlab/chimera_to_matlab_velx","w");
  vely = fopen("src/chimera/to_matlab/chimera_to_matlab/chimera_to_matlab_vely","w");
  pres = fopen("src/chimera/to_matlab/chimera_to_matlab/chimera_to_matlab_pres","w");

  /* slow, check... */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of chimera_to_matlab */



void chimera_write_soln_to_matlab()
{
  INT          i,j;
  NODE        *actnode;

#ifdef DEBUG
  dstrc_enter("chimera_write");
#endif
/*----------------------------------------------------------------------*/

  /* write solution to file */
  for (i=0; i<actfield->ndis; i++)
  {
    for (j=0; j<actfield->dis[i].numnp; j++)
    {
      actnode = &(actfield->dis[i].node[j]);
      fprintf(velx,"%10.5E\n",actnode->sol_increment.a.da[3][0]);
      fprintf(vely,"%10.5E\n",actnode->sol_increment.a.da[3][1]);
      fprintf(pres,"%10.5E\n",actnode->sol_increment.a.da[3][2]);
    }
    fprintf(velx,"\n");
    fprintf(vely,"\n");
    fprintf(pres,"\n");
  }

/*----------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of chimera_write */
#endif /* D_CHIMERA */
