/*----------------------------------------------------------------------*/
/*!
\file
\brief TSI - quasi static solution of thermal field

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 03/06
*/


/*----------------------------------------------------------------------*/
#ifdef D_TSI

/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#ifdef BINIO
#include "../io/io.h"
#endif
#include "tsi_prototypes.h"

/*----------------------------------------------------------------------*/
/*!
\brief File pointers
      This structure is defined in input_control_global.c
      and the type is in standardtypes.h
      It holds all file pointers and some variables needed for the FRSYSTEM
\auther bborn
\date 03/06
*/
extern FILES allfiles;


/*----------------------------------------------------------------------*/
/*!
\brief General problem data
       struct _GENPROB       genprob; defined in global_control.c 
\author bborn
\date 03/06
*/
extern GENPROB genprob;


/*----------------------------------------------------------------------*/
/*!
\brief Fields
       vector of numfld FIELDs, defined in global_control.c
\author bborn
\date 03/06
*/
extern FIELD *field;


/*----------------------------------------------------------------------*/
/*!
\brief Global nodal solution vectors in solver-specific format
       global variable *solv, vector of lenght numfld of structures SOLVAR
       defined in solver_control.c
\author bborn
\date 03/06
*/
extern SOLVAR *solv;


/*----------------------------------------------------------------------*/
/*!
\brief One proc's info about his partition
       - the partition of one proc (all discretizations)
       - the type is in partition.h
\author bborn
\date 03/06
*/
extern PARTITION *partition;


/*----------------------------------------------------------------------*/
/*!
\brief Input/output control flags
       structure of flags to control output defined in out_global.c
\author bborn
\date 03/06
*/
extern IO_FLAGS ioflags;


/*----------------------------------------------------------------------*/
/*!
\brief Rank and communicators (Parallelism!)
       This structure struct _PAR par; is defined in main_ccarat.c
       and the type is in partition.h
\author born
\date 03/06
*/
extern PAR par;


/*----------------------------------------------------------------------*/
/*!
\brief Dynamic control
       pointer to allocate dynamic variables if needed
       dedfined in global_control.c
\author bborn
\date 03/06
*/
extern ALLDYNA* alldyn;


/*----------------------------------------------------------------------*/
/*!
\brief Load curve, load factor function
       number of load curves numcurve
       vector of structures of curves
       defined in input_curves.c
\author bborn
\date 03/06
*/
extern INT numcurve;
extern CURVE* curve;


/*----------------------------------------------------------------------*/
/*!
\brief CALC_ACTIONs
       command passed from control routine to the element level
       to tell element routines what to do
       defined globally in global_calelm.c
\author bborn
\date 03/06
*/
extern CALC_ACTION calc_action[MAXFIELD];


/*---------------------------------------------------------------------*/
/*!
\brief Actual (current) time globally given
\author bborn
\date 03/06
*/
DOUBLE acttime;


/*======================================================================*/
/*!
\brief Initialise thermal field
\param   actpart      PARTITION*   (i)   partition
\param   actintra     INTRA*       (i)   intra-communicator
\param   actfield     FIELD        (i)   thermal field
\param   disnum       INT          (i)   discretisation index
\param   actsolv      SOLVAR*      (io)  solution variables
\param   numeq        INT*         (o)   number of equations
\param   numeq_total  INT*         (o)   total number of equations
\param   actsysarray  INT*         (o)   tangent array index
\param   container    CONTAINER*   (i)   container
\param   dirich_a     ARRAY*       (o)   Dirich. RHS contribution
\author bborn
\date 05/07
*/
void tsi_th_ost_init(PARTITION* actpart,
                     INTRA* actintra,
                     FIELD* actfield,
                     INT disnum,
                     ARRAY_POSITION* ipos,
                     ARRAY_POSITION_SOL* isol,
                     SOLVAR* actsolv,
                     INT* numeq,
                     INT* numeq_total,
                     INT* stiff_array,
                     INT* mass_array,
                     THERM_DYNAMIC* actdyn,
                     CONTAINER* container,
                     DIST_VECTOR** tem,
                     DIST_VECTOR** fext,
                     DIST_VECTOR** fextn,
                     DIST_VECTOR** fint,
                     ARRAY* intforce_a,
                     ARRAY* dirich_a)
{
  INT numtf = genprob.numtf;
  DOUBLE* dirich;
  CALC_ACTION* action = &(calc_action[numtf]);  /* thermal cal_action */
  INT i;  /* index */
  INT init;  /* init solver flag */

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_th_ost_init");
#endif

  /*--------------------------------------------------------------------*/
  /* distributed system matrix, which is used for solving */
  if (actsolv->nsysarray == 1)
  {
    *stiff_array = 0;
    *mass_array = 1;
    actsolv->nsysarray = 2;
  }
  else
  {
    dserror("More than 1 system arrays (actsolv->nsysarray)!");
  }
  /* allocate mass matrix */
  actsolv->sysarray_typ 
    = (SPARSE_TYP*) CCAREALLOC(actsolv->sysarray_typ,
                               actsolv->nsysarray*sizeof(SPARSE_TYP));
  if (!actsolv->sysarray_typ) 
  {
    dserror("Allocation of memory failed");
  }
  actsolv->sysarray
    = (SPARSE_ARRAY*) CCAREALLOC(actsolv->sysarray,
                                 actsolv->nsysarray*sizeof(SPARSE_ARRAY));
  if (!actsolv->sysarray) 
  {
    dserror("Allocation of memory failed");
  }
  /* copy the matrices sparsity mask from stiff_array to mass_array */
  solserv_alloc_cp_sparsemask(actintra,
                              &(actsolv->sysarray_typ[*stiff_array]),
                              &(actsolv->sysarray[*stiff_array]),
                              &(actsolv->sysarray_typ[*mass_array]),
                              &(actsolv->sysarray[*mass_array]));
  /* init the dist sparse matrices to zero */
  for (i=0; i<actsolv->nsysarray; i++)
  {
    solserv_zero_mat(actintra,
                     &(actsolv->sysarray[i]),
                     &(actsolv->sysarray_typ[i]));
  }

  /*--------------------------------------------------------------------*/
  /* get global and local number of equations */
  solserv_getmatdims(&(actsolv->sysarray[*stiff_array]),
                     actsolv->sysarray_typ[*stiff_array],
                     numeq, numeq_total);

  /*--------------------------------------------------------------------*/
  /* number of rhs and solution vectors */
  actsolv->nrhs = 2;  /* WHY 2? */
  solserv_create_vec(&(actsolv->rhs), actsolv->nrhs,
                     *numeq_total, *numeq, "DV");
  /* init the created dist. vectors to zero */
  for (i=0; i<actsolv->nrhs; i++)
  {
    solserv_zero_vec(&(actsolv->rhs[i]));
  }

  /*--------------------------------------------------------------------*/
  /* number of rhs and solution vectors */
  actsolv->nsol = 2;  /* WHY 2? */
  solserv_create_vec(&(actsolv->sol), actsolv->nsol, 
                     *numeq_total, *numeq, "DV");
  /* init the created dist. vectors to zero */
  for (i=0; i<actsolv->nsol; i++)
  {
    solserv_zero_vec(&(actsolv->sol[i]));
  }
  
  /*--------------------------------------------------------------------*/
  /* create 1 vector of full length for Dirichlet part of RHS */
  dirich = amdef("dirich_t", dirich_a, *numeq_total, 1, "DV");
  amzero(dirich_a);

  /*--------------------------------------------------------------------*/
  /* global temperature vectors */
  solserv_create_vec(tem, 1, *numeq_total, *numeq, "DV");
  /* solserv_zero_vec(tem[0]); */
  /* set initial temperature in field */
  {
    INT idof;
    for (idof=0; idof<*numeq; idof++)
    {
      tem[0]->vec.a.dv[idof] = actdyn->initmpr;
    }
  }
  /* spawn initial temperature to nodes */
  solserv_sol_zero(actfield, disnum, node_array_sol, isol->tem0);

  /*--------------------------------------------------------------------*/
  /* external heat loads */
  solserv_create_vec(fext, 1, *numeq_total, *numeq, "DV");
  solserv_zero_vec(fext[0]);
  solserv_create_vec(fextn, 1, *numeq_total, *numeq, "DV");
  solserv_zero_vec(fextn[0]);

  /*--------------------------------------------------------------------*/
  /* internal heat force */
  amdef("intforce_t", intforce_a, *numeq_total, 1, "DV");
  amzero(intforce_a);
  solserv_create_vec(fint, 1, *numeq_total, *numeq, "DV");
  solserv_zero_vec(fint[0]);

  /*-------------------------------------------------------------------*/
  /* initialize solver */
  /* NOTE: solver init phase has to be called with each matrix one wants to
   *       solve with. Solver init phase has to be called with all matrices
   *       one wants to do matrix-vector products and matrix scalar products.
   *       This is not needed by all solver libraries, but the solver-init 
   *       phase is cheap in computation (can be costly in memory)
   *       There will be no solver call on mass or damping array. */
  init = 1;
  solver_control(actfield, disnum, actsolv, actintra,
                 &(actsolv->sysarray_typ[*stiff_array]),
                 &(actsolv->sysarray[*stiff_array]),
                 &(actsolv->sol[*stiff_array]),
                 &(actsolv->rhs[*stiff_array]),
                 init);
  solver_control(actfield, disnum, actsolv, actintra,
                 &(actsolv->sysarray_typ[*mass_array]),
                 &(actsolv->sysarray[*mass_array]),
                 &(actsolv->sol[*mass_array]),
                 &(actsolv->rhs[*mass_array]),
                 init);

  /*--------------------------------------------------------------------*/
  /* init the dist sparse matrix to zero */
  /* NOTE: Has to be called after solver_control(init=1) */
/*   solserv_zero_mat(actintra, */
/*                    &(actsolv->sysarray[*stiff_array]), */
/*                    &(actsolv->sysarray_typ[*stiff_array])); */

  /*--------------------------------------------------------------------*/
  /* init the assembly for ONE sparse matrix */
  init_assembly(actpart, actsolv, actintra, actfield, *stiff_array, disnum);
  init_assembly(actpart, actsolv, actintra, actfield, *mass_array, disnum);

  /*--------------------------------------------------------------------*/
  /* init the element calculating routines */
  *action = calc_therm_init;
  calinit(actfield, actpart, action, container);

  /*--------------------------------------------------------------------*/
  /* put a zero to the place ipos->num=12 in node->sol to init the 
   * velocities and accels of prescribed displacements */
  /* HINT: This actually redefines/reallocates/enlarges the sol
   *       array of each structure node to dimension 12x2 (or 12x3)
   *       from originally 1x2 (or 1x3) */
  solserv_sol_zero(actfield, disnum, node_array_sol, ipos->num-1);

  /*--------------------------------------------------------------------*/
#ifdef BINIO
  /* initialize binary output */
  if (ioflags.output_bin == 1)
  {
    dserror("BINIO is not available!");
  }
#endif

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void tsi_th_ost_init() */


/*======================================================================*/
/*!
\brief Equilibrium
\param   actpart      PARTITION*   (i)   partition
\param   actintra     INTRA*       (i)   intra-communicator
\param   actfield     FIELD        (i)   thermal field
\param   disnum       INT          (i)   discretisation index
\param   actsolv      SOLVAR*      (io)  solution variables
\param   isol         ARRAY_POSITION_SOL* (i) indices in NODE sol arrays
\param   numeq        INT          (i)   number of equations
\param   numeq_total  INT          (i)   total number of equations
\param   actsysarray  INT          (i)   tangent array index
\param   actdyn       THERM_DYNAMIC* (i) thermal control
\param   container    CONTAINER*   (i)   container
\param   dirich_a     ARRAY*       (o)   Dirich. RHS contribution
\author bborn
date 05/07
*/
void tsi_th_ost_equi(PARTITION* actpart,
                     INTRA* actintra,
                     FIELD* actfield,
                     INT disnum,
                     ARRAY_POSITION_SOL* isol,
                     SOLVAR* actsolv,
                     INT numeq,
                     INT numeq_total,
                     INT stiff_array,
                     INT mass_array,
                     THERM_DYNAMIC* actdyn,
                     CONTAINER* container,
                     DIST_VECTOR* tem,
                     DIST_VECTOR* fint,
                     DIST_VECTOR* fext,
                     DIST_VECTOR* fextn,
                     ARRAY* intforce_a,
                     ARRAY* dirich_a)
{
  const INT numtf = genprob.numtf;  /* index of thermal field */
  CALC_ACTION* action = &(calc_action[numtf]);  /* thermal cal_action */
  DOUBLE timcur = actdyn->time;  /* current time */
  INT init;  /* init flag for solver */
  DOUBLE* intforce = intforce_a->a.dv;
  DOUBLE dirichfacs[2];  /* factors needed for Dirichlet-part of RHS */
  DOUBLE* dirich = dirich_a->a.dv;
  DOUBLE gamma = actdyn->gamma;
  DOUBLE dt = actdyn->dt;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_th_ost_equi");
#endif

  /*--------------------------------------------------------------------*/
  /* put the scaled prescribed temperatures to the nodes
   * in field sol (1st 0) at place 0 (2nd 0) together with 
   * free temperatures */
  /* HINT: time curve is called indirectly */
  solserv_putdirich_to_dof(actfield, disnum, 
                           node_array_sol, isol->temdn, timcur);

  /*------------------------------------------------------------------*/
  /* initialise tangent and Dirichlet-RHS loads */
  solserv_zero_mat(actintra, 
                   &(actsolv->sysarray[stiff_array]),
                   &(actsolv->sysarray_typ[stiff_array]));
  solserv_zero_mat(actintra, 
                   &(actsolv->sysarray[mass_array]),
                   &(actsolv->sysarray_typ[mass_array]));
  amzero(dirich_a);
  amzero(intforce_a);  /* initialse internal force vector */

  /*--------------------------------------------------------------------*/
  /* call element routines to calculate & assemble tangent matrices */
  *action = calc_therm_tang_instat;
  container->dvec = intforce;
  container->dirich = dirich;  /* Dirichlet forces */
  container->global_numeq = numeq_total;
  container->kstep = 0;  /* WORK */
  container->isdyn = 1;
  container->isoltemdn = isol->temdn;
  dirichfacs[0] = 1.0/dt;  /* dtinv */
  dirichfacs[1] = gamma;  /* 'theta' in one-step-theta */
  container->dirichfacs = dirichfacs;
  calelm(actfield, actsolv, actpart, actintra, stiff_array,
         mass_array, container, action);

  /*--------------------------------------------------------------------*/
  /* store positive internal forces on fint[0]
   * this is the load at t_{n} 'coz temperature has not been updated
   * at the very first call it is zero */
  solserv_zero_vec(&(fint[0]));
  assemble_vec(actintra, 
               &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(fint[0]), intforce, 1.0);

  /*--------------------------------------------------------------------*/
  /* build effective tangent
   *    Keff = 1/dt*M + gamma*K  (is stored on K) */
  solserv_scal_mat(&(actsolv->sysarray_typ[stiff_array]),
                   &(actsolv->sysarray[stiff_array]),
                   gamma);
  solserv_add_mat(actintra,
                  &(actsolv->sysarray_typ[stiff_array]),
                  &(actsolv->sysarray[stiff_array]),
                  &(actsolv->sysarray_typ[mass_array]),
                  &(actsolv->sysarray[mass_array]),
                  1.0/dt);
  
  /*--------------------------------------------------------------------*/
  /* call RHS-routines to assemble RHS */
  solserv_zero_vec(&(fextn[0]));
  *action = calc_therm_heatload;  /* set action before call of calrhs */
  container->kstep = 0;  /* ? */
  container->inherit = 1;  /* ? */
  container->point_neum = 1;  /* ? */
  calrhs(actfield, actsolv, actpart, actintra, stiff_array,
         &(fextn[0]), action, container);
  
  /*--------------------------------------------------------------------*/
  /* build effective RHS
   *    Feff = 1/dt * M . Th_n
   *         - (1-gamma) * Fint_n
   *         + gamma * Fext_{n+1}
   *         + (1-gamma) * Fext_n 
   *         - Fdirichlet */
  solserv_zero_vec(&(actsolv->rhs[stiff_array]));
#if 0
  printf("load %g %g\n",
         actsolv->rhs[stiff_array].vec.a.dv[0],
         actsolv->rhs[stiff_array].vec.a.dv[0]);
#endif
  solserv_sparsematvec(actintra,
                       &(actsolv->rhs[stiff_array]),
                       &(actsolv->sysarray[mass_array]),
                       &(actsolv->sysarray_typ[mass_array]),
                       &(tem[0]));
#if 0
  printf("load %g %g\n",
         actsolv->rhs[stiff_array].vec.a.dv[0],
         actsolv->rhs[stiff_array].vec.a.dv[1]);
#endif
  solserv_scalarprod_vec(&(actsolv->rhs[stiff_array]), 1.0/dt);
#if 0
  printf("load %g %g\n",
         actsolv->rhs[stiff_array].vec.a.dv[0],
         actsolv->rhs[stiff_array].vec.a.dv[1]);
#endif
  solserv_add_vec(&(fint[0]), &(actsolv->rhs[stiff_array]), -(1.0-gamma));
#if 0
  printf("load %g %g\n",
         actsolv->rhs[stiff_array].vec.a.dv[0],
         actsolv->rhs[stiff_array].vec.a.dv[1]);
#endif
  solserv_add_vec(&(fextn[0]), &(actsolv->rhs[stiff_array]), gamma);
#if 0
  printf("load %g %g\n",
         actsolv->rhs[stiff_array].vec.a.dv[0],
         actsolv->rhs[stiff_array].vec.a.dv[1]);
#endif
  solserv_add_vec(&(fext[0]), &(actsolv->rhs[stiff_array]), 1.0-gamma);
#if 0
  printf("load %g %g\n",
         actsolv->rhs[stiff_array].vec.a.dv[0],
         actsolv->rhs[stiff_array].vec.a.dv[1]);
#endif
  /* add Dirichlet forces */
  assemble_vec(actintra, &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]),
               &(actsolv->rhs[stiff_array]), dirich, -1.0);
#if 0
  printf("load %g %g\n",
         actsolv->rhs[stiff_array].vec.a.dv[0],
         actsolv->rhs[stiff_array].vec.a.dv[1]);
#endif
  
  /*--------------------------------------------------------------------*/
  /* call solver */
  init = 0;
  solserv_zero_vec(&(actsolv->sol[stiff_array]));
  solver_control(actfield, disnum, actsolv, actintra,
                 &(actsolv->sysarray_typ[stiff_array]),
                 &(actsolv->sysarray[stiff_array]),
                 &(actsolv->sol[stiff_array]),
                 &(actsolv->rhs[stiff_array]),
                 init);

  /*--------------------------------------------------------------------*/
  /* store solution */
  solserv_copy_vec(&(actsolv->sol[stiff_array]), &(tem[0]));

  /*--------------------------------------------------------------------*/
  /* put the scaled prescribed temperatures to the nodes
   * in field sol (1st 0) at place 0 (2nd 0) together with 
   * free temperatures */
  /* HINT: time curve is called indirectly */
  solserv_putdirich_to_dof(actfield, disnum, 
                           node_array_sol, isol->temn, timcur);

  /*--------------------------------------------------------------------*/
  /* allreduce the result and put it to the node sol arrays */
  solserv_result_total(actfield,
                       disnum,
                       actintra,
                       &(actsolv->sol[stiff_array]),
                       isol->temn,
                       &(actsolv->sysarray[stiff_array]),
                       &(actsolv->sysarray_typ[stiff_array]));

  /*--------------------------------------------------------------------*/
  /* update temperatures at nodes */
  solserv_sol_copy(actfield, disnum,
                   node_array_sol, node_array_sol,
                   isol->temn, isol->tem);

  /*--------------------------------------------------------------------*/
  /* update of external force */
  solserv_copy_vec(&(fextn[0]), &(fext[0]));

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void tsi_th_ost_equi() */

/*======================================================================*/
/*!
\brief Output
\param   actpart      PARTITION*   (io)  partition
\param   actintra     INTRA*       (i)   intra-communicator
\param   actfield     FIELD        (i)   thermal field
\param   disnum       INT          (i)   discretisation index
\param   actsolv      SOLVAR*      (i)  solution variables
\param   isol         ARRAY_POSITION_SOL* (i) indices in NODE sol arrays
\param   actsysarray  INT          (i)   tangent array index
\param   actdyn       THERM_DYNAMIC* (i) thermal control
\param   container    CONTAINER*   (i)   container
\author bborn
\date 05/07
*/
void tsi_th_ost_out(PARTITION* actpart,
                     INTRA* actintra,
                     FIELD* actfield,
                     INT disnum,
                     ARRAY_POSITION_SOL* isol,
                     SOLVAR* actsolv,
                     INT actsysarray,
                     THERM_DYNAMIC* actdyn,
                     CONTAINER* container)
{
  INT numtf = genprob.numtf;  /* index of thermal field */
  CALC_ACTION *action = &(calc_action[numtf]);  /* action */
  INT iplace;
  INT i2ndsysmat;
  DOUBLE timcur = actdyn->time;  /* current time */
  INT timstp = actdyn->step;  /* current time step */
  INT mod_res = actdyn->step % actdyn->out_res_ev;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_th_ost_out");
#endif

  if (mod_res == 0)
  {

  /*--------------------------------------------------------------------*/
  /* perform heat flux calculation */
  if (ioflags.therm_heatflux == 1)
  {
    *action = calc_therm_heatflux;
    container->dvec = NULL;
    container->dirich = NULL;
    container->global_numeq = 0;
    container->kstep = 0;  /* WORK */
    i2ndsysmat = -1;  /* we do not have a second system matrix */
    calelm(actfield, actsolv, actpart, actintra, actsysarray,
           i2ndsysmat, container, action);
  }

  /*--------------------------------------------------------------------*/
  /* printout results to out */
  if ( (ioflags.output_out == 1) && (ioflags.therm_temper == 1) )
  {
    iplace = 0;  /* must be zero! used at heatflux arrays,
                  * e.g. dynamic elementwise array hflux_gp is
                  * allocated like hflux_gp[1][3*NUMHFLX][MAXGAUSS]
                  * and iplace==0 is used to access the first (left-most)
                  * index */
    /* ATTENTION : iplace is also used for 
     *             actnode->sol.a.da[iplace][nodeindex]
     *             while printing the temperature */
    out_sol(actfield, actpart, disnum, actintra, timstp, isol->temn);
  }

  /*--------------------------------------------------------------------*/
  /* printout results to binary file */
#ifdef BINIO
  if (ioflags.output_bin == 1)
  {
    dserror("BINIO is not available!");
  }
#endif

  /*--------------------------------------------------------------------*/
  /* printout results to Gid */
  if ( (ioflags.output_gid == 1)
       && (ioflags.therm_temper == 1) 
       && (par.myrank == 0) )
  {
    iplace = 0;
    out_gid_sol("temperature", actfield, disnum, actintra, 
                timstp, iplace, timcur);
    /* out_gid_domains(actfield, disnum); */
  }

  /*--------------------------------------------------------------------*/
  /* printout heat flux to Gid */
  if ( (ioflags.output_gid == 1)
       && (ioflags.therm_heatflux == 1)
       && (par.myrank == 0) )
  {
    iplace = 0;
    out_gid_sol("heatflux", actfield, disnum, actintra, 
                timstp, iplace, timcur);
  }

  }

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void tsi_th_ost_out() */

/*======================================================================*/
/*!
\brief Finalise
\param   actsolv      SOLVAR*        (io)  solution variables
\param   out_context  BIN_OUT_FIELD* (io)  BINIO context
\param   dirich_a     ARRAY*         (io)  Dirichlet RHS contribution
\author bborn
\date 05/07
*/
void tsi_th_ost_final(SOLVAR* actsolv,
#ifdef BINIO
                      BIN_OUT_FIELD* out_context,
#endif
                      DIST_VECTOR** tem,
                      DIST_VECTOR** fext,
                      DIST_VECTOR** fextn,
                      DIST_VECTOR** fint,
                      ARRAY* intforce_a,
                      ARRAY* dirich_a)
{
  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_th_ost_final");
#endif

  /*--------------------------------------------------------------------*/
  /* deallocate stuff */
  solserv_del_vec(tem, 1);
  solserv_del_vec(fint, 1);
  solserv_del_vec(fext, 1);
  solserv_del_vec(fextn, 1);
  solserv_del_vec(&(actsolv->rhs), actsolv->nrhs);
  solserv_del_vec(&(actsolv->sol), actsolv->nsol);
  amdel(dirich_a);
  amdel(intforce_a);

  /*--------------------------------------------------------------------*/
#ifdef BINIO
  if (ioflags.output_bin == 1)
  {
    destroy_bin_out_field(out_context);
  }
#endif

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}  /* end void tsi_th_ost_final() */


/*======================================================================*/
/*!
\brief Static iso-thermal solution of thermal field
       Modularised
\param   disnum_s     INT    (i)   index of structural discretisation
\param   disnum_t     INT    (i)   index of thermal discretisation
\author bborn
\date 03/06
*/
void tsi_th_ost_sub(INT disnum_s,
                     INT disnum_t)
{
  INT numtf = genprob.numtf;  /* number (index) of thermal field */

  PARTITION* actpart = &(partition[numtf]);  /* thermal partition */
  INTRA* actintra;  /* active intra-communicator */

  /* field */
  FIELD* actfield = &(field[numtf]);  /* pointer to the thermal FIELD */

  /* discretisation */
  /* disnum_t == index (number) of discretisation */
  /* only a single discretisation
   * this variable name is a little confusing
   * it sets an _actual_ number (better index?)
   * of one of the actfield->ndis discretisations.
   * Simply said: disnum != n(um)dis */
  /* named positions of NODE sol etc. arrays */
  ARRAY_POSITION* ipos = &(actfield->dis[disnum_t].ipos);
  /* named positions (indices) of NODE sol array */
  ARRAY_POSITION_SOL* isol = &(ipos->isol);

  /* solution variable */
  SOLVAR* actsolv = &(solv[numtf]);  /* pointer to field SOLVAR */
  INT numeq;  /* number of equations on this proc */
  INT numeq_total;  /* total number of equations on all procs */
  INT stiff_array;  /* index of stiffness sparse system matrix in 
                     * actsolv->sysarray[] */
  INT mass_array;  /* index of mass matrix in sparse system matrix 
                    * actsolv->sysarray[] */
  
  /* dynamic control (allright, static sontrol) */
  THERM_DYNAMIC* actdyn = alldyn[numtf].tdyn;  /* thermal dynamic control */

  /* output */
#ifdef BINIO
  BIN_OUT_FIELD out_context;
#endif

  /* container */
  CONTAINER container;  /* contains variables defined in container.h */
  
  /* variables */
  DOUBLE timcur;  /* current time */
  DOUBLE timstp;  /* current time step */

  ARRAY dirich_a;   /* redund. full length vect. for Dirichlet-part of RHS */
  ARRAY intforce_a;  /* redund. full length vect. for internal heat */

  DIST_VECTOR* tem;  /* global temperature (1) vector */
  DIST_VECTOR* fint;  /* internal heat vector */
  DIST_VECTOR* fext;  /* external heat vector */
  DIST_VECTOR* fextn;  /* new external heat vector */



  /*====================================================================*/
#ifdef DEBUG
  dstrc_enter("tsi_th_ost_sub");
#endif

  /*--------------------------------------------------------------------*/
  /* a word to the user */
  printf("==============================================================="
         "=========\n");
  printf("TSI quasi static thermal solution routine reached.\n");
  printf("SUB TYPE\n");
  printf("---------------------------------------------------------------"
         "---------\n");

  /*--------------------------------------------------------------------*/
  /* set up container */
  container.fieldtyp = actfield->fieldtyp;  /* thermal field */
  container.isdyn = 0;  /* static calculation */  /* ? */
  container.kintyp = 0;  /* kintyp  = 0: geo_lin */  /* ? */
  container.disnum = disnum_t;  /* current thermal discretisation index */
  container.disnum_t = disnum_t;  /* thermal discretisation index */
  container.disnum_s = disnum_s;  /* structural discretisation index */

  /*--------------------------------------------------------------------*/
  /* intra communicator */
#ifdef PARALLEL
  actintra = &(par.intra[numtf]);
#else
  actintra = (INTRA*) CCACALLOC(1, sizeof(INTRA));
  if (!actintra)
  {
    dserror("Allocation of INTRA failed");
  }
  actintra->intra_fieldtyp = thermal;
  actintra->intra_rank = 0;
  actintra->intra_nprocs = 1;
#endif
  /* there are only procs allowed in here, that belong to the thermal
   * intracommunicator (in case of linear statics, this should be all) */
  if (actintra->intra_fieldtyp != thermal)
  {
    goto end;
  }

  /*--------------------------------------------------------------------*/
  /* one-step-theta (the theta is gamma) time integration scheme */
  actdyn->gamma = 0.5;
  actdyn->initmpr = 0.0;

  /*====================================================================*/
  /* initialise */
  tsi_th_ost_init(actpart,
                  actintra,
                  actfield,
                  disnum_t,
                  ipos,
                  isol,
                  actsolv,
                  &(numeq),
                  &(numeq_total),
                  &(stiff_array),
                  &(mass_array),
                  actdyn,
                  &(container),
                  &(tem),
                  &(fint),
                  &(fext),
                  &(fextn),
                  &(intforce_a),
                  &(dirich_a));

  /*------------------------------------------------------------------*/
  /* set time controls */
  timcur = actdyn->time;
  timstp = actdyn->step;

  /*====================================================================*/
  /* solution */
  tsi_th_ost_equi(actpart,
                  actintra,
                  actfield,
                  disnum_s,
                  isol,
                  actsolv,
                  numeq,
                  numeq_total,
                  stiff_array,
                  mass_array,
                  actdyn,
                  &(container),
                  tem,
                  fint,
                  fext,
                  fextn,
                  &(intforce_a),
                  &(dirich_a));

  /*====================================================================*/
  /* output */
  tsi_th_ost_out(actpart,
                 actintra,
                 actfield,
                 disnum_s,
                 isol,
                 actsolv,
                 stiff_array,
                 actdyn,
                 &(container));

  /*--------------------------------------------------------------------*/
  /* a sag-beim-abschied-leise-servus to the user */
  printf("---------------------------------------------------------------"
         "---------\n");
  printf("TSI quasi static thermal solution routine finished.\n");
  printf("==============================================================="
         "=========\n\n");

  /*--------------------------------------------------------------------*/
  /* the end */
  end:

  /*====================================================================*/
  /* finalise */
  tsi_th_ost_final(actsolv,
#ifdef BINIO
                   &(out_context),
#endif
                   &(tem),
                   &(fint),
                   &(fext),
                   &(fextn),
                   &(intforce_a),
                   &(dirich_a));

#ifndef PARALLEL
  CCAFREE(actintra);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif

  return;

}  /* end of tsi_th_ost_sub() */


/*----------------------------------------------------------------------*/
#endif  /* end of #ifdef D_TSI */
