/*======================================================================*/
/*!
\file
\brief TSI - time integration of structure field

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
/* header files */
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../io/io.h"


/*----------------------------------------------------------------------*/
/*!
\brief File pointers

This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
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
/*----------------------------------------------------------------------*/
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
ALLDYNA               *alldyn;

\author bborn
\date 03/06
*/
extern ALLDYNA *alldyn;


/*----------------------------------------------------------------------*/
/*!
\brief Load curve, load factor function

number of load curves numcurve
vector of structures of curves
defined in input_curves.c
INT                   numcurve;
struct _CURVE        *curve;

\author bborn
\date 03/06
*/
extern INT numcurve;
extern CURVE *curve;


/*----------------------------------------------------------------------*/
/*!
\brief CALC_ACTIONs

enum _CALC_ACTION
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
DOUBLE deltat;


/*======================================================================*/
/*!
\brief 

\author bborn
\date 03/06
*/
void tsi_st_genalp(INT disnum_s,
                   INT disnum_t)
{

  INT i;  /* simply a counter */
  INT numeq;  /* number of equations on this proc */
  INT numeq_total;  /* total number of equations */
  INT init;  /* flag for solver_control call */
  INT mod_disp, mod_stress;
  INT mod_res_write;
  INT disnum = 0;
  INT timeadapt;  /* flag to switch time adaption on/off */

  DOUBLE deltaepot=0.0;

  INT convergence;  /* convergence flag */
  INT itnum;  /* iterator */
  DOUBLE dmax;  /* infinity norm of residual displacements */

  DOUBLE t0, t1;
  
  INT stiff_array;  /* indice of the active system sparse matrix */
  INT mass_array;  /* indice of the active system sparse matrix */
  INT damp_array;  /* indice of the active system sparse matrix */
  INT actcurve;  /* indice of active time curve */
  
  SOLVAR *actsolv;  /* pointer to active solution structure */
  PARTITION *actpart;  /* pointer to active partition */
  INT numsf;  /* number (index) of structure field */
  FIELD *actfield;  /* pointer to active field */
  INTRA *actintra;  /* pointer to active intra-communicator */
  CALC_ACTION *action;  /* pointer to the structure cal_action enum */
  STRUCT_DYNAMIC *actdyn;  /* pointer to structural dynamic input data */
  INT actsysarray;  /* index of actual system array */
  
  DIST_VECTOR *vel;  /* total velocities */
  DIST_VECTOR *acc;  /* total accelerations */
  DIST_VECTOR *fie;  /* internal forces and working array */
  DIST_VECTOR *dispi;  /* distributed vector to hold increm. displacments */
  DIST_VECTOR *work;  /* working vectors */
  
  ARRAY intforce_a;  /* redundant vect. of full length for internal forces */
  DOUBLE *intforce;
  ARRAY dirich_a;  /* redund. vect. of full length for dirichlet-part of rhs */
  DOUBLE *dirich;
  DOUBLE dirichfacs[10];  /* factors needed for dirichlet-part of rhs */
  
  STRUCT_DYN_CALC dynvar;  /* variables to perform dynamic structural simulation */

  ARRAY_POSITION *ipos;   /* named positions of NODE sol etc. arrays */
  ARRAY_POSITION_SOL *isol;  /* named positions (indices) of NODE sol array */
  ARRAY_POSITION_SOLINC *isolinc;  /* named positions (indices) of 
                                    * NODE sol_increment array */
  ARRAY_POSITION_SOLRES *isolres;  /* named positions (indices) of 
                                    * NODE sol_residual array */

#ifdef BINIO
  BIN_OUT_FIELD out_context;
#endif

  CONTAINER container;  /* transfers variables given
                         * in (this) solution technique
                         * to element level */

  /*====================================================================*/
  /* begin body */

#ifdef DEBUG
  dstrc_enter("tsi_dyn_struct");
#endif

  /*--------------------------------------------------------------------*/
  /* a word to the user */
  printf("==============================================================="
         "=========\n");
  printf("TSI structural time integration with generalised-alpha\n");
  printf("---------------------------------------------------------------"
         "---------\n");

  /*--------------------------------------------------------------------*/
  /* set up pointers and container */
  numsf = genprob.numsf;  /* index of structure field */
  actfield = &(field[numsf]);  /* pointer to structure field */
  actdyn = alldyn[numsf].sdyn;  /* dynamic control */
  timeadapt = actdyn->timeadapt;  /* flag for time adaptivity */
  acttime = actdyn->time;  /* initial time */
  container.fieldtyp = actfield->fieldtyp;  /* field type : structure */
  container.isdyn = 1;  /* dynamic calculation */
  container.kintyp = 2;  /* kintyp  = 2: total Lagrangean */  /* ??? */
  disnum = disnum_s;  /* structure discretisation index:
                       * disnum_s == 0,
                       * as only a single discretisation
                       * this variable name is a little confusing
                       * it sets an _actual_ number (better index?)
                       * of one of the actfield->ndis discretisations.
                       * Simply said: disnum != n(um)dis */
  container.disnum = disnum;
  container.disnum_s = disnum_s;  /* structure discretisation index ( ==0 ) */
  container.disnum_t = disnum_t;  /* thermo-discretisation index ( ==0 ) */
  actsysarray = numsf;  /* ? */
  actsolv = &(solv[numsf]);
  if (actsolv->nsysarray == 1)
  {
    actsysarray = 0;
  }
  else
  {
    dserror("More than 1 system arrays (actsolv->nsysarray)!");
  }
  actpart = &(partition[numsf]);  /* of structure field */
  action = &(calc_action[numsf]);  /* ? */

  /*--------------------------------------------------------------------*/
  /* intra communicator */
#ifdef PARALLEL
  actintra = &(par.intra[0]);  /* ???? => numsf ??? */
#else
  actintra = (INTRA*) CCACALLOC(1, sizeof(INTRA));
  if (!actintra)
  {
    dserror("Allocation of INTRA failed");
  }
  actintra->intra_fieldtyp = structure;
  actintra->intra_rank = 0;
  actintra->intra_nprocs = 1;
#endif
  /* there are only procs allowed in here, that belong to the thermal
   * intracommunicator (in case of linear statics, this should be all) */
  if (actintra->intra_fieldtyp != structure)
  {
    goto end;
  }

  /*--------------------------------------------------------------------*/
  /* init the variables in dynvar to zero */
  /* Set all variables to zero. No matter what changes in future. */
  memset(&dynvar, 0, sizeof(STRUCT_DYN_CALC));

  /*--------------------------------------------------------------------*/
  /* damping */
  if (actdyn->damp == 1)
  {
    stiff_array = 0;
    mass_array = 1;
    damp_array = 2;
    actsolv->nsysarray = 3;
  }
  else
  {
    stiff_array = 0;
    mass_array = 1;
    damp_array = -1;
    actsolv->nsysarray = 2;
  }

  /*--------------------------------------------------------------------*/
  /* allocate sparse mass (and damping) matrix */
  /* reallocate the vector of sparse matrices and the vector of 
   * their types : formerly lenght 1, now lenght 2 or 3 dependent on 
   * presence of damp_array */
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
                              &(actsolv->sysarray_typ[stiff_array]),
                              &(actsolv->sysarray[stiff_array]),
                              &(actsolv->sysarray_typ[mass_array]),
                              &(actsolv->sysarray[mass_array]));
  
  if (damp_array > 0)
  {
    solserv_alloc_cp_sparsemask(actintra,
                                &(actsolv->sysarray_typ[stiff_array]),
                                &(actsolv->sysarray[stiff_array]),
                                &(actsolv->sysarray_typ[damp_array]),
                                &(actsolv->sysarray[damp_array]));
  }

  /*--------------------------------------------------------------------*/
  /* init the dist sparse matrices to zero */
  for (i=0; i<actsolv->nsysarray; i++)
  {
    solserv_zero_mat(actintra,
                     &(actsolv->sysarray[i]),
                     &(actsolv->sysarray_typ[i]));
  }

  /*--------------------------------------------------------------------*/
  /* get global and local number of equations */
  solserv_getmatdims(&(actsolv->sysarray[stiff_array]),
                     actsolv->sysarray_typ[stiff_array],
                     &numeq, &numeq_total);

  /*--------------------------------------------------------------------*/
  /* allocate 4 distributed vectors for RHS */
  /* these hold original load vector, 
   *            load vector at time t 
   *            load vctor at time at t-dt 
   *            and interpolated load vector */
  actsolv->nrhs = 4;
  solserv_create_vec(&(actsolv->rhs), 
                     actsolv->nrhs, numeq_total, numeq, "DV");
  for (i=0; i<actsolv->nrhs; i++)
  {
    solserv_zero_vec(&(actsolv->rhs[i]));
  }

  /*--------------------------------------------------------------------*/
  /* allocate 2 dist. solution/displacement vectors */
  /* displacement vector at t
   * displacement vector at t-dt */
  actsolv->nsol = 2;
  solserv_create_vec(&(actsolv->sol),
                     actsolv->nsol, numeq_total, numeq, "DV");
  for (i=0; i<actsolv->nsol; i++) 
  {
    solserv_zero_vec(&(actsolv->sol[i]));
  }

  /*--------------------------------------------------------------------*/
  /* allocate 1 dist vector for iterative displacements increments */
  solserv_create_vec(&dispi, 1, numeq_total, numeq, "DV");
  solserv_zero_vec(&(dispi[0]));

  /*--------------------------------------------------------------------*/
  /* allocate 1 dist vector for velocities */
  solserv_create_vec(&vel, 1, numeq_total, numeq, "DV");
  solserv_zero_vec(&(vel[0]));

  /*--------------------------------------------------------------------*/
  /* allocate 1 dist vector for accelerations */
  solserv_create_vec(&acc, 1, numeq_total, numeq, "DV");
  solserv_zero_vec(&(acc[0]));

  /*--------------------------------------------------------------------*/
  /* create 1 redundant full-length vector for internal forces */
  intforce = amdef("intforce", &intforce_a, numeq_total, 1, "DV");

  /*--------------------------------------------------------------------*/
  /* create 1 vector of full length for Dirichlet part of RHS */
  dirich = amdef("dirich", &dirich_a, numeq_total, 1, "DV");

  /*--------------------------------------------------------------------*/
  /* allocate 3 dist. vectors for internal forces */
  /* internal force at t
   * internal force at t-dt
   * mid-internal forece at t-dt/2 */
  solserv_create_vec(&fie, 3, numeq_total, numeq, "DV");
  for (i=0; i<3; i++)
  {
    solserv_zero_vec(&(fie[i]));
  }

  /*--------------------------------------------------------------------*/
  /* allocate 3 dist. working vectors */
  solserv_create_vec(&work, 3, numeq_total, numeq, "DV");
  for (i=0; i<3; i++)
  {
    solserv_zero_vec(&(work[i]));
  }

  /*--------------------------------------------------------------------*/
  /* initialize solver on all matrices */
  /* NOTE: solver init phase has to be called with each matrix one wants to
   *       solve with. Solver init phase has to be called with all matrices
   *       one wants to do matrix-vector products and matrix scalar products.
   *       This is not needed by all solver libraries, but the solver-init 
   *       phase is cheap in computation (can be costly in memory)
   *       There will be no solver call on mass or damping array. */

  /*--------------------------------------------------------------------*/
  /* initialize solver */
  init = 1;
  solver_control(actfield, disnum, actsolv, actintra, 
		 &(actsolv->sysarray_typ[stiff_array]),
                 &(actsolv->sysarray[stiff_array]),
                 &(dispi[0]), &(actsolv->rhs[0]), init);

  solver_control(actfield, disnum, actsolv, actintra, 
		  &(actsolv->sysarray_typ[mass_array]),
                 &(actsolv->sysarray[mass_array]),
                 &work[0], &work[1], init);

  if (damp_array > 0)
  {
    solver_control(actfield, disnum, actsolv, actintra, 
		   &(actsolv->sysarray_typ[damp_array]),
                   &(actsolv->sysarray[damp_array]),
                   &work[0], &work[1], init);
  }

  /*--------------------------------------------------------------------*/
  /* init the assembly for stiffness and for mass matrix */
  /* (damping is not assembled) */
  init_assembly(actpart, actsolv, actintra, actfield, stiff_array, disnum);
  init_assembly(actpart, actsolv, actintra, actfield, mass_array, disnum);

  /*--------------------------------------------------------------------*/
  /* init the element calculating routines */
  *action = calc_struct_init;
  calinit(actfield, actpart, action, &container);
  
  /*--------------------------------------------------------------------*/
  /* write output of mesh to gid */
  /* if ( (par.myrank == 0) && (ioflags.output_gid == 1) ) */
  /*   { */
  /*     out_gid_msh(); */
  /*   } */

  /*--------------------------------------------------------------------*/
  /* call elements to calculate stiffness and mass */
  *action = calc_struct_nlnstiffmass;
  container.dvec          = NULL;
  container.dirich        = NULL;
  container.global_numeq  = 0;
  container.dirichfacs    = NULL;
  container.kstep         = 0;
  deltat                  = actdyn->dt;
  calelm(actfield, actsolv, actpart, actintra,
         stiff_array, mass_array, &container, action);

  /*--------------------------------------------------------------------*/
  /* calculate damping matrix */
  if (damp_array > 0)
  {
    /* stiffness proportional contribution */
    solserv_add_mat(actintra,
                    &(actsolv->sysarray_typ[damp_array]),
                    &(actsolv->sysarray[damp_array]),
                    &(actsolv->sysarray_typ[stiff_array]),
                    &(actsolv->sysarray[stiff_array]),
                    actdyn->k_damp);
    /* mass proportional contribution */
    solserv_add_mat(actintra,
                    &(actsolv->sysarray_typ[damp_array]),
                    &(actsolv->sysarray[damp_array]),
                    &(actsolv->sysarray_typ[mass_array]),
                    &(actsolv->sysarray[mass_array]),
                    actdyn->m_damp);
    solserv_close_mat(actintra,
		      &(actsolv->sysarray_typ[damp_array]),
		      &(actsolv->sysarray[damp_array]));
  }

  /*--------------------------------------------------------------------*/
  /* set initial step and time */
  actdyn->step = -1;
  actdyn->time = 0.0;

  /*--------------------------------------------------------------------*/
  /* init all applied time curves -*/
  for (actcurve=0; actcurve<numcurve; actcurve++)
  {
    /* dyn_init_curve(actcurve, actdyn->nstep, actdyn->dt, actdyn->maxtime); */
  }

  /*--------------------------------------------------------------------*/
  /* sol indices */
  ipos = &(actfield->dis[disnum].ipos);  /* position array */
  isol = &(ipos->isol);
  isolinc = &(ipos->isolinc);
  isolres = &(ipos->isolres);

  /*--------------------------------------------------------------------*/
  /* put a zero to the place ipos->num=12 in node->sol to init the 
   * velocities and accels of prescribed displacements */
  /* HINT: This actually redefines/reallocates/enlarges the sol
   *       array of each structure node to dimension 12x2 (or 12x3)
   *       from originally 1x2 (or 1x3) */
  solserv_sol_zero(actfield, disnum, node_array_sol, ipos->num-1);

  /*--------------------------------------------------------------------*/
  /* put a zero to the place ipos->numincr=2 in sol_increment of NODEs */
  /* later this will hold internal forces at t and t-dt */
  /* HINT: This actually redefines/reallocates/enlarges the sol_increment
   *       array at each structure node to dimenion 2x2 (or 2x3)
   *       from originally 1x2 (or 1x3) */
  solserv_sol_zero(actfield, disnum, node_array_sol_increment, 
                   ipos->numincr-1);
  /* initialise internal forces f_{int;n} to zero */
  solserv_sol_zero(actfield, disnum, node_array_sol_increment, 
                   isolinc->fint);
  /* initialise internal forces f_{int;n+1} to zero */
  solserv_sol_zero(actfield, disnum, node_array_sol_increment, 
                   isolinc->fintn);


  /*--------------------------------------------------------------------*/
  /* WARNING : BINIO is not available --- work needs to be done */
#ifdef BINIO
  /* initialize binary output
   * It's important to do this only after all the node arrays are set
   * up because their sizes are used to allocate internal memory. */
  init_bin_out_field(&out_context,
                     &(actsolv->sysarray_typ[stiff_array]), 
                     &(actsolv->sysarray[stiff_array]),
                     actfield, actpart, actintra, 0);
#endif

  /*--------------------------------------------------------------------*/
  /*  output to GID postprozessor */
  if ( (par.myrank == 0) && (ioflags.output_gid == 1) )
  {
    /* out_gid_domains(actfield, disnum); */
  }

  /*--------------------------------------------------------------------*/
  /* printout head */
  if (par.myrank == 0) 
  {
    dyn_nlnstruct_outhead(&dynvar, actdyn);
  }

  /*====================================================================*/
  /* START LOOP OVER ALL TIME STEPS */
  /*====================================================================*/
  /*
   * rhs[3]    original load vector
   * rhs[2]             load vector at time t-dt
   * rhs[1]             load vector at time t
   * rhs[0]    interpolated load vector and working array
   *
   * fie[2]    internal forces at step t
   * fie[1]    internal forces at step t-dt
   * fie[0]    interpolated internal forces and working array
   *
   * dispi[0]  displacement increment from t-dt to t
   *
   * sol[0]    total displacements at time t-dt
   * sol[1]    total displacements at time t
   *
   * vel[0]    velocities    at t-dt
   * acc[0]    accelerations at t-dt
   *
   * work[2]   working vector for sums and matrix-vector products
   * work[1]   working vector for sums and matrix-vector products
   * work[0]   working vector for sums and matrix-vector products
   * work[0]   is used to hold residual displacements in corrector
   *           iteration
   *
   * in the nodes, displacements are kept in node[].sol[0][0..numdf-1]
   *               velocities    are kept in node[].sol[1][0..numdf-1]
   *               accelerations are kept in node[].sol[2][0..numdf-1]
   *
   * Values of the different vectors from above in one loop:
   *    /   ...   no change in this step
   *    =   ...   evaluation in this step
   *    +=  ...   evaluation in this step
   *
   */    
  while ( (actdyn->step < actdyn->nstep-1) 
          && (actdyn->time <= actdyn->maxtime) )
  {
    
    /*------------------------------------------------------------------*/
    /* wall clock time ?????? */
    t0 = ds_cputime();

    /*------------------------------------------------------------------*/
    /* increment step */
    actdyn->step++;
    /* repeatcount = 0; */

    /*------------------------------------------------------------------*/
    /* set new time */
    actdyn->time = actdyn->time + actdyn->dt;
    /* put time to global variable for time-dependent load distributions */
    acttime = actdyn->time;

    /*------------------------------------------------------------------*/
    /* set some constants */
    dyn_setconstants(&dynvar, actdyn, actdyn->dt);

    /*------------------------------------------------------------------*/
    /* set incremental displacements dispi[0] to zero */
    solserv_zero_vec(&dispi[0]);

    /*------------------------------------------------------------------*/
    /* set residual displacements in nodes to zero */
    solserv_result_resid(actfield, disnum, actintra, &dispi[0], 
                         isolres->disres,
                         &(actsolv->sysarray[stiff_array]),
                         &(actsolv->sysarray_typ[stiff_array]));


    /*==================================================================*/
    /* PREDICTOR */
    /*==================================================================*/
    
    /*------------------------------------------------------------------*/
    /*  this vector holds loads due to external forces */
    solserv_zero_vec(&(actsolv->rhs[1]));
    container.kstep = 0;
    container.inherit = 1;  /* ??? */
    container.point_neum = 1;  /* ??? */
    *action = calc_struct_eleload;
    calrhs(actfield, actsolv, actpart, actintra, stiff_array,
	   &(actsolv->rhs[1]), action, &container);
    
    /*------------------------------------------------------------------*/
    /* multiply rhs[1] by load factor based on factor rldfac of curve 0 */
    /* WARNING: This control routine at the moment always uses curve 0 
     *          for the RHS */
    actcurve = 0;
    /* Get factor at new time t */
    dyn_facfromcurve(actcurve, actdyn->time, &(dynvar.rldfac));
    dynvar.rldfac = 1.0;
    solserv_scalarprod_vec(&(actsolv->rhs[1]), dynvar.rldfac);

    /*------------------------------------------------------------------*/
    /* put the scaled prescribed displacements to the nodes in field sol
     * at place 4 separate of the free DOFs
     * These are used to calculate the RHS due to the Dirichlet conditions */
    solserv_putdirich_to_dof(actfield, disnum, 0, 4, actdyn->time);

    /*------------------------------------------------------------------*/
    /* put presdisplacements(t) - presdisplacements(t-dt) in place 5 */
    solserv_adddirich(actfield, disnum, 0, 3, 4, 5, -1.0, 1.0);

    /*------------------------------------------------------------------*/
    /* set factors needed for prescribed displacement terms on eff RHS */
    /* dirichfacs[0] = -(1.0-alpham)*(1.0/beta)/(DSQR(dt))
     * dirichfacs[1] =  (1.0-alpham)*(1.0/beta)/dt
     * dirichfacs[2] =  (1.0-alpham)/(2*beta) - 1
     * dirichfacs[3] = -(1.0-alphaf)*(gamma/beta)/dt
     * dirichfacs[4] =  (1.0-alphaf)*gamma/beta - 1
     * dirichfacs[5] =  (gamma/(2*beta)-1)*(1.0-alphaf)
     * dirichfacs[6] = -(1.0-alphaf) or 0
     * dirichfacs[7] =  raleigh damping factor for mass
     * dirichfacs[8] =  raleigh damping factor for stiffness
     * dirichfacs[9] =  dt
     * see PhD theses Mok page 165: Generalized-alpha time integration 
     *                              with prescribed displ. */
    dirichfacs[0] = -dynvar.constants[0];
    dirichfacs[1] =  dynvar.constants[1];
    dirichfacs[2] =  dynvar.constants[2];
    dirichfacs[3] = -dynvar.constants[3];
    dirichfacs[4] =  dynvar.constants[4];
    dirichfacs[5] =  dynvar.constants[5];
    dirichfacs[6] = -dynvar.constants[6];
    if (damp_array > 0) 
    {
      dirichfacs[7] = actdyn->m_damp;
      dirichfacs[8] = actdyn->k_damp;
    }
    else 
    {
      dirichfacs[7] = 0.0;
      dirichfacs[8] = 0.0;
    }
    dirichfacs[9] =  actdyn->dt;

    /*------------------------------------------------------------------*/
    /* calculate tangential stiffness/mass 
     * and internal forces at time t-dt */
    solserv_zero_mat(actintra, 
                     &(actsolv->sysarray[stiff_array]),
                     &(actsolv->sysarray_typ[stiff_array]));
    solserv_zero_mat(actintra, &(actsolv->sysarray[mass_array]),
                     &(actsolv->sysarray_typ[mass_array]));
    amzero(&dirich_a);
    amzero(&intforce_a);

    /*------------------------------------------------------------------*/
    /*  call elements */
    *action = calc_struct_nlnstiffmass;
    container.dvec          = intforce;
    container.dirich        = dirich;
    container.global_numeq  = numeq_total;
    container.dirichfacs    = dirichfacs;
    container.kstep         = 0;
    calelm(actfield, actsolv, actpart, actintra, 
           stiff_array, mass_array, &container, action);

    /*------------------------------------------------------------------*/
    /* store positive internal forces on fie[1] */
    solserv_zero_vec(&fie[1]);
    assemble_vec(actintra, 
                 &(actsolv->sysarray_typ[stiff_array]),
                 &(actsolv->sysarray[stiff_array]),
                 &(fie[1]), intforce, 1.0);

    /*------------------------------------------------------------------*/
    /* determine external mid-force vector by interpolating */
    /* forces rhs[0] = (1-alphaf)rhs[1] + alphaf*rhs[2] */
    solserv_copy_vec(&(actsolv->rhs[2]), &(actsolv->rhs[0]));
    solserv_scalarprod_vec(&(actsolv->rhs[0]), actdyn->alpha_f);
    solserv_add_vec(&(actsolv->rhs[1]), &(actsolv->rhs[0]), 
                    (1.0-actdyn->alpha_f));

    /*------------------------------------------------------------------*/
    /* subtract internal forces from interpolated external forces */
    solserv_add_vec(&(fie[1]), &(actsolv->rhs[0]), -1.0);

    /*------------------------------------------------------------------*/
    /* add rhs from prescribed displacements to RHS */
    assemble_vec(actintra, &(actsolv->sysarray_typ[stiff_array]),
                 &(actsolv->sysarray[stiff_array]),
                 &(actsolv->rhs[0]), dirich, 1.0);

    /*------------------------------------------------------------------*/
    /* create effective load vector (rhs[0]-fie[2])eff */
    /* Peff = rhs[0] - fie[0]
     *   + M*(-a1*dispi[0]+a2*vel[0]+a3*acc[0])
     *   + D*(-a4*dispi[0]+a5*vel[0]+a6*acc[0]) (if present)
     *
     *   a1 = dynvar.constants[0] 
     *      =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
     *   a2 =                        
     *      = ((1.0-alpham) * (1.0/beta)/(DSQR(dt)))*dt
     *   a3 = dynvar.constants[2] 
     *      =  (1.0-alpham) / (2.0*beta) - 1.0
     *   a4 = dynvar.constants[3] 
     *      =  (1.0-alphaf) * ((gamma/beta)/dt)
     *   a5 = dynvar.constants[4] 
     *      =  ((1.0-alphaf) * ((gamma/beta)/dt))*dt - 1.0
     *   a6 =                        
     *      = (gamma/beta)/2.0 - 1.0) * dt * (1.0-alphaf) */
    pefnln_struct(&dynvar, actdyn, actfield, actsolv, actintra,
                  dispi, vel, acc, work, mass_array, damp_array);

    /*------------------------------------------------------------------*/
    /* create effective stiffness matrix */
    /* keff = constants[6] * K + constants[0] * M + constants[3] * D
     *   constants[6] =  (1.0-alphaf)
     *   constants[0] =  (1.0-alpham) * (1.0/beta)/(DSQR(dt))
     *   constants[3] =  (1.0-alphaf) * ((gamma/beta)/dt) */
    kefnln_struct(&dynvar, actdyn, actfield, actsolv, actintra,
                  work, stiff_array, mass_array, damp_array);

    /*------------------------------------------------------------------*/
    /* call for solution of system dispi[0] = Keff^-1 * rhs[0] */
    init = 0;
    solver_control(actfield, disnum, actsolv, actintra,
                   &(actsolv->sysarray_typ[stiff_array]),
                   &(actsolv->sysarray[stiff_array]),
                   &(dispi[0]), &(actsolv->rhs[0]), init);

    /*==================================================================*/
    /* update */

    /*------------------------------------------------------------------*/
    /* update displacements sol[1] = sol[0] + dispi[0] */
    solserv_copy_vec(&(actsolv->sol[0]), &(actsolv->sol[1]));
    solserv_add_vec(&dispi[0], &(actsolv->sol[1]), 1.0);

    /*------------------------------------------------------------------*/
    /* put the scaled prescribed displacements to the nodes */
    /* in field sol (0) at place 0 together with free displacements
     * these are used to calculate the stiffness matrix */
    solserv_putdirich_to_dof(actfield, disnum, 
                             node_array_sol, isol->disn, actdyn->time);

    /*------------------------------------------------------------------*/
    /* return total displacements to the nodes */
    solserv_result_total(actfield, disnum, actintra,
                         &(actsolv->sol[1]),
                         isol->disn, 
                         &(actsolv->sysarray[stiff_array]),
                         &(actsolv->sysarray_typ[stiff_array]));

    /*------------------------------------------------------------------*/
    /* return incremental displacements to the nodes */
    solserv_result_incre(actfield, disnum, actintra,
                         &dispi[0],
                         isolinc->disinc,
                         &(actsolv->sysarray[stiff_array]),
                         &(actsolv->sysarray_typ[stiff_array]));

    /*------------------------------------------------------------------*/
    /* check convergence of predictor */
    convergence = 0;
    dmax        = 0.0;
    solserv_vecnorm_euclid(actintra, &(dispi[0]), &(dynvar.dinorm));
    solserv_vecnorm_euclid(actintra, &(dispi[0]), &(dynvar.dnorm));
    solserv_vecnorm_Linf(actintra, &(dispi[0]), &dmax);
    if (par.myrank == 0) 
    {
      printf("                                                   "
             "Residual %10.5E\n", dynvar.dinorm);
    }
    fflush(stdout);
    if ( (dynvar.dinorm < actdyn->toldisp)
         || (dynvar.dnorm < EPS14)
         || ( (dynvar.dinorm < EPS14) && (dmax < EPS12) ) )
    {
      convergence = 1;  /* inefficient ... otherwise residuals */
    }

    /*==================================================================*/
    /* PERFORM EQUILIBRIUM ITERATION */
    /*==================================================================*/
    itnum = 0;
    while ( (convergence != 1) && (itnum <= actdyn->maxiter) )
    {

      /*----------------------------------------------------------------*/
      /* check if maximally permitted iterations reached */
      if ( (itnum == actdyn->maxiter) && !timeadapt )
      {
        dserror("No convergence in maxiter steps");
      }

      /*----------------------------------------------------------------*/
      /* set factors needed for prescribed displacement terms on eff RHS */
      dirichfacs[0] = -dynvar.constants[0];
      dirichfacs[1] =  dynvar.constants[1];
      dirichfacs[2] =  dynvar.constants[2];
      dirichfacs[3] = -dynvar.constants[3];
      dirichfacs[4] =  dynvar.constants[4];
      dirichfacs[5] =  dynvar.constants[5];
      dirichfacs[6] =  0.0;
      if (damp_array > 0) 
      {
        dirichfacs[7] =  actdyn->m_damp;
        dirichfacs[8] =  actdyn->k_damp;
      } 
      else 
      {
        dirichfacs[7] =  0.0;
        dirichfacs[8] =  0.0;
      }
      dirichfacs[9] =  actdyn->dt;

      /*----------------------------------------------------------------*/
      /* zero the stiffness matrix 
       * and vector for internal forces and dirichlet forces */
      solserv_zero_mat(actintra, 
                       &(actsolv->sysarray[stiff_array]),
                       &(actsolv->sysarray_typ[stiff_array]));
      solserv_zero_mat(actintra,
                       &(actsolv->sysarray[mass_array]),
                       &(actsolv->sysarray_typ[mass_array]));
      amzero(&intforce_a);
      amzero(&dirich_a);

      /*----------------------------------------------------------------*/
      /* call element routines for calculation of 
       * tangential stiffness and intforce */
      *action = calc_struct_nlnstiffmass;
      solserv_sol_zero(actfield, disnum, node_array_sol_increment,
                       isolinc->fintn);
      container.dvec          = intforce;
      container.dirich        = dirich;
      container.global_numeq  = numeq_total;
      container.dirichfacs    = dirichfacs;
      container.kstep         = 0;
      calelm(actfield, actsolv, actpart, actintra, 
             stiff_array, mass_array, &container, action);

      /*----------------------------------------------------------------*/
      /*  store positive internal forces on fie[2] */
      solserv_zero_vec(&fie[2]);
      assemble_vec(actintra, &(actsolv->sysarray_typ[stiff_array]),
                   &(actsolv->sysarray[stiff_array]),
                   &(fie[2]), intforce, 1.0);

      /*----------------------------------------------------------------*/
      /* mid external force by interpolating */
      /* rhs[0] = (1-alphaf)rhs[1] + alphaf*rhs[2] */
      solserv_copy_vec(&(actsolv->rhs[2]), &(actsolv->rhs[0]));
      solserv_scalarprod_vec(&(actsolv->rhs[0]), actdyn->alpha_f);
      solserv_add_vec(&(actsolv->rhs[1]),
                      &(actsolv->rhs[0]),
                      (1.0 - actdyn->alpha_f));

      /*----------------------------------------------------------------*/
      /* mid internal force by interpolating */
      /*  fie[0] = (1-alfaf)fie[2] + alphaf*fie[1] */
      solserv_copy_vec(&fie[2],&fie[0]);
      solserv_scalarprod_vec(&fie[0], (1.0-actdyn->alpha_f));
      solserv_add_vec(&fie[1], &fie[0], actdyn->alpha_f);

      /*----------------------------------------------------------------*/
      /* subtract mid internal forces from mid external forces */
      solserv_add_vec(&fie[0], &(actsolv->rhs[0]), -1.0);

      /*----------------------------------------------------------------*/
      /* add Dirichlet forces from prescribed displacements */
      /* ===> GENERALLY THIS SHOULD BE WRONG!!! --- HOWEVER, CCARAT MAY
       *      NEED IT ???? */
      assemble_vec(actintra,
                   &(actsolv->sysarray_typ[stiff_array]),
                   &(actsolv->sysarray[stiff_array]),
                   &(actsolv->rhs[0]), dirich, 1.0);

      /*----------------------------------------------------------------*/
      /*  create effective load vector (rhs[0]-fie[0])eff */
      pefnln_struct(&dynvar, actdyn, actfield, actsolv, actintra,
                    dispi, vel, acc, work, mass_array, damp_array);

      /*----------------------------------------------------------------*/
      /* create effective stiffness matrix */
      kefnln_struct(&dynvar, actdyn, actfield, actsolv, actintra,
                    work, stiff_array, mass_array, damp_array);

      /*----------------------------------------------------------------*/
      /* solve keff * rsd[0] = rhs[0] */
      /* solve for residual displacements 
       * to correct iterative incremental displacements*/
      init = 0;
      solver_control(actfield, disnum, actsolv, actintra,
                     &(actsolv->sysarray_typ[stiff_array]),
                     &(actsolv->sysarray[stiff_array]),
                     &(work[0]), &(actsolv->rhs[0]), init);

      /*----------------------------------------------------------------*/
      /* return residual displacements to the nodes */
      solserv_result_resid(actfield, disnum, actintra, 
                           &(work[0]),
                           isolres->disres,
                           &(actsolv->sysarray[stiff_array]),
                           &(actsolv->sysarray_typ[stiff_array]));

      /*================================================================*/
      /* update */

      /*----------------------------------------------------------------*/
      /* update the incremental displacements by the residual 
       * displacements */
      solserv_add_vec(&(work[0]), &(dispi[0]), 1.0);

      /*----------------------------------------------------------------*/
      /* update displacements : sol[1] = sol[0] + dispi[0] */
      solserv_copy_vec(&(actsolv->sol[0]), &(actsolv->sol[1]));
      solserv_add_vec(&dispi[0], &(actsolv->sol[1]), 1.0);

      /*----------------------------------------------------------------*/
      /* return total displacements to the nodes */
      solserv_result_total(actfield, disnum, actintra, 
                           &(actsolv->sol[1]),
                           isol->disn,
                           &(actsolv->sysarray[stiff_array]),
                           &(actsolv->sysarray_typ[stiff_array]));

      /*----------------------------------------------------------------*/
      /* return incremental displacements to the nodes */
      solserv_result_incre(actfield, disnum, actintra,
                           &(dispi[0]),
                           isolinc->disinc,
                           &(actsolv->sysarray[stiff_array]),
                           &(actsolv->sysarray_typ[stiff_array]));

      /*================================================================*/
      /* CHECK CONVERGENCE */
      /* convergence = 0; */
      dmax        = 0.0;
      solserv_vecnorm_euclid(actintra, &(work[0]), &(dynvar.dinorm));
      solserv_vecnorm_euclid(actintra, &(dispi[0]), &(dynvar.dnorm));
      solserv_vecnorm_Linf(actintra, &(work[0]), &dmax);
      if (par.myrank == 0) 
      {
        printf("                                                   "
               "Residual %10.5E\n", dynvar.dinorm);
      }
      fflush(stdout);
      if ( (dynvar.dinorm < actdyn->toldisp)
           || (dynvar.dnorm < EPS14)
           || ( (dynvar.dinorm < EPS14) && (dmax < EPS12) ) )
      {
        printf("Convergence reached\n");
        convergence = 1;
      }

      /*----------------------------------------------------------------*/
      /* increase iteration counter */
      itnum = itnum + 1;

    }  /* end of equilibrium iteration */
    /*==================================================================*/
    /* END OF EQUILIBRIUM ITERATION */
    /*==================================================================*/

    /*------------------------------------------------------------------*/
    /* make temporary copy of actsolv->rhs[2] to actsolv->rhs[0] */
    /* (load at t-dt) because in  dyn_nlnstructupd actsolv->rhs[2] is
     *  overwritten but is  still needed to compute energies */
    solserv_copy_vec(&(actsolv->rhs[2]), &(actsolv->rhs[0]));

    /*------------------------------------------------------------------*/
    /*  copy disp from sol place 0 to place 10 */
    solserv_sol_copy(actfield, disnum, 
                     node_array_sol, node_array_sol, 
                     isol->disn, isol->dis);

    /*------------------------------------------------------------------*/
    /* copy vels from sol place 1 to place 11 */
    solserv_sol_copy(actfield, disnum,
                     node_array_sol, node_array_sol, 
                     isol->veln, isol->vel);

    /*------------------------------------------------------------------*/
    /* copy accs from sol place 2 to place 12 */
    solserv_sol_copy(actfield, disnum,
                     node_array_sol, node_array_sol, 
                     isol->accn, isol->acc);
    
    /*------------------------------------------------------------------*/
    /* update displacements, velocities and accelerations */
    dyn_nlnstructupd(actfield,
                     disnum,
                     &dynvar, actdyn, actsolv,
                     &(actsolv->sol[0]),/* total displ. at time t-dt */
                     &(actsolv->sol[1]),/* total displ. at time t    */
                     &(actsolv->rhs[1]),/* load vector  at time t    */
                     &(actsolv->rhs[2]),/* load vector  at time t-dt */
                     &(vel[0]),         /* velocities   at time t    */
                     &(acc[0]),         /* accelerat.   at time t    */
                     &(work[0]),        /* working arrays            */
                     &(work[1]),        /* working arrays            */
                     &(work[2]));       /* working arrays            */
    
    /*------------------------------------------------------------------*/
    /* return velocities to the nodes */
    solserv_result_total(actfield, disnum, actintra, 
                         &(vel[0]), isol->veln, 
                         &(actsolv->sysarray[stiff_array]),
                         &(actsolv->sysarray_typ[stiff_array]));

    /*------------------------------------------------------------------*/
    /* velocities for prescribed dofs to velocities */
    solserv_adddirich(actfield, disnum, 
                      node_array_sol, isol->veldn, isol->disn, isol->veln, 
                      1.0, 0.0);

    /*------------------------------------------------------------------*/
    /* return accel. to the nodes */
    solserv_result_total(actfield, disnum, actintra, 
                         &(acc[0]), isol->accn, 
                         &(actsolv->sysarray[stiff_array]),
                         &(actsolv->sysarray_typ[stiff_array]));
    
    /*------------------------------------------------------------------*/
    /* accel. for prescribed dofs */
    solserv_adddirich(actfield, disnum, 
                      node_array_sol, isol->accdn, isol->disn, isol->accn, 
                      1.0, 0.0);

    /*------------------------------------------------------------------*/
    /* It is a bit messed up, but anyway:
     * in the nodes the results are stored the following way:
     * 
     * in ARRAY sol.a.da[place][0..numdf-1]:
     * place 0  holds total displacements  time t      (free/prescr)
     * place 1  holds velocities           time t      (free/prescr)
     * place 2  holds accels               time t      (free/prescr)
     * place 3  holds displacements        time t-dt   (prescr only)
     * place 4  holds displacements        time t      (prescr only)
     * place 5  holds place 4 - place 3
     * place 6  holds velocities           time t      (prescr only)
     * place 7  holds accels               time t      (prescr only)
     * place 8  is working space
     * place 9  holds contact forces       time t      (free only)
     * place 10 holds total displacements  time t-dt   (free/prescr)
     * place 11 holds velocities           time t-dt   (free/prescr)
     * place 12 holds accels               time t-dt   (free/prescr)
     * 
     * in ARRAY sol_increment.a.da[place][0..numdf-1]
     * place 0 holds converged incremental displacements 
     *         (without prescribed dofs)
     * place 1 holds converged internal forces at time t-dt
     * place 2 holds converged internal forces at time t
     * 
     * in ARRAY sol_residual
     * place 0 holds residual displacements during iteration 
     *         (without prescribed dofs) */

    /*------------------------------------------------------------------*/
    /* make incremental potential energy at the nodes */
    dyn_epot(actfield, disnum, actintra, &dynvar, &deltaepot);
    dynvar.epot += deltaepot;

    /*------------------------------------------------------------------*/
    /* make kinetic energy at element level */
    dyn_ekin(actfield, actsolv, actpart, actintra, action, &container,
             stiff_array, mass_array);
    dynvar.ekin = container.ekin;

    /*------------------------------------------------------------------*/
    /* make external energy */
    dyn_eout(&dynvar, actdyn, actintra, actsolv, 
             &(dispi[0]), &(actsolv->rhs[1]), 
             &(actsolv->rhs[0]), &(work[0]));

    /*------------------------------------------------------------------*/
    /* make total energy */
    dynvar.etot = dynvar.epot + dynvar.ekin;

    /*------------------------------------------------------------------*/
    /* update the internal forces in sol_increment */
    /* copy from sol_increment.a.da[2][i] to sol_increment.a.da[1][i] */
    solserv_sol_copy(actfield, disnum,
                     node_array_sol_increment, node_array_sol_increment,
                     isolinc->fintn, isolinc->fint);

    /*------------------------------------------------------------------*/
    /* check whether to write results or not */
    mod_disp      = actdyn->step % actdyn->updevry_disp;
    mod_stress    = actdyn->step % actdyn->updevry_stress;

    /*------------------------------------------------------------------*/
    /* check whether to write restart or not */
    mod_res_write = actdyn->step % actdyn->res_write_evry;

    /*------------------------------------------------------------------*/
    /* perform stress calculation */
    if ( (mod_stress == 0) || (mod_disp == 0) )
    {
      if (ioflags.struct_stress == 1)
      {
        *action = calc_struct_stress;
        container.dvec          = NULL;
        container.dirich        = NULL;
        container.global_numeq  = 0;
        container.dirichfacs    = NULL;
        container.kstep         = 0;
        calelm(actfield, actsolv, actpart, actintra,
               stiff_array, -1, &container, action);
        /* reduce stresses, so they can be written */
        *action = calc_struct_stressreduce;
        container.kstep = 0;
        calreduce(actfield, actpart, disnum, actintra, action, &container);
      }  /* end of if (ioflags.struct_stress == 1) */
    }  /* end of  if ( (mod_stress == 0) || (mod_disp == 0) ) */

    /*------------------------------------------------------------------*/
    /* print out results to out */
    if ( (mod_stress == 0) || (mod_disp == 0) )
    {
      if ( (ioflags.struct_stress == 1) 
           && (ioflags.struct_disp == 1) 
           && (ioflags.output_out == 1) )
      {
        out_sol(actfield, actpart, disnum, actintra, actdyn->step, 0);
      }  /* end of if */
    }  /* end of if */

    /*------------------------------------------------------------------*/
    /* printout results to gid no time adaptivity */
    if ( (timeadapt == 0) 
         && (par.myrank == 0) 
         && (ioflags.output_gid == 1) )
    {
      if ( (mod_disp == 0) && (ioflags.struct_disp == 1) )
      {
        out_gid_soldyn("displacement", actfield, disnum_s, actintra, 
                       actdyn->step, 0, actdyn->time);
        /*out_gid_soldyn("velocity",actfield,disnum_s,actintra,actdyn->step,1,actdyn->time);*/
        /*out_gid_soldyn("accelerations",actfield,disnum_s,actintra,actdyn->step,2,actdyn->time);*/
      }  /* end of if */
      if ( (mod_stress == 0) && (ioflags.struct_stress == 1) )
      {
        out_gid_soldyn("stress", actfield, disnum_s, actintra, 
                       actdyn->step, 0, actdyn->time);
      }  /* end of if */
    }  /* end of if  */

    /*------------------------------------------------------------------*/
    /* write restart data to pss file */
    if (mod_res_write == 0)
    {
#ifdef BINIO
      restart_write_bin_nlnstructdyn(&out_context,
                                     actdyn, &dynvar,
                                     actsolv->nrhs, actsolv->rhs,
                                     actsolv->nsol, actsolv->sol,
                                     1            , dispi       ,
                                     1            , vel         ,
                                     1            , acc         ,
                                     3            , fie         ,
                                     3            , work);
#else
      restart_write_nlnstructdyn(actdyn, &dynvar, actfield, actpart,
                                 actintra, action,
                                 actsolv->nrhs, actsolv->rhs,
                                 actsolv->nsol, actsolv->sol,
                                 1            , dispi       ,
                                 1            , vel         ,
                                 1            , acc         ,
                                 3            , fie         ,
                                 3            , work        ,
                                 &intforce_a,
                                 &dirich_a,
                                 &container);
#endif
    }

    /*------------------------------------------------------------------*/
    /* print time step */
    if (par.myrank==0 && !timeadapt)
    {
      dyn_nlnstruct_outstep(&dynvar, actdyn, itnum, actdyn->dt);
    }

    /*------------------------------------------------------------------*/
    /*  measure wall clock time for this step */
    t1 = ds_cputime();
    fprintf(allfiles.out_err, "TIME for step %d is %f sec\n", 
            actdyn->step, t1-t0);

    }  /* end of time loop */
  /*====================================================================*/
  /* END OF TIME STEP LOOP */
  /*====================================================================*/
  
  /*--------------------------------------------------------------------*/
  end:
  /*--------------------------------------------------------------------*/
  /* cleaning up phase */
  amdel(&intforce_a);
  amdel(&dirich_a);
  solserv_del_vec(&(actsolv->rhs), actsolv->nrhs);
  solserv_del_vec(&(actsolv->sol), actsolv->nsol);
  solserv_del_vec(&dispi, 1);
  solserv_del_vec(&vel, 1);
  solserv_del_vec(&acc, 1);
  solserv_del_vec(&fie, 3);
  solserv_del_vec(&work, 3);
  /* clean BINIO */
#ifdef BINIO
  destroy_bin_out_field(&out_context);
#endif
  /* clean PARALLEL */
#ifndef PARALLEL
  CCAFREE(actintra);
#endif

  /*--------------------------------------------------------------------*/
  /* a last word to the nervously waiting user */
  printf("---------------------------------------------------------------"
         "---------\n");
  printf("TSI structural time integration generalised-alpha finished.\n");
  printf("==============================================================="
         "=========\n");

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}  /* end of tsi_st_genalp */
