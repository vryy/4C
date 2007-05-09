/*----------------------------------------------------------------------*/
/*!
\file
\brief Fehlberg4/5 time integration

Features:
   - explicit
   - 4th order accurate

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 02/07
*/

/*----------------------------------------------------------------------*/
/* header files */
#include "../headers/standardtypes.h"
#include "../solver/solver.h"
#include "../io/io.h"
#include "tsi_prototypes.h"
#include "tsi_fehlbg4.h"



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

\author bborn
\date 03/06
*/
extern INT numcurve;
extern CURVE *curve;

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
DOUBLE deltat;


/*---------------------------------------------------------------------*/
/*!
\brief Fehlberg4 parameters

\author bborn
\date 03/07
*/
TSI_FEHLBG4 tsi_fehlbg4;


/*======================================================================*/
/*!
\brief Central differences time integration of structural field

\author bborn
\date 02/07
*/
void tsi_st_fehlbg(INT disnum_s,
                   INT disnum_t)

{
  INT i; /* simply a counter */
  INT istg, jstg;  /* stage indices */
  INT numeq;  /* number of equations on this proc */
  INT numeq_total;  /* total number of equations */
  INT init;  /* flag for solver_control call */
  INT mod_disp, mod_stress;
  INT mod_res_write;
  INT restart;
  DOUBLE t0_res, t1_res;
  INT timeadapt;  /* flag to switch time adaption on/off */

  /* Fehlberg4 parameters */
  const INT nstg = tsi_fehlbg4_stages();

  DOUBLE dt;
  INT nstep;

  DOUBLE tim[5] = { 0., 0., 0., 0., 0. };

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
  STRUCT_DYNAMIC *sdyn;  /* pointer to structural dynamic input data */
  TSI_DYNAMIC *tsidyn;  /* pointer to TSI dynamics data */

  DIST_VECTOR *vel;  /* total velocities: V_n and V_{n+1} */
  DIST_VECTOR *acc;  /* total accelerations: A_n and A_{n+1} */
  DIST_VECTOR *fint;  /* internal force: 
                       * F_{Int;n}, F_{Int;n+1} and F_{Int;n+c_i} */
  DIST_VECTOR *fvsc;  /* viscous force vector */
  DIST_VECTOR *fext;  /* external force vector */
  DIST_VECTOR *dispi;   /* distributed vector to hold incremental 
                         * displacments */
  DIST_VECTOR *work;  /* working vectors */

  DIST_VECTOR *accn;  /* acceleration at stages 
                       * (in RK terms, these are generally called 'k') */
  DIST_VECTOR *veln;  /* velocity stages */
  DIST_VECTOR *disn;  /* displacement stages */

  ARRAY intforce_a;  /* redundant vector of full length for internal forces */
  DOUBLE *intforce;
  ARRAY dirich_a;  /* redundant vector of full length for dirichlet-part 
                    * of rhs */
  DOUBLE *dirich;

  STRUCT_DYN_CALC dynvar;  /* variables to perform dynamic structural 
                            * simulation */

  ARRAY_POSITION *ipos;   /* named positions of NODE sol etc. arrays */
  ARRAY_POSITION_SOL *isol;  /* named positions (indices) of NODE sol array */
  ARRAY_POSITION_SOLINC *isolinc;  /* named positions (indices) of 
                                    * NODE sol_increment array */
  ARRAY_POSITION_SOLRES *isolres;  /* named positions (indices) of 
                                    * NODE sol_residual array */

#ifdef BINIO
  BIN_OUT_FIELD out_context;
#endif

  INT disnum = 0;

  CONTAINER container;  /* contains variables defined in container.h */

  /*--------------------------------------------------------------------*/
  container.isdyn = 1;   /* dynamic calculation */

  /*--------------------------------------------------------------------*/
  /* set Fehlberg4 parameters */
  tsi_fehlbg4_setpara(&tsi_fehlbg4);

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_st_fehlbg");
#endif
  
  /*--------------------------------------------------------------------*/
  /* set up */
  restart = genprob.restart;
  /* set some pointers */
  numsf = genprob.numsf;  /* index of structure field */
  actfield = &(field[numsf]);  /* pointer to structure field */
  sdyn = alldyn[numsf].sdyn;  /* dynamic control */
  timeadapt = sdyn->timeadapt;  /* flag for time adaptivity */
  container.fieldtyp = actfield->fieldtyp;
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
  tsidyn = alldyn[genprob.numfld].tsidyn;  /* TSI dynamic control */
  tsidyn->dt = sdyn->dt;
/*   actsysarray = numsf;  /\* ? *\/ */
  actsolv = &(solv[numsf]);
  actpart = &(partition[numsf]);
  action = &(calc_action[numsf]);
  /* this control routine at the moment always uses load curve 0 */
  actcurve = 0;

  /*--------------------------------------------------------------------*/
  /* intra communicator */
#ifdef PARALLEL
  actintra = &(par.intra[0]);  /* ???? => numsf ??? */
#else
  /* if we are not parallel, we have to allocate an alibi 
   * intra-communicator structure */
  actintra = (INTRA*) CCACALLOC(1, sizeof(INTRA));
  if (!actintra)
  {
    dserror("Allocation of INTRA failed");
  }
  actintra->intra_fieldtyp = structure;
  actintra->intra_rank = 0;
  actintra->intra_nprocs = 1;
#endif
  /* there are only procs allowed in here, that belong to the structural
   * intracommunicator (in case of nonlinear struct. dyn., this should 
   * be all) */
  if (actintra->intra_fieldtyp != structure) 
  {
    goto end;
  }

  /*--------------------------------------------------------------------*/
  /* init the variables in dynvar to zero */
  memset(&dynvar, 0, sizeof(STRUCT_DYN_CALC));

  /*--------------------------------------------------------------------*/
  /* initial time */
  acttime = 0.0;  /* This is typically an input parameter, isn't it? */
  sdyn->time = acttime;

  /*--------------------------------------------------------------------*/
  /* check presence of damping matrix
   * and set indices of stiffness and mass sparse matrices */
  if (sdyn->damp == 1)  /* damped */
  {
    actsolv->nsysarray = 3;
    stiff_array = 0;
    mass_array = 1;
    damp_array = 2;
  }
  else  /* undamped */
  {
    actsolv->nsysarray = 2;
    stiff_array = 0;
    mass_array = 1;
    damp_array = -1;
  }

  /*--------------------------------------------------------------------*/
  /* stiff_array already exists, so copy the mask of it to */
  /* mass_array (and damp_array if needed) */
  /* reallocate the vector of sparse matrices and the vector of there types */
  /* formerly lenght 1, now lenght 2 or 3 dependent on presence of 
   * damp_array */
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
  /* copy the matrices sparsity mask from stiff_array to damp_array */
  if (damp_array>0)
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
                     &numeq,
                     &numeq_total);

  /*--------------------------------------------------------------------*/
  /* allocate 4 dist. vectors for RHS */
  /* these hold original load vector, 
   *            load vector at time t,
   *            load vector at t-dt and
   *            interpolated load vector */
  actsolv->nrhs = 4;
  solserv_create_vec(&(actsolv->rhs), 
                     actsolv->nrhs, numeq_total, numeq, "DV");
  for (i=0; i<actsolv->nrhs; i++) 
  {
    solserv_zero_vec(&(actsolv->rhs[i]));
  }

  /*--------------------------------------------------------------------*/
  /* there are 2 solution vectors to hold total displ.
   * one at time t_{n+1} and one at time t_{n} */
  actsolv->nsol = 2;
  solserv_create_vec(&(actsolv->sol), 
                     actsolv->nsol, numeq_total, numeq, "DV");
  for (i=0; i<actsolv->nsol; i++)
  {
    solserv_zero_vec(&(actsolv->sol[i]));
  }

  /*--------------------------------------------------------------------*/
  /* there is one vector to hold incremental displacements */
  solserv_create_vec(&dispi, 1, numeq_total, numeq, "DV");
  for (i=0; i<1; i++) 
  {
    solserv_zero_vec(&(dispi[i]));
  }

  /*--------------------------------------------------------------------*/
  /* allocate two vector vel */
  solserv_create_vec(&vel, 2, numeq_total, numeq, "DV");
  solserv_zero_vec(&(vel[0]));
  solserv_zero_vec(&(vel[1]));

  /*--------------------------------------------------------------------*/
  /* allocate one vector acc */
  solserv_create_vec(&acc, 2, numeq_total, numeq, "DV");
  solserv_zero_vec(&(acc[0]));
  solserv_zero_vec(&(acc[1]));

  /*--------------------------------------------------------------------*/
  /* allocate 5 displacement, velocity and acceleration stage vectors */
  solserv_create_vec(&disn, nstg, numeq_total, numeq, "DV");
  solserv_create_vec(&veln, nstg, numeq_total, numeq, "DV");
  solserv_create_vec(&accn, nstg, numeq_total, numeq, "DV");
  for (istg=0; istg<nstg; istg++)
  {
    solserv_zero_vec(&(disn[istg]));
    solserv_zero_vec(&(veln[istg]));
    solserv_zero_vec(&(accn[istg]));
  }

  /*--------------------------------------------------------------------*/
  /* allocate viscous force vector */
  solserv_create_vec(&fvsc, 1, numeq_total, numeq, "DV");
  solserv_zero_vec(&(fvsc[0]));

  /*--------------------------------------------------------------------*/
  /* allocate external force vector */
  solserv_create_vec(&fext, 3, numeq_total, numeq, "DV");
  solserv_zero_vec(&(fext[0]));
  solserv_zero_vec(&(fext[1]));
  solserv_zero_vec(&(fext[2]));

  /*--------------------------------------------------------------------*/
  /* allocate one redundant vector intforce of full lenght
   * this is used by the element routines to assemble the 
   * internal forces*/
  intforce = amdef("intforce", &intforce_a, numeq_total, 1, "DV");
  /* allocate 3 DIST_VECTOR fie
   * to hold internal forces at t_n, t_{n+1} and in between */
  solserv_create_vec(&fint, 3, numeq_total, numeq, "DV");
  for (i=0; i<3; i++)
  {
    solserv_zero_vec(&(fint[i]));
  }

  /*--------------------------------------------------------------------*/
  /* create a vector of full length for dirichlet part of rhs */
  dirich = amdef("dirich", &dirich_a, numeq_total, 1, "DV");

  /*--------------------------------------------------------------------*/
  /* allocate three working vectors */
  /* By optimizing this routine one could live with one or two working
   * vectors, I needed three to make things straight-forward and easy */
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
   *   There will be no solver call on mass or damping array. */
  /* initialize solver */
  init = 1;
  solver_control(actfield, disnum, actsolv, actintra,
                 &(actsolv->sysarray_typ[stiff_array]),
                 &(actsolv->sysarray[stiff_array]),
                 &(dispi[0]),
                 &(actsolv->rhs[0]),
                 init);
  solver_control(actfield, disnum, actsolv, actintra,
                 &(actsolv->sysarray_typ[mass_array]),
                 &(actsolv->sysarray[mass_array]),
                 &work[0],
                 &work[1],
                 init);
  if (damp_array > 0)
  {
    solver_control(actfield, disnum, actsolv, actintra,
                   &(actsolv->sysarray_typ[damp_array]),
                   &(actsolv->sysarray[damp_array]),
                   &work[0],
                   &work[1],
                   init);
  }

  /*--------------------------------------------------------------------*/
  /* init the assembly for stiffness and for mass matrix
   * (damping is not assembled) */
  init_assembly(actpart, actsolv, actintra, actfield, stiff_array, disnum);
  init_assembly(actpart, actsolv, actintra, actfield, mass_array, disnum);

  /*--------------------------------------------------------------------*/
  /* init the element calculating routines */
  *action = calc_struct_init;
  calinit(actfield, actpart, action, &container);

  /*--------------------------------------------------------------------*/
  /* write output of mesh to gid */
  /* This has ALREADY been done */
/*   if (par.myrank==0 && ioflags.output_gid==1) */
/*   { */
/*     out_gid_msh(); */
/*   } */

  /*--------------------------------------------------------------------*/
  /* call elements to calculate stiffness and mass */
  /* REMARK: stiffness array may be unusable: not all materials 
   *         are linearised. In this case Rayleigh's damping matrix
   *         will be screwed up. */
  *action = calc_struct_nlnstiffmass;
  container.dvec = NULL;
  container.dirich = NULL;
  container.global_numeq = 0;
  container.dirichfacs = NULL;
  container.kstep = 0;
  calelm(actfield, actsolv, actpart, actintra, 
         stiff_array, mass_array, &container, action);

  /*--------------------------------------------------------------------*/
  /* calculate mass matrix and Rayleigh damping matrix */
  /*   C = k_damp * K_{T}(D_0) + m_damp * M */
  if (damp_array > 0)
  {
    /* stiffness-proportional contribution */
    solserv_add_mat(actintra,
                    &(actsolv->sysarray_typ[damp_array]),
                    &(actsolv->sysarray[damp_array]),
                    &(actsolv->sysarray_typ[stiff_array]),
                    &(actsolv->sysarray[stiff_array]),
                    sdyn->k_damp);
    /* mass-proportional contribution */
    solserv_add_mat(actintra,
                    &(actsolv->sysarray_typ[damp_array]),
                    &(actsolv->sysarray[damp_array]),
                    &(actsolv->sysarray_typ[mass_array]),
                    &(actsolv->sysarray[mass_array]),
                    sdyn->m_damp);
  }

  /*--------------------------------------------------------------------*/
  /* NODE's sol indices */
  ipos = &(actfield->dis[disnum].ipos);  /* position array */
  isol = &(ipos->isol);
  isolinc = &(ipos->isolinc);
  isolres = &(ipos->isolres);

  /*--------------------------------------------------------------------*/
  /* compute initial kinetic energy */
  dyne(&dynvar, actintra, actsolv, mass_array, &vel[0], &work[0]);

  /*--------------------------------------------------------------------*/
  /* set initial step and time */
  sdyn->step = 0;

  /*--------------------------------------------------------------------*/
  /* binary IO */
#ifdef BINIO
  /* initialize binary output
   * It's important to do this only after all the node arrays are set
   * up because their sizes are used to allocate internal memory. */
  init_bin_out_field(&out_context,
                     &(actsolv->sysarray_typ[stiff_array]), 
                     &(actsolv->sysarray[stiff_array]),
                     actfield, actpart, actintra, disnum);
#endif

  /*--------------------------------------------------------------------*/
  /* output to GID post-processor */
  if ( (par.myrank == 0) && (ioflags.output_gid == 1) )
  {
    /* out_gid_domains(actfield, disnum); */
  }

  /*--------------------------------------------------------------------*/
  /* printout head */
  if (par.myrank == 0)
  {
    /* a word to the user */
    printf("============================================================="
           "===========\n");
    printf("TSI structural time integration with Fehlberg4\n");
    printf("-------------------------------------------------------------"
           "-----------\n"); 
  }

  /*--------------------------------------------------------------------*/
  /* set some constants */
  /* dyn_setconstants_expl(&dynvar, sdyn, sdyn->dt); */
  dt = sdyn->dt;

  /*--------------------------------------------------------------------*/
  /* form effective left hand side */
  /*    K_eff = M 
   * The effective tangent is constant in time */
  solserv_zero_mat(actintra,
                   &(actsolv->sysarray[stiff_array]),
                   &(actsolv->sysarray_typ[stiff_array]));
  /* add mass contribution  M  */
  solserv_add_mat(actintra,
                  &(actsolv->sysarray_typ[stiff_array]),
                  &(actsolv->sysarray[stiff_array]),
                  &(actsolv->sysarray_typ[mass_array]),
                  &(actsolv->sysarray[mass_array]),
                  1.0);
  /*--------------------------------------------------------------------*/
  /* make triangulation of left hand side,
   * ie factorise/decompose "stiff_array" matrix */
  init = 0;
  solver_control(actfield, disnum, actsolv, actintra,
                 &(actsolv->sysarray_typ[stiff_array]),
                 &(actsolv->sysarray[stiff_array]),
                 &(dispi[0]),
                 &(actsolv->rhs[2]),
                 init);
  /* solver UMFPACK:  matrix format CCF */
  if (actsolv->sysarray_typ[stiff_array] == ccf)
  {
    if (actsolv->sysarray[stiff_array].ccf->is_factored == 1)
    {
      /* the factorised mass matrix will be readily used from now on */
      actsolv->sysarray[stiff_array].ccf->reuse = 1;
    }
    else
    {
      dserror("Cannot reuse factorised mass matrix: Factorisation failed");
    }
  }


  /*====================================================================*/
  /* determine initial accelerations
   *    M.A_0 + C.V_0 + F_{Int;0} = F_{Ext;0}  ==>  A_0 */

  /*--------------------------------------------------------------------*/
  /* set incremental displacements dispi[0] to zero */
  solserv_zero_vec(&dispi[0]);

  /*--------------------------------------------------------------------*/
  /* return total displacements to the nodes */
  solserv_result_total(actfield, disnum, actintra, 
                       &(dispi[0]), isol->disn,
                       &(actsolv->sysarray[stiff_array]),
                       &(actsolv->sysarray_typ[stiff_array]));

  /*--------------------------------------------------------------------*/
  /* return incremental displacements to the nodes */
  solserv_result_incre(actfield, disnum, actintra, 
                       &(dispi[0]), isolinc->disinc,
                       &(actsolv->sysarray[stiff_array]),
                       &(actsolv->sysarray_typ[stiff_array]));

  /*--------------------------------------------------------------------*/
  /* create initial external force vector F_{Ext;0}  (actsolv->rhs[1])*/
  /* the appropriate action is set inside calrhs */
  /* this vector holds loads due to external forces */
  container.kstep = 0;
  container.inherit = 1;
  container.point_neum = 1;
  *action = calc_struct_eleload;
  calrhs(actfield, actsolv, actpart, actintra, stiff_array,
         &(actsolv->rhs[1]), action, &container);
  /* copy the rhs vector */
  solserv_copy_vec(&(actsolv->rhs[1]), &(actsolv->rhs[3]));
  /*  dyn_init_curve(actcurve, sdyn->nstep, sdyn->dt, sdyn->maxtime); */
  /* get factor at a certain time t=0.0 */
  dyn_facfromcurve(actcurve, acttime, &(dynvar.rldfac));
  /* multiply global load vector by load factor */
  solserv_scalarprod_vec(&(actsolv->rhs[1]), dynvar.rldfac);

  /*--------------------------------------------------------------------*/
  /* initial internal forces F_{Int;0} (fie[2])*/
  amzero(&intforce_a);
  *action = calc_struct_internalforce;
  container.dvec = intforce;
  container.dirich = NULL;
  container.global_numeq = numeq_total;
  container.dirichfacs = NULL;
  container.kstep = 0;
  calelm(actfield, actsolv, actpart, actintra, stiff_array, -1,
         &container, action);
  /* store positive internal forces on fie[2] */
  solserv_zero_vec(&fint[2]);
  assemble_vec(actintra, &(actsolv->sysarray_typ[stiff_array]),
               &(actsolv->sysarray[stiff_array]), &(fint[2]), intforce, 1.0);

  /*--------------------------------------------------------------------*/
  /* initial inertial force (actsolv->rhs[0])
   *    F_{Inertial;0} = M . A_0
   *                   = F_{Ext;0} - F_{Int;0} - C . V_0 */
  /* copy rhs[0] <- rhs[1] */
  solserv_copy_vec(&(actsolv->rhs[1]), &(actsolv->rhs[0]));
  /* rhs[0] <- rhs[0] - fie[2] */
  solserv_add_vec(&(fint[2]), &(actsolv->rhs[0]), -1.0);
  /* subtract global viscous forces C . V_0 */
  if (sdyn->damp == 1)
  {
    /* copy work[0] <- vel[0] */
    solserv_copy_vec(&vel[0], &work[0]);
    /* work[1] := damp_array * work[0] */
    solserv_sparsematvec(actintra,
                         &work[1],
                         &(actsolv->sysarray[damp_array]),
                         &(actsolv->sysarray_typ[damp_array]),
                         &work[0]);
    /* subtract viscous forces of inertial forces
     * rhs[0] <- rhs[0] - work[1] */
    solserv_add_vec(&work[1], &(actsolv->rhs[0]), -1.0);
  }

  /*--------------------------------------------------------------------*/
  /* new accelerations A_{0} = M^{-1} . F_{Inertial;0} */
  /* solve for system */
  solserv_zero_vec(&work[0]);
  init = 0;
  solver_control(actfield, disnum, actsolv, actintra,
                 &(actsolv->sysarray_typ[mass_array]),
                 &(actsolv->sysarray[mass_array]),
                 &(work[0]),
                 &(actsolv->rhs[0]),
                 init);
  /* acc[0] <- work[0] */
  solserv_copy_vec(&(work[0]), &(acc[0]));
  /* debug: */ /*solserv_zero_vec(&acc[0]);*/

  /*====================================================================*/
  /* initialise loop */

  /*--------------------------------------------------------------------*/
  /* make norm of initial RHS */
  solserv_vecnorm_euclid(actintra, &(actsolv->rhs[1]), &(dynvar.rnorm));

  /*--------------------------------------------------------------------*/
  /* initialise wall clock time */
  t0 = ds_cputime();
  /* write memory report */
  if (par.myrank == 0) 
  {
    dsmemreport();
  }

  /*====================================================================*/
  /* START LOOP OVER ALL TIME STEPS */
  /*====================================================================*/
  while( (sdyn->step < sdyn->nstep) && (sdyn->time <= sdyn->maxtime) )
  {
    if (sdyn->step == 0)
    {
      printf("Entered time loop\n");
    }

    /* check for restart */
    if (restart)
    {
      t0_res = ds_cputime();
      /* save the stepsize as it will be overwritten in sdyn */
      dt = sdyn->dt;
      tsidyn->dt = dt;
      /* save the number of steps, as it will be overwritten in sdyn */
      nstep = sdyn->nstep;
      /* save the restart interval, as it will be overwritten */
      mod_res_write = sdyn->res_write_evry;
      /* the step to read in is restart */
#ifdef BINIO
      restart_read_bin_nlnstructdyn(sdyn, &dynvar,
                                    &(actsolv->sysarray_typ[stiff_array]),
                                    &(actsolv->sysarray[stiff_array]),
                                    actfield, actpart, 0, actintra,
                                    actsolv->nrhs, actsolv->rhs,
                                    actsolv->nsol, actsolv->sol,
                                    1            , dispi       ,
                                    1            , vel         ,
                                    1            , acc         ,
                                    3            , fint         ,
                                    3            , work        ,
                                    restart);
#else
      restart_read_nlnstructdyn(restart,
                                sdyn,
                                &dynvar,
                                actfield,
                                actpart,
                                actintra,
                                action,
                                actsolv->nrhs, actsolv->rhs,
                                actsolv->nsol, actsolv->sol,
                                1            , dispi       ,
                                1            , vel         ,
                                1            , acc         ,
                                3            , fint         ,
                                3            , work        ,
                                &intforce_a,
                                &dirich_a,
                                &container); /* contains variables defined 
                                              in container.h */
#endif
      /* put the dt to the structure */
      sdyn->dt = dt;
      /* put nstep to the structure */
      sdyn->nstep = nstep;
      /* put restart interval to structure */
      sdyn->res_write_evry = mod_res_write;
      /* switch the restart off */
      restart = 0;
      /* measure time */
      t1_res = ds_cputime();
      fprintf(allfiles.out_err,
              "TIME for restart reading of step %d is %f sec\n",
              sdyn->step,t1_res-t0_res);
    }  /* end if(restart) */


    /*==================================================================*/
    /* determine new displacements and velocities */
    /* initialise new displacements (sol[1])
     *    D_{n+1} = D_n + dt*V_n */
    solserv_copy_vec(&(actsolv->sol[0]), &(actsolv->sol[1]));
    solserv_add_vec(&(vel[0]), &(actsolv->sol[1]), dt);
    /* initialise new velocities (vel[1])
     *    V_{n+1} = V_n */
    solserv_copy_vec(&(vel[0]), &(vel[1]));
    /* loop over Runge-Kutta stages */
    for (istg=0; istg<nstg; istg++)
    {
      /*----------------------------------------------------------------*/
      /* time stage 
       *    t_{n+1} = t_n + dt*c_i */
      tim[istg] = acttime + dt*tsi_fehlbg4.c[istg];
      /*----------------------------------------------------------------*/
      /* velocity at stage
       *    V_{n+c_i} += dt_n * \sum{j=1}{s} ( a_{ij} * k_j ) */
      solserv_copy_vec(&(vel[0]), &(veln[istg]));
      for (jstg=0; jstg<istg; jstg++)
      {
        solserv_add_vec(&(accn[jstg]), &(veln[istg]), 
                        dt*tsi_fehlbg4.a[istg][jstg]);
      }
      /*----------------------------------------------------------------*/
      /* displacement at stage
       *    D_{n+c_i} += c_i * dt * V_n 
       *               + dt_n^2 * \sum{j=1}{i-1} ( aa_{ij} * k_j )*/
      solserv_copy_vec(&(actsolv->sol[0]), &(disn[istg]));
      solserv_add_vec(&(vel[0]), &(disn[istg]), dt*tsi_fehlbg4.c[istg]);
      for (jstg=0; jstg<istg-1; jstg++)
      {
        solserv_add_vec(&(accn[jstg]), &(disn[istg]), 
                        dt*dt*tsi_fehlbg4.aa[istg][jstg]);
      }
      /* put the displacements to the nodes */
      solserv_result_total(actfield, disnum, actintra, &(disn[istg]), 
                           isol->disn,
                           &(actsolv->sysarray[stiff_array]),
                           &(actsolv->sysarray_typ[stiff_array]));
      /* put the prescribed scaled displacements to the nodes
       * in field sol==node_array_sol at 0==isol->disn
       * this overwrites the calculated displacements on the Dirichlet BC */
      solserv_putdirich_to_dof(actfield, disnum, node_array_sol, 
                               isol->disn, tim[istg]);
      /*----------------------------------------------------------------*/
      /* right-hand-side at stage */
      /*----------------------------------------------------------------*/
      /* external force at t_{n+c_i} */
      /* new global load vector F_{ext;n+c_i} (fext[0])*/
      /* make load at time t_{n+1} in rhs[1] */
      solserv_zero_vec(&(fext[0]));
      container.kstep = 0;
      container.inherit = 1;
      container.point_neum = 1;
      *action = calc_struct_eleload;
      calrhs(actfield, actsolv, actpart, actintra, stiff_array,
             &(fext[0]), action, &container);
      dyn_facfromcurve(actcurve, tim[istg], &(dynvar.rldfac));
      solserv_scalarprod_vec(&(fext[0]), dynvar.rldfac);
      /*----------------------------------------------------------------*/
      /* internal force at t_{n+c_i} */
      /* new global internal force vector (fint[2])
       *    F_{Int;n+c_i} = F_{Int}(D_{n+c_i})
       * This is actually at t_{n+c_i} due to putting the stage 
       * displacemens to the nodes (see above) */
      amzero(&intforce_a);
      *action = calc_struct_internalforce;
      container.dvec = intforce;
      container.dirich = NULL;
      container.global_numeq = numeq_total;
      container.dirichfacs = NULL;
      container.kstep = 0;
      tsidyn->actstg = istg;
      calelm(actfield, actsolv, actpart, actintra, stiff_array, -1,
             &container, action);
      /* store positive internal forces on fint[2] */
      solserv_zero_vec(&fint[2]);
      assemble_vec(actintra, &(actsolv->sysarray_typ[stiff_array]),
                   &(actsolv->sysarray[stiff_array]), 
                   &(fint[2]), intforce, 1.0);
      /*----------------------------------------------------------------*/
      /* viscous force at t_{n+c_i} */
      if (sdyn->damp == 1)
      {
        solserv_zero_vec(&fvsc[0]);
        solserv_sparsematvec(actintra,
                             &fvsc[0],
                             &(actsolv->sysarray[damp_array]),
                             &(actsolv->sysarray_typ[damp_array]),
                             &veln[istg]);
      }
      /*----------------------------------------------------------------*/
      /* inertial force at stage (actsolv->rhs[0])
       *    F_{Inertial;n+c_i} = M . A_{n+c_i} 
       *                       = F_{Ext;n+c_i} - F_{Int;n+c_i} 
       *                         - C . V_{n+c_i} */
      solserv_copy_vec(&(fext[0]), &(actsolv->rhs[0]));
      solserv_add_vec(&(fint[2]), &(actsolv->rhs[0]), -1.0);
      solserv_add_vec(&(fvsc[0]), &(actsolv->rhs[0]), -1.0);
      /*----------------------------------------------------------------*/
      /* acceleration at stage A_{n+c_i} */
      /* new accelerations A_{n+1} = M^{-1} . F_{Inertial;n+c_i} */
      /* solve for system --- stiff_array is the mass matrix */
      solserv_zero_vec(&work[0]);
      init = 0;
      solver_control(actfield, disnum, actsolv, actintra,
                     &(actsolv->sysarray_typ[stiff_array]),
                     &(actsolv->sysarray[stiff_array]),
                     &(work[0]),  /* solution */
                     &(actsolv->rhs[0]),  /* RHS */
                     init);
      solserv_copy_vec(&(work[0]), &(accn[istg]));
      /*----------------------------------------------------------------*/
      /* add contribution of current stage to new displacement */
      solserv_add_vec(&(accn[istg]), &(actsolv->sol[1]), 
                      dt*dt*tsi_fehlbg4.bb[istg]);
      /* add contribution of current stage to new velocity */
      solserv_add_vec(&(accn[istg]), &(vel[1]), dt*tsi_fehlbg4.b[istg]);
    }  /* end of loop Runge-Kutta stages */

    /*------------------------------------------------------------------*/
    /* set new external force, last stage load is c_nstg==1 */
    solserv_copy_vec(&(fext[0]), &(fext[2]));

    /*------------------------------------------------------------------*/
    /* make all types of energies */
    dynnle(&dynvar, sdyn, actintra, actsolv, 
           &dispi[0], &fint[1], &fint[2],
           &(fext[1]), &(fext[2]), &work[0]);
    dyne(&dynvar, actintra, actsolv, mass_array, &vel[1], &work[0]);
    dynvar.etot = dynvar.epot + dynvar.ekin;  /* total */

    /*==================================================================*/
    /* update */
    /* update displacements D_{n} := D_{n+1} */
    solserv_copy_vec(&(actsolv->sol[1]), &(actsolv->sol[0]));
    /* update velocities V_{n} := V_{n+1} */
    solserv_copy_vec(&(vel[1]), &(vel[0]));
    /* update accelerations --- purely for output */
    solserv_copy_vec(&(accn[nstg-1]), &(acc[0]));
    /* update external forces --- purely for output */
    solserv_copy_vec(&(fext[2]), &(fext[1]));
    /* update time */
    sdyn->time += sdyn->dt;
    acttime = sdyn->time;
    /* increment step */
    (sdyn->step)++;
    /*------------------------------------------------------------------*/
    /* put new quantities to nodes */
    /* put new displacements sol[0] to nodes */
    solserv_result_total(actfield, disnum, actintra, 
                         &(actsolv->sol[0]), isol->disn,
                         &(actsolv->sysarray[stiff_array]),
                         &(actsolv->sysarray_typ[stiff_array]));
    /* return velocities to the nodes */
    solserv_result_total(actfield,disnum,actintra, 
                         &vel[0], isol->veln,
                         &(actsolv->sysarray[stiff_array]),
                         &(actsolv->sysarray_typ[stiff_array]));
    /* return accel. to the nodes */
    solserv_result_total(actfield,disnum,actintra, 
                         &acc[0], isol->accn,
                         &(actsolv->sysarray[stiff_array]),
                         &(actsolv->sysarray_typ[stiff_array]));
    /*------------------------------------------------------------------*/
    /* update material internal variables */
    *action = calc_struct_update_istep;
    container.dvec = NULL;
    container.dirich = NULL;
    container.global_numeq = 0;
    container.kstep = 0;
    calelm(actfield, actsolv, actpart, actintra,
           stiff_array, -1, &container, action);

    /*------------------------------------------------------------------*/
    /* check whether to write results or not */
    mod_disp = sdyn->step % sdyn->updevry_disp;
    /* check whether print stress or not */
    if (sdyn->updevry_stress > 0)
    {
      mod_stress = sdyn->step % sdyn->updevry_stress;
    }
    else
    {
      mod_stress = -1;
    }
    /* check whether to write restart or not */
    if (sdyn->res_write_evry > 0)
    {
      /* if mod_res_write becomes 0, i.e. current time step sdyn->step
       * is a integer multiple of sdyn->res_write_evry, the restart
       * will be written */
      mod_res_write = sdyn->step % sdyn->res_write_evry;
    }
    else
    {
      /* prevent the attempt to write a restart file */
      mod_res_write = -1;
    }
    /* perform stress calculation */
    if ( (mod_stress == 0) || (mod_disp == 0) )
    {
      if (ioflags.struct_stress == 1)
      {
        *action = calc_struct_stress;
        container.dvec = NULL;
        container.dirich = NULL;
        container.global_numeq = 0;
        container.dirichfacs = NULL;
        container.kstep = sdyn->step;
        calelm(actfield, actsolv, actpart, actintra,
               stiff_array, -1, &container, action);
        /* reduce stresses, so they can be written */
        *action = calc_struct_stressreduce;
        container.kstep = 0;
        calreduce(actfield, actpart, disnum, actintra, action, &container);
      }
    }
    /* print out results to out */
    if ( (mod_stress == 0) || (mod_disp == 0) )
    {
      if ( (ioflags.struct_stress == 1)
           && (ioflags.struct_disp == 1) 
           && (ioflags.output_out == 1) )
      {
        out_sol(actfield, actpart, disnum, actintra, sdyn->step, 0);
      }
    }
    /*------------------------------------------------------------------*/
    /* printout results to gid */
    /* binary output */
#ifdef BINIO
    if (ioflags.output_bin == 1)
    {
      if (mod_disp == 0)
      {
        if (ioflags.struct_disp == 1) 
        {
          out_results(&out_context, sdyn->time, sdyn->step, 0, 
                      OUTPUT_DISPLACEMENT);
          out_results(&out_context, sdyn->time, sdyn->step, 1, 
                      OUTPUT_VELOCITY);
          out_results(&out_context, sdyn->time, sdyn->step, 2, 
                      OUTPUT_ACCELERATION);
        }
      }
      if (mod_stress == 0)
      {
        if (ioflags.struct_stress == 1)
        {
          out_results(&out_context, sdyn->time, sdyn->step, 0, 
                      OUTPUT_STRESS);
        }
      }
    }
#endif
    /* ascii output */
    if ( (par.myrank == 0) && (ioflags.output_gid == 1) )
    {
      if (mod_disp == 0)
      {
        if (ioflags.struct_disp == 1)
        {
          out_gid_sol("displacement", actfield, disnum, actintra, 
                      sdyn->step, 0, ZERO);
          out_gid_sol("velocities", actfield, disnum, actintra, 
                      sdyn->step, 1, ZERO);
          out_gid_sol("accelerations", actfield, disnum, actintra,
                      sdyn->step, 2, ZERO);
        }
      }
      if (mod_stress == 0)
      {
        if (ioflags.struct_stress == 1)
        {
          out_gid_sol("stress", actfield, disnum, actintra, 
                      sdyn->step, 0, ZERO);
        }
      }
    }
    /* write restart data to pss file */
    if (mod_res_write == 0)
    {
#ifdef BINIO
      restart_write_bin_nlnstructdyn(&out_context,
                                     sdyn, &dynvar,
                                     actsolv->nrhs, actsolv->rhs,
                                     actsolv->nsol, actsolv->sol,
                                     1            , dispi       ,
                                     1            , vel         ,
                                     1            , acc         ,
                                     3            , fint        ,
                                     3            , work);
#else
      restart_write_nlnstructdyn(sdyn,
                                 &dynvar,
                                 actfield,
                                 actpart,
                                 actintra,
                                 action,
                                 actsolv->nrhs, actsolv->rhs,
                                 actsolv->nsol, actsolv->sol,
                                 1            , dispi       ,
                                 1            , vel         ,
                                 1            , acc         ,
                                 3            , fint        ,
                                 3            , work        ,
                                 &intforce_a,
                                 &dirich_a,
                                 &container);   /* contains variables 
                                                 * defined in container.h */
#endif
    }
    /*------------------------------------------------------------------*/
    /* print time step */
    if ( (par.myrank == 0) && (mod_disp == 0) )
    {
      /* dyn_nlnstruct_outstep(&dynvar, sdyn, 0, sdyn->dt); */
      /* to STDOUT */
      printf("STEP=%6d | NSTEP=%6d | TIME=%-14.8E | DT=%-14.8E | "
             "ETOT=%-14.8E | \n",
             sdyn->step, sdyn->nstep, sdyn->time, sdyn->dt, dynvar.etot);
      /* to ERR file */
      fprintf(allfiles.out_err,
              "STEP=%6d | NSTEP=%6d | TIME=%-14.8E | DT=%-14.8E | "
              "ETOT=%-14.8E | EPOT=%-14.8E | EKIN=%-14.8E | EOUT=%-14.8E\n",
              sdyn->step, sdyn->nstep, sdyn->time, sdyn->dt, 
              dynvar.etot, dynvar.epot, dynvar.ekin, dynvar.eout);
      fflush(allfiles.out_err);
    }
  }  /* end time loop */

  /*--------------------------------------------------------------------*/
  /* measure time for this step */
  t1 = ds_cputime();
  fprintf(allfiles.out_err, "TIME for step %d is %f sec\n", 
          sdyn->step, t1-t0);

  /*--------------------------------------------------------------------*/
  /* this is the end my lonely friend the end */
  end:

  /*--------------------------------------------------------------------*/
  /* cleaning up phase */
  amdel(&intforce_a);
  solserv_del_vec(&(actsolv->rhs), actsolv->nrhs);
  solserv_del_vec(&(actsolv->sol), actsolv->nsol);
  solserv_del_vec(&dispi, 1);
  solserv_del_vec(&vel, 1);
  solserv_del_vec(&acc, 1);
  solserv_del_vec(&fint, 3);
  solserv_del_vec(&work, 3);

  /*--------------------------------------------------------------------*/
  /* finalise the element calculating routines */
  *action = calc_struct_final;
  /* WARNING: This call might cause problems,
   *          however, a new (e.g. `calfinal') is postponed */
  calinit(actfield, actpart, action, &container);

  /*--------------------------------------------------------------------*/
  /* printout footer */
  if (par.myrank == 0)
  {
    /* a word to the user */
    printf("-------------------------------------------------------------"
           "-----------\n"); 
    printf("End: TSI structural time integration with Fehlberg4\n");
    printf("============================================================="
           "===========\n");
  }

#ifdef BINIO
  destroy_bin_out_field(&out_context);
#endif

  /*--------------------------------------------------------------------*/
#ifndef PARALLEL
  CCAFREE(actintra);
#endif
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
} /* end of dyn_nln_stru_expl */
