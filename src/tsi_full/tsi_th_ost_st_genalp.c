/*======================================================================*/
/*!
\file
\brief Thermal-structure interaction
       - semi-coupled & staggered
       - instationary temperature field solved by one-step-thete
       - dynamic displacement field solved by generalised-alpha

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>

\author bborn
\date 05/07
*/


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
       This structure struct _FILES allfiles is defined in 
       input_control_global.c and the type is in standardtypes.h
       It holds all file pointers and some variables needed for the FRSYSTEM
\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern FILES allfiles;


/*----------------------------------------------------------------------*/
/*!
\brief General problem data
       defined in global_control.c 
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


/*======================================================================*/
/*!
\brief Staggered, semi-coupled thermo-structure interaction
       Static solution of thermal field (redone in every time step)
       Dynamic solution of structural field (with generalised-alpha)
\param   disnum_s     INT    (i)   index of structural discretisation
\param   disnum_t     INT    (i)   index of thermal discretisation
\author bborn
\date 03/06
*/
void tsi_th_ost_st_genalp(INT disnum_s,
                          INT disnum_t)
{
  /* field */
  INT numfld = genprob.numfld;  /* number of fields, i.e. ==2 */
  INT numsf = genprob.numsf;  /* number (index) of structure field */
  INT numtf = genprob.numtf;  /* number (index) of thermal field */
  FIELD* field_s = &(field[numsf]);  /* active structural field */
  FIELD* field_t = &(field[numtf]);  /* active thermal field */
  
  /* partition */
  PARTITION* part_s = &(partition[numsf]);  /* structure partition */
  PARTITION* part_t = &(partition[numtf]);  /* thermal partition */

  /* intra-communicator */
  INTRA* intra;  /* active intra-communicator */
  INTRA* intra_s;  /* active structure intra-communicator */
  INTRA* intra_t;  /* active thermal intra-communicator */

  /* discretisation */
  /* disnum_s is curr. discretisation index */
  /* structure discretisation index:
   * disnum_s == 0,
   * as only a single discretisation
   * this variable name is a little confusing
   * it sets an _actual_ number (better index?)
   * of one of the actfield->ndis discretisations.
   * Simply said: disnum != n(um)dis */
  /* named positions of NODE sol etc. arrays */
  ARRAY_POSITION* ipos_s = &(field_s->dis[disnum_s].ipos);
  /* named positions (indices) of NODE sol array */
  ARRAY_POSITION_SOL* isol_s = &(ipos_s->isol);
  /* named positions (indices) of NODE sol_increment array */
  ARRAY_POSITION_SOLINC* isolinc_s = &(ipos_s->isolinc);
  /* named positions (indices) of NODE sol_residual array */
  ARRAY_POSITION_SOLRES* isolres_s = &(ipos_s->isolres);
  /* disnum_t == index (number) of discretisation */
  /* only a single discretisation
   * this variable name is a little confusing
   * it sets an _actual_ number (better index?)
   * of one of the actfield->ndis discretisations.
   * Simply said: disnum != n(um)dis */
  /* named positions of NODE sol etc. arrays */
  ARRAY_POSITION* ipos_t = &(field_t->dis[disnum_t].ipos);
  /* named positions (indices) of NODE sol array */
  ARRAY_POSITION_SOL* isol_t = &(ipos_t->isol);

  /* solution variable */
  SOLVAR* solv_s = &(solv[numsf]);  /* active solution structure */
  INT numeq_s;  /* number of equations on this proc */
  INT numeq_total_s;  /* total number of equations */
  INT stiff_array_s;  /* index of the active system sparse matrix */
  INT mass_array_s;  /* index of the active system sparse matrix */
  INT damp_array_s;  /* index of the active system sparse matrix */
  SOLVAR* solv_t = &(solv[numtf]);  /* pointer to field SOLVAR */
  INT numeq_t;  /* number of equations on this proc */
  INT numeq_total_t;  /* total number of equations on all procs */
  INT stiff_array_t;  /* active sparse system matrix in 
                       * actsolv->sysarray[] */
  INT mass_array_t;  /* index of the active system sparse matrix */
  
  /* dynamic control */
  STRUCT_DYNAMIC* sdyn = alldyn[numsf].sdyn;  /* structural dynamics */
  THERM_DYNAMIC* tdyn = alldyn[numtf].tdyn;  /* thermal dynamic control */
  TSI_DYNAMIC* tsidyn = alldyn[numfld].tsidyn;  /* TSI dynamics */
  STRUCT_DYN_CALC sdynvar;  /* variables for dynamic structural simulation */
  INT timeadapt = sdyn->timeadapt;  /* flag to switch time adaption on/off */

  /* output */
  INT mod_stdout = 0;  /* indicates whether to write to STDOUT */
  /* INT mod_res_write; */
#ifdef BINIO
  BIN_OUT_FIELD out_context_s;
  BIN_OUT_FIELD out_context_t;
#endif
  
  /* container */
  CONTAINER container_s;  /* transfers variables given
                         * in (this) solution technique
                         * to element level */
  CONTAINER container_t;  /* contains variables defined in container.h */

  /* iteration */
  INT converged;  /* convergence flag */
  INT itnum;  /* iterator */

  /* timing */
  DOUBLE t0, t1;  /* wall clock times for measuring */
  
  /* global vectors of structural field */
  const INT vel_num = 1;  /* number of veloc. vectors */
  DIST_VECTOR* vel;  /* total velocities */
  const INT acc_num = 1;  /* number of accel. vectors */
  DIST_VECTOR* acc;  /* total accelerations */
  const INT fie_num = 3;  /* number of force vectors */
  DIST_VECTOR* fie;  /* internal forces and working array */
  const INT dispi_num = 1;  /* number of increm. displ. vectors */
  DIST_VECTOR* dispi;  /* distributed vector to hold increm. displacments */
  const INT work_num = 3;  /* number of work vectors */
  DIST_VECTOR* work;  /* working vectors */
  ARRAY intforce_a;  /* redundant vect. of full length for internal forces */
  ARRAY dirforce_a;  /* redund. full length vector for Dirichlet-part of RHS */

  /* global vectors of thermal field */
  DIST_VECTOR* tem;  /* global temperature (1) vector */
  DIST_VECTOR* hint;  /* internal heat vector */
  DIST_VECTOR* hext;  /* external heat vector */
  DIST_VECTOR* hextn;  /* new external heat vector */
  ARRAY intheat_a;  /* redund. full length vect. for internal heat */
  ARRAY dirheat_a;   /* redund. full length vect. for Dirichlet-part of RHS */

  INT i;
  
  /*====================================================================*/
  /* begin body */

#ifdef DEBUG
  dstrc_enter("tsi_th_ost_st_genalp");
#endif

  /*====================================================================*/
  /* a word to the user */
  if (par.myrank == 0)
  { 
    printf("============================================================="
           "=============\n");
    printf("Thermo-structure interaction: staggered, semi-coupled\n");
    printf("Thermal in-stationary solution with one-step-theta\n");
    printf("Structural dynamic solution with generalised-alpha\n");
    printf("-------------------------------------------------------------"
           "-------------\n");
  }

  /*--------------------------------------------------------------------*/
  /* set up structure container */
  container_s.fieldtyp = field_s->fieldtyp;  /* field type : structure */
  container_s.isdyn = 1;  /* dynamic calculation */
  container_s.kintyp = 2;  /* kintyp  = 2: total Lagrangean */  /* ??? */
  container_s.disnum = disnum_s;
  container_s.disnum_s = disnum_s;  /* structure discretisation index (==0) */
  container_s.disnum_t = disnum_t;  /* thermo-discretisation index (==0) */

  /*--------------------------------------------------------------------*/
  /* set up thermal container */
  container_t.fieldtyp = field_t->fieldtyp;  /* thermal field */
  container_t.isdyn = 0;  /* static calculation */  /* ? */
  container_t.kintyp = 0;  /* kintyp  = 0: geo_lin */  /* ? */
  container_t.disnum = disnum_t;  /* current thermal discretisation index */
  container_t.disnum_t = disnum_t;  /* thermal discretisation index */
  container_t.disnum_s = disnum_s;  /* structural discretisation index */

  /*--------------------------------------------------------------------*/
  /* synchronise structural dynamic control */
  sdyn->nstep = tsidyn->nstep;
  sdyn->dt = tsidyn->dt;
  sdyn->updevry_disp = tsidyn->out_res_ev;
  sdyn->updevry_stress = tsidyn->out_res_ev;

  /*--------------------------------------------------------------------*/
  /* synchronise thermal dynamic control */
  tdyn->nstep = tsidyn->nstep;
  tdyn->dt = tsidyn->dt;
  tdyn->out_res_ev = tsidyn->out_res_ev;
  tdyn->gamma = tsidyn->th_gamma;

  /*--------------------------------------------------------------------*/
  /* check time step adaptivity */
  if (timeadapt)
  {
    dserror("Time step size adaptivity is not available!");
  }

  /*--------------------------------------------------------------------*/
  /* intra communicator */
#ifdef PARALLEL
  intra_s = &(par.intra[numsf]);
  intra_t = &(par.intra[numtf]);
#else
  intra = (INTRA*) CCACALLOC(numfld, sizeof(INTRA));
  intra_s = &(intra[numsf]);
  intra_t = &(intra[numtf]);
  if ( (!intra) || (!intra_s) || (!intra_t) )
  {
    dserror("Allocation of intra-communicator failed");
  }
  intra_s->intra_fieldtyp = structure;
  intra_s->intra_rank = 0;
  intra_s->intra_nprocs = 1;
  intra_t->intra_fieldtyp = thermal;
  intra_t->intra_rank = 0;
  intra_t->intra_nprocs = 1;
#endif
  /* there are only procs allowed in here, that belong to the structural
   * intracommunicator (in case of linear statics, this should be all) */
  if ( (intra_s->intra_fieldtyp != structure)
       || (intra_t->intra_fieldtyp != thermal) )
  {
    goto end;
  }

  /*====================================================================*/
  /* initialise thermal field */
  tsi_th_ost_init(part_t,
                  intra_t,
                  field_t,
                  disnum_t,
                  ipos_t,
                  solv_t,
                  &(numeq_t),
                  &(numeq_total_t),
                  &(stiff_array_t),
                  &(mass_array_t),
                  &(container_t),
                  &(tem),
                  &(hint),
                  &(hext),
                  &(hextn),
                  &(intheat_a),
                  &(dirheat_a));

  /*====================================================================*/
  /* initialise structural field */
  tsi_st_genalp_init(part_s,
                     intra_s,
                     field_s,
                     disnum_s,
                     ipos_s, isol_s, isolinc_s,
                     solv_s,
                     &(numeq_total_s), &(numeq_s),
                     &(stiff_array_s), &(mass_array_s), &(damp_array_s),
                     sdyn, &(sdynvar),
                     &(container_s),
#ifdef BINIO
                     &(out_context_s),
#endif
                     vel_num, &(vel),
                     acc_num, &(acc),
                     fie_num, &(fie),
                     dispi_num, &(dispi),
                     work_num, &(work),
                     &(intforce_a),
                     &(dirforce_a));
#if 0
  solserv_create_vec(&fie, fie_num, numeq_total_s, numeq_s, "DV");
  for (i=0; i<fie_num; i++)
  {
    solserv_zero_vec(&(fie[i]));
  }
  solserv_create_vec(&work, work_num, numeq_total_s, numeq_s, "DV");
  for (i=0; i<work_num; i++)
  {
    solserv_zero_vec(&(work[i]));
  }
#endif

  /*--------------------------------------------------------------------*/
  /* set initial step and time */
/*   acttime = tsidyn->time;  /\* initial time *\/ */
  tsidyn->step = -1;

  /*--------------------------------------------------------------------*/
  /* printout head */
/*   if (par.myrank == 0)  */
/*   { */
/*     dyn_nlnstruct_outhead(&sdynvar, sdyn); */
/*   } */


  /*====================================================================*/
  /* START LOOP OVER ALL TIME STEPS */
  /*====================================================================*/
  /*
   * rhs[3]    original load vector
   * rhs[2]             load vector at time t_{n}
   * rhs[1]             load vector at time t_{n+1}
   * rhs[0]    interpolated load vector and working array
   *
   * fie[2]    internal forces at step t_{n+1}
   * fie[1]    internal forces at step t_{n}
   * fie[0]    interpolated internal forces and working array
   *
   * dispi[0]  displacement increment \inc\D_{n+1} from t_{n} to t_{n+1}
   *
   * sol[0]    total displacements \D_{n} at time t_{n+1}
   * sol[1]    total displacements \D_{n} at time t_n
   *
   * vel[0]    velocities \V_{n} at t_{n}
   * acc[0]    accelerations \V_{n} at t_{n}
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
  while ( (tsidyn->step < tsidyn->nstep-1) 
          && (tsidyn->time <= tsidyn->maxtime) )
  {
    /*------------------------------------------------------------------*/
    /* wall clock time in the beginning of current time step*/
    t0 = ds_cputime();

    /*------------------------------------------------------------------*/
    /* increment step */
    tsidyn->step++;
    sdyn->step = tsidyn->step;
    tdyn->step = tsidyn->step;

    /*------------------------------------------------------------------*/
    /* set new time t_{n+1} */
    tsidyn->time += tsidyn->dt;
    sdyn->time = tsidyn->time;
    tdyn->time = tsidyn->time;

    /*------------------------------------------------------------------*/
    /* output to STDOUT */
    mod_stdout = tsidyn->step % tsidyn->out_std_ev;  

    /*==================================================================*/
    /* THERMAL FIELD  :  SOLUTION */
    /*==================================================================*/
    /* inform user */
    if ( (par.myrank == 0) && (mod_stdout == 0) )
    {
      printf("\nStep %d: Solve thermal field...\n", tdyn->step);
    }
#if 0
    {
      INT i;
      for (i=0; i<2; i++)
      {
        printf("temp %d: %24.10e\n", i, tem->vec.a.dv[i]);
      }
    }
#endif
    /* solve thermal field */
    tsi_th_ost_equi(part_t,
                    intra_t,
                    field_t,
                    disnum_t,
                    isol_t,
                    solv_t,
                    numeq_t, numeq_total_t,
                    stiff_array_t, mass_array_t,
                    tdyn,
                    &(container_t),
                    tem,
                    hint,
                    hext,
                    hextn,
                    &(intheat_a),
                    &(dirheat_a));
#if 0
    {
      INT i;
      for (i=0; i<2; i++)
      {
        printf("temp %d: %24.10e\n", i, tem->vec.a.dv[i]);
      }
    }
#endif


    /*==================================================================*/
    /* STRUCTURAL FIELD  :  PREDICTOR */
    /*==================================================================*/
    /* inform user */
    if ( (par.myrank == 0) && (mod_stdout == 0) )
    {
      printf("Step %d: Solve structural field...\n", sdyn->step);
    }
    /*------------------------------------------------------------------*/
    /* predictor */
    tsi_st_genalp_pred(part_s,
                       intra_s,
                       field_s,
                       disnum_s,
                       isol_s, isolinc_s, isolres_s,
                       solv_s,
                       numeq_total_s,
                       stiff_array_s, mass_array_s, damp_array_s,
                       sdyn,
                       &(sdynvar),
                       &(container_s),
                       vel,
                       acc,
                       fie,
                       dispi,
                       work,
                       &(dirforce_a),
                       &(intforce_a));
#if 0
    {
      INT i;
      for (i=0; i<numeq_t; i++)
      {
        printf("temp %d: %24.10e\n", i, dispi->vec.a.dv[i]);
      }
    }
#endif

    /*------------------------------------------------------------------*/
    /* convergence check */
    converged = 0;
    tsi_st_genalp_chkcnv(intra_s,
                         sdyn,
                         &(sdynvar),
                         dispi,
                         dispi,
                         mod_stdout,
                         &(converged));

    /*==================================================================*/
    /* STRUCTURAL FIELD  :  EQUILIBRIUM ITERATION */
    /*==================================================================*/
    itnum = 0;
    while ( (converged != 1) && (itnum <= sdyn->maxiter) )
    {

      /*----------------------------------------------------------------*/
      /* check if maximally permitted iterations reached */
      if ( (itnum == sdyn->maxiter) && !timeadapt )
      {
        printf("No convergence in maxiter steps\n");
      }

      /*----------------------------------------------------------------*/
      /* make the equilibrium iteration */
      tsi_st_genalp_equi(part_s,
                         intra_s,
                         field_s,
                         disnum_s,
                         isol_s, isolinc_s, isolres_s,
                         solv_s,
                         numeq_total_s,
                         stiff_array_s, mass_array_s, damp_array_s,
                         sdyn, &(sdynvar),
                         &(container_s),
                         vel,
                         acc,
                         fie,
                         dispi,
                         work,
                         &(intforce_a),
                         &(dirforce_a));

      /*----------------------------------------------------------------*/
      /* convergence check */
      tsi_st_genalp_chkcnv(intra_s,
                           sdyn,
                           &(sdynvar),
                           work,
                           dispi,
                           mod_stdout,
                           &(converged));

      /*----------------------------------------------------------------*/
      /* increase iteration counter */
      ++itnum;

    }  /* end of equilibrium iteration */
    /*==================================================================*/
    /* END OF EQUILIBRIUM ITERATION */
    /*==================================================================*/


    /*------------------------------------------------------------------*/
    /* incremental update of structural field */
    tsi_st_genalp_updincr(part_s,
                          intra_s,
                          field_s,
                          disnum_s,
                          isol_s, isolinc_s,
                          solv_s,
                          mass_array_s, stiff_array_s,
                          sdyn,
                          &(sdynvar),
                          &(container_s),
                          vel,
                          acc,
                          dispi,
                          work);

    /*==================================================================*/
    /* OUTPUT */
    /*==================================================================*/
    /*------------------------------------------------------------------*/
    /* output thermal field */
    tsi_th_ost_out(part_t,
                   intra_t,
                   field_t,
                   disnum_s,
                   isol_t,
                   solv_t,
                   stiff_array_t,
                   tdyn,
                   &(container_t));

    /*------------------------------------------------------------------*/
    /* output structural field */
    tsi_st_genalp_out(part_s,
                      intra_s,
                      field_s,
                      disnum_s,
                      isol_s,
                      solv_s,
                      stiff_array_s,
                      sdyn, &(sdynvar),
                      &(container_s),
#ifdef BINIO
                      &(out_context_s),
#endif
                      dispi_num, dispi,
                      vel_num, vel,
                      acc_num, acc,
                      fie_num, fie,
                      work_num, work,
                      &(intforce_a),
                      &(dirforce_a));

    /*------------------------------------------------------------------*/
    /* print time step to STDOUT*/
    if ( (par.myrank==0) && (!timeadapt) && (mod_stdout==0) )
    {
      dyn_nlnstruct_outstep(&(sdynvar), sdyn, itnum, sdyn->dt);
    }

    /*------------------------------------------------------------------*/
    /*  measure wall clock time for this step */
    t1 = ds_cputime();
    fprintf(allfiles.out_err, "TIME for step %d is %f sec\n", 
            tsidyn->step, t1-t0);

  }  /* end of time loop */
  /*====================================================================*/
  /* END OF TIME STEP LOOP */
  /*====================================================================*/
  
  /*--------------------------------------------------------------------*/
  /* the end, my friend */
  end:

  /*====================================================================*/
  /* deallocate thermo stuff */
  tsi_th_ost_final(solv_t,
#ifdef BINIO
                   &(out_context_t),
#endif
                   &(tem),
                   &(hint),
                   &(hext),
                   &(hextn),
                   &(intheat_a),
                   &(dirheat_a));

  /*====================================================================*/
  /* deallocate structure stuff */
  tsi_st_genalp_final(solv_s,
#ifdef BINIO
                      &(out_context_s),
#endif
                      dispi_num, &(dispi),
                      vel_num, &(vel),
                      acc_num, &(acc),
                      fie_num, &(fie),
                      work_num, &(work),
                      &(intforce_a),
                      &(dirforce_a));

  /* clean PARALLEL */
#ifndef PARALLEL
  CCAFREE(intra);
#endif

  /*====================================================================*/
  /* a last word to the nervously waiting user */
  if (par.myrank == 0)
  {
    printf("-------------------------------------------------------------"
           "-----------\n");
    printf("Finished: Thermo-structure interaction\n");
    /* printf("Thank you for your attention.\n"); */
    printf("============================================================="
           "===========\n");
  }

#ifdef DEBUG
  dstrc_exit();
#endif

  return;
}  /* end of tsi_st_genalp */
