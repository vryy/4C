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

#ifndef CCADISCRET

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

This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM

\auther bborn
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
/*----------------------------------------------------------------------*/
extern GENPROB genprob;


/*----------------------------------------------------------------------*/
/*!
\brief Fields

vector of numfld FIELDs, defined in global_control.c

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
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
/*----------------------------------------------------------------------*/
extern PARTITION *partition;


/*----------------------------------------------------------------------*/
/*!
\brief Input/output control flags

structure of flags to control output defined in out_global.c

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
extern IO_FLAGS ioflags;


/*----------------------------------------------------------------------*/
/*!
\brief Rank and communicators (Parallelism!)

This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h

\author born
\date 03/06
*/
/*----------------------------------------------------------------------*/
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
/*----------------------------------------------------------------------*/
extern ALLDYNA* alldyn;


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
/*----------------------------------------------------------------------*/
extern INT numcurve;
extern CURVE* curve;


/*----------------------------------------------------------------------*/
/*!
\brief CALC_ACTIONs

enum _CALC_ACTION
command passed from control routine to the element level
to tell element routines what to do
defined globally in global_calelm.c

\author bborn
\date 03/06
/*----------------------------------------------------------------------*/
extern CALC_ACTION calc_action[MAXFIELD];


/*---------------------------------------------------------------------*/
/*!
\brief Actual (current) time globally given

\author bborn
\date 03/06
*/
/*---------------------------------------------------------------------*/
DOUBLE acttime;


/*----------------------------------------------------------------------*/
/*!
\brief

\author bborn
\date 03/06
*/
/*----------------------------------------------------------------------*/
void tsi_th_stat(INT disnum_s,
                 INT disnum_t)
{
  INT numtf = genprob.numtf;  /* number (index) of thermal field */
  FIELD* actfield = &(field[numtf]);  /* pointer to the thermal FIELD */
  PARTITION* actpart = &(partition[numtf]);  /* pointer to the fields PARTITION structure */

  INT i;  /* a counter */
  INT numeq;  /* number of equations on this proc */
  INT numeq_total;  /* total number of equations on all procs */
  INT init;  /* init flag for solver */
  INT actsysarray;  /* active sparse system matrix in
                     * actsolv->sysarray[] */
  INT disnum;  /* index (number) of discretisation */
  INT iplace;
  INT i2ndsysmat;  /* flag 2nd system matrix or not */
  DOUBLE rhsfact;  /* factor to multiply RHS */

  ARRAY dirich_a;
  DOUBLE* dirich;

  THERM_DYNAMIC* actdyn = alldyn[numtf].tdyn;  /* pointer to dynamic control
                                                       * --- allright, static control */
  TSI_DYNAMIC* tsidyn = alldyn[genprob.numfld].tsidyn;

  SOLVAR* actsolv = &(solv[numtf]);  /* pointer to field SOLVAR */


  INTRA* actintra;  /* pointer to the fields intra-communicator structure */
  CALC_ACTION *action = &(calc_action[numtf]);  /* pointer to the structures cal_action enum */

  DOUBLE timcur;  /* current time */
  DOUBLE timstp;  /* current time step */

  CONTAINER container;  /* contains variables defined in container.h */

  ARRAY_POSITION *ipos;  /* named positions of NODE sol etc. arrays */
  ARRAY_POSITION_SOL *isol;  /* named positions (indices) of NODE sol array */

  SPARSE_TYP array_typ;  /* type of psarse system matrix */

  DIST_VECTOR* initem;  /* initial temperature */

#ifdef BINIO
  BIN_OUT_FIELD out_context;
#endif



  /*====================================================================*/
  /* begin body */

#ifdef DEBUG
  dstrc_enter("tsi_th_stat");
#endif

  /*--------------------------------------------------------------------*/
  /* a word to the user */
  printf("==============================================================="
         "=========\n");
  printf("TSI quasi static thermal solution routine reached.\n");
  printf("---------------------------------------------------------------"
         "---------\n");

  /*--------------------------------------------------------------------*/
  /* link to global variables & a settings */
  container.fieldtyp = actfield->fieldtyp;
  container.isdyn = 0;  /* static calculation */  /* ? */
  container.kintyp = 0;  /* kintyp  = 0: geo_lin */  /* ? */
  disnum = disnum_t;  /* only a single discretisation
                       * this variable name is a little confusing
                       * it sets an _actual_ number (better index?)
                       * of one of the actfield->ndis discretisations.
                       * Simply said: disnum != n(um)dis */
  container.disnum = disnum;
  container.disnum_t = disnum_t;
  container.disnum_s = disnum_s;
  /* distributed system matrix, which is used for solving */
  actsysarray = numtf;  /* ? */
  if (actsolv->nsysarray == 1)
  {
    actsysarray = 0;
  }
  else
  {
    dserror("More than 1 system arrays (actsolv->nsysarray)!");
  }

  /*--------------------------------------------------------------------*/
  /* intra communicator */
#ifdef PARALLEL
  actintra = &(par.intra[0]);  /* ? */
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
  /* typ of global matrix */
  array_typ = actsolv->sysarray_typ[actsysarray];

  /* get global and local number of equations */
  solserv_getmatdims(&(actsolv->sysarray[actsysarray]),
                     actsolv->sysarray_typ[actsysarray],
                     &numeq, &numeq_total);

  /* number of rhs and solution vectors */
  actsolv->nrhs = 2;  /* WHY 2? */
  actsolv->nsol = 2;  /* WHY 2? */
  solserv_create_vec(&(actsolv->rhs), actsolv->nrhs,
                     numeq_total, numeq, "DV");
  solserv_create_vec(&(actsolv->sol), actsolv->nsol,
                     numeq_total, numeq, "DV");
  dirich = amdef("dirich", &dirich_a, numeq_total, 1, "DV");
  amzero(&dirich_a);

  /* init the created dist. vectors to zero */
  for (i=0; i<actsolv->nrhs; i++)
  {
    solserv_zero_vec(&(actsolv->rhs[i]));
  }
  for (i=0; i<actsolv->nsol; i++)
  {
    solserv_zero_vec(&(actsolv->sol[i]));
  }

  /*------------------------------------------------------------------*/
  /* set time controls */
  timcur = actdyn->time;
  timstp = actdyn->step;

  /*-------------------------------------------------------------------*/
  /* initialize solver */
  init = 1;
  solver_control(actfield, disnum, actsolv, actintra,
                 &(actsolv->sysarray_typ[actsysarray]),
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sol[actsysarray]),
                 &(actsolv->rhs[actsysarray]),
                 init);

  /*--------------------------------------------------------------------*/
  /* init the dist sparse matrix to zero */
  /* NOTE: Has to be called after solver_control(init=1) */
  solserv_zero_mat(actintra,
                   &(actsolv->sysarray[actsysarray]),
                   &(actsolv->sysarray_typ[actsysarray]));

  /*--------------------------------------------------------------------*/
  /* init the assembly for ONE sparse matrix */
  init_assembly(actpart, actsolv, actintra, actfield, actsysarray, disnum);

  /*--------------------------------------------------------------------*/
  /* init the element calculating routines */
  *action = calc_therm_init;
  calinit(actfield, actpart, action, &container);

  /*--------------------------------------------------------------------*/
  /* sol indices */
  ipos = &(actfield->dis[disnum].ipos);  /* position array */
  isol = &(ipos->isol);  /* position array of sol */

  /*--------------------------------------------------------------------*/
#ifdef BINIO
  /* initialize binary output */
  if (ioflags.output_bin == 1)
  {
    dserror("BINIO is not available!");
  }
#endif

  /*--------------------------------------------------------------------*/
  /* write output of mesh to gid */
  /* THIS IS NOT NECESSARY, AS ntacal ALREADY CALLS IT!
   * HOWEVER, IT IS DONE for prb_structure AND prb_fsi
   * SEPARATELY. POSSIBLY, IT WILL BE REQUIRED IN THE FUTURE */
  /*   if ( (par.myrank == 0) && (ioflags.output_gid == 1) ) */
  /*   { */
  /*     out_gid_msh(); */
  /*   } */

  /*--------------------------------------------------------------------*/
  /* put a zero to the place ipos->num=12 in node->sol to init the
   * velocities and accels of prescribed displacements */
  /* HINT: This actually redefines/reallocates/enlarges the sol
   *       array of each structure node to dimension 12x2 (or 12x3)
   *       from originally 1x2 (or 1x3) */
  solserv_sol_zero(actfield, disnum, node_array_sol, ipos->num-1);

  /*--------------------------------------------------------------------*/
  /* initial temperature in field */
  solserv_create_vec(&initem, 1, numeq_total, numeq, "DV");
  /* set initial temperature in field */
  actdyn->initmpr = tsidyn->th_initmpr;
  {
    INT idof;
    for (idof=0; idof<numeq; idof++)
    {
      initem->vec.a.dv[idof] = actdyn->initmpr;
    }
  }
  /* spawn initial temperature to nodes */
  solserv_result_total(actfield,
                       disnum,
                       actintra,
                       initem,
                       isol->tem0,
                       &(actsolv->sysarray[actsysarray]),
                       &(actsolv->sysarray_typ[actsysarray]));

  /*--------------------------------------------------------------------*/
  /* put the scaled prescribed temperatures to the nodes
   * in field sol (1st 0) at place 0 (2nd 0) together with
   * free temperatures */
  /* HINT: time curve is called indirectly */
  solserv_putdirich_to_dof(actfield, disnum,
                           node_array_sol, isol->temdn, timcur);

  /*--------------------------------------------------------------------*/
  /* call element routines to calculate & assemble tangent matrix */
  *action = calc_therm_tang_stat;
  container.dvec = NULL;
  container.dirich = dirich;
  container.global_numeq = numeq_total;
  container.kstep = 0;  /* WORK */
  container.isoltemdn = isol->temdn;
  i2ndsysmat = -1;  /* we do not have a second system matrix */
  calelm(actfield, actsolv, actpart, actintra, actsysarray,
         i2ndsysmat, &container, action);

  /*--------------------------------------------------------------------*/
  /* call rhs-routines to assemble rhs */
  container.kstep = 0;  /* WORK */
  container.inherit = 1;  /* ? */
  container.point_neum = 1;  /* ? */
  *action = calc_therm_heatload;  /* set action before call of calrhs */
  calrhs(actfield, actsolv, actpart, actintra, actsysarray,
         &(actsolv->rhs[actsysarray]), action, &container);
  rhsfact = -1.0;  /* multiply RHS by this factor */
  assemble_vec(actintra, &(actsolv->sysarray_typ[actsysarray]),
               &(actsolv->sysarray[actsysarray]),
               &(actsolv->rhs[actsysarray]), dirich, rhsfact);

  /*--------------------------------------------------------------------*/
  /* call solver */
  init = 0;
  solver_control(actfield, disnum, actsolv, actintra,
                 &(actsolv->sysarray_typ[actsysarray]),
                 &(actsolv->sysarray[actsysarray]),
                 &(actsolv->sol[actsysarray]),
                 &(actsolv->rhs[actsysarray]),
                 init);

  /*--------------------------------------------------------------------*/
  /* put the scaled prescribed temperatures to the nodes
   * in field sol (node_array_sol==0) at place isol->tem==0 together with
   * free temperatures */
  /* HINT: time curve is called indirectly */
  solserv_putdirich_to_dof(actfield, disnum,
                           node_array_sol, isol->temn, timcur);

  /*--------------------------------------------------------------------*/
  /* allreduce the result and put it to the node sol arrays */
  solserv_result_total(actfield,
                       disnum,
                       actintra,
                       &(actsolv->sol[actsysarray]),
                       isol->temn,
                       &(actsolv->sysarray[actsysarray]),
                       &(actsolv->sysarray_typ[actsysarray]));

  /*--------------------------------------------------------------------*/
  /* perform heat flux calculation */
  if (ioflags.therm_heatflux == 1)
  {
    *action = calc_therm_heatflux;
    container.dvec = NULL;
    container.dirich = NULL;
    container.global_numeq = 0;
    container.kstep = 0;  /* WORK */
    i2ndsysmat = -1;  /* we do not have a second system matrix */
    calelm(actfield, actsolv, actpart, actintra, actsysarray,
           i2ndsysmat, &container, action);
    /* reduce stresses, so they can be written */
    /* this only has an effect for SHELL8/9 */
    /* *action = calc_heatflux_stressreduce; */
    /*     container.kstep = 0;  /\* WORK *\/ */
    /*     calreduce(actfield, actpart, disnum, actintra, action, &container); */
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

  /*--------------------------------------------------------------------*/
  /* a sag-beim-abschied-leise-servus to the user */
  printf("---------------------------------------------------------------"
         "---------\n");
  printf("TSI quasi static thermal solution routine finished.\n");
  printf("==============================================================="
         "=========\n\n");

  end:

#ifdef BINIO
  if (ioflags.output_bin == 1)
  {
    destroy_bin_out_field(&out_context);
  }
#endif

#ifndef PARALLEL
  CCAFREE(actintra);
#endif

#ifdef DEBUG
  dstrc_exit();
#endif

  return;

}  /* end of tsi_th_stat() */


/*----------------------------------------------------------------------*/
#endif  /* end of #ifdef D_TSI */
#endif
