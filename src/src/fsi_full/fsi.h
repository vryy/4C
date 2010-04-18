#ifndef FSI_H
#define FSI_H

#include "../headers/standardtypes.h"
#include "../io/io.h"


typedef struct _FSI_FLUID_WORK
{
  INT             outstep;            /* counter for output control       */
  INT             restartstep;        /* counter for restart control      */

  DOUBLE          tes;
  DOUBLE          tss;

  ARRAY           frhs_a;

#ifdef PARALLEL
  ARRAY           fcouple_a;
  ARRAY           recvfcouple_a;
#endif

#ifndef PARALLEL
  INTRA           dummy_intra;
#endif

  ARRAY           totarea_a;

  FIELD* fluid_field;

#if defined(TRILINOS_PACKAGE) && defined(PM_TRILINOS)
  ARRAY           fgradprhs_a;
  TRILINOSMATRIX  grad;
  TRILINOSMATRIX  lmass;

  DIST_VECTOR   *press_rhs;
  DIST_VECTOR   *press_sol;

  INT numpdof;
#endif

#ifdef BINIO
  BIN_OUT_FIELD out_context;
  BIN_OUT_FIELD restart_context;
#endif
} FSI_FLUID_WORK;



typedef struct _FSI_STRUCT_WORK
{
  INT           outstep;        /* counter for output control        */
  INT           restartstep;    /* counter for restart control       */
  DOUBLE        t0;
  DOUBLE        dmax;           /* infinity norm of residual displacements    */

  INT           stiff_array;    /* indice of the active system sparse matrix  */
  INT           mass_array;     /* indice of the active system sparse matrix  */
  INT           damp_array;     /* indice of the active system sparse matrix  */

  DIST_VECTOR  *vel;            /* total velocities                           */
  DIST_VECTOR  *acc;            /* total accelerations                        */
  DIST_VECTOR  *fie;            /* internal forces and working array          */
  DIST_VECTOR  *dispi;          /* distributed vector to hold incr. disp.     */
  DIST_VECTOR  *work;           /* working vectors                            */

  ARRAY         intforce_a;     /* red. vector of full length for int forces  */
  ARRAY         fsiforce_a;     /* red. vector of full length for int forces  */
  ARRAY         dirich_a;       /* red. vector of full length for dirich rhs  */
  STRUCT_DYN_CALC dynvar;       /* variables to perform dynamic struct sim    */

  FIELD* struct_field;

#ifndef PARALLEL
  INTRA           dummy_intra;
#endif

#ifdef BINIO
  BIN_OUT_FIELD out_context;
  BIN_OUT_FIELD restart_context;
#endif
} FSI_STRUCT_WORK;



typedef struct _FSI_ALE_WORK
{
  INT      outstep;          /* counter for output to .out                         */
  INT      restartstep;      /* counter for output of restart data                 */
  ARRAY    dirich_a;	     /* red. vector of full length for dirich-part of rhs  */

#ifndef PARALLEL
  INTRA           dummy_intra;
#endif

  /* 2step */
  DOUBLE         min, max;     /*<! scaling parameters for ale two_step */
  DOUBLE         min_stiff;
  DOUBLE         max_stiff;

  /* LAS */
  ARRAY     index_a;

  FIELD* ale_field;

#ifdef BINIO
  BIN_OUT_FIELD out_context;
  BIN_OUT_FIELD restart_context;
#endif
} FSI_ALE_WORK;


#endif
