/*!----------------------------------------------------------------------
\file
\brief ls_prototypes.h

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_LS

/*----------------------------------------------------------------------*
 | structures for level set                               irhan 04/04   |
 *----------------------------------------------------------------------*/
#include "ls.h"
/*----------------------------------------------------------------------*
 | structures for xfem                                    irhan 04/04   |
 *----------------------------------------------------------------------*/
#include "../xfem/xfem.h"

/* RULE HOW TO ADD NEW FILES AND FUNCTIONS:
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE !!!
*/
/*!
\addtogroup LEVELSET
*//*! @{ (documentation module open)*/


/* ls2_calfuncderiv.c */
void ls2_funct(
  DOUBLE*,
  DOUBLE**,
  DOUBLE,
  DOUBLE,
  DIS_TYP
  );


void ls2_jaco(
  DOUBLE**,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE*,
  INT,
  ELEMENT*
  );


void ls2_gder(
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE,
  INT
  );


/* ls2_calset.c */
void ls2_calset(
  ELEMENT*,
  DOUBLE **,
  DOUBLE*,
  DOUBLE*
  );


void ls2_calset1(
  ELEMENT*,
  INT,
  DOUBLE*
  );


/* ls2_inpele.c */
void ls2inp(
  ELEMENT*
  );


void fluid_to_ls(
  const FIELD*,
  const FIELD*
  );


void ls_find_compatible_ele(
  const ELEMENT*,
  const ELEMENT*,
  INT*
  );


/* ls2_intg.c */
void ls2_intg(
  LS2_INTG_DATA*
  );


/* ls2_main.c */
void ls2(
  PARTITION*,
  INTRA*,
  ELEMENT*,
  ARRAY*,
  ARRAY*,
  ARRAY*,
  ARRAY*,
  CALC_ACTION*,
  CONTAINER*
  );


void ls2_init(void);


void ls2_calc(
  ELEMENT*
  );


void ls2_calc_advection(
  ELEMENT*
  );


void ls2_calc_reinitialization(
  ELEMENT*
  );


void ls2_calc_gradient(
  ELEMENT*
  );


void ls2_calc_curvature(
  ELEMENT*
  );


void ls2_calc_localized(
  ELEMENT*
  );


void ls2_setfluidvel(
  ELEMENT*
  );


void ls2_setfluidvel_byuser(
  INT,
  INT
  );


void ls2_setfluidvel_byfluid(
  ELEMENT*
  );


/* ls_algout.c */
void ls_algout();


/* ls_convcheck.c */
INT ls_convcheck(
  LS_DYNAMIC*,
  DOUBLE,
  INT,
  DOUBLE,
  DOUBLE
  );


/* ls_dirich.c */
void ls_initdirich(
  FIELD*
  );


void ls_setdirich(
  FIELD*,
  INT
  );


/* ls_dyn.c */
void ls_dyn(void);


/* ls_fluid.c */
void ls_fluid(
  INT
  );


void ls_fluid_init_data();


void ls_fluid_init_fluid();


void ls_fluid_init_solver();


void ls_fluid_reinit_solver();


void ls_fluid_solv(void);


void ls_fluid_fina(void);


void ls_fluid_clea(void);


/* ls_init.c */
void ls_init(
  FIELD*,
  LS_DYNAMIC*,
  INT
  );


void ls_init_single_circle(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_double_circle(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_triple_circle(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_single_circle_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_double_circle_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_triple_circle_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_single_line(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_double_line(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_triple_line(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_single_line_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_double_line_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


void ls_init_triple_line_sharp(
  FIELD*,
  LS_DYNAMIC*
  );


/* ls_levelset.c */
void ls_levelset(
  INT
  );


void ls_levelset_init_data();


void ls_levelset_init_levelset();


void ls_levelset_init_solver();


void ls_levelset_solv(void);


void ls_levelset_fina(void);


void ls_levelset_clea(void);


void ls_levelset_write_restart();


void ls_levelset_read_restart();


void ls_levelset_reset_sol_increment();


/* ls_solserv_sol_copy_reinit.c */
void solserv_sol_copy_reinit(
  FIELD*,
  INT,
  INT,
  INT,
  INT,
  INT
  );


/* ls_restart.c */
void ls2_write_restart(
  ELEMENT*,
  INT,
  long int*
  );


void ls2_read_restart(
  ELEMENT *actele,
  INT nhandle,
  long int *handles
  );


/* ls_update.c */
void ls_main(
  LSFLAG
  );


void ls_initialize(void);


void ls_updt(void);


void ls_construct(
  ELEMENT*
  );


void ls_updatelement(
  ELEMENT*,
  DOUBLE*
  );


void is_tricut(
  DOUBLE*,
  INT*,
  INT*,
  INT*,
  INT*
  );


void ls_compintersection(
  DOUBLE*,
  DOUBLE*,
  INT,
  INT
  );


void ls_comppoint(
  DOUBLE*,
  DOUBLE*,
  INT,
  INT,
  INT
  );


void ls_updateneighbor(
  ELEMENT*,
  INT,
  INT
  );


void ls_makepatch(
  ELEMENT*,
  ELEMENT *[],
  INT*
  );


void ls_activate(void);


void ls_write(void);


void ls_write_fluid_soln();


void ls_reset(void);


void ls_resetintdata(
  LS_INT_DATA*
  );


void ls_resetpolydata(
  LS_POLY_DATA*
  );


void ls_finalize(void);


void ls_printelinfo(
  ELEMENT*
  );


void ls_printelinfo_to_file(
  ELEMENT*
  );


void ls_printnodeinfo(
  NODE*
  );


void ls_printnodeinfo_to_file(
  NODE*
  );


void ls_setdata(void);


void ls_printelements();


void ls_printelement(
  ELEMENT*
  );


void ls_localize(void);


void ls_to_matlab();


void ls_polygon_con(
  INT*,
  INT,
  INT
  );


void ls_polygonize(void);


void ls_init_material(void);


void ls_updt_material(void);


void ls_check_profile(void);


void ls_init_pres_bd();


ELEMENT* ls_find_last_in_list();


void ls_mark_tip_elements();


void ls_reset_fld_sol_increment();


void ls_result_incre(
  FIELD             *actfield,
  INTRA             *actintra,
  DIST_VECTOR       *sol,
  INT                place,
  SPARSE_ARRAY      *sysarray,
  SPARSE_TYP        *sysarray_typ,
  DOUBLE            *lrat,
  LS_DYNAMIC       *lsdyn
  );


/* twophase_inpctr_data.c */
void twophase_inpctr_data();
/*! @} (documentation module close)*/
#endif
