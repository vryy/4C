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
  CALC_ACTION*
  );


void ls2_init();


void ls2_calc(
  ELEMENT*
  );


void ls2_calc_nonlocalized(
  ELEMENT*
  );


void ls2_calc_reinitialized(
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
void ls_algout(
  LS_DYNAMIC*
  );


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
  FIELD*,
  LS_DYNAMIC*
  );


void ls_setdirich(
  FIELD*,
  LS_DYNAMIC*,
  INT
  );


/* ls_dyn.c */
void ls_dyn();


/* ls_fluid.c */
void ls_fluid(
  INT
  );


void ls_fluid_init();


void ls_fluid_solv();


void ls_fluid_fina();


void ls_fluid_clea();


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


void ls_levelset_init();


void ls_levelset_solv();


void ls_levelset_fina();


void ls_levelset_clea();


/* ls_update.c */
void ls_main(
  LSFLAG
  );


void ls_update(
  FRONTLSFLAG
  );


void ls_initialize();


void ls_updt();


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


void ls_activate();


void ls_write();


void ls_write_soln();


void ls_reset();


void ls_resetintdata(
  LS_INT_DATA*
  );


void ls_resetpolydata(
  LS_POLY_DATA*
  );


void ls_finalize();


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


void ls_setdata();


void ls_printelements();				 


void ls_printelement(
  ELEMENT*
  );


void ls_localize();


void ls_to_matlab();									  


void ls_polygon_con(
  INT*,
  INT,
  INT
  );


void ls_polygonize();


void ls_init_material();


void ls_updt_material();


void ls_check_profile();

void ls_init_pres_bd();
/*! @} (documentation module close)*/	    
#endif
