/*!----------------------------------------------------------------------
\file
\brief xfem_prototypes.h

<pre>
Maintainer: Baris Irhan
            irhan@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/irhan/
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_XFEM

/*----------------------------------------------------------------------*
 | structures for level set                               irhan 04/04   |
 *----------------------------------------------------------------------*/
#include "../ls/ls.h"
/*----------------------------------------------------------------------*
 | structures for xfem                                    irhan 04/04   |
 *----------------------------------------------------------------------*/
#include "xfem.h"

/* RULE HOW TO ADD NEW FILES AND FUNCTIONS:
   1.) THE FILENAMES ARE IN ALPHABETICAL ORDER !!!
   2.) FUNCTIONS ARE IN THE SAME ORDER LIKE IN THE FILE !!!
*/
/*!
\addtogroup XFEM
*//*! @{ (documentation module open)*/

/* xfem_calfuncderiv.c */
void xfem_f2_funct(
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE,
  DOUBLE,
  DIS_TYP,
  DOUBLE*,
  INT,
  INT
  );


void xfem_f2_funct1(
  DOUBLE*,
  DOUBLE,
  DOUBLE,
  DIS_TYP,
  INT,
  DOUBLE*,
  INT
  );


void xfem_f2_derxy(
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE,
  INT,
  DOUBLE*,
  DOUBLE*,
  INT
  );


void xfem_f2_derxy2(
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  INT,
  DOUBLE*,
  DOUBLE*,
  INT
  );


void xfem_f2_veli(
  DOUBLE*,
  DOUBLE*,
  DOUBLE**,
  INT
  );


void xfem_f2_vder(
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  INT
  );


void xfem_f2_vder2(
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  INT
  );


/* xfem_calgalmat.c */
void xfem_f2_calkvv(
  ELEMENT*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE**,
  DOUBLE,
  DOUBLE,
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calkvp(
  DOUBLE**,
  DOUBLE*,
  DOUBLE**,
  DOUBLE,
  INT,
  INT*
  );


void xfem_f2_calmvv(
  DOUBLE**,
  DOUBLE*,
  DOUBLE,
  INT,
  INT*,
  DOUBLE
  );


/* xfem_calstabmat.c */
void xfem_f2_calstabkvv(
  ELEMENT*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE,
  DOUBLE,
  INT,
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabkvp(
  ELEMENT*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE,
  DOUBLE,
  INT,
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabmvv(
  ELEMENT*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE,
  DOUBLE,
  INT,
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabkpv(
  ELEMENT*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE,
  DOUBLE,
  INT,
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabmpv(
  DOUBLE**,
  DOUBLE*,
  DOUBLE**,
  DOUBLE,
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabkpp(
  DOUBLE**,
  DOUBLE**,
  DOUBLE,
  INT
  );


/* xfem_f2_calele.c */
void xfem_f2_calele(
  FLUID_DATA*,
  ELEMENT*,
  ARRAY*,
  ARRAY*,
  ARRAY*,
  ARRAY*,
  ARRAY*,
  INT*,
  INT*,
  INT,
  INT,
  INT
  );


void xfem_f2_loc_con(void);


void xfem_f2_loc_ass_tangent(void);


void xfem_f2_loc_ass_intforce(
  DOUBLE*,
  DOUBLE*
  );


void xfem_f2_array_init(void);


void xfem_f2_init(
  ELEMENT*
  );


void xfem_f2_iand(void);


/* xfem_f2_calelesize.c */
void xfem_f2_calelesize(
  ELEMENT*,
  FLUID_DATA*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE*,
  DOUBLE**
  );


/* xfem_f2_calextrhs.c */
void xfem_f2_calgalexfv(
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE,
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabexfv(
  ELEMENT*,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE*,
  DOUBLE*,
  DOUBLE,
  DOUBLE,
  INT,
  INT,
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabexfp(
  DOUBLE*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE,
  INT,
  INT,
  DOUBLE
  );


/* xfem_f2_calint.c */
void xfem_f2_calint(
  FLUID_DATA*,
  ELEMENT*,
  INT*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE*,
  DOUBLE*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE**,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**
  );


/* xfem_f2_calset.c */
void xfem_f2_calset(
  ELEMENT*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  INT*
  );


/* xfem_f2_caltimerhs.c */
void xfem_f2_calgaltfv(
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE,
  DOUBLE,
  DOUBLE,
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabtfv(
  ELEMENT*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE**,
  DOUBLE*,
  DOUBLE,
  DOUBLE,
  INT,
  INT,
  INT*,
  DOUBLE
  );


void xfem_f2_calstabtfp(
  DOUBLE*,
  DOUBLE**,
  DOUBLE**,
  DOUBLE*,
  DOUBLE*,
  DOUBLE*,
  DOUBLE,
  DOUBLE,
  INT,
  INT,
  DOUBLE
  );


/* xfem_f2_intg.c */
void xfem_f2_intg(
  FLUID_DATA*
  );


/* xfem_f2_main.c */
void xfem_fluid2(
  PARTITION*,
  INTRA*,
  ELEMENT*,
  ARRAY*,
  ARRAY*,
  ARRAY*,
  ARRAY*,
  ARRAY*,
  CALC_ACTION*,
  INT*,
  INT*,
  CONTAINER*
  );


/* xfem_polygon.c */
void xfem_polygon(
  XFEMPOLYFLAG,
  ELEMENT*
  );


void xfem_polygon_init(
  ELEMENT*
  );


void xfem_polygon_cons(
  ELEMENT*
  );


void xfem_polygon_GP(
  ELEMENT*
  );


void xfem_polygon_getsubp(
  INT,
  INT
  );


void xfem_polygon_compGP(
  DIS_TYP
  );


void xfem_polygon_funct(void);


void xfem_polygon_deriv(void);


void xfem_polygon_resNewton(void);


void xfem_polygon_tanNewton(void);


void xfem_polygon_write(
  ELEMENT*
  );


void xfem_polygon_open(void);


void xfem_polygon_close(void);


void xfem_polygon_area_rect(void);


DOUBLE xfem_polygon_area_subtri(
  INT
  );


DOUBLE xfem_polygon_area_tri(void);


void xfem_polygon_target_tri(void);


void xfem_polygon_target_subtri(
  INT
  );


void xfem_polygon_mat(
  ELEMENT*
  );
/*! @} (documentation module close)*/
#endif
