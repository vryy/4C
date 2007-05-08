/*----------------------------------------------------------------------*/
/*!
\file
\brief Prototypes (declarations) of thermal structure interaction
functions

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089-289-15237
</pre>
*/


/*----------------------------------------------------------------------*/
/* header files */
#ifdef BINIO
#include "../io/io.h"
#endif

/*======================================================================*/
/* The entries are ordered in the following way:
 * (1)  Filenames are in alphabetical order.
 * (2)  Functions are in the same order as in the file name.
 */


/*----------------------------------------------------------------------*/
/* file tsi_coupling.c */
void tsi_coupling(FIELD *structfield,
                  INT disnum_s,
                  FIELD *thermfield,
                  INT disnum_t);


/*----------------------------------------------------------------------*/
/* file tsi_dyn.c */
void tsi_dyn();


/*----------------------------------------------------------------------*/
/* file tsi_init.c */
void tsi_init_chkfld(INT *numfld, /* total number of fields */
                     INT *numsf,  /* number (index) of structure field */
                     INT *numtf,  /* number (index) of thermal field */
                     FIELD **structfield,  /* structure field */
                     FIELD **thermfield);  /* thermal field */
void tsi_init_chkdis(FIELD *actfield, 
                     INT *disnum);
void tsi_init_alldyn(INT numsf,
                     INT numtf,
                     STRUCT_DYNAMIC **sdyn,
                     THERM_DYNAMIC **tdyn,
                     TSI_DYNAMIC **tsidyn);
void tsi_init_nodsol(FIELD *structfield,
                     INT disnum_s,
                     FIELD *thermfield, 
                     INT disnum_t);
void tsi_init_curve();

/*----------------------------------------------------------------------*/
/* file tsi_st_cendif.c */
void tsi_st_cendif(INT disnum_s,
                   INT disnum_t);

/*----------------------------------------------------------------------*/
/* file tsi_st_fehlbg.c */
void tsi_st_fehlbg(INT disnum_s,
                   INT disnum_t);

/*----------------------------------------------------------------------*/
/* file tsi_st_genalp.c */
void tsi_st_genalp(INT disnum_s,
                   INT disnum_t);

/*----------------------------------------------------------------------*/
/* file tsi_st_genalp_sub.c */
void tsi_st_genalp_init(PARTITION* actpart,
                        INTRA* actintra,
                        FIELD* actfield,
                        INT disnum,
                        ARRAY_POSITION* ipos,
                        ARRAY_POSITION_SOL* isol,
                        ARRAY_POSITION_SOLINC* isolinc,
                        SOLVAR* actsolv,
                        INT* numeq_total,
                        INT* numeq,
                        INT* stiff_array,
                        INT* mass_array,
                        INT* damp_array,
                        STRUCT_DYNAMIC* actdyn,
                        STRUCT_DYN_CALC* dynvar,
                        CONTAINER* container,
#ifdef BINIO
                        BIN_OUT_FIELD* out_context,
#endif
                        INT vel_num,
                        DIST_VECTOR** vel,
                        INT acc_num,
                        DIST_VECTOR** acc,
                        INT fie_num,
                        DIST_VECTOR** fie,
                        INT dispi_num,
                        DIST_VECTOR** dispi,
                        INT work_num,
                        DIST_VECTOR** work,
                        ARRAY* intforce_a,
                        ARRAY* dirich_a);
void tsi_st_genalp_pred(PARTITION* actpart,
                        INTRA* actintra,
                        FIELD* actfield,
                        INT disnum,
                        ARRAY_POSITION_SOL* isol,
                        ARRAY_POSITION_SOLINC* isolinc,
                        ARRAY_POSITION_SOLRES* isolres,
                        SOLVAR* actsolv,
                        INT numeq_total,
                        INT stiff_array,
                        INT mass_array,
                        INT damp_array,
                        STRUCT_DYNAMIC* actdyn,
                        STRUCT_DYN_CALC* dynvar,
                        CONTAINER* container,
                        DIST_VECTOR* vel,
                        DIST_VECTOR* acc,
                        DIST_VECTOR* fie,
                        DIST_VECTOR* dispi,
                        DIST_VECTOR* work,
                        ARRAY* dirich_a,
                        ARRAY* intforce_a);
void tsi_st_genalp_equi(PARTITION* actpart,
                        INTRA* actintra,
                        FIELD* actfield,
                        INT disnum,
                        ARRAY_POSITION_SOL* isol,
                        ARRAY_POSITION_SOLINC* isolinc,
                        ARRAY_POSITION_SOLRES* isolres,
                        SOLVAR* actsolv,
                        INT numeq_total,
                        INT stiff_array,
                        INT mass_array,
                        INT damp_array,
                        STRUCT_DYNAMIC* actdyn,
                        STRUCT_DYN_CALC* dynvar,
                        CONTAINER* container,
                        DIST_VECTOR* vel,
                        DIST_VECTOR* acc,
                        DIST_VECTOR* fie,
                        DIST_VECTOR* dispi,
                        DIST_VECTOR* work,
                        ARRAY* intforce_a,
                        ARRAY* dirich_a);
void tsi_st_genalp_chkcnv(INTRA* actintra,
                          STRUCT_DYNAMIC* actdyn,
                          STRUCT_DYN_CALC* dynvar,
                          DIST_VECTOR* work,
                          DIST_VECTOR* dispi,
                          INT mod_stdout,
                          INT* converged);
void tsi_st_genalp_updincr(PARTITION* actpart,
                           INTRA* actintra,
                           FIELD* actfield,
                           INT disnum,
                           ARRAY_POSITION_SOL* isol,
                           ARRAY_POSITION_SOLINC* isolinc,
                           SOLVAR* actsolv,
                           INT mass_array,
                           INT stiff_array,
                           STRUCT_DYNAMIC* actdyn,
                           STRUCT_DYN_CALC* dynvar,
                           CONTAINER* container,
                           DIST_VECTOR* vel,
                           DIST_VECTOR* acc,
                           DIST_VECTOR* dispi,
                           DIST_VECTOR* work);
void tsi_st_genalp_out(PARTITION* actpart,
                       INTRA* actintra,
                       FIELD* actfield,
                       INT disnum_s,
                       ARRAY_POSITION_SOL* isol,
                       SOLVAR* actsolv,
                       INT stiff_array,
                       STRUCT_DYNAMIC* actdyn,
                       STRUCT_DYN_CALC* dynvar,
                       CONTAINER* container,
#ifdef BINIO
                       BIN_OUT_FIELD* out_context,
#endif
                       INT dispi_num,
                       DIST_VECTOR* dispi,
                       INT vel_num,
                       DIST_VECTOR* vel,
                       INT acc_num,
                       DIST_VECTOR* acc,
                       INT fie_num,
                       DIST_VECTOR* fie,
                       INT work_num,
                       DIST_VECTOR* work,
                       ARRAY* intforce_a,
                       ARRAY* dirich_a);
void tsi_st_genalp_final(SOLVAR* actsolv,
			 const INT dispi_num,
			 DIST_VECTOR** dispi,
			 const INT vel_num,
			 DIST_VECTOR** vel,
			 const INT acc_num,
			 DIST_VECTOR** acc,
			 const INT fie_num,
			 DIST_VECTOR** fie,
			 const INT work_num,
			 DIST_VECTOR** work,
			 ARRAY* intforce_a,
			 ARRAY* dirich_a);
void tsi_st_genalp_sub(INT disnum_s,
                       INT disnum_t);

/*----------------------------------------------------------------------*/
/* file tsi_th_stat.c */
void tsi_th_stat(INT disnum_s,
                 INT disnum_t);
