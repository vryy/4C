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
/* file tsi_dyn_struct.c */
void tsi_dyn_struct(INT disnum_s,
                    INT disnum_t);


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
/* file tsi_stat_therm.c */
void tsi_stat_therm(INT disnum_s,
                    INT disnum_t);
