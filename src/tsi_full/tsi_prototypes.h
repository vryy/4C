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


/*----------------------------------------------------------------------*/
/* The entries are ordered in the following way:
 * (1)  Filenames are in alphabetical order.
 * (2)  Functions are in the same order as in the file name.
 */


/*----------------------------------------------------------------------*/
/* file tsi_coupling.c */
void tsi_initcoupling(
  FIELD *structfield,
  INT disnum_s,
  FIELD *thermfield,
  INT disnum_t
);


/*----------------------------------------------------------------------*/
/* file tsi_dyn.c */
void tsi_dyn();


/*----------------------------------------------------------------------*/
/* file tsi_dyn_struct.c */
void tsi_dyn_struct();


/*----------------------------------------------------------------------*/
/* file tsi_stat_therm.c */
void tsi_stat_therm();
