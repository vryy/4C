/*----------------------------------------------------------------------*/
/*!
\file
\brief Fehlberg4 prototypes to access its parameters

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
\date 03/07
*/

/*----------------------------------------------------------------------*/
/* define constants */
#define TSI_FEHLBG4_STG (5)  /* number of stages */

/*----------------------------------------------------------------------*/
/* Fehlberg4 parameters */
typedef struct _TSI_FEHLBG4
{
  /* stages */
  DOUBLE stg;
  /* A matrix */
  DOUBLE a[TSI_FEHLBG4_STG][TSI_FEHLBG4_STG];
  /* b vector */
  DOUBLE b[TSI_FEHLBG4_STG];
  /* c vector */
  DOUBLE c[TSI_FEHLBG4_STG];
  /* A.A matrix */
  DOUBLE aa[TSI_FEHLBG4_STG][TSI_FEHLBG4_STG];
  /* b.A vector */
  DOUBLE bb[TSI_FEHLBG4_STG];
} TSI_FEHLBG4;

/*----------------------------------------------------------------------*/
/* file tsi_fehlbg.c */
void tsi_fehlbg4_setpara(TSI_FEHLBG4* para);

INT tsi_fehlbg4_stages();
