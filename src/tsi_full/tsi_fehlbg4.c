/*----------------------------------------------------------------------*/
/*!
\file
\brief Fehlberg4 parameters

Features:
   - explicit
   - 4th order accurate
   - with embedded 5th order associate

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
/* header files */
#include "../headers/standardtypes.h"
#include "tsi_fehlbg4.h"

/*======================================================================*/
/*!
\brief get number of stages

\return number of stages

\author bborn
\date 03/07
*/
void tsi_fehlbg4_setpara(TSI_FEHLBG4* para)
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_fehlbg4_setpara");
#endif

  /*--------------------------------------------------------------------*/
  /* set parameters */

  /* number of stages */
  para->stg = TSI_FEHLBG4_STG;

  /* A matrix */
  para->a[0][0] = 0.;
  para->a[0][1] = 0.;
  para->a[0][2] = 0.;
  para->a[0][3] = 0.;
  para->a[0][4] = 0.;
  para->a[1][0] = 1./4.;
  para->a[1][1] = 0.;
  para->a[1][2] = 0.;
  para->a[1][3] = 0.;
  para->a[1][4] = 0.;
  para->a[2][0] = 3./32.;
  para->a[2][1] = 9./32.;
  para->a[2][2] = 0.;
  para->a[2][3] = 0.;
  para->a[2][4] = 0.;
  para->a[3][0] = 1932./2197.;
  para->a[3][1] = -7200./2197.;
  para->a[3][2] = 7296./2197.;
  para->a[3][3] = 0.;
  para->a[3][4] = 0.;
  para->a[4][0] = 439./216.;
  para->a[4][1] = -8.;
  para->a[4][2] = 3680./513.;
  para->a[4][3] = -845./4104.;
  para->a[4][4] = 0.;

  /* b vector */
  para->b[0] = 25./216.;
  para->b[1] = 0.;
  para->b[2] = 1408./2565.;
  para->b[3] = 2197./4104;
  para->b[4] = -1./5.;

  /* c vector */
  para->c[0] = 0.;
  para->c[1] = 1./4.;
  para->c[2] = 3./8.;
  para->c[3] = 12./13.;
  para->c[4] = 1.;

  /* A.A matrix */
  para->aa[0][0] = 0.;
  para->aa[0][1] = 0.;
  para->aa[0][2] = 0.;
  para->aa[0][3] = 0.;
  para->aa[0][4] = 0.;
  para->aa[1][0] = 0.;
  para->aa[1][1] = 0.;
  para->aa[1][2] = 0.;
  para->aa[1][3] = 0.;
  para->aa[1][4] = 0.;
  para->aa[2][0] = 9./128.;
  para->aa[2][1] = 0.;
  para->aa[2][2] = 0.;
  para->aa[2][3] = 0.;
  para->aa[2][4] = 0.;
  para->aa[3][0] = -1116./2197.;
  para->aa[3][1] = 2052./2197.;
  para->aa[3][2] = 0.;
  para->aa[3][3] = 0.;
  para->aa[3][4] = 0.;
  para->aa[4][0] = -353./234.;
  para->aa[4][1] = 35./13.;
  para->aa[4][2] = -80./117.;
  para->aa[4][3] = 0.;
  para->aa[4][4] = 0.;

  /* b.A vector */
  para->bb[0] = 25./216.;
  para->bb[1] = 0.;
  para->bb[2] = 176./513.;
  para->bb[3] = 169./4104.;
  para->bb[4] = 0.;

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}

/*======================================================================*/
/*!
\brief get number of stages

\return number of stages

\author bborn
\date 03/07
*/
INT tsi_fehlbg4_stages()
{

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_fehlbg4_stages");
#endif

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return TSI_FEHLBG4_STG;
}


