/*----------------------------------------------------------------------*/
/*!
\file
\brief Fehlberg4 parameters

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
/* header files */
#include "../headers/standardtypes.h"
#include "tsi_fehlbg4.h"


/*======================================================================*/
/*!
\brief get compononent of A array

\author bborn
\date 03/07
*/
DOUBLE tsi_fehlbg4_a(const INT i,  /*!< 1st index */
                     const INT j)  /*!< 2nd index */
{
  DOUBLE fb4a[5][5] = { { 0., 0., 0., 0., 0. },
                        { 1./4., 0., 0., 0., 0. },
                        { 3./32., 9./32., 0., 0., 0. },
                        { 1932./2197., -7200./2197., 7296./2197., 0., 0. },
                        { 439./216., -8., 3680./513., -845./4104., 0. } };

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_fehlbg4_a");
#endif

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return fb4a[i][j];
}

/*======================================================================*/
/*!
\brief get compononent of b vector

\author bborn
\date 03/07
*/
DOUBLE tsi_fehlbg4_b(const INT i)
{
  DOUBLE fb4b[5] = { 25./216., 0., 1408./2565., 2197./4104, -1./5. };

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_fehlbg4_b");
#endif

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return fb4b[i];
}

/*======================================================================*/
/*!
\brief get compononent of c vector

\author bborn
\date 03/07
*/
DOUBLE tsi_fehlbg4_c(const INT i)
{
  DOUBLE fb4c[5] = { 0., 1./4., 3./8., 12./13., 1. };

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_fehlbg4_c");
#endif

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return fb4c[i];
}

/*======================================================================*/
/*!
\brief get compononent of A.A array

\author bborn
\date 03/07
*/
DOUBLE tsi_fehlbg4_aa(const INT i,
                      const INT j)
{
  DOUBLE fb4aa[5][5] = { { 0., 0., 0., 0., 0. },
                         { 0., 0., 0., 0., 0. },
                         { 9./128., 0., 0., 0., 0. },
                         { -1116./2197., 2052./2197., 0., 0., 0. },
                         { -353./234., 35./13., -80./117., 0., 0. } };

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_fehlbg4_aa");
#endif

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return fb4aa[i][j];
}

/*======================================================================*/
/*!
\brief get compononent of b vector

\author bborn
\date 03/07
*/
DOUBLE tsi_fehlbg4_bb(const INT i)
{
  DOUBLE fb4bb[5] = { 25./216., 0., 176./513., 169./4104., 0. };

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_enter("tsi_fehlbg4_bb");
#endif

  /*--------------------------------------------------------------------*/
#ifdef DEBUG
  dstrc_exit();
#endif
  return fb4bb[i];
}
