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
/* file tsi_fehlbg.c */
DOUBLE tsi_fehlbg4_a(const INT i,
                     const INT j);

DOUBLE tsi_fehlbg4_b(const INT i);

DOUBLE tsi_fehlbg4_c(const INT i);

DOUBLE tsi_fehlbg4_aa(const INT i,
                      const INT j);

DOUBLE tsi_fehlbg4_bb(const INT i);
