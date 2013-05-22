/*----------------------------------------------------------------------*/
/*!
\file scatra_utils.H

\brief utility functions for scalar transport problems

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_UTILS_H
#define SCATRA_UTILS_H

#include "../drt_inpar/inpar_scatra.H"


namespace SCATRA
{

  /// does the given scatratype belong to the ELCH group of problems?
  bool IsElch(const enum INPAR::SCATRA::ScaTraType scatratype);

} // namespace SCATRA

#endif
