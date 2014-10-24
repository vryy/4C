/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_crack.cpp

 <pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            089 - 289-15257
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_crack.H"
#include "../drt_lib/drt_dserror.H"


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleCrackWrapper::AleCrackWrapper(Teuchos::RCP<Ale> ale)
  : AleWrapper(ale)
{
  // we just have an empty constructor since there is no interface defined for crack
  // we have BC only on the crack tip nodes which varies at each time step
  return;
}

