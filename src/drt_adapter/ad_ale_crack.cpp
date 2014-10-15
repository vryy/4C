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
  dserror("ALE Crack not implemented yet.");

  return;
}

