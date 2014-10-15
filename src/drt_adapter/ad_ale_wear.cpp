/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_wear.cpp

 \brief Wrapper for the ALE time integration

 <pre>
       Maintainer: Philip Farah
       farah@lnm.mw.tum.de
       089 - 289-15257
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_wear.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleWearWrapper::AleWearWrapper(Teuchos::RCP<Ale> ale)
  : AleWrapper(ale)
{
  dserror("ALE wear not implemented yet!");
  return;
}
