/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_fpsi.cpp

 <pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de
            089 - 289-15240
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_fpsi.H"
#include "../drt_lib/drt_dserror.H"

#include "../drt_ale_new/ale_utils_mapextractor.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleFpsiWrapper::AleFpsiWrapper(Teuchos::RCP<Ale> ale)
  : AleFsiWrapper(ale)
{
  dserror("FPSI adapter not implemented yet!");

  return;
}
