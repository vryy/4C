/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_biofilm_fsi.cpp

 \brief Wrapper for the ALE time integration

 <pre>
       Maintainer: Christoph Ager
       ager@lnm.mw.tum.de
       089 - 289-15249
 </pre>
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_biofilm_fsi.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleBiofilmFsiWrapper::AleBiofilmFsiWrapper(Teuchos::RCP<Ale> ale)
  : AleFsiWrapper(ale)
{
  //Biofilm FSI not implemented yet!");
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleBiofilmFsiWrapper::SolveBioGr()
{
  //Todo (ager): do we need this??
}
