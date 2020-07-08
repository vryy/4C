/*----------------------------------------------------------------------------*/
/*! \file

 \brief Wrapper for the ALE time integration


 \level 2
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_wear.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_ale/ale_utils_mapextractor.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleWearWrapper::AleWearWrapper(Teuchos::RCP<Ale> ale) : AleWrapper(ale)
{
  // create the Wear interface
  interface_ = Teuchos::rcp(new ALE::UTILS::MapExtractor);
  interface_->Setup(*Discretization());
  SetupDBCMapEx(ALE::UTILS::MapExtractor::dbc_set_wear, interface_);

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALE::UTILS::MapExtractor> ADAPTER::AleWearWrapper::Interface() const
{
  return interface_;
}
