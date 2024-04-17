/*----------------------------------------------------------------------------*/
/*! \file

 \brief Wrapper for the ALE time integration


 \level 2
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "baci_adapter_ale_wear.hpp"

#include "baci_ale_utils_mapextractor.hpp"
#include "baci_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

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

FOUR_C_NAMESPACE_CLOSE
