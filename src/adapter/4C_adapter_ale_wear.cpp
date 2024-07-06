/*----------------------------------------------------------------------------*/
/*! \file

 \brief Wrapper for the ALE time integration


 \level 2
 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_adapter_ale_wear.hpp"

#include "4C_ale_utils_mapextractor.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Adapter::AleWearWrapper::AleWearWrapper(Teuchos::RCP<Ale> ale) : AleWrapper(ale)
{
  // create the Wear interface
  interface_ = Teuchos::rcp(new ALE::UTILS::MapExtractor);
  interface_->setup(*discretization());
  setup_dbc_map_ex(ALE::UTILS::MapExtractor::dbc_set_wear, interface_);

  return;
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALE::UTILS::MapExtractor> Adapter::AleWearWrapper::interface() const
{
  return interface_;
}

FOUR_C_NAMESPACE_CLOSE
