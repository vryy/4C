/*----------------------------------------------------------------------------*/
/*! \file

 \brief FSI Wrapper for the ALE time integration

 \level 1


 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_adapter_ale_fsi.hpp"

#include "4C_ale_utils_mapextractor.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Adapter::AleFsiWrapper::AleFsiWrapper(Teuchos::RCP<Ale> ale) : AleWrapper(ale)
{
  // create the FSI interface
  interface_ = Teuchos::rcp(new ALE::UTILS::MapExtractor);
  interface_->Setup(*discretization());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALE::UTILS::MapExtractor> Adapter::AleFsiWrapper::Interface() const
{
  return interface_;
}

FOUR_C_NAMESPACE_CLOSE
