/*----------------------------------------------------------------------------*/
/*!
 \file ad_ale_fsi.cpp

 \brief FSI Wrapper for the ALE time integration

 \level 1

 \maintainer Matthias Mayr

 */
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_fsi.H"

#include "../drt_ale/ale_utils_mapextractor.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
ADAPTER::AleFsiWrapper::AleFsiWrapper(Teuchos::RCP<Ale> ale) : AleWrapper(ale)
{
  // create the FSI interface
  interface_ = Teuchos::rcp(new ALE::UTILS::MapExtractor);
  interface_->Setup(*Discretization());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const ALE::UTILS::MapExtractor> ADAPTER::AleFsiWrapper::Interface() const
{
  return interface_;
}
