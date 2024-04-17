/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create the appropriate geometry coupler strategy for fluid-beam interaction


\level 1
*/
/*----------------------------------------------------------------------*/


#include "baci_fbi_immersedcoupler_factory.hpp"

#include "baci_fbi_immersed_geometry_coupler.hpp"
#include "baci_fbi_immersed_geometry_coupler_binning.hpp"
#include "baci_inpar_fbi.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<FBI::FBIGeometryCoupler> FBI::GeometryCouplerFactory::CreateGeometryCoupler(
    const Teuchos::ParameterList& fbidyn)
{
  INPAR::FBI::BeamToFluidPreSortStrategy presort_strategy =
      Teuchos::getIntegralValue<INPAR::FBI::BeamToFluidPreSortStrategy>(fbidyn, "PRESORT_STRATEGY");

  Teuchos::RCP<FBI::FBIGeometryCoupler> coupler;

  if (presort_strategy == INPAR::FBI::BeamToFluidPreSortStrategy::bruteforce)
  {
    coupler = Teuchos::rcp(new FBI::FBIGeometryCoupler());
  }
  else if (presort_strategy == INPAR::FBI::BeamToFluidPreSortStrategy::binning)
  {
    coupler = Teuchos::rcp(new FBI::FBIBinningGeometryCoupler());
  }
  else
    dserror("Unknown Beam to Fluid PreSort Strategy");

  return coupler;
}

FOUR_C_NAMESPACE_CLOSE
