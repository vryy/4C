/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create the appropriate constraint enforcement strategy for fluid-beam interaction


\level 1
*/
/*----------------------------------------------------------------------*/


#include "baci_fbi_constraintenforcer_factory.hpp"

#include "baci_fbi_adapter_constraintbridge.hpp"
#include "baci_fbi_adapter_constraintbridge_penalty.hpp"
#include "baci_fbi_constraintenforcer_penalty.hpp"
#include "baci_fbi_immersed_geometry_coupler.hpp"
#include "baci_fbi_immersedcoupler_factory.hpp"
#include "baci_inpar_fbi.hpp"
#include "baci_inpar_fsi.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::FBIConstraintenforcer> ADAPTER::ConstraintEnforcerFactory::CreateEnforcer(
    const Teuchos::ParameterList& fsidyn, const Teuchos::ParameterList& fbidyn)
{
  Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge =
      Teuchos::rcp(new ADAPTER::FBIConstraintBridgePenalty());

  Teuchos::RCP<FBI::FBIGeometryCoupler> coupler =
      FBI::GeometryCouplerFactory::CreateGeometryCoupler(fbidyn);

  return Teuchos::rcp(new ADAPTER::FBIPenaltyConstraintenforcer(bridge, coupler));
}

FOUR_C_NAMESPACE_CLOSE
