/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create the appropriate constraint enforcement strategy for fluid-beam interaction


\level 1
*/
/*----------------------------------------------------------------------*/


#include "fbi_constraintenforcer_factory.H"
#include "fbi_constraintenforcer_penalty.H"
#include "fbi_adapter_constraintbridge.H"
#include "fbi_adapter_constraintbridge_penalty.H"
#include "fbi_immersed_geometry_coupler.H"
#include "fbi_immersedcoupler_factory.H"

#include "inpar_fbi.H"
#include "inpar_fsi.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
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
