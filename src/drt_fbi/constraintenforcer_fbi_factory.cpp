/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create the appropriate constraint enforcement strategy for fluid-beam interaction


\level 1
*/
/*----------------------------------------------------------------------*/


#include "constraintenforcer_fbi_factory.H"
#include "constraintenforcer_fbi_penalty.H"
#include "ad_fbi_constraintbridge.H"
#include "ad_fbi_constraintbridge_penalty.H"
#include "immersed_geometry_coupler_fbi.H"

#include "../drt_inpar/inpar_fsi.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::FBIConstraintenforcer> ADAPTER::ConstraintEnforcerFactory::CreateEnforcer(
    const Teuchos::ParameterList& fsidyn)
{
  Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge =
      Teuchos::rcp(new ADAPTER::FBIConstraintBridgePenalty());
  Teuchos::RCP<FBI::FBIGeometryCoupler> coupler = Teuchos::rcp(new FBI::FBIGeometryCoupler());
  return Teuchos::rcp(new ADAPTER::FBIPenaltyConstraintenforcer(bridge, coupler));
}
