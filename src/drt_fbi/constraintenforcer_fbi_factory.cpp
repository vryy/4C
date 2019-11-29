/*----------------------------------------------------------------------*/
/*! \file
\file constraintenforcer_fbi_factory.cpp

\brief Factory to create the appropriate constraint enforcement strategy for fluid-beam interaction

\maintainer Nora Hagmeyer

\level 1
*/
/*----------------------------------------------------------------------*/


#include "constraintenforcer_fbi_factory.H"
#include "constraintenforcer_fbi_penalty.H"
#include "ad_fbi_constraintbridge.H"
#include "ad_fbi_constraintbridge_penalty.H"

#include "../drt_inpar/inpar_fsi.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<ADAPTER::FBIConstraintenforcer> ADAPTER::ConstraintEnforcerFactory::CreateEnforcer(
    const Teuchos::ParameterList& fsidyn)
{
  Teuchos::RCP<ADAPTER::FBIConstraintBridge> bridge =
      Teuchos::rcp(new ADAPTER::FBIConstraintBridgePenalty());
  return Teuchos::rcp(new ADAPTER::FBIPenaltyConstraintenforcer(bridge));
}
