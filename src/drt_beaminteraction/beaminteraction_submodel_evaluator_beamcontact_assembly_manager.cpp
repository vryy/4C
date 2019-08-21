/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices.

\maintainer Ivo Steinbrecher

\level 3

*/
/*-----------------------------------------------------------*/


#include "beaminteraction_submodel_evaluator_beamcontact_assembly_manager.H"

#include "beam_contact_pair.H"


/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManager::BeamContactAssemblyManager(
    std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs)
    : assembly_contact_elepairs_(assembly_contact_elepairs)
{
}
