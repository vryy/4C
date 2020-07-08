/*-----------------------------------------------------------*/
/*! \file

\brief Class to assemble pair based contributions into global matrices.


\level 1

*/
/*-----------------------------------------------------------*/


#include "partitioned_penaltycoupling_assembly_manager.H"

#include "../drt_beaminteraction/beam_contact_pair.H"


/**
 *
 */
BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManager::
    PartitionedBeamInteractionAssemblyManager(
        std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> assembly_contact_elepairs)
    : assembly_contact_elepairs_(assembly_contact_elepairs)
{
}
