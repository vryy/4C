/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create appropriate beam to fluid meshtying assembly managers for the desired
constraint discretization approach

\level 2
\maintainer Nora Hagmeyer
*/

#include "beam_to_fluid_assembly_manager_factory.H"
#include "partitioned_penaltycoupling_assembly_manager_direct.H"
#include "partitioned_penaltycoupling_assembly_manager_indirect.H"
#include "beam_to_fluid_meshtying_params.H"
#include "../drt_beaminteraction/beam_contact_pair.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_fbi.H"
/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManager>
BEAMINTERACTION::BeamToFluidAssemblyManagerFactory::CreateAssemblyManager(
    Teuchos::RCP<const DRT::Discretization> discretization1,
    Teuchos::RCP<const DRT::Discretization> discretization2,
    std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> interaction_pairs,
    const Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> params_ptr)
{
  // Get the meshtying discretization method.
  INPAR::FBI::BeamToFluidDiscretization meshtying_discretization =
      params_ptr->GetContactDiscretization();

  switch (meshtying_discretization)
  {
    case INPAR::FBI::BeamToFluidDiscretization::mortar:
      return Teuchos::rcp(
          new BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerIndirect(
              interaction_pairs, discretization1, discretization2, params_ptr));
      break;
    case INPAR::FBI::BeamToFluidDiscretization::gauss_point_to_segment:
      return Teuchos::rcp(
          new BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerDirect(
              interaction_pairs));
      break;
    default:
      dserror("Beam To Fluid Meshtying Discretization Type not supported!");
      return Teuchos::null;
  }
}
