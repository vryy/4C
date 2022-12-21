/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create appropriate beam to fluid meshtying assembly managers for the desired
constraint discretization approach

\level 2
*/

#include "beam_to_fluid_assembly_manager_factory.H"
#include "fluid_assembly_strategy.H"
#include "partitioned_penaltycoupling_assembly_manager_direct.H"
#include "partitioned_penaltycoupling_assembly_manager_indirect.H"
#include "beam_to_fluid_meshtying_params.H"
#include "beam_contact_pair.H"
#include "discret.H"
#include "inpar_fbi.H"
/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManager>
BEAMINTERACTION::BeamToFluidAssemblyManagerFactory::CreateAssemblyManager(
    Teuchos::RCP<const DRT::Discretization> discretization1,
    Teuchos::RCP<const DRT::Discretization> discretization2,
    std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> interaction_pairs,
    const Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> params_ptr,
    Teuchos::RCP<FBI::UTILS::FBIAssemblyStrategy> assemblystrategy)
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
              interaction_pairs, assemblystrategy));
      break;
    default:
      dserror("Beam To Fluid Meshtying Discretization Type not supported!");
      return Teuchos::null;
  }
}
