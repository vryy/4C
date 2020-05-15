/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create appropriate beam to fluid meshtying assembly managers for the desired
constraint discretization approach

\maintainer Nora Hagmeyer

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef BEAM_TO_FLUID_ASSEMBLY_MANAGER_FACTORY_H
#define BEAM_TO_FLUID_ASSEMBLY_MANAGER_FACTORY_H

#include <Teuchos_RCP.hpp>
#include <vector>

namespace DRT
{
  class Element;
  class Discretization;
}  // namespace DRT

namespace FBI
{
  class BeamToFluidMeshtyingParams;
}
namespace BEAMINTERACTION
{
  class BeamContactPair;

  namespace SUBMODELEVALUATOR
  {
    class PartitionedBeamInteractionAssemblyManager;
  }
  // namespace BEAMINTERACTION
  /**
   *  \brief Factory that creates the appropriate beam to fluid meshtying assembly manager for the
   * desired discretization
   *
   */
  class BeamToFluidAssemblyManagerFactory
  {
   private:
    /// constructor
    BeamToFluidAssemblyManagerFactory(){};

   public:
    /// destructor
    ~BeamToFluidAssemblyManagerFactory(){};

    /**
     *  \brief Creates the appropriate beam to fluid meshtying assembly manager for the desired
     * discretization
     *
     * This function is static so that it can be called without creating a factory object first.
     * It can be called directly.
     *
     * \param[in] params_ptr Container containing the Fluid beam interaction parameters
     * \param[in] pairs Vector of possible fluid beam interaction pairs
     *
     * \return beam interaction assembly manager
     */
    static Teuchos::RCP<
        BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManager>
    CreateAssemblyManager(Teuchos::RCP<const DRT::Discretization> discretization1,
        Teuchos::RCP<const DRT::Discretization> discretization2,
        std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> interaction_pairs,
        const Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> params_ptr);
  };
}  // namespace BEAMINTERACTION

#endif
