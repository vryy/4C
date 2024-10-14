/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create appropriate beam to fluid meshtying assembly managers for the desired
constraint discretization approach


\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FBI_BEAM_TO_FLUID_ASSEMBLY_MANAGER_FACTORY_HPP
#define FOUR_C_FBI_BEAM_TO_FLUID_ASSEMBLY_MANAGER_FACTORY_HPP

#include "4C_config.hpp"

#include <Teuchos_RCP.hpp>

#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

namespace FBI
{
  class BeamToFluidMeshtyingParams;
  namespace Utils
  {
    class FBIAssemblyStrategy;
  }
}  // namespace FBI
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
    BeamToFluidAssemblyManagerFactory() = delete;

   public:
    /**
     *  \brief Creates the appropriate beam to fluid meshtying assembly manager for the desired
     * discretizations
     *
     * This function is static so that it can be called without creating a factory object first.
     * It can be called directly.
     *
     * \param[in] params_ptr Container containing the Fluid beam interaction parameters
     * \param[in] interaction_pairs Vector of possible fluid beam interaction pairs
     * \params[in] assemblystrategy object handling the assembly into the global fluid matrix
     *
     * \return beam interaction assembly manager
     */
    static Teuchos::RCP<
        BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManager>
    create_assembly_manager(Teuchos::RCP<const Core::FE::Discretization> discretization1,
        Teuchos::RCP<const Core::FE::Discretization> discretization2,
        std::vector<Teuchos::RCP<BEAMINTERACTION::BeamContactPair>> interaction_pairs,
        const Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> params_ptr,
        Teuchos::RCP<FBI::Utils::FBIAssemblyStrategy> assemblystrategy);
  };
}  // namespace BEAMINTERACTION

FOUR_C_NAMESPACE_CLOSE

#endif
