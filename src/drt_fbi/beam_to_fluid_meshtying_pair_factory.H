/*----------------------------------------------------------------------*/
/*! \file

\brief Factory to create appropriate beam to fluid meshtying pair for the desired discretization
approach

\maintainer Nora Hagmeyer

\level 2
*/
/*----------------------------------------------------------------------*/

#ifndef BEAM_TO_FLUID_MESHTYING_PAIR_FACTORY_H
#define BEAM_TO_FLUID_MESHTYING_PAIR_FACTORY_H

#include <vector>

#include <Teuchos_RCP.hpp>

namespace DRT
{
  class Element;
}
namespace BEAMINTERACTION
{
  class BeamContactPair;
}  // namespace BEAMINTERACTION
namespace FBI
{
  class BeamToFluidMeshtyingParams;
  /**
   *  \brief Factory that creates the appropriate beam to fluid meshtying pair for the desired
   * discretization
   *
   */
  class PairFactory
  {
   private:
    /// constructor
    PairFactory(){};

   public:
    /// destructor
    ~PairFactory(){};

    /**
     *  \brief Creates the appropriate beam to fluid meshtying pair for the desired
     * discretization
     *
     * This function is static so that it can be called without creating a factory object first.
     * It can be called directly.
     *
     * \param[in] params_ptr Container containing the Fluid beam interaction parameters
     * \param[in] ele_ptrs Pointer containing the elements which make up a pair
     *
     * \return Beam contact pair
     */
    static Teuchos::RCP<BEAMINTERACTION::BeamContactPair> CreatePair(
        std::vector<DRT::Element const*> const& ele_ptrs,
        const Teuchos::RCP<FBI::BeamToFluidMeshtyingParams> params_ptr);
  };
}  // namespace FBI

#endif
