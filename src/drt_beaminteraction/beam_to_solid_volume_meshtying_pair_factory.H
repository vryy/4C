/*!
\file beam_to_solid_volume_meshtying_pair_factory.H

\brief Create the pairs for beam to solid volume meshtying depending on the input parameters..

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef BEAM_TO_SOLID_VOLUME_MESHTYING_FACTORY_H_
#define BEAM_TO_SOLID_VOLUME_MESHTYING_FACTORY_H_


#include "beam_contact_params.H"
#include "beam_contact_pair.H"


namespace BEAMINTERACTION
{
  /**
   * Create the correct type of beam to solid volume meshtying element, depending on the simulation
   * parameters.
   * @param ele_ptrs Vector containing pointer to the elements of the pair.
   * @param params_ptr Beam contact parameters.
   * @return Teuchos pointer to the created pair.
   */
  Teuchos::RCP<BeamContactPair> BeamToSolidVolumeMeshtyingPairFactory(
      std::vector<DRT::Element const*> const& ele_ptrs,
      const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr);
}  // namespace BEAMINTERACTION

#endif
