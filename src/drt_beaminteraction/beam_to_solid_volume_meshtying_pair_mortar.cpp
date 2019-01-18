/*!
\file beam_to_solid_volume_meshtying_pair_mortar.cpp

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element using mortar shape
functions for the traction.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_volume_meshtying_pair_mortar.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "beam_contact_params.H"
#include "beam_to_solid_volume_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume.H"


/**
 *
 */
template <typename beam, typename solid, typename mortar>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam, solid,
    mortar>::BeamToSolidVolumeMeshtyingPairMortar()
    : BeamToSolidVolumeMeshtyingPairBase<beam, solid>()
{
  // Empty constructor.
}


/**
 * Explicit template initialization of template class.
 */
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8, BEAMINTERACTION::t_mortar_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20, BEAMINTERACTION::t_mortar_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27, BEAMINTERACTION::t_mortar_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4, BEAMINTERACTION::t_mortar_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10, BEAMINTERACTION::t_mortar_line2>;

template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8, BEAMINTERACTION::t_mortar_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20, BEAMINTERACTION::t_mortar_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27, BEAMINTERACTION::t_mortar_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4, BEAMINTERACTION::t_mortar_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10, BEAMINTERACTION::t_mortar_line3>;
