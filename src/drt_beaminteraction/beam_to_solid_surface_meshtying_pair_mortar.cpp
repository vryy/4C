/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_mortar.H"

#include "beam_contact_params.H"
#include "beam_to_solid_surface_meshtying_params.H"
#include "beaminteraction_calc_utils.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"

#include "Epetra_FEVector.h"


/**
 *
 */
template <typename beam, typename surface, typename mortar>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairMortar<beam, surface,
    mortar>::BeamToSolidSurfaceMeshtyingPairMortar()
    : base_class()
{
  // Empty constructor.
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri3, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri6, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad4, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad8, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad9, t_line2>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_nurbs9, t_line2>;

  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri3, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri6, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad4, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad8, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad9, t_line3>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_nurbs9, t_line3>;

  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri3, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_tri6, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad4, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad8, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_quad9, t_line4>;
  template class BeamToSolidSurfaceMeshtyingPairMortar<t_hermite, t_nurbs9, t_line4>;
}  // namespace BEAMINTERACTION
