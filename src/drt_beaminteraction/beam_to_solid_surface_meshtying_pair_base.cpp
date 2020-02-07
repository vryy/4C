/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_base.H"

#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"


/**
 *
 */
template <typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<beam,
    surface>::BeamToSolidSurfaceMeshtyingPairBase()
    : base_class(), meshtying_is_evaluated_(false)
{
  // Empty constructor.
}

/**
 *
 */
template <typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<beam, surface>::CreateGeometryPair(
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>& geometry_evaluation_data_ptr)
{
  // Call the method of the base class.
  BeamContactPair::CreateGeometryPair(geometry_evaluation_data_ptr);

  // Set up the geometry pair, it will be initialized in the Init call of the base class.
  this->geometry_pair_ = GEOMETRYPAIR::GeometryPairLineToSurfaceFactory<double, beam, surface>(
      geometry_evaluation_data_ptr);
}

/**
 *
 */
template <typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<beam, surface>::SetFaceElement(
    Teuchos::RCP<GEOMETRYPAIR::FaceElement>& face_element)
{
  face_element_ =
      Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::FaceElementTempalte<surface>>(face_element, true);
}

/**
 *
 */
template <typename beam, typename surface>
Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToSurface<double, beam, surface>>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<beam, surface>::CastGeometryPair() const
{
  return Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::GeometryPairLineToSurface<double, beam, surface>>(
      this->geometry_pair_, true);
};


/**
 * Explicit template initialization of template class.
 */
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri3>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri6>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad4>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad8>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad9>;
