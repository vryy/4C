/*----------------------------------------------------------------------*/
/*! \file

\brief Base meshtying element for meshtying between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_base.H"

#include "beam_to_solid_vtu_output_writer_base.H"
#include "beam_to_solid_vtu_output_writer_visualization.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"


/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::BeamToSolidSurfaceMeshtyingPairBase()
    : base_class(), meshtying_is_evaluated_(false)
{
  // Empty constructor.
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  if (!meshtying_is_evaluated_)
  {
    CastGeometryPair()->PreEvaluate(this->ele1posref_,
        this->face_element_->GetFaceReferencePosition(), this->line_to_3D_segments_,
        this->face_element_->GetReferenceNormals());
  }
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::GetPairVisualization(Teuchos::RCP<BeamToSolidVtuOutputWriterBase>
                                       visualization_writer,
    const Teuchos::ParameterList& visualization_params) const
{
  // Get visualization of base class.
  base_class::GetPairVisualization(visualization_writer, visualization_params);
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::CreateGeometryPair(const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>&
        geometry_evaluation_data_ptr)
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
template <typename scalar_type, typename beam, typename surface>
void BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam,
    surface>::SetFaceElement(Teuchos::RCP<GEOMETRYPAIR::FaceElement>& face_element)
{
  face_element_ =
      Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::FaceElementTemplate<surface>>(face_element, true);

  // The second element in the pair has to be the face element.
  CastGeometryPair()->SetElement2(face_element_->GetDrtFaceElement());
}

/**
 *
 */
template <typename scalar_type, typename beam, typename surface>
Teuchos::RCP<GEOMETRYPAIR::GeometryPairLineToSurface<double, beam, surface>>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<scalar_type, beam, surface>::CastGeometryPair()
    const
{
  return Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::GeometryPairLineToSurface<double, beam, surface>>(
      this->geometry_pair_, true);
};


/**
 * Explicit template initialization of template class.
 */
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_tri3::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri3>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_tri6::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tri6>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_quad4::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad4>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_quad8::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad8>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_quad9::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_quad9>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairBase<
    Sacado::ELRFad::SLFad<double, GEOMETRYPAIR::t_hermite::n_dof_ + GEOMETRYPAIR::t_nurbs9::n_dof_>,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_nurbs9>;
