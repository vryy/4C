/*----------------------------------------------------------------------*/
/*! \file

\brief Gauss point to segment mesh tying element for between a 3D beam and a surface element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_surface_meshtying_pair_gauss_point.H"

#include "beam_contact_params.H"
#include "beam_to_solid_surface_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_line_to_surface.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_factory.H"
#include "../drt_geometry_pair/geometry_pair_element_faces.H"


/**
 *
 */
template <typename beam, typename surface>
BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<beam,
    surface>::BeamToSolidSurfaceMeshtyingPairGaussPoint()
    : base_class()
{
  // Empty constructor.
}

/**
 *
 */
template <typename beam, typename surface>
bool BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<beam, surface>::Evaluate(
    LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  // Call Evaluate on the geometry Pair. Only do this once for mesh tying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(this->ele1posref_,
        this->face_element_->GetFaceReferencePosition(), this->line_to_3D_segments_,
        &(this->face_element_->GetReferenceNormals()));
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, return no contact status.
  if (this->line_to_3D_segments_.size() == 0) return false;
#if 0
  // Initialize variables for position and force vectors.
  LINALG::Matrix<3, 1, double> dr_beam_ref;
  LINALG::Matrix<3, 1, scalar_type_fad> r_beam;
  LINALG::Matrix<3, 1, scalar_type_fad> r_surface;
  LINALG::Matrix<3, 1, scalar_type_fad> force;
  LINALG::Matrix<beam::n_dof_, 1, scalar_type_fad> force_element_1(true);
  LINALG::Matrix<surface::n_dof_, 1, scalar_type_fad> force_element_2(true);

  // Initialize scalar variables.
  double segment_jacobian, beam_segmentation_factor;
  double penalty_parameter =
      this->Params()->BeamToSolidSurfaceMeshtyingParams()->GetPenaltyParameter();

  // Calculate the meshtying forces.
  // Loop over segments.
  for (unsigned int i_segment = 0; i_segment < this->line_to_3D_segments_.size(); i_segment++)
  {
    // Factor to account for a segment length not from -1 to 1.
    beam_segmentation_factor = 0.5 * this->line_to_3D_segments_[i_segment].GetSegmentLength();

    // Gauss point loop.
    for (unsigned int i_gp = 0;
         i_gp < this->line_to_3D_segments_[i_segment].GetProjectionPoints().size(); i_gp++)
    {
      // Get the current Gauss point.
      const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& projected_gauss_point =
          this->line_to_3D_segments_[i_segment].GetProjectionPoints()[i_gp];

      // Get the jacobian in the reference configuration.
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Get the current positions on beam and solid.
      GEOMETRYPAIR::EvaluatePosition<beam>(
          projected_gauss_point.GetEta(), this->ele1pos_, r_beam, this->Element1());



      LINALG::Matrix<3, 1, double> xi2;
      LINALG::Matrix<3, 1, double> r2;
      GEOMETRYPAIR::EvaluateSurfacePosition<surface>(
          this->face_element_->GetFaceReferencePosition(), xi2, r2, this->Element2(),
          &(this->face_element_->GetReferenceNormals()));



      // Calculate the force in this Gauss point. The sign of the force calculated here is the one
      // that acts on the beam.
      force = r_surface;
      force -= r_beam;
      force.Scale(penalty_parameter);

      // The force vector is in R3, we need to calculate the equivalent nodal forces on the element
      // dof. This is done with the virtual work equation $F \delta r = f \delta q$.
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_element_1(i_dof) += force(i_dir) * r_beam(i_dir).dx(i_dof) *
                                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
      for (unsigned int i_dof = 0; i_dof < surface::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_element_2(i_dof) -= force(i_dir) * r_surface(i_dir).dx(i_dof + beam::n_dof_) *
                                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
    }
  }


  // Fill in the entries for the local matrices and vectors.
  {
    // Resize and initialize the return variables.
    if (forcevec1 != NULL) forcevec1->Size(beam::n_dof_);
    if (forcevec2 != NULL) forcevec2->Size(surface::n_dof_);
    if (stiffmat11 != NULL) stiffmat11->Shape(beam::n_dof_, beam::n_dof_);
    if (stiffmat12 != NULL) stiffmat12->Shape(beam::n_dof_, surface::n_dof_);
    if (stiffmat21 != NULL) stiffmat21->Shape(surface::n_dof_, beam::n_dof_);
    if (stiffmat22 != NULL) stiffmat22->Shape(surface::n_dof_, surface::n_dof_);

    if (forcevec1 != NULL && forcevec2 != NULL)
    {
      // $f_1$
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        (*forcevec1)(i_dof) = force_element_1(i_dof).val();
      // $f_2$
      for (unsigned int i_dof = 0; i_dof < surface::n_dof_; i_dof++)
        (*forcevec2)(i_dof) = force_element_2(i_dof).val();
    }

    if (stiffmat11 != NULL && stiffmat12 != NULL && stiffmat21 != NULL && stiffmat22 != NULL)
    {
      // $k_{11}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < beam::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < beam::n_dof_; i_dof_2++)
          (*stiffmat11)(i_dof_1, i_dof_2) = -force_element_1(i_dof_1).dx(i_dof_2);

      // $k_{12}, k_{21}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < beam::n_dof_; i_dof_1++)
      {
        for (unsigned int i_dof_2 = 0; i_dof_2 < surface::n_dof_; i_dof_2++)
        {
          (*stiffmat12)(i_dof_1, i_dof_2) = -force_element_1(i_dof_1).dx(beam::n_dof_ + i_dof_2);
          (*stiffmat21)(i_dof_2, i_dof_1) = -force_element_2(i_dof_2).dx(i_dof_1);
        }
      }

      // $k_{22}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < surface::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < surface::n_dof_; i_dof_2++)
          (*stiffmat22)(i_dof_1, i_dof_2) = -force_element_2(i_dof_1).dx(beam::n_dof_ + i_dof_2);
    }
  }
#endif
  // Return true as there are meshtying contributions.
  return false;
}


/**
 * Explicit template initialization of template class.
 */
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri3>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri6>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad4>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad8>;
template class BEAMINTERACTION::BeamToSolidSurfaceMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad9>;
