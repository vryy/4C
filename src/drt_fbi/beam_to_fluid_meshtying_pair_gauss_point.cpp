/*!
\file beam_to_fluid_meshtying_pair_gauss_point.cpp

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_fluid_meshtying_pair_gauss_point.H"
#include "beam_to_fluid_meshtying_pair_base.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_beaminteraction/beam_contact_params.H"
#include "beam_to_fluid_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_element_types.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume.H"


/**
 *
 */
template <typename beam, typename fluid>
BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<beam,
    fluid>::BeamToFluidMeshtyingPairGaussPoint()
    : BeamToFluidMeshtyingPairBase<beam, fluid>()
{
  // Empty constructor.
}


/**
 *
 */
template <typename beam, typename fluid>
bool BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<beam, fluid>::Evaluate(
    LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(
        this->ele1poscur_, this->ele2poscur_, this->line_to_volume_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, return no contact status.
  if (this->line_to_volume_segments_.size() == 0) return false;

  // Initialize variables for position and force vectors.
  LINALG::TMatrix<double, 3, 1> dr_beam_ref;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> r_beam;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> r_fluid;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> v_beam;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> v_fluid;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> force;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> force2;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, beam::n_dof_, 1> force_element_1(true);
  LINALG::TMatrix<TYPE_BTS_VMT_AD, fluid::n_dof_, 1> force_element_2(true);
  LINALG::TMatrix<TYPE_BTS_VMT_AD, fluid::n_dof_, 1> force_element_f(true);

  // Initialize scalar variables.
  double segment_jacobian, beam_segmentation_factor;

  // Calculate the meshtying forces.
  // Loop over segments.
  for (unsigned int i_segment = 0; i_segment < this->line_to_volume_segments_.size(); i_segment++)
  {
    // Factor to account for a segment length not from -1 to 1.
    beam_segmentation_factor = 0.5 * this->line_to_volume_segments_[i_segment].GetSegmentLength();

    // Gauss point loop.
    for (unsigned int i_gp = 0;
         i_gp < this->line_to_volume_segments_[i_segment].GetProjectionPoints().size(); i_gp++)
    {
      // Get the current Gauss point.
      const GEOMETRYPAIR::ProjectionPointLineToVolume<double>& projected_gauss_point =
          this->line_to_volume_segments_[i_segment].GetProjectionPoints()[i_gp];

      // Get the jacobian in the reference configuration.
      GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
          projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());

      // Jacobian including the segment length.
      segment_jacobian = dr_beam_ref.Norm2() * beam_segmentation_factor;

      // Get the current positions on beam and fluid.
      GEOMETRYPAIR::EvaluatePosition<beam>(projected_gauss_point.GetEta(), this->ele1pos_, r_beam,
          this->Element1());  // need beam velocity here -> hand in from outside? No..
      GEOMETRYPAIR::EvaluatePosition<fluid>(projected_gauss_point.GetXi(), this->ele2pos_,
          r_fluid);  // need fluid velocity here -> hand in from outside? No.. we do not compute the
                     // force here

      // extrapolate the current velocity on beam and fluid.
      GEOMETRYPAIR::EvaluateVelocity<beam>(projected_gauss_point.GetEta(), this->ele1vel_, v_beam,
          this->Element1());  // need beam velocity here -> hand in from outside -> ele1vel_?
      GEOMETRYPAIR::EvaluateVelocity<fluid>(projected_gauss_point.GetXi(), this->ele2vel_,
          v_fluid);  // need fluid velocity here -> hand in from outside -> ele2vel_?

      // Calculate the force in this Gauss point. The sign of the force calculated here is the one
      // that acts on the fluid.  todo check consistency
      force += v_fluid;
      force -= v_beam;
      //      force.Scale(penalty_parameter); we want to scale this later on in order to seperate
      //      discretization approach from constraitn enforcement technique

      // The force vector is in R3, we need to calculate the equivalent nodal forces on the
      // element dof. This is done by a basis transformation. We do not need this for fixed grid
      // fluid, but for beams an ALE todo Check/implement update routine for elepos_ and thus
      // r_fluid
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_element_1(i_dof) += force(i_dir) * r_beam(i_dir).dx(i_dof) *
                                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
      for (unsigned int i_dof = 0; i_dof < fluid::n_dof_; i_dof++)
        for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
          force_element_2(i_dof) -= force(i_dir) * r_fluid(i_dir).dx(i_dof + beam::n_dof_) *
                                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
      if (Teuchos::rcp_dynamic_cast<FBI::BeamToFluidMeshtyingParams>(this->Params(), true)
              ->GetWeakDirichletFlag())
      {
        force2 -= v_beam;
        for (unsigned int i_dof = 0; i_dof < fluid::n_dof_; i_dof++)
          for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
            force_element_f(i_dof) -= force2(i_dir) * r_fluid(i_dir).dx(i_dof + beam::n_dof_) *
                                      projected_gauss_point.GetGaussWeight() * segment_jacobian;
      }
    }
  }

  // Fill in the entries for the local matrices and vectors.
  {
    // Resize and initialize the return variables.
    if (forcevec1 != NULL) forcevec1->Size(beam::n_dof_);
    if (forcevec2 != NULL) forcevec2->Size(fluid::n_dof_);
    if (stiffmat11 != NULL) stiffmat11->Shape(beam::n_dof_, beam::n_dof_);
    if (stiffmat12 != NULL) stiffmat12->Shape(beam::n_dof_, fluid::n_dof_);
    if (stiffmat21 != NULL) stiffmat21->Shape(fluid::n_dof_, beam::n_dof_);
    if (stiffmat22 != NULL) stiffmat22->Shape(fluid::n_dof_, fluid::n_dof_);

    if (forcevec1 != NULL)
    {
      // $f_1$
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        (*forcevec1)(i_dof) = force_element_1(i_dof).val();
    }
    if (forcevec2 != NULL)
    {
      // $f_2$
      if (Teuchos::rcp_dynamic_cast<FBI::BeamToFluidMeshtyingParams>(this->Params(), true)
              ->GetWeakDirichletFlag())
      {
        for (unsigned int i_dof = 0; i_dof < fluid::n_dof_; i_dof++)
          (*forcevec2)(i_dof) = force_element_f(i_dof).val();
      }
      else
      {
        for (unsigned int i_dof = 0; i_dof < fluid::n_dof_; i_dof++)
          (*forcevec2)(i_dof) = force_element_2(i_dof).val();
      }
    }

    if (stiffmat11 != NULL)
    {
      // $k_{11}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < beam::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < beam::n_dof_; i_dof_2++)
          (*stiffmat11)(i_dof_1, i_dof_2) = -force_element_1(i_dof_1).dx(i_dof_2);
    }
    if (stiffmat12 != NULL || stiffmat21 != NULL)
      // $k_{12}, k_{21}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < beam::n_dof_; i_dof_1++)
      {
        for (unsigned int i_dof_2 = 0; i_dof_2 < fluid::n_dof_; i_dof_2++)
        {
          (*stiffmat12)(i_dof_1, i_dof_2) = -force_element_1(i_dof_1).dx(beam::n_dof_ + i_dof_2);
          (*stiffmat21)(i_dof_2, i_dof_1) = -force_element_2(i_dof_2).dx(i_dof_1);
        }
      }

    if (stiffmat22 != NULL)
    {
      // $k_{22}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < fluid::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < fluid::n_dof_; i_dof_2++)
          (*stiffmat22)(i_dof_1, i_dof_2) = -force_element_2(i_dof_1).dx(beam::n_dof_ + i_dof_2);
    }
  }
  // Return true as there are meshtying contributions.
  return true;
}


/**
 * Explicit template initialization of template class.
 */
// Hermite beam element, hex8 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8>;
// Hermite beam element, hex20 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20>;
// Hermite beam element, hex27 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27>;
// Hermite beam element, tet4 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4>;
// Hermite beam element, tet10 solid element.
template class BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10>;
