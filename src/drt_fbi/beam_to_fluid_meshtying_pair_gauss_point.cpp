/*----------------------------------------------------------------------*/
/*! \file

\brief Meshtying element for meshtying between a 1D beam and a 3D fluid element.

\level 1
\maintainer Nora Hagmeyer
*/


#include "beam_to_fluid_meshtying_pair_gauss_point.H"
#include "beam_to_fluid_meshtying_pair_base.H"

#include "../linalg/linalg_fixedsizematrix.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_beaminteraction/beam_contact_params.H"
#include "beam_to_fluid_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
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

/*------------------------------------------------------------------------------------------------*/

template <typename beam, typename fluid>
bool BEAMINTERACTION::BeamToFluidMeshtyingPairGaussPoint<beam, fluid>::Evaluate(
    LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  // Call Evaluate on the geometry Pair
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(
        this->ele1poscur_, this->ele2poscur_, this->line_to_3D_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, return no contact status.
  if (this->line_to_3D_segments_.size() == 0) return false;

  // Initialize variables for position and force vectors.
  LINALG::Matrix<3, 1, double> dr_beam_ref;
  LINALG::Matrix<3, 1, scalar_type> r_beam;
  LINALG::Matrix<3, 1, scalar_type> r_fluid;
  LINALG::Matrix<3, 1, scalar_type> v_beam;
  LINALG::Matrix<3, 1, scalar_type> v_fluid;
  LINALG::Matrix<3, 1, scalar_type> force;
  LINALG::Matrix<3, 1, scalar_type> force2;
  LINALG::Matrix<beam::n_dof_, 1, scalar_type> force_element_1(true);
  LINALG::Matrix<fluid::n_dof_, 1, scalar_type> force_element_2(true);
  LINALG::Matrix<fluid::n_dof_, 1, scalar_type> force_element_f(true);
  LINALG::Matrix<1, beam::n_nodes_ * beam::n_val_, double> N_beam(true);
  LINALG::Matrix<1, fluid::n_nodes_ * fluid::n_val_, double> N_fluid(true);

  // Resize and initialize the return variables.
  if (forcevec1 != NULL) forcevec1->Size(beam::n_dof_);
  if (forcevec2 != NULL) forcevec2->Size(fluid::n_dof_);
  if (stiffmat11 != NULL) stiffmat11->Shape(beam::n_dof_, beam::n_dof_);
  if (stiffmat12 != NULL) stiffmat12->Shape(beam::n_dof_, fluid::n_dof_);
  if (stiffmat21 != NULL) stiffmat21->Shape(fluid::n_dof_, beam::n_dof_);
  if (stiffmat22 != NULL) stiffmat22->Shape(fluid::n_dof_, fluid::n_dof_);

  // Initialize scalar variables.
  double segment_jacobian, beam_segmentation_factor;

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

      // Get the current positions on beam and fluid.
      GEOMETRYPAIR::EvaluatePosition<beam>(
          projected_gauss_point.GetEta(), this->ele1pos_, r_beam, this->Element1());
      GEOMETRYPAIR::EvaluatePosition<fluid>(projected_gauss_point.GetXi(), this->ele2pos_, r_fluid);

      N_beam.Clear();
      N_fluid.Clear();

      // Evaluate the chapefunctions at the current gauss point
      beam::EvaluateShapeFunction(N_beam, projected_gauss_point.GetEta(),
          std::integral_constant<unsigned int, 1>{}, this->Element1());
      fluid::EvaluateShapeFunction(N_fluid, projected_gauss_point.GetXi(),
          std::integral_constant<unsigned int, 3>{}, this->Element2());

      // assemble fluid mass matrix
      if (stiffmat22 != NULL)
      {
        for (unsigned int i_fluid_node1 = 0; i_fluid_node1 < fluid::n_nodes_; i_fluid_node1++)
          for (unsigned int i_fluid_val1 = 0; i_fluid_val1 < fluid::n_val_; i_fluid_val1++)
            for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
              for (unsigned int i_fluid_node2 = 0; i_fluid_node2 < fluid::n_nodes_; i_fluid_node2++)
                for (unsigned int i_fluid_val2 = 0; i_fluid_val2 < fluid::n_val_; i_fluid_val2++)
                  (*stiffmat22)(i_fluid_node1 * fluid::n_val_ * 3 + i_fluid_val1 * 3 + i_dir,
                      i_fluid_node2 * fluid::n_val_ * 3 + i_fluid_val2 * 3 + i_dir) +=
                      N_fluid(i_fluid_node1) * N_fluid(i_fluid_node2) *
                      projected_gauss_point.GetGaussWeight() * segment_jacobian;
      }

      // assemble beam mass matrix
      if (stiffmat11 != NULL)
      {
        for (unsigned int i_beam_node1 = 0; i_beam_node1 < beam::n_nodes_; i_beam_node1++)
          for (unsigned int i_beam_val1 = 0; i_beam_val1 < beam::n_val_; i_beam_val1++)
            for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
              for (unsigned int i_beam_node2 = 0; i_beam_node2 < beam::n_nodes_; i_beam_node2++)
                for (unsigned int i_beam_val2 = 0; i_beam_val2 < beam::n_val_; i_beam_val2++)
                  (*stiffmat11)(i_beam_node1 * beam::n_val_ * 3 + 3 * i_beam_val1 + i_dir,
                      i_beam_node2 * beam::n_val_ * 3 + 3 * i_beam_val2 + i_dir) +=
                      N_beam(i_beam_node1 * beam::n_val_ + i_beam_val1) *
                      N_beam(i_beam_node2 * beam::n_val_ + i_beam_val2) *
                      projected_gauss_point.GetGaussWeight() * segment_jacobian;
      }

      // assemble fluid beam coupling matrix
      if (stiffmat21 != NULL)
      {
        for (unsigned int i_fluid_node1 = 0; i_fluid_node1 < fluid::n_nodes_; i_fluid_node1++)
          for (unsigned int i_fluid_val1 = 0; i_fluid_val1 < fluid::n_val_; i_fluid_val1++)
            for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
              for (unsigned int i_beam_node2 = 0; i_beam_node2 < beam::n_nodes_; i_beam_node2++)
                for (unsigned int i_beam_val2 = 0; i_beam_val2 < beam::n_val_; i_beam_val2++)
                  (*stiffmat21)(i_fluid_node1 * fluid::n_val_ * 3 + i_fluid_val1 * 3 + i_dir,
                      i_beam_node2 * beam::n_val_ * 3 + 3 * i_beam_val2 + i_dir) +=
                      N_fluid(i_fluid_node1) * N_beam(i_beam_node2 * beam::n_val_ + i_beam_val2) *
                      projected_gauss_point.GetGaussWeight() * segment_jacobian;
      }

      // assemble fluid beam coupling matrix
      if (stiffmat12 != NULL)
      {
        for (unsigned int i_beam_node1 = 0; i_beam_node1 < beam::n_nodes_; i_beam_node1++)
          for (unsigned int i_beam_val1 = 0; i_beam_val1 < beam::n_val_; i_beam_val1++)
            for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
              for (unsigned int i_fluid_node2 = 0; i_fluid_node2 < fluid::n_nodes_; i_fluid_node2++)
                for (unsigned int i_fluid_val2 = 0; i_fluid_val2 < fluid::n_val_; i_fluid_val2++)
                  (*stiffmat12)(i_beam_node1 * beam::n_val_ * 3 + 3 * i_beam_val1 + i_dir,
                      i_fluid_node2 * fluid::n_val_ * 3 + 3 * i_fluid_val2 + i_dir) +=
                      N_fluid(i_fluid_node2) * N_beam(i_beam_node1 * beam::n_val_ + i_beam_val1) *
                      projected_gauss_point.GetGaussWeight() * segment_jacobian;
      }
    }
  }

  // assemble structure force vector
  {
    if (forcevec1 != NULL)
    {
      for (unsigned int i_dof1 = 0; i_dof1 < beam::n_dof_; i_dof1++)
      {
        for (unsigned int i_dof2 = 0; i_dof2 < beam::n_dof_; i_dof2++)
          (*forcevec1)(i_dof1) +=
              (*stiffmat11)(i_dof1, i_dof2) * FADUTILS::CastToDouble(this->ele1vel_(i_dof2));
        for (unsigned int i_dof2 = 0; i_dof2 < fluid::n_dof_; i_dof2++)
          (*forcevec1)(i_dof1) -=
              (*stiffmat12)(i_dof1, i_dof2) * FADUTILS::CastToDouble(this->ele2vel_(i_dof2));
      }
    }

    // assemble fluid force vector
    if (forcevec2 != NULL)
    {
      if (!Teuchos::rcp_dynamic_cast<FBI::BeamToFluidMeshtyingParams>(this->Params(), true)
               ->GetWeakDirichletFlag())

      {
        for (unsigned int i_dof1 = 0; i_dof1 < fluid::n_dof_; i_dof1++)
        {
          for (unsigned int i_dof2 = 0; i_dof2 < fluid::n_dof_; i_dof2++)
            (*forcevec2)(i_dof1) +=
                (*stiffmat22)(i_dof1, i_dof2) * FADUTILS::CastToDouble(this->ele2vel_(i_dof2));
          for (unsigned int i_dof2 = 0; i_dof2 < beam::n_dof_; i_dof2++)
            (*forcevec2)(i_dof1) -=
                (*stiffmat21)(i_dof1, i_dof2) * FADUTILS::CastToDouble(this->ele1vel_(i_dof2));
        }
      }
      else

        for (unsigned int i_dof1 = 0; i_dof1 < fluid::n_dof_; i_dof1++)
        {
          for (unsigned int i_dof2 = 0; i_dof2 < beam::n_dof_; i_dof2++)
            (*forcevec2)(i_dof1) -=
                (*stiffmat21)(i_dof1, i_dof2) * FADUTILS::CastToDouble(this->ele1vel_(i_dof2));
        }
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
