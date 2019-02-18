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

#include "../drt_geometry_pair/geometry_pair_element_types.H"
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
 *
 */
template <typename beam, typename solid, typename mortar>
bool BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<beam, solid, mortar>::EvaluateDM(
    LINALG::SerialDenseMatrix& mortar_D, LINALG::SerialDenseMatrix& mortar_M)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(
        this->ele1posref_, this->ele2posref_, this->line_to_volume_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no intersection segments, return no contact status.
  if (this->line_to_volume_segments_.size() == 0) return false;

  // Initialize variables for local mortar matrices.
  LINALG::TMatrix<double, mortar::n_dof_, beam::n_dof_> D(true);
  LINALG::TMatrix<double, mortar::n_dof_, solid::n_dof_> M(true);

  // Initialize variables for shape function values.
  LINALG::TMatrix<double, 1, mortar::n_nodes_ * mortar::n_val_> N_mortar(true);
  LINALG::TMatrix<double, 1, beam::n_nodes_ * beam::n_val_> N_beam(true);
  LINALG::TMatrix<double, 1, solid::n_nodes_ * solid::n_val_> N_solid(true);

  // Initialize variable for beam position derivative.
  LINALG::TMatrix<double, 3, 1> dr_beam_ref(true);

  // Initialize scalar variables.Clear
  double segment_jacobian, beam_segmentation_factor;

  // Calculate the mortar matrices.
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

      // Get the shape function matrices.
      N_mortar.Clear();
      N_beam.Clear();
      N_solid.Clear();
      mortar::EvaluateShapeFunction(
          N_mortar, projected_gauss_point.GetEta(), std::integral_constant<unsigned int, 1>{});
      beam::EvaluateShapeFunction(N_beam, projected_gauss_point.GetEta(),
          std::integral_constant<unsigned int, 1>{}, this->Element1());
      solid::EvaluateShapeFunction(
          N_solid, projected_gauss_point.GetXi(), std::integral_constant<unsigned int, 3>{});

      // Fill in the local mortar matrix D.
      for (unsigned int i_mortar_node = 0; i_mortar_node < mortar::n_nodes_; i_mortar_node++)
        for (unsigned int i_mortar_val = 0; i_mortar_val < mortar::n_val_; i_mortar_val++)
          for (unsigned int i_beam_node = 0; i_beam_node < beam::n_nodes_; i_beam_node++)
            for (unsigned int i_beam_val = 0; i_beam_val < beam::n_val_; i_beam_val++)
              for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
                D(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim,
                    i_beam_node * beam::n_val_ * 3 + i_beam_val * 3 + i_dim) +=
                    N_mortar(i_mortar_node * mortar::n_val_ + i_mortar_val) *
                    N_beam(i_beam_node * beam::n_val_ + i_beam_val) *
                    projected_gauss_point.GetGaussWeight() * segment_jacobian;

      // Fill in the local mortar matrix M.
      for (unsigned int i_mortar_node = 0; i_mortar_node < mortar::n_nodes_; i_mortar_node++)
        for (unsigned int i_mortar_val = 0; i_mortar_val < mortar::n_val_; i_mortar_val++)
          for (unsigned int i_solid_node = 0; i_solid_node < solid::n_nodes_; i_solid_node++)
            for (unsigned int i_solid_val = 0; i_solid_val < solid::n_val_; i_solid_val++)
              for (unsigned int i_dim = 0; i_dim < 3; i_dim++)
                M(i_mortar_node * mortar::n_val_ * 3 + i_mortar_val * 3 + i_dim,
                    i_solid_node * solid::n_val_ * 3 + i_solid_val * 3 + i_dim) +=
                    N_mortar(i_mortar_node * mortar::n_val_ + i_mortar_val) *
                    N_solid(i_solid_node * solid::n_val_ + i_solid_val) *
                    projected_gauss_point.GetGaussWeight() * segment_jacobian;
    }
  }

  // If we get to this point, the pair has a mortar contribution.
  return true;
}


/**
 * Explicit template initialization of template class.
 */
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line2>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line2>;

template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line3>;
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line3>;
