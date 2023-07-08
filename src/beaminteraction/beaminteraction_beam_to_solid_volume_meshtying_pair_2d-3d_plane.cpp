/*----------------------------------------------------------------------*/
/*! \file

\brief Class for 2D-3D beam-to-solid volume mesh tying based on a plane beam element. This
simplifies the triad construction and torsion free beam elements can be used.

\level 3
*/


#include "beaminteraction_beam_to_solid_volume_meshtying_pair_2d-3d_plane.H"

#include "linalg_utils_densematrix_inverse.H"
#include "linalg_serialdensematrix.H"
#include "linalg_serialdensevector.H"

#include "beaminteraction_contact_params.H"
#include "beaminteraction_beam_to_solid_volume_meshtying_params.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_utility_classes.H"
#include "geometry_pair_line_to_3D_evaluation_data.H"
#include "geometry_pair_line_to_volume_gauss_point_projection_cross_section.H"
#include "beam3_triad_interpolation_local_rotation_vectors.H"


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair2D3DPlane<beam, solid>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  if (!this->meshtying_is_evaluated_)
  {
    CORE::LINALG::Matrix<beam::n_dof_, 1, double> beam_coupling_ref;
    CORE::LINALG::Matrix<solid::n_dof_, 1, double> solid_coupling_ref;
    this->GetCouplingReferencePosition(beam_coupling_ref, solid_coupling_ref);
    this->CastGeometryPair()->PreEvaluate(
        beam_coupling_ref, solid_coupling_ref, this->line_to_3D_segments_);
  }
}

/**
 *
 */
template <typename beam, typename solid>
bool BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair2D3DPlane<beam, solid>::Evaluate(
    CORE::LINALG::SerialDenseVector* forcevec1, CORE::LINALG::SerialDenseVector* forcevec2,
    CORE::LINALG::SerialDenseMatrix* stiffmat11, CORE::LINALG::SerialDenseMatrix* stiffmat12,
    CORE::LINALG::SerialDenseMatrix* stiffmat21, CORE::LINALG::SerialDenseMatrix* stiffmat22)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    CORE::LINALG::Matrix<beam::n_dof_, 1, double> beam_coupling_ref;
    CORE::LINALG::Matrix<solid::n_dof_, 1, double> solid_coupling_ref;
    this->GetCouplingReferencePosition(beam_coupling_ref, solid_coupling_ref);
    this->CastGeometryPair()->Evaluate(
        beam_coupling_ref, solid_coupling_ref, this->line_to_3D_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no segments, this pair has no contribution. Also there can be no more than one
  // segment.
  if (this->line_to_3D_segments_.size() == 0)
    return false;
  else if (this->line_to_3D_segments_.size() > 1)
    dserror("There can be a maximum of one segment!");

  // Get the vector with the projection points for this pair.
  const std::vector<GEOMETRYPAIR::ProjectionPoint1DTo3D<double>>& projection_points =
      this->line_to_3D_segments_[0].GetProjectionPoints();

  // If there are no projection points, return no contact status.
  if (projection_points.size() == 0) return false;

  // Initialize variables for position and force vectors.
  CORE::LINALG::Matrix<3, 1, double> dr_beam_ref;
  CORE::LINALG::Matrix<3, 1, scalar_type> r_beam;
  CORE::LINALG::Matrix<3, 3, scalar_type> triad;
  CORE::LINALG::Matrix<3, 1, scalar_type> r_cross_section_ref;
  CORE::LINALG::Matrix<3, 1, scalar_type> r_cross_section;
  CORE::LINALG::Matrix<3, 1, scalar_type> r_solid;
  CORE::LINALG::Matrix<3, 1, scalar_type> force;
  CORE::LINALG::Matrix<beam::n_dof_, 1, scalar_type> force_element_1(true);
  CORE::LINALG::Matrix<solid::n_dof_, 1, scalar_type> force_element_2(true);

  // Initialize scalar variables.
  double beam_jacobian;
  double penalty_parameter =
      this->Params()->BeamToSolidVolumeMeshtyingParams()->GetPenaltyParameter();

  // Calculate the meshtying forces.
  // Loop over segments.
  for (unsigned int i_integration_point = 0; i_integration_point < projection_points.size();
       i_integration_point++)
  {
    // Get the current Gauss point.
    const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>& projected_gauss_point =
        projection_points[i_integration_point];

    // Get the jacobian in the reference configuration.
    GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
        projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());
    beam_jacobian = 0.5 * dr_beam_ref.Norm2();

    // Get the current positions on beam and solid.
    GEOMETRYPAIR::EvaluatePosition<beam>(
        projected_gauss_point.GetEta(), this->ele1pos_, r_beam, this->Element1());
    GEOMETRYPAIR::EvaluateTriadAtPlaneCurve<beam>(
        projected_gauss_point.GetEta(), this->ele1pos_, triad, this->Element1());
    r_cross_section_ref(0) = 0.0;
    r_cross_section_ref(1) = projected_gauss_point.GetEtaCrossSection()(0);
    r_cross_section_ref(2) = projected_gauss_point.GetEtaCrossSection()(1);
    r_cross_section.Multiply(triad, r_cross_section_ref);
    r_beam += r_cross_section;
    GEOMETRYPAIR::EvaluatePosition<solid>(projected_gauss_point.GetXi(), this->ele2pos_, r_solid);

    // Calculate the force in this Gauss point. The sign of the force calculated here is the one
    // that acts on the beam.
    force = r_solid;
    force -= r_beam;
    force.Scale(penalty_parameter);

    // The force vector is in R3, we need to calculate the equivalent nodal forces on the element
    // dof. This is done with the virtual work equation $F \delta r = f \delta q$.
    for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        force_element_1(i_dof) += force(i_dir) * r_beam(i_dir).dx(i_dof) *
                                  projected_gauss_point.GetGaussWeight() * beam_jacobian;
    for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        force_element_2(i_dof) -= force(i_dir) * r_solid(i_dir).dx(i_dof + beam::n_dof_) *
                                  projected_gauss_point.GetGaussWeight() * beam_jacobian;
  }


  // Fill in the entries for the local matrices and vectors.
  {
    // Resize and initialize the return variables.
    if (forcevec1 != nullptr) forcevec1->Size(beam::n_dof_);
    if (forcevec2 != nullptr) forcevec2->Size(solid::n_dof_);
    if (stiffmat11 != nullptr) stiffmat11->Shape(beam::n_dof_, beam::n_dof_);
    if (stiffmat12 != nullptr) stiffmat12->Shape(beam::n_dof_, solid::n_dof_);
    if (stiffmat21 != nullptr) stiffmat21->Shape(solid::n_dof_, beam::n_dof_);
    if (stiffmat22 != nullptr) stiffmat22->Shape(solid::n_dof_, solid::n_dof_);

    if (forcevec1 != nullptr && forcevec2 != nullptr)
    {
      // $f_1$
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        (*forcevec1)(i_dof) = CORE::FADUTILS::CastToDouble(force_element_1(i_dof));
      // $f_2$
      for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
        (*forcevec2)(i_dof) = CORE::FADUTILS::CastToDouble(force_element_2(i_dof));
    }

    if (stiffmat11 != nullptr && stiffmat12 != nullptr && stiffmat21 != nullptr &&
        stiffmat22 != nullptr)
    {
      // $k_{11}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < beam::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < beam::n_dof_; i_dof_2++)
          (*stiffmat11)(i_dof_1, i_dof_2) =
              -CORE::FADUTILS::CastToDouble(force_element_1(i_dof_1).dx(i_dof_2));

      // $k_{12}, k_{21}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < beam::n_dof_; i_dof_1++)
      {
        for (unsigned int i_dof_2 = 0; i_dof_2 < solid::n_dof_; i_dof_2++)
        {
          (*stiffmat12)(i_dof_1, i_dof_2) =
              -CORE::FADUTILS::CastToDouble(force_element_1(i_dof_1).dx(beam::n_dof_ + i_dof_2));
          (*stiffmat21)(i_dof_2, i_dof_1) =
              -CORE::FADUTILS::CastToDouble(force_element_2(i_dof_2).dx(i_dof_1));
        }
      }

      // $k_{22}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < solid::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < solid::n_dof_; i_dof_2++)
          (*stiffmat22)(i_dof_1, i_dof_2) =
              -CORE::FADUTILS::CastToDouble(force_element_2(i_dof_1).dx(beam::n_dof_ + i_dof_2));
    }
  }

  // Return true as there are mesh tying contributions.
  return true;
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPair2D3DPlane<beam, solid>::GetTriadAtXiDouble(
    const double xi, CORE::LINALG::Matrix<3, 3, double>& triad, const bool reference) const
{
  if (reference)
  {
    CORE::LINALG::Matrix<beam::n_dof_, 1, double> beam_coupling_ref;
    CORE::LINALG::Matrix<solid::n_dof_, 1, double> dummy;
    this->GetCouplingReferencePosition(beam_coupling_ref, dummy);
    GEOMETRYPAIR::EvaluateTriadAtPlaneCurve<beam>(xi, beam_coupling_ref, triad, this->Element1());
  }
  else
  {
    GEOMETRYPAIR::EvaluateTriadAtPlaneCurve<beam>(
        xi, CORE::FADUTILS::CastToDouble(this->ele1pos_), triad, this->Element1());
  }
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidVolumeMeshtyingPair2D3DPlane<t_hermite, t_hex8>;
  template class BeamToSolidVolumeMeshtyingPair2D3DPlane<t_hermite, t_hex20>;
  template class BeamToSolidVolumeMeshtyingPair2D3DPlane<t_hermite, t_hex27>;
  template class BeamToSolidVolumeMeshtyingPair2D3DPlane<t_hermite, t_tet4>;
  template class BeamToSolidVolumeMeshtyingPair2D3DPlane<t_hermite, t_tet10>;
}  // namespace BEAMINTERACTION
