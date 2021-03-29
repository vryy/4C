/*----------------------------------------------------------------------*/
/*! \file

\brief Mesh tying for 2D rotational coupling examples with integration over the beam surface.

\level 3
*/


#include "beam_to_solid_volume_meshtying_pair_gauss_point_cross_section_rotation.H"

#include "../linalg/linalg_utils_densematrix_inverse.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "beam_contact_params.H"
#include "beam_to_solid_volume_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_element_functions.H"
#include "../drt_geometry_pair/geometry_pair_utility_classes.H"
#include "../drt_geometry_pair/geometry_pair_line_to_3D_evaluation_data.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume_gauss_point_projection_cross_section.H"



#include "beaminteraction_calc_utils.H"
#include "../drt_beam3/beam3r.H"
#include <Epetra_FEVector.h>


/**
 *
 */
template <typename beam, typename solid>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation<beam,
    solid>::BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation()
    : BeamToSolidVolumeMeshtyingPairBase<beam, solid>()
{
  // Empty constructor.
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation<beam,
    solid>::EvaluateAndAssemble(const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<Epetra_FEVector>& force_vector,
    const Teuchos::RCP<LINALG::SparseMatrix>& stiffness_matrix,
    const Teuchos::RCP<const Epetra_Vector>& displacement_vector)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(
        this->ele1posref_, this->ele2posref_, this->line_to_3D_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no segments, this pair has no contribution. Also there can be no more than one
  // segment.
  if (this->line_to_3D_segments_.size() == 0)
    return;
  else if (this->line_to_3D_segments_.size() > 1)
    dserror("There can be a maximum of one segment!");

  // Get the vector with the projection points for this pair.
  const std::vector<GEOMETRYPAIR::ProjectionPoint1DTo3D<double>>& projection_points =
      this->line_to_3D_segments_[0].GetProjectionPoints();

  // If there are no projection points, return no contact status.
  if (projection_points.size() == 0) return;

  // Set the FAD variables for the beam DOFs.
  LINALG::Matrix<beam::n_dof_, 1, scalar_type_rot_2nd> q_beam(true);
  for (unsigned int i_beam = 0; i_beam < beam::n_dof_; i_beam++)
    q_beam(i_beam) = FADUTILS::HigherOrderFadValue<scalar_type_rot_2nd>::apply(
        n_dof_pair_, i_beam, FADUTILS::CastToDouble(this->ele1pos_(i_beam)));

  // Check that the beam element is a SR beam.
  auto beam_ele = dynamic_cast<const DRT::ELEMENTS::Beam3r*>(this->Element1());
  if (beam_ele == nullptr)
    dserror("GetBeamTriadInterpolationScheme is only implemented for SR beams.");

  // Get the rotations of the beam rotation nodes.
  std::vector<double> beam_displacement_vector_full_double;
  BEAMINTERACTION::UTILS::GetCurrentElementDis(
      *discret, beam_ele, displacement_vector, beam_displacement_vector_full_double);
  unsigned int rot_dof_indices[n_dof_rot_] = {3, 12, 18};
  LINALG::Matrix<n_dof_rot_, 1, scalar_type_rot_2nd> q_rot(true);
  for (unsigned int i_rot = 0; i_rot < n_dof_rot_; i_rot++)
    q_rot(i_rot) = FADUTILS::HigherOrderFadValue<scalar_type_rot_2nd>::apply(n_dof_pair_,
        beam::n_dof_ + i_rot, beam_displacement_vector_full_double[rot_dof_indices[i_rot]]);

  {
    // Other rotations have to be 0.
    unsigned int other_rot_dof_indices[] = {4, 5, 13, 14, 19, 20};
    const double tol = 1e-10;
    double other_values = 0.0;
    for (unsigned int i_dim = 1; i_dim < 6; i_dim++)
      other_values += pow(beam_displacement_vector_full_double[other_rot_dof_indices[i_dim]], 2.0);
    if (sqrt(other_values) > tol)
      dserror(
          "The rotation values for y and z rotations habe to be 0. Other values: %f, tolerance: "
          "%f",
          FADUTILS::sqrt(other_values), tol);
  }

  // Set the FAD variables for the solid DOFs.
  LINALG::Matrix<solid::n_dof_, 1, scalar_type_rot_2nd> q_solid(true);
  for (unsigned int i_solid = 0; i_solid < solid::n_dof_; i_solid++)
    q_solid(i_solid) = FADUTILS::HigherOrderFadValue<scalar_type_rot_2nd>::apply(n_dof_pair_,
        n_dof_rot_ + beam::n_dof_ + i_solid, FADUTILS::CastToDouble(this->ele2pos_(i_solid)));

  // Initialize local matrices.
  LINALG::Matrix<n_dof_pair_, 1, double> local_force(true);
  LINALG::Matrix<n_dof_pair_, n_dof_pair_, double> local_stiff(true);

  // Initialize variables for position and force vectors.
  LINALG::Matrix<3, 1, double> dr_beam_ref;
  LINALG::Matrix<3, 1, scalar_type_rot_2nd> dr_beam;
  LINALG::Matrix<3, 3, scalar_type_rot_2nd> beam_triad(true);
  LINALG::Matrix<1, 1, scalar_type_rot_2nd> phi_beam;
  LINALG::Matrix<3, 1, scalar_type_rot_2nd> r_beam;
  LINALG::Matrix<3, 1, scalar_type_rot_2nd> r_solid;
  LINALG::Matrix<3, 1, scalar_type_rot_2nd> coupling_vector;
  scalar_type_rot_2nd penalty_potential = 0.0;

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

    // Reference tangent has to be in x direction.
    const double tol = 1e-10;
    double out_of_plane_values = 0.0;
    for (unsigned int i_dim = 1; i_dim < 3; i_dim++)
      out_of_plane_values += pow(dr_beam_ref(i_dim), 2.0);
    if (sqrt(out_of_plane_values) > tol)
      dserror(
          "The reference beam tangent has to be in x-direction. Out of plane value: %f, tolerance: "
          "%f",
          FADUTILS::sqrt(out_of_plane_values), tol);

    // Get the beam triad.
    GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
        projected_gauss_point.GetEta(), q_beam, dr_beam, this->Element1());
    dr_beam.Scale(1.0 / FADUTILS::VectorNorm(dr_beam));
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++) beam_triad(i_dim, 0) = dr_beam(i_dim);
    GEOMETRYPAIR::EvaluatePosition<GEOMETRYPAIR::t_line3>(
        projected_gauss_point.GetEta(), q_rot, phi_beam);
    beam_triad(1, 1) = cos(phi_beam(0));
    beam_triad(2, 1) = sin(phi_beam(0));
    beam_triad(1, 2) = -sin(phi_beam(0));
    beam_triad(2, 2) = cos(phi_beam(0));

    // Get the current positions on beam and solid as well as the beam rotation.
    GEOMETRYPAIR::EvaluatePositionLineCrossSection<beam>(projected_gauss_point.GetEta(),
        projected_gauss_point.GetEtaCrossSection(), q_beam, r_beam, this->Element1(), &beam_triad);
    GEOMETRYPAIR::EvaluatePosition<solid>(projected_gauss_point.GetXi(), q_solid, r_solid);

    // Calculate the force in this Gauss point. The sign of the force calculated here is the one
    // that acts on the beam.
    coupling_vector = r_solid;
    coupling_vector -= r_beam;

    // Add to penalty potential.
    penalty_potential += projected_gauss_point.GetGaussWeight() * beam_jacobian *
                         coupling_vector.Dot(coupling_vector) * penalty_parameter * 0.5;
  }

  // Get the GIDs of this pair.
  std::vector<int> lm_beam, lm_solid, lmowner, lmstride;
  this->Element1()->LocationVector(*discret, lm_beam, lmowner, lmstride);
  this->Element2()->LocationVector(*discret, lm_solid, lmowner, lmstride);

  int pos_dof_indices[] = {0, 1, 2, 6, 7, 8, 9, 10, 11, 15, 16, 17};
  LINALG::Matrix<n_dof_pair_, 1, int> gid_pair;
  for (unsigned int i = 0; i < beam::n_dof_; i++) gid_pair(i) = lm_beam[pos_dof_indices[i]];
  for (unsigned int i = 0; i < n_dof_rot_; i++)
    gid_pair(beam::n_dof_ + i) = lm_beam[rot_dof_indices[i]];
  for (unsigned int i = 0; i < solid::n_dof_; i++)
    gid_pair(beam::n_dof_ + n_dof_rot_ + i) = lm_solid[i];

  // If given, assemble force terms into the global force vector.
  if (force_vector != Teuchos::null)
  {
    LINALG::Matrix<n_dof_pair_, 1> force_vector_double;
    for (unsigned int i_dof = 0; i_dof < n_dof_pair_; i_dof++)
      force_vector_double(i_dof) = FADUTILS::CastToDouble(penalty_potential.dx(i_dof));
    force_vector->SumIntoGlobalValues(n_dof_pair_, gid_pair.A(), force_vector_double.A());
  }

  // If given, assemble force terms into the global stiffness matrix.
  if (stiffness_matrix != Teuchos::null)
    for (unsigned int i_dof = 0; i_dof < n_dof_pair_; i_dof++)
      for (unsigned int j_dof = 0; j_dof < n_dof_pair_; j_dof++)
        stiffness_matrix->FEAssemble(FADUTILS::CastToDouble(penalty_potential.dx(i_dof).dx(j_dof)),
            gid_pair(i_dof), gid_pair(j_dof));
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation<beam,
    solid>::CreateGeometryPair(const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataBase>&
        geometry_evaluation_data_ptr)
{
  // Call the method of the base class.
  BeamContactPair::CreateGeometryPair(geometry_evaluation_data_ptr);

  // Cast the geometry evaluation data to the correct format.
  auto line_to_3d_evaluation_data = Teuchos::rcp_dynamic_cast<GEOMETRYPAIR::LineTo3DEvaluationData>(
      geometry_evaluation_data_ptr, true);

  // Check that the cylinder strategy is given in the input file.
  INPAR::GEOMETRYPAIR::LineTo3DStrategy strategy = line_to_3d_evaluation_data->GetStrategy();
  if (strategy != INPAR::GEOMETRYPAIR::LineTo3DStrategy::gauss_point_projection_cross_section)
    dserror(
        "The cross section projection only works with cross section projection in the geometry "
        "pairs.");

  // Explicitly create the cylinder pair here, as this contact pair only works with this kind of
  // geometry pair.
  this->geometry_pair_ = Teuchos::rcp(
      new GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double, beam,
          solid>(line_to_3d_evaluation_data));
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation<beam,
    solid>::EvaluateBeamPosition(const GEOMETRYPAIR::ProjectionPoint1DTo3D<double>&
                                     integration_point,
    LINALG::Matrix<3, 1, scalar_type>& r_beam, bool reference) const
{
  if (reference)
    GEOMETRYPAIR::EvaluatePositionLineCrossSection<beam>(integration_point.GetEta(),
        integration_point.GetEtaCrossSection(), this->ele1posref_, r_beam, this->Element1());
  else
    GEOMETRYPAIR::EvaluatePositionLineCrossSection<beam>(integration_point.GetEta(),
        integration_point.GetEtaCrossSection(), this->ele1pos_, r_beam, this->Element1());
}


/**
 * Explicit template initialization of template class.
 */
namespace BEAMINTERACTION
{
  using namespace GEOMETRYPAIR;

  template class BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation<t_hermite, t_hex8>;
  template class BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation<t_hermite, t_hex20>;
  template class BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation<t_hermite, t_hex27>;
  template class BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation<t_hermite, t_tet4>;
  template class BeamToSolidVolumeMeshtyingPairGaussPointCrossSectionRotation<t_hermite, t_tet10>;
}  // namespace BEAMINTERACTION
