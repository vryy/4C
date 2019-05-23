/*!

\brief Meshtying element for meshtying between a beam and a 3D solid element using Gauss points
on the surface of the (circular) beam cross section.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "beam_contact_params.H"
#include "beam_to_solid_volume_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_element_types.H"
#include "../drt_geometry_pair/geometry_pair_utility_functions.H"
#include "../drt_geometry_pair/geometry_pair_evaluation_data_global.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume_evaluation_data.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume_gauss_point_projection_cross_section.H"
#include "beam_to_solid_volume_meshtying_pair_gauss_point_cross_section.H"


/**
 *
 */
template <typename beam, typename solid>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<beam,
    solid>::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection()
    : BeamToSolidVolumeMeshtyingPairBase<beam, solid>()
{
  // Empty constructor.
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<beam, solid>::Init(
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr,
    std::vector<DRT::Element const*> elements)
{
  // Check that the correct geometry pair is given.
  if (INPAR::GEOMETRYPAIR::LineToVolumeStrategy::gauss_point_projection_cross_section !=
      geometry_evaluation_data_ptr->LineToVolumeEvaluationData()->GetStrategy())
    dserror(
        "The class BeamToSolidVolumeMeshtyingPairGaussPointCylinder can only be used with the "
        "geometry pair GeometryPairLineToVolumeGaussPointProjectionCylinder set by the input "
        "parameter STRATEGY gauss_point_projection_cross_section");

  // Call Init of base class, the geometry pair will be created and initialized there.
  BeamToSolidVolumeMeshtyingPairBase<beam, solid>::Init(
      params_ptr, geometry_evaluation_data_ptr, elements);
}


/**
 *
 */
template <typename beam, typename solid>
bool BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<beam, solid>::Evaluate(
    LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPair()->Evaluate(
        this->ele1posref_, this->ele2posref_, this->line_to_volume_segments_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no segments, this pair has no contribution. Also there can be no more than one
  // segment.
  if (this->line_to_volume_segments_.size() == 0)
    return false;
  else if (this->line_to_volume_segments_.size() > 1)
    dserror("There can be a maximum of one segment!");

  // Get the vector with the projection points for this pair.
  const std::vector<GEOMETRYPAIR::ProjectionPointLineToVolume<double>>& projection_points =
      this->line_to_volume_segments_[0].GetProjectionPoints();

  // If there are no projection points, return no contact status.
  if (projection_points.size() == 0) return false;

  // Initialize variables for position and force vectors.
  LINALG::TMatrix<double, 3, 1> dr_beam_ref;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> r_beam;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> r_solid;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1> force;
  LINALG::TMatrix<TYPE_BTS_VMT_AD, beam::n_dof_, 1> force_element_1(true);
  LINALG::TMatrix<TYPE_BTS_VMT_AD, solid::n_dof_, 1> force_element_2(true);

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
    const GEOMETRYPAIR::ProjectionPointLineToVolume<double>& projected_gauss_point =
        projection_points[i_integration_point];

    // Get the jacobian in the reference configuration.
    GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
        projected_gauss_point.GetEta(), this->ele1posref_, dr_beam_ref, this->Element1());
    beam_jacobian = 0.5 * dr_beam_ref.Norm2();

    // Get the current positions on beam and solid.
    GEOMETRYPAIR::EvaluatePositionLineCrossSection<beam>(projected_gauss_point.GetEta(),
        projected_gauss_point.GetEtaCrossSection(), this->ele1pos_, r_beam, this->Element1());
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
    if (forcevec1 != NULL) forcevec1->Size(beam::n_dof_);
    if (forcevec2 != NULL) forcevec2->Size(solid::n_dof_);
    if (stiffmat11 != NULL) stiffmat11->Shape(beam::n_dof_, beam::n_dof_);
    if (stiffmat12 != NULL) stiffmat12->Shape(beam::n_dof_, solid::n_dof_);
    if (stiffmat21 != NULL) stiffmat21->Shape(solid::n_dof_, beam::n_dof_);
    if (stiffmat22 != NULL) stiffmat22->Shape(solid::n_dof_, solid::n_dof_);

    if (forcevec1 != NULL && forcevec2 != NULL)
    {
      // $f_1$
      for (unsigned int i_dof = 0; i_dof < beam::n_dof_; i_dof++)
        (*forcevec1)(i_dof) = force_element_1(i_dof).val();
      // $f_2$
      for (unsigned int i_dof = 0; i_dof < solid::n_dof_; i_dof++)
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
        for (unsigned int i_dof_2 = 0; i_dof_2 < solid::n_dof_; i_dof_2++)
        {
          (*stiffmat12)(i_dof_1, i_dof_2) = -force_element_1(i_dof_1).dx(beam::n_dof_ + i_dof_2);
          (*stiffmat21)(i_dof_2, i_dof_1) = -force_element_2(i_dof_2).dx(i_dof_1);
        }
      }

      // $k_{22}$
      for (unsigned int i_dof_1 = 0; i_dof_1 < solid::n_dof_; i_dof_1++)
        for (unsigned int i_dof_2 = 0; i_dof_2 < solid::n_dof_; i_dof_2++)
          (*stiffmat22)(i_dof_1, i_dof_2) = -force_element_2(i_dof_1).dx(beam::n_dof_ + i_dof_2);
    }
  }

  // Return true as there are meshtying contributions.
  return true;
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<beam,
    solid>::CreateGeometryPair(const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal>
        geometry_evaluation_data_ptr)
{
  // Check that the cylinder strategy is given in the input file.
  INPAR::GEOMETRYPAIR::LineToVolumeStrategy strategy =
      geometry_evaluation_data_ptr->LineToVolumeEvaluationData()->GetStrategy();
  if (strategy != INPAR::GEOMETRYPAIR::LineToVolumeStrategy::gauss_point_projection_cross_section)
    dserror(
        "The cross section projection only works with cross section projection in the geometry "
        "pairs.");

  // Explicitly create the cylinder pair here, as this contact pair only works with this kind of
  // geometry pair.
  this->geometry_pair_ = Teuchos::rcp(
      new GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double, beam,
          solid>());
}

/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<beam,
    solid>::EvaluateBeamPosition(const GEOMETRYPAIR::ProjectionPointLineToVolume<double>&
                                     integration_point,
    LINALG::TMatrix<TYPE_BTS_VMT_AD, 3, 1>& r_beam, bool reference) const
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
// Hermite beam element, hex8 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>;
// Hermite beam element, hex20 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>;
// Hermite beam element, hex27 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>;
// Hermite beam element, tet4 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>;
// Hermite beam element, tet10 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>;
