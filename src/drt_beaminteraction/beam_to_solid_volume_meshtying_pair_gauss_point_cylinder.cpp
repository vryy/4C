/*!
\file beam_to_solid_volume_meshtying_pair_gauss_point_cylinder.cpp

\brief Meshtying element for meshtying between a 3D beam and a 3D solid element using Gauss points
on the surface of the (circular) beam cross section.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_volume_meshtying_pair_gauss_point_cylinder.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"

#include "beam_contact_params.H"
#include "beam_to_solid_volume_meshtying_params.H"
#include "../drt_geometry_pair/geometry_pair_element_types.H"
#include "../drt_geometry_pair/geometry_pair_utility_functions.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume_gauss_point_projection_cylinder.H"
#include "../drt_geometry_pair/geometry_pair_evaluation_data_global.H"
#include "../drt_geometry_pair/geometry_pair_line_to_volume_evaluation_data.H"


/**
 *
 */
template <typename beam, typename solid>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCylinder<beam,
    solid>::BeamToSolidVolumeMeshtyingPairGaussPointCylinder()
    : BeamToSolidVolumeMeshtyingPairBase<beam, solid>()
{
  // Empty constructor.
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCylinder<beam, solid>::Init(
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams> params_ptr,
    const Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> geometry_evaluation_data_ptr,
    std::vector<DRT::Element const*> elements)
{
  // Check that the correct geometry pair is given.
  if (INPAR::GEOMETRYPAIR::LineToVolumeStrategy::gauss_point_projection_cylinder !=
      geometry_evaluation_data_ptr->LineToVolumeEvaluationData()->GetStrategy())
    dserror(
        "The class BeamToSolidVolumeMeshtyingPairGaussPointCylinder can only be used with the "
        "geometry pair GeometryPairLineToVolumeGaussPointProjectionCylinder set by the input "
        "parameter STRATEGY gauss_point_projection_cylinder");

  // Call Init of base class, the geometry pair will be created and initialized there.
  BeamToSolidVolumeMeshtyingPairBase<beam, solid>::Init(
      params_ptr, geometry_evaluation_data_ptr, elements);
}


/**
 *
 */
template <typename beam, typename solid>
void BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCylinder<beam, solid>::PreEvaluate()
{
  // Call PreEvaluate on the geometry Pair.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPairCylinder()->PreEvaluateCylinder(
        this->ele1posref_, this->ele2posref_, cylinder_to_volume_points_);
  }
}


/**
 *
 */
template <typename beam, typename solid>
bool BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCylinder<beam, solid>::Evaluate(
    LINALG::SerialDenseVector* forcevec1, LINALG::SerialDenseVector* forcevec2,
    LINALG::SerialDenseMatrix* stiffmat11, LINALG::SerialDenseMatrix* stiffmat12,
    LINALG::SerialDenseMatrix* stiffmat21, LINALG::SerialDenseMatrix* stiffmat22)
{
  // Call Evaluate on the geometry Pair. Only do this once for meshtying.
  if (!this->meshtying_is_evaluated_)
  {
    this->CastGeometryPairCylinder()->EvaluateCylinder(
        this->ele1posref_, this->ele2posref_, cylinder_to_volume_points_);
    this->meshtying_is_evaluated_ = true;
  }

  // If there are no projection points, return no contact status.
  if (this->cylinder_to_volume_points_.size() == 0) return false;

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
  for (unsigned int i_integration_point = 0;
       i_integration_point < this->cylinder_to_volume_points_.size(); i_integration_point++)
  {
    // Get the current Gauss point.
    const GEOMETRYPAIR::ProjectionPointVolumeToVolume<double>& projected_gauss_point =
        this->cylinder_to_volume_points_[i_integration_point];

    // Get the jacobian in the reference configuration.
    GEOMETRYPAIR::EvaluatePositionDerivative1<beam>(
        projected_gauss_point.GetXi1()(0), this->ele1posref_, dr_beam_ref, this->Element1());
    beam_jacobian = 0.5 * dr_beam_ref.Norm2();

    // Get the current positions on beam and solid.
    GEOMETRYPAIR::EvaluatePositionLineVolume<beam>(
        projected_gauss_point.GetXi1(), this->ele1pos_, r_beam, this->Element1());
    GEOMETRYPAIR::EvaluatePosition<solid>(projected_gauss_point.GetXi2(), this->ele2pos_, r_solid);

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
Teuchos::RCP<
    GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<double, beam, solid>>
BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCylinder<beam,
    solid>::CastGeometryPairCylinder() const
{
  return Teuchos::rcp_dynamic_cast<
      GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<double, beam, solid>>(
      this->geometry_pair_, true);
};


/**
 * Explicit template initialization of template class.
 */
// Hermite beam element, hex8 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCylinder<
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>;
// Hermite beam element, hex20 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCylinder<
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>;
// Hermite beam element, hex27 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCylinder<
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>;
// Hermite beam element, tet4 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCylinder<
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>;
// Hermite beam element, tet10 solid element.
template class BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCylinder<
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>;
