/*!incomplete
\file geometry_pair_line_to_volume_gauss_point_projection_cylinder.cpp

\brief Line to volume interaction with simple Gauss point projection and boundary segmentation.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_volume_gauss_point_projection_cylinder.H"
#include "geometry_pair_element_types.H"
#include "geometry_pair_evaluation_data_global.H"
#include "geometry_pair_line_to_volume_evaluation_data.H"
#include "geometry_pair_utility_classes.H"

#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include <math.h>


#define radius 0.1


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<scalar_type, line,
    volume>::Setup()
{
  // Call Setup on the base class.
  GeometryPairLineToVolume<scalar_type, line, volume>::Setup();

  // Check if a projection tracking vector exists for this line element. If not a new one is
  // created.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->EvaluationData()->LineToVolumeEvaluationData()->GetGaussPointProjectionTrackerMutable();

  if (projection_tracker.find(line_element_id) == projection_tracker.end())
  {
    int n_gauss_points =
        this->EvaluationData()->LineToVolumeEvaluationData()->GetNumberOfGaussPoints() *
        this->EvaluationData()->LineToVolumeEvaluationData()->GetGaussPointsCircumfence();
    std::vector<bool> new_tracking_vector;
    new_tracking_vector.resize(n_gauss_points, false);
    projection_tracker[line_element_id] = new_tracking_vector;
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<scalar_type, line,
    volume>::PreEvaluateCylinder(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
    const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
    std::vector<GEOMETRYPAIR::ProjectionPointVolumeToVolume<scalar_type>>&
        cylinder_to_volume_points) const
{
  // Check if the element is initialized.
  this->CheckInitSetup();

  // Get the Gauss point projection tracker for this line element.
  std::vector<bool>& line_projection_tracker = GetLineProjectionVectorMutable();

  // Gauss rule.
  DRT::UTILS::IntegrationPoints1D gauss_points_axis =
      this->EvaluationData()->LineToVolumeEvaluationData()->GetGaussPoints();
  unsigned int n_gauss_points_axis =
      this->EvaluationData()->LineToVolumeEvaluationData()->GetNumberOfGaussPoints();
  unsigned int n_gauss_points_circ =
      this->EvaluationData()->LineToVolumeEvaluationData()->GetGaussPointsCircumfence();

  // Initilaize variables for the projection.
  scalar_type eta;
  double alpha;
  LINALG::TMatrix<scalar_type, 3, 1> r_beam;
  LINALG::TMatrix<scalar_type, 3, 1> xi_beam;
  LINALG::TMatrix<scalar_type, 3, 1> xi_solid;
  ProjectionResult projection_result;
  cylinder_to_volume_points.clear();

  // Loop over Gauss points and check if they project to this volume.
  for (unsigned int index_gp_axis = 0; index_gp_axis < n_gauss_points_axis; index_gp_axis++)
  {
    for (unsigned int index_gp_circ = 0; index_gp_circ < n_gauss_points_circ; index_gp_circ++)
    {
      // Index of the current Gauss point in the tracking vector.
      unsigned int index_gp = index_gp_axis * n_gauss_points_circ + index_gp_circ;

      // Only check points that do not already have a valid projection.
      if (line_projection_tracker[index_gp] == false)
      {
        // Centerline coordinate.
        eta = gauss_points_axis.qxg[index_gp_axis][0];

        // Get the spatial position of the beam centerline.
        GEOMETRYPAIR::EvaluatePosition<line>(eta, q_line, r_beam, this->Element1());

        // Add the in crossection position.
        alpha = 2. * M_PI / double(n_gauss_points_circ) * index_gp_circ;
        r_beam(1) += radius * cos(alpha);
        r_beam(2) += radius * sin(alpha);

        // Parameter coordinates on the line.
        xi_beam(0) = eta;
        xi_beam(1) = cos(alpha);
        xi_beam(2) = sin(alpha);

        // Project point to the volume.
        this->ProjectPointToVolume(r_beam, q_volume, xi_solid, projection_result);
        if (projection_result == ProjectionResult::projection_found_valid)
        {
          // Valid Gauss point was found, add to this segment and set tracking point to true.
          cylinder_to_volume_points.push_back(ProjectionPointVolumeToVolume<scalar_type>(xi_beam,
              xi_solid, gauss_points_axis.qwgt[index_gp_axis] * 2. / double(n_gauss_points_circ)));
          line_projection_tracker[index_gp] = true;
        }
      }
    }
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<scalar_type, line,
    volume>::EvaluateCylinder(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
    const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
    std::vector<GEOMETRYPAIR::ProjectionPointVolumeToVolume<scalar_type>>&
        cylinder_to_volume_points) const
{
  // Check if the element is initialized.
  this->CheckInitSetup();

  // Check if one point projected in PreEvaluate.
  if (cylinder_to_volume_points.size() > 1)
  {
    // Check if all points of this beam projected.
    const std::vector<bool>& projection_vector = GetLineProjectionVectorMutable();
    bool all_projected =
        std::all_of(projection_vector.begin(), projection_vector.end(), [](bool v) { return v; });
    if (!all_projected)
      dserror("The cylinder projection currently only works if all points on a beam projected!");
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
std::vector<bool>& GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<scalar_type,
    line, volume>::GetLineProjectionVectorMutable() const
{
  // Get the Gauss point projection tracker for this line element.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->EvaluationData()->LineToVolumeEvaluationData()->GetGaussPointProjectionTrackerMutable();
  return projection_tracker[line_element_id];
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>;
