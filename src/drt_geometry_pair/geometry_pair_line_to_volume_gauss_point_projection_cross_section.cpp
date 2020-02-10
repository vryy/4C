/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with Gauss point projection on the cylinder surface along the
line.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_volume_gauss_point_projection_cross_section.H"

#include "geometry_pair_element_functions.H"
#include "geometry_pair_line_to_3D_evaluation_data.H"
#include "geometry_pair_utility_classes.H"

#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_integration.H"

#include <math.h>


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<scalar_type, line,
    volume>::Setup()
{
  // Call Setup on the base class.
  GeometryPairLineToVolume<scalar_type, line, volume>::Setup();

  // Check if a projection tracking vector exists for this line element. If not a new one is
  // created.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->line_to_3d_evaluation_data_->GetGaussPointProjectionTrackerMutable();

  if (projection_tracker.find(line_element_id) == projection_tracker.end())
  {
    int n_gauss_points =
        this->line_to_3d_evaluation_data_->GetNumberOfGaussPoints() *
        this->line_to_3d_evaluation_data_->GetNumberOfIntegrationPointsCircumfence();
    std::vector<bool> new_tracking_vector;
    new_tracking_vector.resize(n_gauss_points, false);
    projection_tracker[line_element_id] = new_tracking_vector;
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<scalar_type, line,
    volume>::PreEvaluate(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Check if the element is initialized.
  this->CheckInitSetup();

  // Get the Gauss point projection tracker for this line element.
  std::vector<bool>& line_projection_tracker = GetLineProjectionVectorMutable();

  // Gauss rule.
  DRT::UTILS::IntegrationPoints1D gauss_points_axis =
      this->line_to_3d_evaluation_data_->GetGaussPoints();
  unsigned int n_gauss_points_axis = this->line_to_3d_evaluation_data_->GetNumberOfGaussPoints();
  unsigned int n_integration_points_circ =
      this->line_to_3d_evaluation_data_->GetNumberOfIntegrationPointsCircumfence();

  // Initilaize variables for the projection.
  scalar_type eta;
  double alpha;
  LINALG::Matrix<3, 1, scalar_type> r_line_cross_section;
  LINALG::Matrix<2, 1, scalar_type> eta_cross_section;
  LINALG::Matrix<3, 1, scalar_type> xi_solid;
  ProjectionResult projection_result;
  segments.clear();
  bool one_projects = false;
  LineSegment<scalar_type> projection_point_segment;

  // Loop over Gauss points and check if they project to this volume.
  for (unsigned int index_gp_axis = 0; index_gp_axis < n_gauss_points_axis; index_gp_axis++)
  {
    for (unsigned int index_gp_circ = 0; index_gp_circ < n_integration_points_circ; index_gp_circ++)
    {
      // Index of the current Gauss point in the tracking vector.
      unsigned int index_gp = index_gp_axis * n_integration_points_circ + index_gp_circ;

      // Only check points that do not already have a valid projection.
      if (line_projection_tracker[index_gp] == false)
      {
        // Centerline coordinate and coodrinates in the cross section.
        eta = gauss_points_axis.qxg[index_gp_axis][0];
        alpha = 2. * M_PI / double(n_integration_points_circ) * index_gp_circ;
        eta_cross_section(0) = cos(alpha);
        eta_cross_section(1) = sin(alpha);

        // Position of point in the cross section.
        GEOMETRYPAIR::EvaluatePositionLineCrossSection<line>(
            eta, eta_cross_section, q_line, r_line_cross_section, this->Element1());

        // Project point to the volume.
        this->ProjectPointToOther(r_line_cross_section, q_volume, xi_solid, projection_result);
        if (projection_result == ProjectionResult::projection_found_valid)
        {
          // Valid Gauss point was found, add to this segment and set tracking point to true.
          ProjectionPoint1DTo3D<scalar_type> new_point(eta, xi_solid,
              gauss_points_axis.qwgt[index_gp_axis] * 2. / double(n_integration_points_circ));
          new_point.SetEtaCrossSection(eta_cross_section);
          projection_point_segment.AddProjectionPoint(new_point);
          line_projection_tracker[index_gp] = true;

          one_projects = true;
        }
      }
    }
  }

  if (one_projects)
  {
    // Clear the segment vector and add the found segment for the current line to volume pair.
    segments.clear();
    segments.push_back(projection_point_segment);
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<scalar_type, line,
    volume>::Evaluate(const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Check if the element is initialized.
  this->CheckInitSetup();

  // Only zero one segments are expected.
  if (segments.size() > 1)
    dserror(
        "There should be zero or one segments for the Gauss point cylinder projection method. The "
        "actual value is %d!",
        segments.size());

  // Check if one point projected in PreEvaluate.
  if (segments.size() == 1 && segments[0].GetNumberOfProjectionPoints() > 0)
  {
    // Check if all points of this line projected.
    const std::vector<bool>& projection_vector = GetLineProjectionVectorMutable();
    bool all_projected =
        std::all_of(projection_vector.begin(), projection_vector.end(), [](bool v) { return v; });
    if (!all_projected)
    {
      unsigned int valid_projection_points = 0;
      for (auto const& value : projection_vector)
        if (value) valid_projection_points += 1;
      dserror(
          "The cross section projection currently only works if all points on a line project! Of "
          "the %d points, only %d projected.",
          projection_vector.size(), valid_projection_points);
    }
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
std::vector<bool>& GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<
    scalar_type, line, volume>::GetLineProjectionVectorMutable() const
{
  // Get the Gauss point projection tracker for this line element.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->line_to_3d_evaluation_data_->GetGaussPointProjectionTrackerMutable();
  return projection_tracker[line_element_id];
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<double,
    GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>;
