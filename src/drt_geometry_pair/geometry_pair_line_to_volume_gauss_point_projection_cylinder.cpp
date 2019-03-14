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

#define n_gauss_points_cylinder_axis 5
#define n_gauss_points_cylinder_circ 5

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
    int n_gauss_points = n_gauss_points_cylinder_axis * n_gauss_points_cylinder_circ;
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
    volume>::PreEvaluate(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
    const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Check if the element is initialized.
  this->CheckInitSetup();

  // Get the Gauss point projection tracker for this line element.
  std::vector<bool>& line_projection_tracker = GetLineProjectionVectorMutable();

  // Gauss rule.
  DRT::UTILS::IntegrationPoints1D gauss_points_axis =
      DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GaussRule1D::intrule_line_5point);
  DRT::UTILS::IntegrationPoints1D gauss_points_circumfence =
      DRT::UTILS::IntegrationPoints1D(DRT::UTILS::GaussRule1D::intrule_line_5point);

  // Initilaize variables for the projection.
  scalar_type eta;
  LINALG::TMatrix<scalar_type, 3, 1> xi;
  ProjectionResult projection_result;
  LineSegment<scalar_type> line_segment;
  bool one_projects = false;

  double angle, eta_axis, eta_circ;
  double radius = 0.1;
  LINALG::TMatrix<scalar_type, 3, 1> pos_on_cylinder;

  // Loop over Gauss points and check if they project to this volume.
  for (unsigned int index_gp_axis = 0; index_gp_axis < n_gauss_points_cylinder_axis;
       index_gp_axis++)
  {
    for (unsigned int index_gp_circ = 0; index_gp_circ < n_gauss_points_cylinder_circ;
         index_gp_circ++)
    {
      // Only check points that do not already have a valid projection.
      if (line_projection_tracker[index_gp_axis * index_gp_axis + index_gp_circ] == false)
      {
        eta_axis = gauss_points_axis.qxg[index_gp_axis][0];
        eta_circ = gauss_points_circumfence.qxg[index_gp_circ][0];

        // Get the point on the circumfence of the beam.
        // This only works if the tangent of the beam in the reference configuration is [1, 0, 0].
        angle = 0.5 * (eta_circ + 1.) * 2. * 3.14159265358979323846;
        EvaluatePosition<line>(eta_axis, q_line, pos_on_cylinder, this->Element1());
        pos_on_cylinder(1) += cos(angle) * radius;
        pos_on_cylinder(2) += sin(angle) * radius;
        this->ProjectPointToVolume(pos_on_cylinder, q_volume, xi, projection_result);

        if (projection_result == ProjectionResult::projection_found_valid)
        {
          // Valid Gauss point was found, add to this segment and set tracking point to true.
          line_segment.AddProjectionPoint(
              ProjectionPointLineToVolume<scalar_type>(eta_axis, xi, 1.));
          line_projection_tracker[index_gp_axis * index_gp_axis + index_gp_circ] = true;

          one_projects = true;
        }
      }
    }
  }

  if (one_projects)
  {
    // Clear the segment vector and add the found segment for the current line to volume pair.
    segments.clear();
    segments.push_back(line_segment);
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCylinder<scalar_type, line,
    volume>::Evaluate(const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
    const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Check if the element is initialized.
  this->CheckInitSetup();

  // Only zero one segments are expected.
  if (segments.size() > 1)
    dserror(
        "There should be zero or one segments for the Gauss point projection method. The actual "
        "value is %d!",
        segments.size());

  // Check if one point projected in PreEvaluate.
  if (segments.size() == 1 && segments[0].GetNumberOfProjectionPoints() > 0)
  {
    // Flag if segmentation is needed.
    bool need_segmentation = false;

    // Check if all Gauss points projected for this line.
    const std::vector<bool>& line_projection_tracker = GetLineProjectionVectorMutable();
    for (auto const& projects : line_projection_tracker)
      if (!projects) need_segmentation = true;

    if (need_segmentation)
    {
      // Segmentation is needed. First get the intersection points with the volume.
      std::vector<ProjectionPointLineToVolume<scalar_type>> intersection_points;
      this->IntersectLineWithVolume(q_line, q_volume, intersection_points);

      // This algorithm only works if one intersection point was found.
      if (intersection_points.size() != 1)
        dserror("In the segmentation case we expect exactly one found intersection point. Got: %d!",
            intersection_points.size());

      // Get the limits of the segmented line.
      scalar_type eta_a, eta_b, eta_intersection_point, eta_first_gauss_point;
      eta_intersection_point = intersection_points[0].GetEta();
      eta_first_gauss_point = segments[0].GetProjectionPoints()[0].GetEta();
      if (eta_intersection_point < eta_first_gauss_point)
      {
        eta_a = eta_intersection_point;
        eta_b = 1.;
      }
      else
      {
        eta_a = -1.;
        eta_b = eta_intersection_point;
      }

      // Reproject the Gauss points on the segmented line.
      segments[0] = LineSegment<scalar_type>(eta_a, eta_b);
      this->ProjectGaussPointsOnSegmentToVolume(q_line, q_volume,
          this->EvaluationData()->LineToVolumeEvaluationData()->GetGaussPoints(), segments[0]);
    }
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
