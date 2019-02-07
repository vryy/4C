/*!
\file geometry_pair_line_to_volume_segmentation.cpp

\brief Line to volume interaction with full segmentation of the line.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_volume_segmentation.H"
#include "geometry_pair_element_types.H"
#include "geometry_pair_evaluation_data_global.H"
#include "geometry_pair_line_to_volume_evaluation_data.H"
#include "geometry_pair_utility_classes.H"

#include "../headers/FAD_utils.H"


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<scalar_type, line, volume>::Setup()
{
  // Call Setup on the base class.
  GeometryPairLineToVolume<scalar_type, line, volume>::Setup();

  // Check if a segment tracker exists for this line element. If not a new one is created.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::set<LineSegment<double>>>& segment_tracker_map =
      this->EvaluationData()->LineToVolumeEvaluationData()->SegmentTrackerMutable();

  if (segment_tracker_map.find(line_element_id) == segment_tracker_map.end())
  {
    std::set<LineSegment<double>> new_tracking_set;
    new_tracking_set.clear();
    segment_tracker_map[line_element_id] = new_tracking_set;
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<scalar_type, line, volume>::Evaluate(
    const LINALG::TMatrix<scalar_type, line::n_dof_, 1>& q_line,
    const LINALG::TMatrix<scalar_type, volume::n_dof_, 1>& q_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Check if the element is initialized.
  this->CheckInitSetup();

  // Only zero segments are expected.
  if (segments.size() > 0)
    dserror("There should be zero segments for the segmentation method. The actual value is %d!",
        segments.size());

  // Number of search points.
  unsigned int n_search_points =
      this->EvaluationData()->LineToVolumeEvaluationData()->GetNumberOfSearchPoints();

  // Set up vector with projection points for the search points.
  std::vector<ProjectionPointLineToVolume<scalar_type>> search_points;
  search_points.reserve(n_search_points);
  LINALG::TMatrix<scalar_type, 3, 1> xi_start;
  this->SetStartValuesElement2(xi_start);
  scalar_type eta;
  for (unsigned int i_search_point = 0; i_search_point < n_search_points; i_search_point++)
  {
    eta = -1. + i_search_point * 2. / (n_search_points - 1.);
    search_points.push_back(ProjectionPointLineToVolume<scalar_type>(eta, xi_start));
  }

  // Project all the search points.
  unsigned int dummy;
  unsigned int n_projections;
  this->ProjectPointsOnLineToVolume(q_line, q_volume, search_points, dummy, n_projections);

  // If no point could be projected return, as we assume that we wont find a surface projection.
  // This usually happens for higher order volume elements, more search points can be a solution to
  // that problem.
  if (n_projections != 0)
  {
    // Now we start from every search point that could be projected (not necessary valid) and look
    // for surface intersections of the line with the volume. The intersection points are stored in
    // a std::set that will keep them unique and in order.
    std::set<ProjectionPointLineToVolume<scalar_type>> intersection_points;

    // If the start and/or end points (nodes of the line element) could be projected valid, add them
    // to the intersections set. This helps for some cases where the line just touches the surface
    // of the volume.
    if (search_points.front().GetProjectionResult() == ProjectionResult::projection_found_valid)
      intersection_points.insert(search_points.front());
    if (search_points.back().GetProjectionResult() == ProjectionResult::projection_found_valid)
      intersection_points.insert(search_points.back());

    // Vector for intersection point search.
    std::vector<ProjectionPointLineToVolume<scalar_type>> search_intersection_points;

    // Starting from each search point, try to project to all surfaces of the volume.
    for (auto const& point : search_points)
    {
      // Only use search points that could be projected.
      if (point.GetProjectionResult() != ProjectionResult::projection_not_found)
      {
        // Get the intersections with the volume.
        this->IntersectLineWithVolume(
            q_line, q_volume, search_intersection_points, point.GetEta(), point.GetXi());

        // Add the found intersection points to the set.
        for (auto& found_point : search_intersection_points)
          intersection_points.insert(found_point);
      }
    }

    // In the case of zero and one intersection points, no segmentation is needed. One point only
    // occurs when the line just touches the volume. This is the reason why the start and/or end
    // points are added to the set of found intersection.
    if (intersection_points.size() > 1)
    {
      // The intersection points in intersection_points are in order. Now it is checked if a
      // point between the intersection points is inside or outside the volume and therefore it can
      // be decided if the segment is part of this pair.. By doing so we can avoid complications in
      // cases when a line is exactly between two volumes.

      // We reuse the vector created for the search points.
      search_points.clear();
      search_points.reserve(intersection_points.size() - 1);

      // Loop over the middle points and check if they are projected valid or not.
      std::set<GEOMETRYPAIR::LineSegment<double>>& segment_tracker = GetSegmentTrackingSetMutable();
      ProjectionResult projection_result;
      bool last_segment_active = false;
      scalar_type eta_start;
      unsigned int counter = 0;
      for (typename std::set<ProjectionPointLineToVolume<scalar_type>>::iterator set_iterator =
               intersection_points.begin();
           set_iterator != intersection_points.end(); ++set_iterator)
      {
        // Reference to this current point.
        const ProjectionPointLineToVolume<scalar_type>& start_point = *set_iterator;

        // Get the next intersection point and calculate the projection. This can only be done if
        // the iterator is not on its last iteration.
        if (counter != intersection_points.size() - 1)
        {
          // Reference to the next point.
          const ProjectionPointLineToVolume<scalar_type>& end_point = *std::next(set_iterator);

          // Get starting points for projection.
          eta = 0.5 * (start_point.GetEta() + end_point.GetEta());
          xi_start = start_point.GetXi();
          xi_start += end_point.GetXi();
          xi_start.Scale(0.5);

          // Project and check result.
          this->ProjectPointOnLineToVolume(q_line, q_volume, eta, xi_start, projection_result);
        }
        else
        {
          // In case the iterator is on the last intersection point, finish up the active segment.
          projection_result = ProjectionResult::projection_found_not_valid;
        }

        // Check if a new segment is found.
        if (projection_result == ProjectionResult::projection_found_valid)
        {
          // The line part between the current start_point and end_point is a segment for this pair.
          // Start its counter if no segment is currently active.
          if (!last_segment_active)
          {
            last_segment_active = true;
            eta_start = start_point.GetEta();
          }
        }
        else
        {
          // This point is outside of the volume -> if current segment exists finish it.
          if (last_segment_active)
          {
            // Create a segment with double as the scalar type.
            LineSegment<scalar_type> new_segment_double(
                FADUTILS::CastToDouble(eta_start), FADUTILS::CastToDouble(start_point.GetEta()));

            // Check if the segment already exists for this line.
            if (segment_tracker.find(new_segment_double) == segment_tracker.end())
            {
              // Add the new segment to this pair and to the evaluation tracker.
              segments.push_back(LineSegment<scalar_type>(eta_start, start_point.GetEta()));
              segment_tracker.insert(new_segment_double);

              // Project the Gauss points on the segment. All points have to project valid.
              this->ProjectGaussPointsOnSegmentToVolume(q_line, q_volume,
                  this->EvaluationData()->LineToVolumeEvaluationData()->GetGaussPoints(),
                  segments.back());
            }

            // Deactivate the current segment.
            last_segment_active = false;
          }
        }

        // Advance the counter.
        counter++;
      }
    }
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
std::set<GEOMETRYPAIR::LineSegment<double>>& GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<
    scalar_type, line, volume>::GetSegmentTrackingSetMutable() const
{
  // Get the segment tracker for this line element.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::set<GEOMETRYPAIR::LineSegment<double>>>& segment_tracker_map =
      this->EvaluationData()->LineToVolumeEvaluationData()->SegmentTrackerMutable();
  return segment_tracker_map[line_element_id];
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4>;
template class GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10>;
