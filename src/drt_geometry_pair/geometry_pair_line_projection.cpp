/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and other geometry types.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_projection.H"

#include "geometry_pair_element.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_line_to_volume_gauss_point_projection.H"
#include "geometry_pair_line_to_volume_segmentation.H"
#include "geometry_pair_line_to_3D_evaluation_data.H"


/**
 *
 */
template <typename pair_type>
void GEOMETRYPAIR::LineTo3DBase<pair_type>::ProjectPointOnLineToOther(const pair_type* pair,
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other, const scalar_type& eta,
    LINALG::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result)
{
  // Get the point on the line.
  LINALG::Matrix<3, 1, scalar_type> r_line;
  GEOMETRYPAIR::EvaluatePosition<line>(eta, q_line, r_line, pair->Element1());

  // Project the point to the solid.
  pair->ProjectPointToOther(r_line, q_other, xi, projection_result);
}

/**
 *
 */
template <typename pair_type>
void GEOMETRYPAIR::LineTo3DBase<pair_type>::ProjectPointsOnLineToOther(const pair_type* pair,
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& projection_points,
    unsigned int& n_projections_valid, unsigned int& n_projections)
{
  // Initialize counters.
  n_projections_valid = 0;
  n_projections = 0;

  // Loop over points and check if they project to this other geomety.
  for (auto& point : projection_points)
  {
    // Project the point.
    ProjectPointOnLineToOther(pair, q_line, q_other, point.GetEta(), point.GetXiMutable(),
        point.GetProjectionResultMutable());

    // Update the counters.
    if (point.GetProjectionResult() == ProjectionResult::projection_found_valid)
    {
      n_projections_valid++;
      n_projections++;
    }
    if (point.GetProjectionResult() == ProjectionResult::projection_found_not_valid)
      n_projections++;
  }
}

/**
 *
 */
template <typename pair_type>
void GEOMETRYPAIR::LineTo3DBase<pair_type>::ProjectPointsOnLineToOther(const pair_type* pair,
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& projection_points,
    unsigned int& n_projections_valid)
{
  // Initialize dummy variable.
  unsigned int n_projections_dummy;

  // Project the points.
  ProjectPointsOnLineToOther(
      pair, q_line, q_other, projection_points, n_projections_valid, n_projections_dummy);
}

/**
 *
 */
template <typename pair_type>
void GEOMETRYPAIR::LineTo3DBase<pair_type>::ProjectGaussPointsOnSegmentToOther(
    const pair_type* pair, const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
    const DRT::UTILS::IntegrationPoints1D& gauss_points, LineSegment<scalar_type>& segment)
{
  // Set up the vector with the projection points.
  std::vector<ProjectionPoint1DTo3D<scalar_type>>& projection_points =
      segment.GetProjectionPointsMutable();
  projection_points.clear();
  projection_points.reserve(gauss_points.nquad);
  LINALG::Matrix<3, 1, scalar_type> xi_start;
  StartValues<other::geometry_type_>::Set(xi_start);
  for (unsigned int i = 0; i < (unsigned int)gauss_points.nquad; i++)
  {
    scalar_type eta = segment.GetEtaA() +
                      (segment.GetEtaB() - segment.GetEtaA()) * 0.5 * (gauss_points.qxg[i][0] + 1.);
    projection_points.push_back(
        ProjectionPoint1DTo3D<scalar_type>(eta, xi_start, gauss_points.qwgt[i]));
  }

  // Project the Gauss points to the other geometry.
  unsigned int n_valid_projections;
  ProjectPointsOnLineToOther(pair, q_line, q_other, projection_points, n_valid_projections);

  // Check if all points could be projected.
  if (n_valid_projections != (unsigned int)gauss_points.nquad)
    dserror(
        "All Gauss points need to have a valid projection. The number of Gauss points is %d, but "
        "the number of valid projections is %d!",
        gauss_points.nquad, n_valid_projections);
}


/**
 *
 */
template <typename pair_type>
void GEOMETRYPAIR::LineTo3DBase<pair_type>::IntersectLineWithOther(const pair_type* pair,
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& intersection_points)
{
  // Set default values for the parameter coordinates.
  scalar_type eta_start;
  LINALG::Matrix<3, 1, scalar_type> xi_start;
  StartValues<line::geometry_type_>::Set(eta_start);
  StartValues<other::geometry_type_>::Set(xi_start);

  // Call the intersect function.
  pair->IntersectLineWithOther(q_line, q_other, intersection_points, eta_start, xi_start);
};


/**
 *
 */
template <typename pair_type>
void GEOMETRYPAIR::LineTo3DGaussPointProjection<pair_type>::PreEvaluate(const pair_type* pair,
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
    std::vector<LineSegment<scalar_type>>& segments)
{
  // Get the Gauss point projection tracker for this line element.
  std::vector<bool>& line_projection_tracker = GetLineProjectionVectorMutable(pair);

  // Gauss rule.
  DRT::UTILS::IntegrationPoints1D gauss_points = pair->GetEvaluationData()->GetGaussPoints();

  // Initilaize variables for the projection.
  scalar_type eta;
  LINALG::Matrix<3, 1, scalar_type> xi;
  ProjectionResult projection_result;
  LineSegment<scalar_type> line_segment;
  bool one_projects = false;

  // Loop over Gauss points and check if they project to this other geometry.
  for (unsigned int index_gp = 0; index_gp < line_projection_tracker.size(); index_gp++)
  {
    // Only check points that do not already have a valid projection.
    if (line_projection_tracker[index_gp] == false)
    {
      eta = gauss_points.qxg[index_gp][0];
      LineTo3DBase<pair_type>::ProjectPointOnLineToOther(
          pair, q_line, q_other, eta, xi, projection_result);
      if (projection_result == ProjectionResult::projection_found_valid)
      {
        // Valid Gauss point was found, add to this segment and set tracking point to true.
        line_segment.AddProjectionPoint(
            ProjectionPoint1DTo3D<scalar_type>(eta, xi, gauss_points.qwgt[index_gp]));
        line_projection_tracker[index_gp] = true;

        one_projects = true;
      }
    }
  }

  if (one_projects)
  {
    // Clear the segment vector and add the found segment for the current line-to-xxx pair.
    segments.clear();
    segments.push_back(line_segment);
  }
}

/**
 *
 */
template <typename pair_type>
void GEOMETRYPAIR::LineTo3DGaussPointProjection<pair_type>::Evaluate(const pair_type* pair,
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
    std::vector<LineSegment<scalar_type>>& segments)
{
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
    const std::vector<bool>& line_projection_tracker = GetLineProjectionVectorMutable(pair);
    for (auto const& projects : line_projection_tracker)
      if (!projects) need_segmentation = true;

    if (need_segmentation)
    {
      // Segmentation is needed. First get the intersection points with the other geometry.
      std::vector<ProjectionPoint1DTo3D<scalar_type>> intersection_points;
      LineTo3DBase<pair_type>::IntersectLineWithOther(pair, q_line, q_other, intersection_points);

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
      LineTo3DBase<pair_type>::ProjectGaussPointsOnSegmentToOther(
          pair, q_line, q_other, pair->GetEvaluationData()->GetGaussPoints(), segments[0]);
    }
  }
}

/**
 *
 */
template <typename pair_type>
std::vector<bool>&
GEOMETRYPAIR::LineTo3DGaussPointProjection<pair_type>::GetLineProjectionVectorMutable(
    const pair_type* pair)
{
  // Get the Gauss point projection tracker for this line element.
  int line_element_id = pair->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      pair->GetEvaluationData()->GetGaussPointProjectionTrackerMutable();
  auto find = projection_tracker.find(line_element_id);
  if (find == projection_tracker.end())
    dserror("Could not find the projection tracker for line id %d.", line_element_id);
  return find->second;
}


/**
 *
 */
template <typename pair_type>
void GEOMETRYPAIR::LineTo3DSegmentation<pair_type>::Evaluate(const pair_type* pair,
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<other::n_dof_, 1, scalar_type>& q_other,
    std::vector<LineSegment<scalar_type>>& segments)
{
  // Only zero segments are expected.
  if (segments.size() > 0)
    dserror("There should be zero segments for the segmentation method. The actual value is %d!",
        segments.size());

  // Number of search points.
  unsigned int n_search_points = pair->GetEvaluationData()->GetNumberOfSearchPoints();

  // Set up vector with projection points for the search points.
  std::vector<ProjectionPoint1DTo3D<scalar_type>> search_points;
  search_points.reserve(n_search_points);
  LINALG::Matrix<3, 1, scalar_type> xi_start;
  StartValues<other::geometry_type_>::Set(xi_start);
  scalar_type eta;
  for (unsigned int i_search_point = 0; i_search_point < n_search_points; i_search_point++)
  {
    eta = -1. + i_search_point * 2. / (n_search_points - 1.);
    search_points.push_back(ProjectionPoint1DTo3D<scalar_type>(eta, xi_start));
  }

  // Project all the search points.
  unsigned int dummy;
  unsigned int n_projections;
  LineTo3DBase<pair_type>::ProjectPointsOnLineToOther(
      pair, q_line, q_other, search_points, dummy, n_projections);

  // If no point could be projected return, as we assume that we wont find a surface projection.
  // This usually happens for higher order elements, more search points can be a solution to
  // that problem.
  if (n_projections != 0)
  {
    // Now we start from every search point that could be projected (not necessary valid) and look
    // for surface intersections of the line with the other geometry. The intersection points are
    // stored in a std::set that will keep them unique and in order.
    std::set<ProjectionPoint1DTo3D<scalar_type>> intersection_points;

    // If the start and/or end points (nodes of the line element) could be projected valid, add them
    // to the intersections set. This helps for some cases where the line just touches the surface
    // of the other geometry.
    if (search_points.front().GetProjectionResult() == ProjectionResult::projection_found_valid)
      intersection_points.insert(search_points.front());
    if (search_points.back().GetProjectionResult() == ProjectionResult::projection_found_valid)
      intersection_points.insert(search_points.back());

    // Vector for intersection point search.
    std::vector<ProjectionPoint1DTo3D<scalar_type>> search_intersection_points;

    // Starting from each search point, try to project to all surfaces of the other geometry.
    for (auto const& point : search_points)
    {
      // Only use search points that could be projected.
      if (point.GetProjectionResult() != ProjectionResult::projection_not_found)
      {
        // Get the intersections with the other geometry.
        pair->IntersectLineWithOther(
            q_line, q_other, search_intersection_points, point.GetEta(), point.GetXi());

        // Add the found intersection points to the set.
        for (auto& found_point : search_intersection_points)
          intersection_points.insert(found_point);
      }
    }

    // In the case of zero and one intersection points, no segmentation is needed. One point only
    // occurs when the line just touches the other geometry. This is the reason why the start and/or
    // end points are added to the set of found intersection.
    if (intersection_points.size() > 1)
    {
      // The intersection points in intersection_points are in order. Now it is checked if a
      // point between the intersection points is inside or outside the other geometry and therefore
      // it can be decided if the segment is part of this pair. By doing so we can avoid
      // complications in cases when a line is exactly between two other geometries.

      // We reuse the vector created for the search points.
      search_points.clear();
      search_points.reserve(intersection_points.size() - 1);

      // Loop over the middle points and check if they are projected valid or not.
      std::set<GEOMETRYPAIR::LineSegment<double>>& segment_tracker =
          pair->GetSegmentTrackingSetMutable();
      ProjectionResult projection_result;
      bool last_segment_active = false;
      scalar_type eta_start;
      unsigned int counter = 0;
      for (typename std::set<ProjectionPoint1DTo3D<scalar_type>>::iterator set_iterator =
               intersection_points.begin();
           set_iterator != intersection_points.end(); ++set_iterator)
      {
        // Reference to this current point.
        const ProjectionPoint1DTo3D<scalar_type>& start_point = *set_iterator;

        // Get the next intersection point and calculate the projection. This can only be done if
        // the iterator is not on its last iteration.
        if (counter != intersection_points.size() - 1)
        {
          // Reference to the next point.
          const ProjectionPoint1DTo3D<scalar_type>& end_point = *std::next(set_iterator);

          // Get starting points for projection.
          eta = 0.5 * (start_point.GetEta() + end_point.GetEta());
          xi_start = start_point.GetXi();
          xi_start += end_point.GetXi();
          xi_start.Scale(0.5);

          // Project and check result.
          LineTo3DBase<pair_type>::ProjectPointOnLineToOther(
              pair, q_line, q_other, eta, xi_start, projection_result);
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
          // This point is outside of the other geometry -> if current segment exists finish it.
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
              LineTo3DBase<pair_type>::ProjectGaussPointsOnSegmentToOther(pair, q_line, q_other,
                  pair->GetEvaluationData()->GetGaussPoints(), segments.back());
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
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::LineTo3DGaussPointProjection<
    GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_hex8>>;
template class GEOMETRYPAIR::LineTo3DGaussPointProjection<
    GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_hex20>>;
template class GEOMETRYPAIR::LineTo3DGaussPointProjection<
    GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_hex27>>;
template class GEOMETRYPAIR::LineTo3DGaussPointProjection<
    GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_tet4>>;
template class GEOMETRYPAIR::LineTo3DGaussPointProjection<
    GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_tet10>>;
template class GEOMETRYPAIR::LineTo3DGaussPointProjection<
    GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjection<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_nurbs27>>;

template class GEOMETRYPAIR::LineTo3DSegmentation<
    GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_hex8>>;
template class GEOMETRYPAIR::LineTo3DSegmentation<
    GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_hex20>>;
template class GEOMETRYPAIR::LineTo3DSegmentation<
    GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_hex27>>;
template class GEOMETRYPAIR::LineTo3DSegmentation<
    GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_tet4>>;
template class GEOMETRYPAIR::LineTo3DSegmentation<
    GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_tet10>>;
template class GEOMETRYPAIR::LineTo3DSegmentation<
    GEOMETRYPAIR::GeometryPairLineToVolumeSegmentation<double, GEOMETRYPAIR::t_hermite,
        GEOMETRYPAIR::t_nurbs27>>;
