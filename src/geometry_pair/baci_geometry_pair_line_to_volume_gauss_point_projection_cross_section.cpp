/*----------------------------------------------------------------------*/
/*! \file

\brief Line to volume interaction with Gauss point projection on the cylinder surface along the
line.

\level 1
*/


#include "baci_geometry_pair_line_to_volume_gauss_point_projection_cross_section.hpp"

#include "baci_beam3_base.hpp"
#include "baci_beam3_triad_interpolation_local_rotation_vectors.hpp"
#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_geometry_pair_element_evaluation_functions.hpp"
#include "baci_geometry_pair_line_to_3D_evaluation_data.hpp"
#include "baci_geometry_pair_utility_classes.hpp"
#include "baci_lib_element.hpp"

#include <math.h>

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename scalar_type, typename line, typename volume>
GEOMETRYPAIR::GeometryPairLineToVolumeGaussPointProjectionCrossSection<scalar_type, line,
    volume>::GeometryPairLineToVolumeGaussPointProjectionCrossSection(const DRT::Element* element1,
    const DRT::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::LineTo3DEvaluationData>& evaluation_data)
    : GeometryPairLineToVolume<scalar_type, line, volume>(element1, element2, evaluation_data)
{
  // Check if a projection tracking vector exists for this line element. If not a new one is
  // created.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->line_to_3d_evaluation_data_->GetGaussPointProjectionTracker();

  if (projection_tracker.find(line_element_id) == projection_tracker.end())
  {
    int n_gauss_points =
        this->line_to_3d_evaluation_data_->GetNumberOfGaussPoints() *
        this->line_to_3d_evaluation_data_->GetNumberOfIntegrationPointsCircumference();
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
    volume>::PreEvaluate(const ElementData<line, scalar_type>& element_data_line,
    const ElementData<volume, scalar_type>& element_data_volume,
    std::vector<LineSegment<scalar_type>>& segments,
    const LARGEROTATIONS::TriadInterpolationLocalRotationVectors<3, double>*
        line_triad_interpolation) const
{
  // Get the Gauss point projection tracker for this line element.
  std::vector<bool>& line_projection_tracker = GetLineProjectionVector();

  // Gauss rule.
  CORE::FE::IntegrationPoints1D gauss_points_axis =
      this->line_to_3d_evaluation_data_->GetGaussPoints();
  unsigned int n_gauss_points_axis = this->line_to_3d_evaluation_data_->GetNumberOfGaussPoints();
  unsigned int n_integration_points_circ =
      this->line_to_3d_evaluation_data_->GetNumberOfIntegrationPointsCircumference();

  // Initilaize variables for the projection.
  scalar_type eta;
  double alpha;
  CORE::LINALG::Matrix<3, 1, scalar_type> r_line;
  CORE::LINALG::Matrix<3, 3, scalar_type> triad;
  CORE::LINALG::Matrix<3, 1, scalar_type> r_cross_section;
  CORE::LINALG::Matrix<3, 1, scalar_type> r_surface;
  CORE::LINALG::Matrix<3, 1, scalar_type> eta_cross_section(true);
  CORE::LINALG::Matrix<2, 1, scalar_type> eta_cross_section_2d;
  CORE::LINALG::Matrix<3, 1, scalar_type> xi_volume;
  ProjectionResult projection_result;
  segments.clear();
  bool one_projects = false;
  LineSegment<scalar_type> projection_point_segment;

  // Get the radius from the beam element.
  const double radius = (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(this->Element1()))
                            ->GetCircularCrossSectionRadiusForInteractions();

  // Loop over Gauss points and check if they project to this volume.
  for (unsigned int index_gp_axis = 0; index_gp_axis < n_gauss_points_axis; index_gp_axis++)
  {
    // Parameter coordinate along the line.
    eta = gauss_points_axis.qxg[index_gp_axis][0];

    // Get the triad on the line.
    if (line_triad_interpolation != nullptr)
      line_triad_interpolation->GetInterpolatedTriadAtXi(triad, eta);
    else
      GEOMETRYPAIR::EvaluateTriadAtPlaneCurve<line>(eta, element_data_line, triad);

    // Get the position on the line.
    GEOMETRYPAIR::EvaluatePosition<line>(eta, element_data_line, r_line);

    for (unsigned int index_gp_circ = 0; index_gp_circ < n_integration_points_circ; index_gp_circ++)
    {
      // Index of the current Gauss point in the tracking vector.
      unsigned int index_gp = index_gp_axis * n_integration_points_circ + index_gp_circ;

      // Only check points that do not already have a valid projection.
      if (line_projection_tracker[index_gp] == false)
      {
        // Coordinates in the cross section.
        alpha = 2.0 * M_PI / double(n_integration_points_circ) * index_gp_circ;
        eta_cross_section(0) = 0;
        eta_cross_section(1) = cos(alpha) * radius;
        eta_cross_section(2) = sin(alpha) * radius;

        // Get the point on the beams surface.
        r_cross_section.Multiply(triad, eta_cross_section);
        r_surface = r_line;
        r_surface += r_cross_section;

        // Project point to the volume.
        this->ProjectPointToOther(r_surface, element_data_volume, xi_volume, projection_result);
        if (projection_result == ProjectionResult::projection_found_valid)
        {
          // Valid Gauss point was found, add to this segment and set tracking point to true.
          ProjectionPoint1DTo3D<scalar_type> new_point(eta, xi_volume,
              gauss_points_axis.qwgt[index_gp_axis] * 2.0 / double(n_integration_points_circ));
          for (unsigned int i_dim = 0; i_dim < 2; i_dim++)
            eta_cross_section_2d(i_dim) = eta_cross_section(i_dim + 1);
          new_point.SetEtaCrossSection(eta_cross_section_2d);
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
    volume>::Evaluate(const ElementData<line, scalar_type>& element_data_line,
    const ElementData<volume, scalar_type>& element_data_volume,
    std::vector<LineSegment<scalar_type>>& segments) const
{
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
    const std::vector<bool>& projection_vector = GetLineProjectionVector();
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
    scalar_type, line, volume>::GetLineProjectionVector() const
{
  // Get the Gauss point projection tracker for this line element.
  int line_element_id = this->Element1()->Id();
  std::map<int, std::vector<bool>>& projection_tracker =
      this->line_to_3d_evaluation_data_->GetGaussPointProjectionTracker();
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

FOUR_C_NAMESPACE_CLOSE
