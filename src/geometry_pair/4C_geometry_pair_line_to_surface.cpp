/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and surfaces.

\level 1
*/


#include "4C_geometry_pair_line_to_surface.hpp"

#include "4C_beam3_base.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"
#include "4C_geometry_pair_constants.hpp"
#include "4C_geometry_pair_element.hpp"
#include "4C_geometry_pair_element_evaluation_functions.hpp"
#include "4C_geometry_pair_scalar_types.hpp"
#include "4C_geometry_pair_utility_classes.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_fad.hpp"

FOUR_C_NAMESPACE_OPEN

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::GeometryPairLineToSurface(
    const Core::Elements::Element* element1, const Core::Elements::Element* element2,
    const Teuchos::RCP<GEOMETRYPAIR::LineToSurfaceEvaluationData>& line_to_surface_evaluation_data)
    : GeometryPair(element1, element2),
      line_to_surface_evaluation_data_(line_to_surface_evaluation_data),
      is_unit_test_(false)
{
  // For the current implementation, the line element has to be on the same processor as the pair
  // object. This is because the tracking vector in LineTo3DEvaluationData is only local and we
  // need this vector for segmentation e.t.c.
  int myrank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (element1->Owner() != myrank)
    FOUR_C_THROW(
        "The GeometryPairLineToSurface pair has to be on the same processor as the line element! "
        "Currently the pair is on rank %d, the line element on %d!",
        myrank, element1->Owner());
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::ProjectPointToOther(
    const Core::LinAlg::Matrix<3, 1, scalar_type>& point,
    const ElementData<surface, scalar_type>& element_data_surface,
    Core::LinAlg::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
    const bool min_one_iteration) const
{
  ProjectPointToSurface(point, element_data_surface, xi, projection_result,
      get_surface_normal_influence_direction(element_data_surface), min_one_iteration);
}


/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::intersect_line_with_other(
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<surface, scalar_type>& element_data_surface,
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& intersection_points,
    const scalar_type& eta_start, const Core::LinAlg::Matrix<3, 1, scalar_type>& xi_start) const
{
  unsigned int n_faces;
  std::vector<unsigned int> face_fixed_parameters;
  std::vector<double> face_fixed_values;
  get_face_fixed_parameters(n_faces, face_fixed_parameters, face_fixed_values);

  // Clear the input vector.
  intersection_points.clear();
  intersection_points.reserve(n_faces);

  // Create variables.
  scalar_type eta;
  Core::LinAlg::Matrix<3, 1, scalar_type> xi;
  ProjectionResult intersection_found;

  // Try to intersect the beam with each face.
  for (unsigned int i = 0; i < n_faces; i++)
  {
    // Set starting values.
    xi = xi_start;
    eta = eta_start;

    // Intersect the line with the surface.
    intersect_line_with_surface_edge(element_data_line, element_data_surface,
        face_fixed_parameters[i], face_fixed_values[i], eta, xi, intersection_found);

    // If a valid intersection is found, add it to the output vector.
    if (intersection_found == ProjectionResult::projection_found_valid)
    {
      intersection_points.push_back(ProjectionPoint1DTo3D<scalar_type>(eta, xi));
      intersection_points.back().SetIntersectionFace(i);
    }
  }
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line,
    surface>::intersect_line_with_surface_edge(const ElementData<line, scalar_type>&
                                                   element_data_line,
    const ElementData<surface, scalar_type>& element_data_surface,
    const unsigned int& fixed_parameter, const double& fixed_value, scalar_type& eta,
    Core::LinAlg::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
    const bool min_one_iteration) const
{
  // Check the input parameters.
  {
    if (surface::geometry_type_ == DiscretizationTypeGeometry::quad && fixed_parameter > 1)
      FOUR_C_THROW(
          "Fixed_parameter in intersect_line_with_surface_edge has to be smaller than 2 with a "
          "quad element.");
    else if (surface::geometry_type_ == DiscretizationTypeGeometry::none)
      FOUR_C_THROW("Wrong DiscretizationTypeGeometry type given.");
    else if (fixed_parameter > 2)
      FOUR_C_THROW("fixed_parameter in intersect_line_with_surface_edge can be 2 at maximum.");
  }

  // Approximated influence size of the surface.
  const double normal_influence_direction =
      get_surface_normal_influence_direction(element_data_surface);

  // Initialize data structures.
  // Point on line.
  Core::LinAlg::Matrix<3, 1, scalar_type> r_line;
  Core::LinAlg::Matrix<3, 1, scalar_type> dr_line;

  // Point on surface.
  Core::LinAlg::Matrix<3, 1, scalar_type> r_surface;
  Core::LinAlg::Matrix<3, 3, scalar_type> dr_surface;

  // Residuum.
  Core::LinAlg::Matrix<4, 1, scalar_type> residuum;
  Core::LinAlg::Matrix<4, 1, scalar_type> delta_xi;
  // Initialize the increment with a value that will not pass the first convergence check.
  delta_xi.put_scalar(10 * Constants::projection_xi_eta_tol);

  // Jacobian / inverse.
  Core::LinAlg::Matrix<4, 4, scalar_type> J_J_inv;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  {
    // Local Newton iteration.
    unsigned int counter = 0;
    while (counter < Constants::local_newton_iter_max)
    {
      // Evaluate the position and its derivative on the line.
      EvaluatePosition<line>(eta, element_data_line, r_line);
      EvaluatePositionDerivative1<line>(eta, element_data_line, dr_line);

      // Evaluate the position and its derivative on the surface.
      EvaluateSurfacePositionAndDerivative(element_data_surface, xi, r_surface, dr_surface);

      // Evaluate the residuum $r_{surface} - r_{line} = R_{pos}$ and $xi(i) - value = R_{edge}$
      J_J_inv.put_scalar(0.);
      residuum.put_scalar(0.);
      for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
        residuum(i_dir) = r_surface(i_dir) - r_line(i_dir);

      if (fixed_parameter < 2)
      {
        residuum(3) = xi(fixed_parameter) - fixed_value;
        J_J_inv(3, fixed_parameter) = 1.;
      }
      else
      {
        for (unsigned int i = 0; i < 2; i++)
        {
          residuum(3) += xi(i);
          J_J_inv(3, i) = 1.;
        }
        residuum(3) -= fixed_value;
      }

      if (counter == 0 and min_one_iteration)
      {
        // if the min_one_iteration flag is set we run at least one iteration, so the dependency on
        // FAD variables is calculated correctly.
      }
      else if (Core::FADUtils::VectorNorm(residuum) < Constants::local_newton_res_tol &&
               Core::FADUtils::VectorNorm(delta_xi) < Constants::projection_xi_eta_tol)
      {
        // System is solved, now check if the parameter coordinates are valid.
        if (ValidParameter1D(eta) &&
            ValidParameterSurface<scalar_type, surface>(xi, normal_influence_direction))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (Core::FADUtils::VectorNorm(residuum) > Constants::local_newton_res_max) break;

      // Fill up the jacobian.
      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < 3; j++)
        {
          J_J_inv(i, j) = dr_surface(i, j);
        }
        J_J_inv(i, 3) = -dr_line(i);
      }

      // Solve the linearized system.
      if (Core::LinAlg::SolveLinearSystemDoNotThrowErrorOnZeroDeterminantScaled(
              J_J_inv, residuum, delta_xi, Constants::local_newton_det_tol))
      {
        // Set the new parameter coordinates.
        eta -= delta_xi(3);
        for (unsigned int i = 0; i < 3; i++) xi(i) -= delta_xi(i);

        // Advance Newton iteration counter.
        counter++;
      }
      else
        break;
    }
  }
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
double GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line,
    surface>::get_surface_normal_influence_direction(const ElementData<surface, scalar_type>&
        element_data_surface) const
{
  if (is_unit_test_)
    return -1.0;
  else
  {
    double surface_size = get_surface_size(element_data_surface);
    double line_tube_size_radius = (dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(Element1()))
                                       ->get_circular_cross_section_radius_for_interactions();
    return std::max(surface_size, 3.0 * line_tube_size_radius);
  }
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
double GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::get_surface_size(
    const ElementData<surface, scalar_type>& element_data_surface) const
{
  // Get the position of the first 3 nodes of the surface.
  Core::LinAlg::Matrix<2, 1, double> xi_corner_node;
  Core::LinAlg::Matrix<3, 1, Core::LinAlg::Matrix<3, 1, double>> corner_nodes;
  Core::LinAlg::SerialDenseMatrix nodal_coordinates =
      Core::FE::getEleNodeNumbering_nodes_paramspace(surface::discretization_);
  const auto element_data_surface_double =
      ElementDataToDouble<surface>::ToDouble(element_data_surface);
  for (unsigned int i_node = 0; i_node < 3; i_node++)
  {
    for (unsigned int i_dim = 0; i_dim < 2; i_dim++)
      xi_corner_node(i_dim) = nodal_coordinates(i_dim, i_node);
    EvaluatePosition<surface>(xi_corner_node, element_data_surface_double, corner_nodes(i_node));
  }

  // Calculate the maximum distance between the three points.
  double max_distance = 0.0;
  double distance = 0.0;
  Core::LinAlg::Matrix<3, 1, double> diff;
  for (unsigned int i_node = 0; i_node < 3; i_node++)
  {
    for (unsigned int j_node = 0; j_node < 3; j_node++)
    {
      if (i_node == j_node) continue;

      diff = corner_nodes(j_node);
      diff -= corner_nodes(i_node);
      distance = Core::FADUtils::VectorNorm(diff);
      if (distance > max_distance) max_distance = distance;
    }
  }

  return max_distance;
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::get_face_fixed_parameters(
    unsigned int& n_faces, std::vector<unsigned int>& face_fixed_parameters,
    std::vector<double>& face_fixed_values) const
{
  if (surface::geometry_type_ == DiscretizationTypeGeometry::quad)
  {
    n_faces = 4;
    face_fixed_parameters = {0, 0, 1, 1};
    face_fixed_values = {-1., 1., -1., 1.};
  }
  else if (surface::geometry_type_ == DiscretizationTypeGeometry::triangle)
  {
    n_faces = 3;
    face_fixed_parameters = {0, 1, 2};
    face_fixed_values = {0., 0., 1.};
  }
  else
  {
    FOUR_C_THROW("Wrong DiscretizationTypeGeometry given!");
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurfaceFADWrapper<scalar_type, line, surface>::pre_evaluate(
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<surface, scalar_type>& element_data_surface,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Call pre_evaluate on the double pair.
  std::vector<LineSegment<double>> segments_double;
  geometry_pair_double_->pre_evaluate(ElementDataToDouble<line>::ToDouble(element_data_line),
      ElementDataToDouble<surface>::ToDouble(element_data_surface), segments_double);

  // Convert the created double segments to a segment of scalar type.
  segments.clear();
  for (auto& segment_double : segments_double)
  {
    // Create the segment with the scalar FAD type.
    segments.push_back(LineSegment<scalar_type>());
    CopySegment(segment_double, segments.back());
  }
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurfaceFADWrapper<scalar_type, line, surface>::evaluate(
    const ElementData<line, scalar_type>& element_data_line,
    const ElementData<surface, scalar_type>& element_data_surface,
    std::vector<LineSegment<scalar_type>>& segments) const
{
  // Convert the input segments to a segment of scalar type double.
  std::vector<LineSegment<double>> segments_double;
  for (auto& segment : segments)
  {
    // Create the segment with the scalar FAD type.
    segments_double.push_back(LineSegment<double>());
    CopySegment(segment, segments_double.back());
  }

  // Call Evaluate on the double pair.
  geometry_pair_double_->evaluate(ElementDataToDouble<line>::ToDouble(element_data_line),
      ElementDataToDouble<surface>::ToDouble(element_data_surface), segments_double);

  // Get the face parameters.
  unsigned int n_faces;
  std::vector<unsigned int> face_fixed_parameters;
  std::vector<double> face_fixed_values;
  this->get_face_fixed_parameters(n_faces, face_fixed_parameters, face_fixed_values);

  // Initialize variables for the projections.
  ProjectionResult projection_result = ProjectionResult::none;
  Core::LinAlg::Matrix<3, 1, scalar_type> point_in_space;

  // If segments are found, convert them to FAD segments.
  segments.clear();
  for (auto& segment_double : segments_double)
  {
    // Create the segment with the scalar FAD type.
    segments.push_back(LineSegment<scalar_type>());
    LineSegment<scalar_type>& new_segment = segments.back();

    // Add the projection point to an array.
    std::array<std::reference_wrapper<ProjectionPoint1DTo3D<scalar_type>>, 2>
        segment_start_end_points = {new_segment.GetStartPoint(), new_segment.GetEndPoint()};
    segment_start_end_points[0].get().set_from_other_point_double(segment_double.GetStartPoint());
    segment_start_end_points[1].get().set_from_other_point_double(segment_double.GetEndPoint());

    // If the start or end points are intersection points, the intersections have to be reevaluated.
    for (auto& point : segment_start_end_points)
    {
      const int intersection_face = point.get().GetIntersectionFace();
      if (intersection_face >= 0)
      {
        this->intersect_line_with_surface_edge(element_data_line, element_data_surface,
            face_fixed_parameters[intersection_face], face_fixed_values[intersection_face],
            point.get().GetEta(), point.get().GetXi(), projection_result, true);
      }
    }

    // Reevaluate the integration points along the segment.
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& projection_points =
        new_segment.GetProjectionPoints();
    projection_points.resize(segment_double.get_number_of_projection_points());
    for (unsigned int i_point = 0; i_point < segment_double.get_number_of_projection_points();
         i_point++)
    {
      // Position of the projection point within the segment.
      auto& projection_point_double = segment_double.GetProjectionPoints()[i_point];
      const double factor = (projection_point_double.GetEta() - segment_double.GetEtadata()) /
                            segment_double.GetSegmentLength();

      // Calculate spatial point.
      auto& projection_point = projection_points[i_point];
      projection_point.set_from_other_point_double(projection_point_double);
      projection_point.SetEta(new_segment.GetEtadata() + new_segment.GetSegmentLength() * factor);

      EvaluatePosition<line>(projection_point.GetEta(), element_data_line, point_in_space);

      // Calculate the projection.
      this->ProjectPointToOther(
          point_in_space, element_data_surface, projection_point.GetXi(), projection_result, true);
    }
  }
}

/**
 *
 */
template <typename scalar_type, typename surface>
void GEOMETRYPAIR::ProjectPointToSurface(const Core::LinAlg::Matrix<3, 1, scalar_type>& point,
    const ElementData<surface, scalar_type>& element_data_surface,
    Core::LinAlg::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
    const double normal_influence_direction, const bool min_one_iteration)
{
  // Vectors in 3D.
  Core::LinAlg::Matrix<3, 1, scalar_type> r_surface;
  Core::LinAlg::Matrix<3, 1, scalar_type> delta_xi;
  // Initialize the increment with a value that will not pass the first convergence check.
  delta_xi.put_scalar(10 * Constants::projection_xi_eta_tol);
  Core::LinAlg::Matrix<3, 1, scalar_type> residuum;

  // Jacobian / inverse.
  Core::LinAlg::Matrix<3, 3, scalar_type> J_J_inv;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  // Local Newton iteration.
  {
    unsigned int counter = 0;
    while (counter < Constants::local_newton_iter_max)
    {
      // Evaluate the position and its derivative on the surface.
      EvaluateSurfacePositionAndDerivative(element_data_surface, xi, r_surface, J_J_inv);

      // Evaluate the residuum $r_{solid} - r_{point} = R_{pos}$.
      residuum = r_surface;
      residuum -= point;

      if (counter == 0 and min_one_iteration)
      {
        // if the min_one_iteration flag is set we run at least one iteration, so the dependency on
        // FAD variables is calculated correctly.
      }
      else if (Core::FADUtils::VectorNorm(residuum) < Constants::local_newton_res_tol &&
               Core::FADUtils::VectorNorm(delta_xi) < Constants::projection_xi_eta_tol)
      {
        if (ValidParameterSurface<scalar_type, surface>(xi, normal_influence_direction))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (Core::FADUtils::VectorNorm(residuum) > Constants::local_newton_res_max) break;

      // Solve the linearized system.
      if (Core::LinAlg::SolveLinearSystemDoNotThrowErrorOnZeroDeterminantScaled(
              J_J_inv, residuum, delta_xi, Constants::local_newton_det_tol))
      {
        // Set the new parameter coordinates.
        xi -= delta_xi;

        // Advance Newton iteration counter.
        counter++;
      }
      else
        break;
    }
  }
}


/**
 * Explicit template initialization of template class.
 */
namespace GEOMETRYPAIR
{
  template class GeometryPairLineToSurface<double, t_line2, t_tri3>;
  template class GeometryPairLineToSurface<double, t_line2, t_tri6>;
  template class GeometryPairLineToSurface<double, t_line2, t_quad4>;
  template class GeometryPairLineToSurface<double, t_line2, t_quad8>;
  template class GeometryPairLineToSurface<double, t_line2, t_quad9>;
  template class GeometryPairLineToSurface<double, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_line2,
      t_tri3>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_line2,
      t_tri6>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_line2,
      t_quad4>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_line2,
      t_quad8>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_line2,
      t_quad9>;
  template class GeometryPairLineToSurface<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_line2, t_tri3>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_line2, t_tri6>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_line2, t_quad4>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_line2, t_quad8>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_line2, t_quad9>;
  template class GeometryPairLineToSurface<
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_tri3>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_tri6>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad4>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad8>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_line2, t_quad9>;
  template class GeometryPairLineToSurfaceFADWrapper<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_line2,
      t_tri3>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_line2,
      t_tri6>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_line2,
      t_quad4>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_line2,
      t_quad8>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_line2,
      t_quad9>;
  template class GeometryPairLineToSurfaceFADWrapper<
      line_to_surface_patch_scalar_type_fixed_size<t_line2, t_nurbs9>, t_line2, t_nurbs9>;

  template class GeometryPairLineToSurface<double, t_hermite, t_tri3>;
  template class GeometryPairLineToSurface<double, t_hermite, t_tri6>;
  template class GeometryPairLineToSurface<double, t_hermite, t_quad4>;
  template class GeometryPairLineToSurface<double, t_hermite, t_quad8>;
  template class GeometryPairLineToSurface<double, t_hermite, t_quad9>;
  template class GeometryPairLineToSurface<double, t_hermite, t_nurbs9>;

  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_hermite,
      t_tri3>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_hermite,
      t_tri6>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_hermite,
      t_quad4>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_hermite,
      t_quad8>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type_1st_order, t_hermite,
      t_quad9>;
  template class GeometryPairLineToSurface<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_hermite, t_tri3>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_hermite, t_tri6>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_hermite, t_quad4>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_hermite, t_quad8>;
  template class GeometryPairLineToSurface<line_to_surface_patch_scalar_type, t_hermite, t_quad9>;
  template class GeometryPairLineToSurface<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;

  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri3>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_tri6>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad4>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad8>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type_1st_order,
      t_hermite, t_quad9>;
  template class GeometryPairLineToSurfaceFADWrapper<
      line_to_surface_patch_scalar_type_fixed_size_1st_order<t_hermite, t_nurbs9>, t_hermite,
      t_nurbs9>;

  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_hermite,
      t_tri3>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_hermite,
      t_tri6>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_hermite,
      t_quad4>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_hermite,
      t_quad8>;
  template class GeometryPairLineToSurfaceFADWrapper<line_to_surface_patch_scalar_type, t_hermite,
      t_quad9>;
  template class GeometryPairLineToSurfaceFADWrapper<
      line_to_surface_patch_scalar_type_fixed_size<t_hermite, t_nurbs9>, t_hermite, t_nurbs9>;
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE
