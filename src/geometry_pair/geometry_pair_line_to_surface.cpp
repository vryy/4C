/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and surfaces.

\level 1
*/


#include "geometry_pair_line_to_surface.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_utility_classes.H"
#include "geometry_pair_constants.H"
#include "geometry_pair_scalar_types.H"

#include "linalg_utils_densematrix_inverse.H"
#include "utils_exceptions.H"
#include "discretization_fem_general_utils_local_connectivity_matrices.H"


/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::Init(
    const DRT::Element* element1, const DRT::Element* element2)
{
  // Call init of base class.
  GeometryPair::Init(element1, element2);

  // For the current implementation, the line element has to be on the same processor as the pair
  // object. This is because the tracking vector in LineTo3DEvaluationData is only local and we
  // need this vector for segmentation e.t.c.
  int myrank = -1;
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  if (element1->Owner() != myrank)
    dserror(
        "The GeometryPairLineToSurface pair has to be on the same processor as the line element! "
        "Currently the pair is on rank %d, the line element on %d!",
        myrank, element1->Owner());
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::ProjectPointToOther(
    const LINALG::Matrix<3, 1, scalar_type>& point,
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
    LINALG::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals,
    const bool min_one_iteration) const
{
  // Initialize data structures

  // Approximated size of surface and beam diameter for valid projection check.
  const scalar_type surface_size = GetSurfaceSize(q_surface);
  const double beam_radius = GetLineRadius();

  // Vectors in 3D.
  LINALG::Matrix<3, 1, scalar_type> r_surface;
  LINALG::Matrix<3, 1, scalar_type> delta_xi;
  // Initialize the increment with a value that will not pass the first convergence check.
  delta_xi.PutScalar(10 * CONSTANTS::projection_xi_eta_tol);
  LINALG::Matrix<3, 1, scalar_type> residuum;

  // Jacobian / inverse.
  LINALG::Matrix<3, 3, scalar_type> J_J_inv;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  // Local Newton iteration.
  {
    unsigned int counter = 0;
    while (counter < CONSTANTS::local_newton_iter_max)
    {
      // Evaluate the position and its derivative on the surface.
      EvaluateSurfacePositionAndDerivative(q_surface, xi, r_surface, J_J_inv, nodal_normals);

      // Evaluate the residuum $r_{solid} - r_{point} = R_{pos}$.
      residuum = r_surface;
      residuum -= point;

      if (counter == 0 and min_one_iteration)
      {
        // if the min_one_iteration flag is set we run at least one iteration, so the dependency on
        // FAD variables is calculated correctly.
      }
      else if (FADUTILS::VectorNorm(residuum) < CONSTANTS::local_newton_res_tol &&
               FADUTILS::VectorNorm(delta_xi) < CONSTANTS::projection_xi_eta_tol)
      {
        if (ValidParameterSurface(xi, surface_size, beam_radius))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (FADUTILS::VectorNorm(residuum) > CONSTANTS::local_newton_res_max) break;

      // Solve the linearized system.
      if (LINALG::SolveLinearSystemDoNotThrowErrorOnZeroDeterminantScaled(
              J_J_inv, residuum, delta_xi, CONSTANTS::local_newton_det_tol))
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
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::IntersectLineWithOther(
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& intersection_points,
    const scalar_type& eta_start, const LINALG::Matrix<3, 1, scalar_type>& xi_start,
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const
{
  unsigned int n_faces;
  std::vector<unsigned int> face_fixed_parameters;
  std::vector<double> face_fixed_values;
  GetFaceFixedParameters(n_faces, face_fixed_parameters, face_fixed_values);

  // Clear the input vector.
  intersection_points.clear();
  intersection_points.reserve(n_faces);

  // Create variables.
  scalar_type eta;
  LINALG::Matrix<3, 1, scalar_type> xi;
  ProjectionResult intersection_found;

  // Try to intersect the beam with each face.
  for (unsigned int i = 0; i < n_faces; i++)
  {
    // Set starting values.
    xi = xi_start;
    eta = eta_start;

    // Intersect the line with the surface.
    IntersectLineWithSurfaceEdge(q_line, q_surface, face_fixed_parameters[i], face_fixed_values[i],
        eta, xi, intersection_found, nodal_normals);

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
    surface>::EvaluateSurfacePositionAndDerivative(const LINALG::Matrix<surface::n_dof_, 1,
                                                       scalar_type>& q_surface,
    const LINALG::Matrix<3, 1, scalar_type>& xi, LINALG::Matrix<3, 1, scalar_type>& r,
    LINALG::Matrix<3, 3, scalar_type>& dr,
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const
{
  // Create a nested FAD type.
  using FAD_outer = Sacado::ELRFad::SLFad<scalar_type, 3>;

  // Setup the AD variables.
  LINALG::Matrix<3, 1, FAD_outer> xi_AD;
  for (unsigned int i_dim = 0; i_dim < 3; i_dim++) xi_AD(i_dim) = FAD_outer(3, i_dim, xi(i_dim));
  LINALG::Matrix<3, 1, FAD_outer> r_AD;

  // Evaluate the position.
  EvaluateSurfacePosition<surface>(xi_AD, q_surface, r_AD, Element2(), nodal_normals);

  // Extract the return values from the AD types.
  for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
  {
    r(i_dir) = r_AD(i_dir).val();
    for (unsigned int i_dim = 0; i_dim < 3; i_dim++) dr(i_dir, i_dim) = r_AD(i_dir).dx(i_dim);
  }
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line,
    surface>::IntersectLineWithSurfaceEdge(const LINALG::Matrix<line::n_dof_, 1, scalar_type>&
                                               q_line,
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
    const unsigned int& fixed_parameter, const double& fixed_value, scalar_type& eta,
    LINALG::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result,
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals,
    const bool min_one_iteration) const
{
  // Check the input parameters.
  {
    if (surface::geometry_type_ == DiscretizationTypeGeometry::quad && fixed_parameter > 1)
      dserror(
          "Fixed_parameter in IntersectLineWithSurfaceEdge has to be smaller than 2 with a "
          "quad element.");
    else if (surface::geometry_type_ == DiscretizationTypeGeometry::none)
      dserror("Wrong DiscretizationTypeGeometry type given.");
    else if (fixed_parameter > 2)
      dserror("fixed_parameter in IntersectLineWithSurfaceEdge can be 2 at maximum.");
  }

  // Approximated size of surface and beam diameter for valid projection check.
  const scalar_type surface_size = GetSurfaceSize(q_surface);
  const double beam_radius = GetLineRadius();

  // Initialize data structures.
  // Point on line.
  LINALG::Matrix<3, 1, scalar_type> r_line;
  LINALG::Matrix<3, 1, scalar_type> dr_line;

  // Point on surface.
  LINALG::Matrix<3, 1, scalar_type> r_surface;
  LINALG::Matrix<3, 3, scalar_type> dr_surface;

  // Residuum.
  LINALG::Matrix<4, 1, scalar_type> residuum;
  LINALG::Matrix<4, 1, scalar_type> delta_xi;
  // Initialize the increment with a value that will not pass the first convergence check.
  delta_xi.PutScalar(10 * CONSTANTS::projection_xi_eta_tol);

  // Jacobian / inverse.
  LINALG::Matrix<4, 4, scalar_type> J_J_inv;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  {
    // Local Newton iteration.
    unsigned int counter = 0;
    while (counter < CONSTANTS::local_newton_iter_max)
    {
      // Evaluate the position and its derivative on the line.
      EvaluatePosition<line>(eta, q_line, r_line, Element1());
      EvaluatePositionDerivative1<line>(eta, q_line, dr_line, Element1());

      // Evaluate the position and its derivative on the surface.
      EvaluateSurfacePositionAndDerivative(q_surface, xi, r_surface, dr_surface, nodal_normals);

      // Evaluate the residuum $r_{surface} - r_{line} = R_{pos}$ and $xi(i) - value = R_{edge}$
      J_J_inv.PutScalar(0.);
      residuum.PutScalar(0.);
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
      else if (FADUTILS::VectorNorm(residuum) < CONSTANTS::local_newton_res_tol &&
               FADUTILS::VectorNorm(delta_xi) < CONSTANTS::projection_xi_eta_tol)
      {
        // System is solved, now check if the parameter coordinates are valid.
        if (ValidParameter1D(eta) && ValidParameterSurface(xi, surface_size, beam_radius))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (FADUTILS::VectorNorm(residuum) > CONSTANTS::local_newton_res_max) break;

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
      if (LINALG::SolveLinearSystemDoNotThrowErrorOnZeroDeterminantScaled(
              J_J_inv, residuum, delta_xi, CONSTANTS::local_newton_det_tol))
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
bool GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::ValidParameterSurface(
    LINALG::Matrix<3, 1, scalar_type>& xi, const scalar_type& surface_size,
    const double beam_radius) const
{
  // We only need to theck the normal distance if the coordinates are within the surface.
  if (!ValidParameter2D<surface>(xi)) return false;

  if ((-surface_size < xi(2) and xi(2) < surface_size) or
      (-3.0 * beam_radius < xi(2) and xi(2) < 3 * beam_radius))
    return true;
  else
    return false;
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
double GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::GetLineRadius() const
{
  if (is_unit_test_)
    return 1e10;
  else
    return (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(Element1()))
        ->GetCircularCrossSectionRadiusForInteractions();
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
scalar_type GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::GetSurfaceSize(
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface) const
{
  // Get the position of the first 3 nodes of the surface.
  LINALG::Matrix<2, 1, double> xi_corner_node;
  LINALG::Matrix<3, 1, LINALG::Matrix<3, 1, scalar_type>> corner_nodes;
  LINALG::SerialDenseMatrix nodal_coordinates =
      DRT::UTILS::getEleNodeNumbering_nodes_paramspace(surface::discretization_);
  for (unsigned int i_node = 0; i_node < 3; i_node++)
  {
    for (unsigned int i_dim = 0; i_dim < 2; i_dim++)
      xi_corner_node(i_dim) = nodal_coordinates(i_dim, i_node);
    EvaluatePosition<surface>(xi_corner_node, q_surface, corner_nodes(i_node), Element2());
  }

  // Calculate the maximum distance between the three points.
  scalar_type max_distance = 0.0;
  scalar_type distance = 0.0;
  LINALG::Matrix<3, 1, scalar_type> diff;
  for (unsigned int i_node = 0; i_node < 3; i_node++)
  {
    for (unsigned int j_node = 0; j_node < 3; j_node++)
    {
      if (i_node == j_node) continue;

      diff = corner_nodes(j_node);
      diff -= corner_nodes(i_node);
      distance = FADUTILS::VectorNorm(diff);
      if (distance > max_distance) max_distance = distance;
    }
  }

  return max_distance;
}

/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::GetFaceFixedParameters(
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
    dserror("Wrong DiscretizationTypeGeometry given!");
  }
}


/**
 *
 */
template <typename scalar_type, typename line, typename surface>
void GEOMETRYPAIR::GeometryPairLineToSurfaceFADWrapper<scalar_type, line, surface>::PreEvaluate(
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
    std::vector<LineSegment<scalar_type>>& segments,
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const
{
  // Call PreEvaluate on the double pair.
  std::vector<LineSegment<double>> segments_double;
  LINALG::Matrix<3 * surface::n_nodes_, 1, double> nodal_normals_double(false);
  geometry_pair_double_->PreEvaluate(FADUTILS::CastToDouble(q_line),
      FADUTILS::CastToDouble(q_surface), segments_double,
      VectorPointerToVectorDouble(nodal_normals, nodal_normals_double));

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
void GEOMETRYPAIR::GeometryPairLineToSurfaceFADWrapper<scalar_type, line, surface>::Evaluate(
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
    std::vector<LineSegment<scalar_type>>& segments,
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const
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
  LINALG::Matrix<3 * surface::n_nodes_, 1, double> nodal_normals_double(false);
  geometry_pair_double_->Evaluate(FADUTILS::CastToDouble(q_line), FADUTILS::CastToDouble(q_surface),
      segments_double, VectorPointerToVectorDouble(nodal_normals, nodal_normals_double));

  // Get the face parameters.
  unsigned int n_faces;
  std::vector<unsigned int> face_fixed_parameters;
  std::vector<double> face_fixed_values;
  this->GetFaceFixedParameters(n_faces, face_fixed_parameters, face_fixed_values);

  // Initialize variables for the projections.
  ProjectionResult projection_result = ProjectionResult::none;
  LINALG::Matrix<3, 1, scalar_type> point_in_space;

  // If segments are found, convert them to FAD segments.
  segments.clear();
  for (auto& segment_double : segments_double)
  {
    // Create the segment with the scalar FAD type.
    segments.push_back(LineSegment<scalar_type>());
    LineSegment<scalar_type>& new_segment = segments.back();

    // Add the projection point to an array.
    std::array<std::reference_wrapper<ProjectionPoint1DTo3D<scalar_type>>, 2>
        segment_start_end_points = {
            new_segment.GetStartPointMutable(), new_segment.GetEndPointMutable()};
    segment_start_end_points[0].get().SetFromOtherPointDouble(segment_double.GetStartPoint());
    segment_start_end_points[1].get().SetFromOtherPointDouble(segment_double.GetEndPoint());

    // If the start or end points are intersection points, the intersections have to be reevaluated.
    for (auto& point : segment_start_end_points)
    {
      const int intersection_face = point.get().GetIntersectionFace();
      if (intersection_face >= 0)
      {
        this->IntersectLineWithSurfaceEdge(q_line, q_surface,
            face_fixed_parameters[intersection_face], face_fixed_values[intersection_face],
            point.get().GetEtaMutable(), point.get().GetXiMutable(), projection_result,
            nodal_normals, true);
      }
    }

    // Reevaluate the integration points along the segment.
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& projection_points =
        new_segment.GetProjectionPointsMutable();
    projection_points.resize(segment_double.GetNumberOfProjectionPoints());
    for (unsigned int i_point = 0; i_point < segment_double.GetNumberOfProjectionPoints();
         i_point++)
    {
      // Position of the projection point within the segment.
      auto& projection_point_double = segment_double.GetProjectionPoints()[i_point];
      const double factor = (projection_point_double.GetEta() - segment_double.GetEtaA()) /
                            segment_double.GetSegmentLength();

      // Calculate spatial point.
      auto& projection_point = projection_points[i_point];
      projection_point.SetFromOtherPointDouble(projection_point_double);
      projection_point.SetEta(new_segment.GetEtaA() + new_segment.GetSegmentLength() * factor);

      // EvaluatePosition<line>(eta, q_line, r_line, this->Element1());
      EvaluatePosition<line>(projection_point.GetEta(), q_line, point_in_space, this->Element1());

      // Calculate the projection.
      this->ProjectPointToOther(point_in_space, q_surface, projection_point.GetXiMutable(),
          projection_result, nodal_normals, true);
    }
  }
}


/**
 * Explicit template initialization of template class.
 */
namespace GEOMETRYPAIR
{
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
