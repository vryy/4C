/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and surfaces.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_surface.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_utility_classes.H"
#include "geometry_pair_constants.H"

#include "../linalg/linalg_utils_densematrix_inverse.H"
#include "../drt_lib/drt_dserror.H"


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
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const
{
  // Initialize data structures

  // Vectors in 3D.
  LINALG::Matrix<3, 1, scalar_type> r_surface;
  LINALG::Matrix<3, 1, scalar_type> delta_xi;
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

      // Check if tolerance is fulfilled.
      if (residuum.Norm2() < CONSTANTS::local_newton_res_tol)
      {
        if (ValidParameter2D<surface>(xi))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (residuum.Norm2() > CONSTANTS::local_newton_res_max) break;

      // Invert the jacobian and check if the system is solvable.
      if (LINALG::Inverse3x3DoNotThrowErrorOnZeroDeterminant(
              J_J_inv, CONSTANTS::local_newton_det_tol))
      {
        // Solve the linearized system.
        delta_xi.Multiply(J_J_inv, residuum);
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
  // Get number of faces for this volume and create a vector with the indices of the faces, so all
  // surfaces of the volume can be checked for an intersection with the line.
  unsigned int n_faces;
  std::vector<unsigned int> face_fixed_parameters;
  std::vector<double> face_fixed_values;
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
  EvaluateSurfacePosition<surface>(q_surface, xi_AD, r_AD, Element2(), nodal_normals);

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
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const
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

  // Initialize data structures.
  // Point on line.
  LINALG::Matrix<3, 1, scalar_type> r_line;
  LINALG::Matrix<3, 1, scalar_type> dr_line;

  // Point on surface.
  LINALG::Matrix<3, 1, scalar_type> r_surface;
  LINALG::Matrix<3, 3, scalar_type> dr_surface;

  // Residuum.
  LINALG::Matrix<4, 1, scalar_type> residuum;
  LINALG::Matrix<4, 1, scalar_type> delta_x;

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

      // Check if tolerance is fulfilled.
      if (residuum.Norm2() < CONSTANTS::local_newton_res_tol)
      {
        // Check if the parameter coordinates are valid.
        if (ValidParameter1D(eta) && ValidParameter2D<surface>(xi))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (residuum.Norm2() > CONSTANTS::local_newton_res_max) break;

      // Fill up the jacobian.
      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < 3; j++)
        {
          J_J_inv(i, j) = dr_surface(i, j);
        }
        J_J_inv(i, 3) = -dr_line(i);
      }

      // Invert the jacobian and check if the determinant is not 0.
      if (LINALG::Inverse4x4DoNotThrowErrorOnZeroDeterminant(
              J_J_inv, CONSTANTS::local_newton_det_tol))
      {
        // Solve the linearized system.
        delta_x.Multiply(J_J_inv, residuum);

        // Set the new parameter coordinates.
        eta -= delta_x(3);
        for (unsigned int i = 0; i < 3; i++) xi(i) -= delta_x(i);

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
template class GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri3>;
template class GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tri6>;
template class GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad4>;
template class GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad8>;
template class GEOMETRYPAIR::GeometryPairLineToSurface<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_quad9>;
