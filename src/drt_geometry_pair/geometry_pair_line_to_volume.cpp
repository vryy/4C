/*----------------------------------------------------------------------*/
/*! \file

\brief Class for interaction of lines and volumes.

\level 1
\maintainer Ivo Steinbrecher
*/


#include "geometry_pair_line_to_volume.H"
#include "geometry_pair_element_functions.H"
#include "geometry_pair_utility_classes.H"
#include "geometry_pair_constants.H"

#include "../linalg/linalg_utils_densematrix_inverse.H"
#include "../drt_lib/drt_dserror.H"


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, line, volume>::Init(
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
        "The GeometryPairLineToVolume pair has to be on the same processor as the line element! "
        "Currently the pair is on rank %d, the line element on %d!",
        myrank, element1->Owner());
}


/**
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, line, volume>::ProjectPointToOther(
    const LINALG::Matrix<3, 1, scalar_type>& point,
    const LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
    LINALG::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result) const
{
  // Initialize data structures

  // Point on volume.
  LINALG::Matrix<3, 1, scalar_type> r_volume;

  // Jacobian / inverse.
  LINALG::Matrix<3, 3, scalar_type> J_J_inv;

  // Increment of xi.
  LINALG::Matrix<3, 1, scalar_type> delta_xi;

  // Residuum.
  LINALG::Matrix<3, 1, scalar_type> residuum;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  // Local Newton iteration.
  {
    unsigned int counter = 0;
    while (counter < CONSTANTS::local_newton_iter_max)
    {
      // Get the point coordinates on the volume.
      GEOMETRYPAIR::EvaluatePosition<volume>(xi, q_volume, r_volume, Element2());

      // Evaluate the residuum $r_{volume} - r_{line} = R_{pos}$
      residuum = r_volume;
      residuum -= point;

      // Check if tolerance is fulfilled.
      if (residuum.Norm2() < CONSTANTS::local_newton_res_tol)
      {
        if (ValidParameter3D<volume>(xi))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (residuum.Norm2() > CONSTANTS::local_newton_res_max) break;

      // Get the jacobian.
      GEOMETRYPAIR::EvaluatePositionDerivative1<volume>(xi, q_volume, J_J_inv, Element2());

      // Invert the jacobian and check if the system is solvable.
      if (LINALG::InverseDoNotThrowErrorOnZeroDeterminant(J_J_inv, CONSTANTS::local_newton_det_tol))
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
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, line, volume>::IntersectLineWithSurface(
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
    const unsigned int& fixed_parameter, const double& fixed_value, scalar_type& eta,
    LINALG::Matrix<3, 1, scalar_type>& xi, ProjectionResult& projection_result) const
{
  // Check the input parameters.
  {
    if (volume::geometry_type_ == DiscretizationTypeGeometry::hexahedron && fixed_parameter > 2)
      dserror(
          "Fixed_parameter in IntersectLineWithVolume has to be smaller than 3 with a hexahedron "
          "element.");
    else if (volume::dim_ != 3)
      dserror("Wrong DiscretizationTypeGeometry type given.");
    else if (fixed_parameter > 3)
      dserror("fixed_parameter in IntersectLineWithVolume can be 3 at maximum.");
  }

  // Initialize data structures
  // Point on line.
  LINALG::Matrix<3, 1, scalar_type> r_line;
  LINALG::Matrix<3, 1, scalar_type> dr_line;

  // Point on volume.
  LINALG::Matrix<3, 1, scalar_type> r_volume;
  LINALG::Matrix<3, 3, scalar_type> dr_volume;

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
      // Get the point coordinates on the line and volume.
      EvaluatePosition<line>(eta, q_line, r_line, Element1());
      EvaluatePosition<volume>(xi, q_volume, r_volume, Element2());

      // Evaluate the residuum $r_{volume} - r_{line} = R_{pos}$ and $xi(i) - value = R_{surf}$
      J_J_inv.PutScalar(0.);
      residuum.PutScalar(0.);
      for (unsigned int i = 0; i < 3; i++)
      {
        residuum(i) = r_volume(i) - r_line(i);
      }
      if (fixed_parameter < 3)
      {
        residuum(3) = xi(fixed_parameter) - fixed_value;
        J_J_inv(3, fixed_parameter) = 1.;
      }
      else
      {
        for (unsigned int i = 0; i < 3; i++)
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
        if (ValidParameter1D(eta) && ValidParameter3D<volume>(xi))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Check if residuum is in a sensible range where we still expect to find a solution.
      if (residuum.Norm2() > CONSTANTS::local_newton_res_max) break;

      // Get the positional derivatives.
      EvaluatePositionDerivative1<line>(eta, q_line, dr_line, Element1());
      EvaluatePositionDerivative1<volume>(xi, q_volume, dr_volume, Element2());

      // Fill up the jacobian.
      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < 3; j++)
        {
          J_J_inv(i, j) = dr_volume(i, j);
        }
        J_J_inv(i, 3) = -dr_line(i);
      }

      // Invert the jacobian and check if the determinant is not 0.
      if (LINALG::InverseDoNotThrowErrorOnZeroDeterminant(J_J_inv, CONSTANTS::local_newton_det_tol))
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
 *
 */
template <typename scalar_type, typename line, typename volume>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, line, volume>::IntersectLineWithOther(
    const LINALG::Matrix<line::n_dof_, 1, scalar_type>& q_line,
    const LINALG::Matrix<volume::n_dof_, 1, scalar_type>& q_volume,
    std::vector<ProjectionPoint1DTo3D<scalar_type>>& intersection_points,
    const scalar_type& eta_start, const LINALG::Matrix<3, 1, scalar_type>& xi_start) const
{
  // Get number of faces for this volume and create a vector with the indices of the faces, so all
  // surfaces of the volume can be checked for an intersection with the line.
  unsigned int n_faces;
  std::vector<unsigned int> face_fixed_parameters;
  std::vector<double> face_fixed_values;
  if (volume::geometry_type_ == DiscretizationTypeGeometry::hexahedron)
  {
    n_faces = 6;
    face_fixed_parameters = {0, 0, 1, 1, 2, 2};
    face_fixed_values = {-1., 1., -1., 1., -1., 1.};
  }
  else if (volume::geometry_type_ == DiscretizationTypeGeometry::tetraeder)
  {
    n_faces = 4;
    face_fixed_parameters = {0, 1, 2, 3};
    face_fixed_values = {0., 0., 0., 1.};
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
    IntersectLineWithSurface(q_line, q_volume, face_fixed_parameters[i], face_fixed_values[i], eta,
        xi, intersection_found);

    // If a valid intersection is found, add it to the output vector.
    if (intersection_found == ProjectionResult::projection_found_valid)
    {
      intersection_points.push_back(ProjectionPoint1DTo3D<scalar_type>(eta, xi));
    }
  }
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex8>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex20>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_hex27>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet4>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_tet10>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, GEOMETRYPAIR::t_hermite,
    GEOMETRYPAIR::t_nurbs27>;
