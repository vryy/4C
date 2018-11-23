/*!
\file geometry_pair_line_to_volume.cpp

\brief Class for interaction of lines and volumes.

<pre>
\level 3
\maintainer Ivo Steinbrecher
            ivo.steinbrecher@unibw.de
            +49 89 6004-4403
</pre>
*/


#include "geometry_pair_line_to_volume.H"
#include "geometry_pair_utility_classes.H"

#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_beam3/beam3.H"


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
template <typename scalar_type_get_pos>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2,
    n_nodal_values_element_2>::GetElement1Position(const scalar_type& xi,
    const LINALG::TMatrix<scalar_type_get_pos, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
        q,
    LINALG::TMatrix<scalar_type_get_pos, 3, 1>& r) const
{
  // Matrix for shape function values.
  LINALG::TMatrix<scalar_type, 1, n_nodes_element_1 * n_nodal_values_element_1> N(true);

  // Get discretization type.
  const DRT::Element::DiscretizationType distype = Element1()->Shape();

  if (n_nodal_values_element_1 == 1)
  {
    dserror("One nodal value for line elements not yet implemented!");
    DRT::UTILS::shape_function_1D(N, xi, distype);
  }
  else if (n_nodal_values_element_1 == 2)
  {
    double length = (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(Element1()))->RefLength();
    const DRT::Element::DiscretizationType distype1herm = DRT::Element::line2;

    // Get values of shape functions.
    DRT::UTILS::shape_function_hermite_1D(N, xi, length, distype1herm);

    // Calculate the position.
    r.Clear();
    for (unsigned int dim = 0; dim < 3; dim++)
    {
      for (unsigned int node = 0; node < n_nodes_element_1; node++)
      {
        for (unsigned int val = 0; val < n_nodal_values_element_1; val++)
        {
          r(dim) += q(3 * n_nodal_values_element_1 * node + 3 * val + dim) *
                    N(n_nodal_values_element_1 * node + val);
        }
      }
    }
  }
  else
    dserror(
        "Only line elements with one (nodal positions) or two "
        "(nodal positions + nodal tangents) values are valid!");
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2,
    n_nodal_values_element_2>::GetElement1PositionDerivative(const scalar_type& xi,
    const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>& q,
    LINALG::TMatrix<scalar_type, 3, 1>& dr) const
{
  // Matrix for shape function values.
  LINALG::TMatrix<scalar_type, 1, 3 * n_nodes_element_1 * n_nodal_values_element_1> dN(true);

  // Get discretization type.
  const DRT::Element::DiscretizationType distype = Element1()->Shape();

  if (n_nodal_values_element_1 == 1)
  {
    dserror("One nodal value for line elements not yet implemented!");
    DRT::UTILS::shape_function_1D_deriv1(dN, xi, distype);
  }
  else if (n_nodal_values_element_1 == 2)
  {
    double length = (dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(Element1()))->RefLength();
    const DRT::Element::DiscretizationType distype1herm = DRT::Element::line2;

    // Get values of shape functions.
    DRT::UTILS::shape_function_hermite_1D_deriv1(dN, xi, length, distype1herm);

    // Calculate the position.
    dr.Clear();
    for (unsigned int dim = 0; dim < 3; dim++)
    {
      for (unsigned int node = 0; node < n_nodes_element_1; node++)
      {
        for (unsigned int val = 0; val < n_nodal_values_element_1; val++)
        {
          dr(dim) += q(3 * n_nodal_values_element_1 * node + 3 * val + dim) *
                     dN(n_nodal_values_element_1 * node + val);
        }
      }
    }
  }
  else
    dserror(
        "Only line elements with one (nodal positions) or two "
        "(nodal positions + nodal tangents) values are valid!");
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
template <typename scalar_type_get_pos>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2,
    n_nodal_values_element_2>::GetElement2Position(const LINALG::TMatrix<scalar_type, 3, 1>& xi,
    const LINALG::TMatrix<scalar_type_get_pos, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
        q,
    LINALG::TMatrix<scalar_type_get_pos, 3, 1>& r) const
{
  // Matrix for shape function values.
  LINALG::TMatrix<scalar_type, 1, n_nodes_element_2 * n_nodal_values_element_2> N(true);

  // Check what type of volume was given.
  if (n_nodal_values_element_2 != 1)
    dserror("Only volume elements with one nodal values are implemented!");

  // Clear shape function matrix.
  N.Clear();

  // Get the shape functions.
  DRT::UTILS::shape_function_3D(N, xi(0), xi(1), xi(2), Element2()->Shape());

  // Calculate the position.
  r.Clear();
  for (unsigned int dim = 0; dim < 3; dim++)
  {
    for (unsigned int node = 0; node < n_nodes_element_2; node++)
    {
      r(dim) += q(3 * node + dim) * N(node);
    }
  }
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>::
    GetElement2PositionDerivative(const LINALG::TMatrix<scalar_type, 3, 1>& xi,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>& q,
        LINALG::TMatrix<scalar_type, 3, 3>& dr) const
{
  // Matrix for shape function values.
  LINALG::TMatrix<scalar_type, 3, n_nodes_element_2 * n_nodal_values_element_2> dN(true);

  // Check what type of volume was given.
  if (n_nodal_values_element_2 != 1)
    dserror("Only volume elements with one nodal values are implemented!");

  // Clear shape function matrix.
  dN.Clear();

  // Get the shape functions.
  DRT::UTILS::shape_function_3D_deriv1(dN, xi(0), xi(1), xi(2), Element2()->Shape());

  // Calculate the position.
  dr.Clear();
  for (unsigned int dim = 0; dim < 3; dim++)
  {
    for (unsigned int direction = 0; direction < 3; direction++)
    {
      for (unsigned int node = 0; node < n_nodes_element_2; node++)
      {
        dr(dim, direction) += q(3 * node + dim) * dN(direction, node);
      }
    }
  }
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>::
    ProjectPointOnLineToVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        const scalar_type& eta, LINALG::TMatrix<scalar_type, 3, 1>& xi,
        ProjectionResult& projection_result) const
{
  // Initialize data structures
  // Point on line.
  LINALG::TMatrix<scalar_type, 3, 1> r_line;

  // Point on volume.
  LINALG::TMatrix<scalar_type, 3, 1> r_volume;

  // Jacobian.
  LINALG::TMatrix<scalar_type, 3, 3> J;
  LINALG::TMatrix<scalar_type, 3, 3> J_inverse;

  // Increment of xi.
  LINALG::TMatrix<scalar_type, 3, 1> delta_xi;

  // Residuum.
  LINALG::TMatrix<scalar_type, 3, 1> residuum;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  // Local Newton iteration.
  {
    // Get the position on the beam that the solid should match.
    GetElement1Position(eta, q_line, r_line);

    unsigned int counter = 0;
    while (counter < LOCALNEWTON_MAX_ITER)
    {
      // Get the point coordinates on the volume.
      GetElement2Position(xi, q_volume, r_volume);

      // Evaluate the residuum $r_{volume} - r_{line} = R_{pos}$
      residuum = r_volume;
      residuum -= r_line;

      // Check if tolerance is fulfilled.
      if (residuum.Norm2() < LOCALNEWTON_RES_TOL)
      {
        // We only check xi, as eta is given by the user and is assumed to be correct.
        if (ValidParameterElement2(xi))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Get the jacobian.
      GetElement2PositionDerivative(xi, q_volume, J);

      // Check the determinant of the jacobian.
      if (J.Determinant() < LOCALNEWTON_DET_TOL) break;

      // Solve the linearized system.
      J_inverse.Invert(J);
      delta_xi.Multiply(J_inverse, residuum);
      xi -= delta_xi;

      // Advance Newton iteration counter.
      counter++;
    }
  }
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>::
    ProjectPointsOnLineToVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        const std::vector<scalar_type>& eta_points,
        std::vector<ProjectionPointLineToVolume<scalar_type>>& projection_points,
        bool allow_not_valid_projections) const
{
  // Initialize variables for the projection.
  LINALG::TMatrix<scalar_type, 3, 1> xi;
  ProjectionResult projection_result;
  bool last_point_valid = false;
  projection_points.clear();

  // Loop over points and check if they project to this volume.
  for (auto const& eta : eta_points)
  {
    // If the last projected point was valid, use that one as start point for this itereation.
    if (!last_point_valid) SetStartValuesElement2(xi);

    // Project the point.
    ProjectPointOnLineToVolume(q_line, q_volume, eta, xi, projection_result);

    // Check the result from projection.
    if (projection_result == ProjectionResult::projection_found_valid ||
        (projection_result == ProjectionResult::projection_found_not_valid &&
            allow_not_valid_projections))
    {
      // Add point to the projection points vector.
      projection_points.push_back(ProjectionPointLineToVolume<scalar_type>(eta, xi));

      last_point_valid = true;
    }
    else
    {
      last_point_valid = false;
    }
  }
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>::
    ProjectPointsOnSegmentToVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        unsigned int n_points, LineSegment<scalar_type>& segment,
        bool allow_not_valid_projections) const
{
  // Create vector with eta values.
  std::vector<scalar_type> eta_values(n_points);
  for (unsigned int i = 0; i < n_points; i++)
    eta_values[i] =
        segment.GetEtaA() + i * (segment.GetEtaB() - segment.GetEtaA()) / (n_points - 1.);

  // Project the points on the segment.
  ProjectPointsOnLineToVolume(
      q_line, q_volume, eta_values, segment.GetGaussPointsMutable(), allow_not_valid_projections);
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>::
    ProjectPointsOnSegmentToVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        const DRT::UTILS::IntegrationPoints1D& gauss_points,
        LineSegment<scalar_type>& segment) const
{
  // Create vector with eta values.
  std::vector<scalar_type> eta_values(gauss_points.nquad);
  for (unsigned int i = 0; i < (unsigned int)gauss_points.nquad; i++)
    eta_values[i] = segment.GetEtaA() +
                    (segment.GetEtaB() - segment.GetEtaA()) * 0.5 * (gauss_points.qxg[i][0] + 1.);

  // Project the points on the segment.
  ProjectPointsOnLineToVolume(q_line, q_volume, eta_values, segment.GetGaussPointsMutable(), false);

  // Check if all points could be projected.
  if (segment.GetGaussPointsMutable().size() != (unsigned int)gauss_points.nquad)
    dserror(
        "All Gauss points need to have a valid projection. The number of Gauss points is %d, but "
        "the number of valid projections is %d!",
        gauss_points.nquad, segment.GetGaussPointsMutable().size());

  // Set the Gauss weights for the found points.
  for (unsigned int i = 0; i < (unsigned int)gauss_points.nquad; i++)
  {
    segment.GetGaussPointsMutable()[i].SetGaussWeight(gauss_points.qwgt[i]);
  }
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>::
    IntersectLineWithSurface(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        const unsigned int& fixed_parameter, const double& fixed_value, scalar_type& eta,
        LINALG::TMatrix<scalar_type, 3, 1>& xi, ProjectionResult& projection_result) const
{
  // Check the input parameters.
  {
    if (GetVolumeType() == DiscretizationTypeVolume::hexaeder && fixed_parameter > 2)
      dserror(
          "Fixed_parameter in IntersectLineWithVolume has to be smaller than 3 with a hexaeder "
          "element.");
    else if (fixed_parameter > 3)
      dserror("fixed_parameter in IntersectLineWithVolume can be 3 at maximum.");
  }

  // Initialize data structures
  // Point on line.
  LINALG::TMatrix<scalar_type, 3, 1> r_line;
  LINALG::TMatrix<scalar_type, 3, 1> dr_line;

  // Point on volume.
  LINALG::TMatrix<scalar_type, 3, 1> r_volume;
  LINALG::TMatrix<scalar_type, 3, 3> dr_volume;

  // Residuum.
  LINALG::TMatrix<scalar_type, 4, 1> residuum;
  LINALG::TMatrix<scalar_type, 4, 1> delta_x;

  // Jacobian.
  LINALG::TMatrix<scalar_type, 4, 4> J;
  LINALG::TMatrix<scalar_type, 4, 4> J_inverse;

  // Solver.
  LINALG::FixedSizeSerialDenseSolver<4, 4> matrix_solver;

  // Reset the projection result flag.
  projection_result = ProjectionResult::projection_not_found;

  {
    // Local Newton iteration.
    unsigned int counter = 0;
    while (counter < LOCALNEWTON_MAX_ITER)
    {
      // Get the point coordinates on the line and volume.
      GetElement1Position(eta, q_line, r_line);
      GetElement2Position(xi, q_volume, r_volume);

      // Evaluate the residuum $r_{volume} - r_{line} = R_{pos}$ and $xi(i) - value = R_{surf}$
      J.PutScalar(0.);
      residuum.PutScalar(0.);
      for (unsigned int i = 0; i < 3; i++)
      {
        residuum(i) = r_volume(i) - r_line(i);
      }
      if (fixed_parameter < 3)
      {
        residuum(3) = xi(fixed_parameter) - fixed_value;
        J(3, fixed_parameter) = 1.;
      }
      else
      {
        for (unsigned int i = 0; i < 3; i++)
        {
          residuum(3) += xi(fixed_parameter);
          J(3, fixed_parameter) = 1.;
        }
        residuum(3) -= fixed_value;
      }

      // Check if tolerance is fulfilled.
      if (residuum.Norm2() < LOCALNEWTON_RES_TOL)
      {
        // Check if the parameter coordinates are valid.
        if (ValidParameterElement1(eta) && ValidParameterElement2(xi))
          projection_result = ProjectionResult::projection_found_valid;
        else
          projection_result = ProjectionResult::projection_found_not_valid;
        break;
      }

      // Get the positional derivatives.
      GetElement1PositionDerivative(eta, q_line, dr_line);
      GetElement2PositionDerivative(xi, q_volume, dr_volume);

      // Fill up the jacobian.
      for (unsigned int i = 0; i < 3; i++)
      {
        for (unsigned int j = 0; j < 3; j++)
        {
          J(i, j) = dr_volume(i, j);
        }
        J(i, 3) = -dr_line(i);
      }

      // Solve the linearized system.
      if (abs(J.Determinant()) < 1e-10) break;
      matrix_solver.SetMatrix(J);
      matrix_solver.SetVectors(delta_x, residuum);
      int err = matrix_solver.Solve();
      if (err != 0) break;

      // Set the new parameter coordinates.
      eta -= delta_x(3);
      for (unsigned int i = 0; i < 3; i++) xi(i) -= delta_x(i);

      // Advance Newton iteration counter.
      counter++;
    }
  }
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>::
    IntersectLineWithVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        std::vector<ProjectionPointLineToVolume<scalar_type>>& intersection_points,
        const scalar_type& eta_start, const LINALG::TMatrix<scalar_type, 3, 1>& xi_start) const
{
  // Get number of faces for this volume and create a vector with the indices of the faces, so all
  // surfaces of the volume can be checked for an intersection with the line.
  unsigned int n_faces;
  std::vector<unsigned int> face_fixed_parameters;
  std::vector<double> face_fixed_values;
  if (GetVolumeType() == DiscretizationTypeVolume::hexaeder)
  {
    n_faces = 6;
    face_fixed_parameters = {0, 0, 1, 1, 2, 2};
    face_fixed_values = {-1., 1., -1., 1., -1., 1.};
  }
  else
  {
    n_faces = 4;
    face_fixed_parameters = {0, 1, 2, 3};
    face_fixed_values = {0., 0., 0., 1.};
  }

  // Clear the input vector.
  intersection_points.clear();
  intersection_points.reserve(n_faces);

  // Create variables.
  scalar_type eta;
  LINALG::TMatrix<scalar_type, 3, 1> xi;
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
      intersection_points.push_back(ProjectionPointLineToVolume<scalar_type>(eta, xi));
    }
  }
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2, n_nodal_values_element_2>::
    IntersectLineWithVolume(
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_1 * n_nodal_values_element_1, 1>&
            q_line,
        const LINALG::TMatrix<scalar_type, 3 * n_nodes_element_2 * n_nodal_values_element_2, 1>&
            q_volume,
        std::vector<ProjectionPointLineToVolume<scalar_type>>& intersection_points) const
{
  // Set default values for the parameter coordinates.
  scalar_type eta_start;
  LINALG::TMatrix<scalar_type, 3, 1> xi_start;
  SetStartValuesElement1(eta_start);
  SetStartValuesElement2(xi_start);

  // Call the intersect function.
  IntersectLineWithVolume(q_line, q_volume, intersection_points, eta_start, xi_start);
};

/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
GEOMETRYPAIR::DiscretizationTypeVolume
GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1, n_nodal_values_element_1,
    n_nodes_element_2, n_nodal_values_element_2>::GetVolumeType() const
{
  if (n_nodes_element_2 == 8 || n_nodes_element_2 == 20 || n_nodes_element_2 == 27)
    return GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder;
  else if (n_nodes_element_2 == 4 || n_nodes_element_2 == 10)
    return GEOMETRYPAIR::DiscretizationTypeVolume::tetraeder;
  else
    dserror("Unknown volume type in GetVolumeType()!");
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
bool GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2,
    n_nodal_values_element_2>::ValidParameterElement1(const scalar_type& eta) const
{
  double xi_limit = 1.0 + PROJECTION_XI_ETA_TOL;
  if (fabs(eta) < xi_limit) return true;

  // Default value.
  return false;
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
bool GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2,
    n_nodal_values_element_2>::ValidParameterElement2(const LINALG::TMatrix<scalar_type, 3, 1>& xi)
    const
{
  double xi_limit = 1.0 + PROJECTION_XI_ETA_TOL;
  if (GetVolumeType() == DiscretizationTypeVolume::hexaeder)
  {
    if (fabs(xi(0)) < xi_limit && fabs(xi(1)) < xi_limit && fabs(xi(2)) < xi_limit) return true;
  }
  else
  {
    if (xi(0) > -PROJECTION_XI_ETA_TOL && xi(1) > -PROJECTION_XI_ETA_TOL &&
        xi(2) > -PROJECTION_XI_ETA_TOL && xi(0) + xi(1) + xi(2) < 1.0 + PROJECTION_XI_ETA_TOL)
      return true;
  }

  // Default value.
  return false;
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2,
    n_nodal_values_element_2>::SetStartValuesElement1(scalar_type& eta) const
{
  eta = 0.;
}


/**
 *
 */
template <typename scalar_type, unsigned int n_nodes_element_1,
    unsigned int n_nodal_values_element_1, unsigned int n_nodes_element_2,
    unsigned int n_nodal_values_element_2>
void GEOMETRYPAIR::GeometryPairLineToVolume<scalar_type, n_nodes_element_1,
    n_nodal_values_element_1, n_nodes_element_2,
    n_nodal_values_element_2>::SetStartValuesElement2(LINALG::TMatrix<scalar_type, 3, 1>& xi) const
{
  if (GetVolumeType() == GEOMETRYPAIR::DiscretizationTypeVolume::hexaeder)
    xi.PutScalar(0.0);
  else
    xi.PutScalar(0.25);
}


/**
 * Explicit template initialization of template class.
 */
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 8, 1>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 20, 1>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 27, 1>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 4, 1>;
template class GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 10, 1>;


/**
 * We need to explicitly initialize the position functions for AD types. For example in case of beam
 * to solid meshtying the geometry interactions are done with the constant reference configuration
 * and therefore doubles, but in the Evaluate function the position needs to be evaluated with AD
 * types to get the difference in the current configuration.
 */
template void GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 8, 1>::GetElement1Position(
    const double&,
    const LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 8>, 3 * 2 * 2, 1>&,
    LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 8>, 3, 1>&) const;
template void GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 20, 1>::GetElement1Position(
    const double&,
    const LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 20>, 3 * 2 * 2, 1>&,
    LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 20>, 3, 1>&) const;
template void GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 27, 1>::GetElement1Position(
    const double&,
    const LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 27>, 3 * 2 * 2, 1>&,
    LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 27>, 3, 1>&) const;
template void GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 4, 1>::GetElement1Position(
    const double&,
    const LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 4>, 3 * 2 * 2, 1>&,
    LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 4>, 3, 1>&) const;
template void GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 10, 1>::GetElement1Position(
    const double&,
    const LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 10>, 3 * 2 * 2, 1>&,
    LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 10>, 3, 1>&) const;

template void GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 8, 1>::GetElement2Position(
    const LINALG::TMatrix<double, 3, 1>&,
    const LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 8>, 3 * 8, 1>&,
    LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 8>, 3, 1>&) const;
template void GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 20, 1>::GetElement2Position(
    const LINALG::TMatrix<double, 3, 1>&,
    const LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 20>, 3 * 20, 1>&,
    LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 20>, 3, 1>&) const;
template void GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 27, 1>::GetElement2Position(
    const LINALG::TMatrix<double, 3, 1>&,
    const LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 27>, 3 * 27, 1>&,
    LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 27>, 3, 1>&) const;
template void GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 4, 1>::GetElement2Position(
    const LINALG::TMatrix<double, 3, 1>&,
    const LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 4>, 3 * 4, 1>&,
    LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 4>, 3, 1>&) const;
template void GEOMETRYPAIR::GeometryPairLineToVolume<double, 2, 2, 10, 1>::GetElement2Position(
    const LINALG::TMatrix<double, 3, 1>&,
    const LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 10>, 3 * 10, 1>&,
    LINALG::TMatrix<Sacado::ELRFad::SLFad<double, 3 * 2 * 2 + 3 * 10>, 3, 1>&) const;
