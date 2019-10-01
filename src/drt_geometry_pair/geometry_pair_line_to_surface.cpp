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
    Teuchos::RCP<GEOMETRYPAIR::GeometryEvaluationDataGlobal> evaluation_data_ptr,
    const DRT::Element* element1, const DRT::Element* element2)
{
  // Call init of base class.
  GeometryPair::Init(evaluation_data_ptr, element1, element2);

  // For the current implementation, the line element has to be on the same processor as the pair
  // object. This is because the tracking vector in LineToVolumeEvaluationData is only local and we
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
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::ProjectPointToSurface(
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
template <typename scalar_type_evaluate>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::EvaluateSurfacePosition(
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
    const LINALG::Matrix<3, 1, scalar_type_evaluate>& xi,
    LINALG::Matrix<3, 1, scalar_type_evaluate>& r,
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const
{
  // Evaluate the normal.
  LINALG::Matrix<3, 1, scalar_type_evaluate> normal;
  EvaluateSurfaceNormal(q_surface, xi, normal, nodal_normals);

  // Evaluate the position on the surface.
  GEOMETRYPAIR::EvaluatePosition<surface>(xi, q_surface, r, Element2());

  // Add the normal part to the position.
  normal.Scale(xi(2));
  r += normal;
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
  EvaluateSurfacePosition(q_surface, xi_AD, r_AD, nodal_normals);

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
template <typename scalar_type_evaluate>
void GEOMETRYPAIR::GeometryPairLineToSurface<scalar_type, line, surface>::EvaluateSurfaceNormal(
    const LINALG::Matrix<surface::n_dof_, 1, scalar_type>& q_surface,
    const LINALG::Matrix<3, 1, scalar_type_evaluate>& xi,
    LINALG::Matrix<3, 1, scalar_type_evaluate>& normal,
    const LINALG::Matrix<3 * surface::n_nodes_, 1, scalar_type>* nodal_normals) const
{
  if (nodal_normals == NULL)
  {
    // Calculate the normal as the geometrical normal on the element.
    LINALG::Matrix<3, 2, scalar_type_evaluate> dr;
    LINALG::Matrix<3, 1, scalar_type_evaluate> dr_0;
    LINALG::Matrix<3, 1, scalar_type_evaluate> dr_1;
    GEOMETRYPAIR::EvaluatePositionDerivative1<surface>(xi, q_surface, dr, Element2());
    for (unsigned int i_dir = 0; i_dir < 3; i_dir++)
    {
      dr_0(i_dir) = dr(i_dir, 0);
      dr_1(i_dir) = dr(i_dir, 1);
    }
    normal.CrossProduct(dr_0, dr_1);
    normal.Scale(1.0 / FADUTILS::VectorNorm(normal));
  }
  else
  {
    // Calculate the normal as a interpolation of nodal normals.
    GEOMETRYPAIR::EvaluatePosition<surface>(xi, *nodal_normals, normal, Element2());
    normal.Scale(1.0 / FADUTILS::VectorNorm(normal));
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
