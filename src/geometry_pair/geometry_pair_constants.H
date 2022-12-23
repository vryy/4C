/*----------------------------------------------------------------------*/
/*! \file

\brief Constants used in the geometry pair namespace.

\level 1
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_CONSTANTS_H_
#define GEOMETRY_PAIR_CONSTANTS_H_


namespace GEOMETRYPAIR
{
  namespace CONSTANTS
  {
    //! Parameter coordinates are considered equal if they are within this tolerance.
    const double projection_xi_eta_tol = 1e-8;

    //! Maximum number of iterations in the local Newton.
    const unsigned int local_newton_iter_max = 10;

    //! Residuum is considered converged if its absolute value is smaller than this tolerance.
    const double local_newton_res_tol = 1e-10;

    //! Newton iteration is considered non solvable if the residuum is larger than this value.
    const double local_newton_res_max = 1e+10;

    //! Jacobian is considered non invertible if it is smaller than this value.
    const double local_newton_det_tol = 1e-8;

    //! Tolerance for positions to be equal.
    const double pos_tol = 1e-10;
  }  // namespace CONSTANTS
}  // namespace GEOMETRYPAIR


#endif
