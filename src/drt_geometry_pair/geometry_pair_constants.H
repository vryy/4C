/*!
\file geometry_pair_constants.H

\brief Class for interaction of lines and volumes.

\level 3
\maintainer Ivo Steinbrecher
*/
// End doxygen header.


#ifndef GEOMETRY_PAIR_CONSTANTS_H_
#define GEOMETRY_PAIR_CONSTANTS_H_


namespace GEOMETRYPAIR
{
  //! Parameter coordinates are considered equal if they are within this tolerance.
  const double PROJECTION_XI_ETA_TOL = 1e-8;


  //! Maximum number of iterations in the local Newton.
  const unsigned int LOCALNEWTON_ITER_MAX = 10;

  //! Residuum is considered converged if its absolute value is smaller than this tolerance.
  const double LOCALNEWTON_RES_TOL = 1e-10;

  //! If the residuum is larger than this value, the local Newton is considered unsolvable.
  const double LOCALNEWTON_RES_MAX = 1e+10;

  //! Jacobian is considered non invertible if it is smaller than this value.
  const double LOCALNEWTON_DET_TOL = 1e-10;

  //! If the determinant is larger than this value, the local Newton is considered unsolvable.
  const double LOCALNEWTON_DET_MAX = 1e+10;
}  // namespace GEOMETRYPAIR


#endif
