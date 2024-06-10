/*---------------------------------------------------------------------*/
/*! \file

\brief Integrates base functions over a line using one-dimensional Gauss quadrature
equations

\level 2


*----------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_LINE_INTEGRATION_HPP
#define FOUR_C_CUT_LINE_INTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_cut_enum.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_linalg_fixedsizematrix.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*!
\brief Performs integration of  base functions along the line using the standard Gaussian rule
*/
class LineIntegration
{
 public:
  LineIntegration(
      Core::LinAlg::Matrix<2, 2> endPts, int inte_num, std::vector<double> alpha, bool bcellInt)
      : end_pts_(endPts), inte_num_(inte_num), alpha_(alpha), bcell_int_(bcellInt)
  {
  }

  /*!
  \brief Integration of a function along the line using standard Gaussian rule
  */
  double integrate_line();

  /*!
  \brief Choose the base function to be integrated
  */
  void set_integ_type(Core::Geo::Cut::ProjectionDirection inttype) { int_type_ = inttype; }

  /*!
  \brief Transform the Gauss integration point available in the limit (-1,1) to the actual line
  coordinates
   */
  void Transform(const Core::LinAlg::Matrix<2, 2> &xyze, const double &eta,
      Core::LinAlg::Matrix<2, 1> &x_gp_lin, Core::LinAlg::Matrix<2, 1> &normal, double &drs);

 private:
  /*!
  \brief Compute the normal vector of the line
  */
  Core::LinAlg::Matrix<2, 1> compute_normal();

  //! end point of the line
  // first index decides the x or y coordinate, second index decides the start point or end point
  Core::LinAlg::Matrix<2, 2> end_pts_;

  //! defines the base function to be integrated
  int inte_num_;

  //! from equation of plane x=alpha0+alpha1*y+alpha2*z
  std::vector<double> alpha_;

  //! whether this is a boundary cell integration or called for volumecell integration
  bool bcell_int_;

  //! over which plane (x, y or z) boundarycell has to be projected when performing boundarycell
  //! integration
  Core::Geo::Cut::ProjectionDirection int_type_;
};

FOUR_C_NAMESPACE_CLOSE

#endif
