/*---------------------------------------------------------------------*/
/*! \file

\brief Integrates base functions over a line using one-dimensional Gauss quadrature
equations

\level 2


*----------------------------------------------------------------------*/

#include "4C_cut_line_integration.hpp"

#include "4C_cut_base.hpp"
#include "4C_cut_base_boundarycell.hpp"
#include "4C_cut_tolerance.hpp"
#include "4C_discretization_fem_general_utils_boundary_integration.hpp"

#include <cmath>
#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------------------------------------------------*
 *      Compute normal vector for the line. if normal(0)==0, the line need not be integrated
 *(Divergence theorem)       *
 *----------------------------------------------------------------------------------------------------------------------*/
Core::LinAlg::Matrix<2, 1> LineIntegration::compute_normal()
{
  Core::LinAlg::Matrix<2, 1> normal;
  double dy = end_pts_(1, 1) - end_pts_(1, 0);
  double dx = -end_pts_(0, 1) + end_pts_(0, 0);
  double modd = sqrt(dx * dx + dy * dy);

  normal(0, 0) = dy / modd;
  normal(1, 0) = dx / modd;

  return normal;
}

/*-----------------------------------------------------------------------------*
 *                performs integration over the given line                      *
 *------------------------------------------------------------------------------*/
double LineIntegration::integrate_line()
{
  Core::LinAlg::Matrix<2, 1> normal;
  normal = compute_normal();

  if (fabs(normal(0, 0)) < TOL_LINE_NORMAL) return 0.0;

  double inte = 0.0;

  // 8 is the order of Gauss integration used in the line integration
  // since we integrate 6th order polynomial in volume, 8th order must be used for line
  Core::FE::GaussIntegration gi(Core::FE::CellType::line2, (DIRECTDIV_GAUSSRULE + 1));

  for (Core::FE::GaussIntegration::iterator iquad = gi.begin(); iquad != gi.end(); ++iquad)
  {
    const Core::LinAlg::Matrix<1, 1> eta(iquad.Point());
    double weight = iquad.Weight();
    Core::LinAlg::Matrix<2, 1> normaltemp, actCoord;
    double drs = 0.0;
    Transform(end_pts_, eta(0, 0), actCoord, normaltemp, drs);

    if (bcell_int_ == false)  // integration over volumecell
    {
      double linein = base_func_line_int(actCoord, inte_num_, alpha_);
      inte = inte + weight * linein * drs;
    }
    else  // integration over boundarycell
    {
      double linein = 0.0;
      if (int_type_ == Core::Geo::Cut::proj_x)
        linein = base_func_surfX(actCoord, inte_num_, alpha_);
      else if (int_type_ == Core::Geo::Cut::proj_y)
        linein = base_func_surfY(actCoord, inte_num_, alpha_);
      else if (int_type_ == Core::Geo::Cut::proj_z)
        linein = base_func_surfZ(actCoord, inte_num_, alpha_);
      else
        FOUR_C_THROW("Integration type unspecified");
      inte = inte + weight * linein * drs;
    }
  }
  inte = inte * normal(0, 0);

  return inte;
}

/*---------------------------------------------------------------------------------------------------------------------*
 *     Transform the Gaussian point coordinates and weight from (-1,1) interval to actual coordinate
 *of the lines       *
 *----------------------------------------------------------------------------------------------------------------------*/
void LineIntegration::Transform(const Core::LinAlg::Matrix<2, 2> &xyze, const double &eta,
    Core::LinAlg::Matrix<2, 1> &x_gp_lin, Core::LinAlg::Matrix<2, 1> &normal, double &drs)
{
  const int numnodes = Core::FE::num_nodes<Core::FE::CellType::line2>;
  Core::LinAlg::Matrix<numnodes, 1> funct;
  Core::LinAlg::Matrix<1, numnodes> deriv;
  Core::LinAlg::Matrix<1, 1> metrictensor;

  Core::FE::shape_function_1D(funct, eta, Core::FE::CellType::line2);
  Core::FE::shape_function_1D_deriv1(deriv, eta, Core::FE::CellType::line2);
  Core::FE::ComputeMetricTensorForBoundaryEle<Core::FE::CellType::line2>(
      xyze, deriv, metrictensor, drs, &normal);

  x_gp_lin.Multiply(xyze, funct);

  return;
}

FOUR_C_NAMESPACE_CLOSE
