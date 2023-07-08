/*---------------------------------------------------------------------*/
/*! \file

\brief Integrates base functions over a line using one-dimensional Gauss quadrature
equations

\level 2


*----------------------------------------------------------------------*/

#include "cut_line_integration.H"
#include "cut_tolerance.H"
#include "cut_base.H"
#include "cut_base_boundarycell.H"
#include "discretization_fem_general_utils_boundary_integration.H"

#include <cmath>
#include <iostream>

/*---------------------------------------------------------------------------------------------------------------------*
 *      Compute normal vector for the line. if normal(0)==0, the line need not be integrated
 *(Divergence theorem)       *
 *----------------------------------------------------------------------------------------------------------------------*/
CORE::LINALG::Matrix<2, 1> LineIntegration::compute_normal()
{
  CORE::LINALG::Matrix<2, 1> normal;
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
  CORE::LINALG::Matrix<2, 1> normal;
  normal = compute_normal();

  if (fabs(normal(0, 0)) < TOL_LINE_NORMAL) return 0.0;

  double inte = 0.0;

  // 8 is the order of Gauss integration used in the line integration
  // since we integrate 6th order polynomial in volume, 8th order must be used for line
  CORE::DRT::UTILS::GaussIntegration gi(DRT::Element::line2, (DIRECTDIV_GAUSSRULE + 1));

  for (CORE::DRT::UTILS::GaussIntegration::iterator iquad = gi.begin(); iquad != gi.end(); ++iquad)
  {
    const CORE::LINALG::Matrix<1, 1> eta(iquad.Point());
    double weight = iquad.Weight();
    CORE::LINALG::Matrix<2, 1> normaltemp, actCoord;
    double drs = 0.0;
    Transform(end_pts_, eta(0, 0), actCoord, normaltemp, drs);

    if (bcellInt_ == false)  // integration over volumecell
    {
      double linein = base_func_line_int(actCoord, inte_num_, alpha_);
      inte = inte + weight * linein * drs;
    }
    else  // integration over boundarycell
    {
      double linein = 0.0;
      if (intType_ == CORE::GEO::CUT::proj_x)
        linein = base_func_surfX(actCoord, inte_num_, alpha_);
      else if (intType_ == CORE::GEO::CUT::proj_y)
        linein = base_func_surfY(actCoord, inte_num_, alpha_);
      else if (intType_ == CORE::GEO::CUT::proj_z)
        linein = base_func_surfZ(actCoord, inte_num_, alpha_);
      else
        dserror("Integration type unspecified");
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
void LineIntegration::Transform(const CORE::LINALG::Matrix<2, 2> &xyze, const double &eta,
    CORE::LINALG::Matrix<2, 1> &x_gp_lin, CORE::LINALG::Matrix<2, 1> &normal, double &drs)
{
  const int numnodes =
      CORE::DRT::UTILS::DisTypeToNumNodePerEle<::DRT::Element::line2>::numNodePerElement;
  CORE::LINALG::Matrix<numnodes, 1> funct;
  CORE::LINALG::Matrix<1, numnodes> deriv;
  CORE::LINALG::Matrix<1, 1> metrictensor;

  CORE::DRT::UTILS::shape_function_1D(funct, eta, ::DRT::Element::line2);
  CORE::DRT::UTILS::shape_function_1D_deriv1(deriv, eta, ::DRT::Element::line2);
  CORE::DRT::UTILS::ComputeMetricTensorForBoundaryEle<::DRT::Element::line2>(
      xyze, deriv, metrictensor, drs, &normal);

  x_gp_lin.Multiply(xyze, funct);

  return;
}
