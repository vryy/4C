/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "baci_bele_vele3.H"
#include "baci_lib_discret.H"
#include "baci_linalg_utils_sparse_algebra_math.H"



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 07/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Vele3Line::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  return 0;
}  // DRT::ELEMENTS::Vele3Line::Evaluate



/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Vele3Line::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  return 0;
}



/*----------------------------------------------------------------------*
 |  Optimal Gauss rule (public)                              u.may 04/09|
 *----------------------------------------------------------------------*/
CORE::DRT::UTILS::GaussRule1D DRT::ELEMENTS::Vele3Line::getOptimalGaussrule(
    const CORE::FE::CellType& distype)
{
  CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
  switch (distype)
  {
    case CORE::FE::CellType::line2:
      rule = CORE::DRT::UTILS::GaussRule1D::line_2point;
      break;
    case CORE::FE::CellType::line3:
      rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}
