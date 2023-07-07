/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "bele_vele3.H"
#include "linalg_utils_sparse_algebra_math.H"
#include "lib_discret.H"



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            g.bau 07/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Vele3Line::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  return 0;
}  // DRT::ELEMENTS::Vele3Line::Evaluate



/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Vele3Line::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}



/*----------------------------------------------------------------------*
 |  Optimal Gauss rule (public)                              u.may 04/09|
 *----------------------------------------------------------------------*/
CORE::DRT::UTILS::GaussRule1D DRT::ELEMENTS::Vele3Line::getOptimalGaussrule(
    const DiscretizationType& distype)
{
  CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
  switch (distype)
  {
    case line2:
      rule = CORE::DRT::UTILS::GaussRule1D::line_2point;
      break;
    case line3:
      rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}
