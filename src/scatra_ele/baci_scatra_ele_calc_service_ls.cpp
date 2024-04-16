/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for level sets

\level 2

*/
/*--------------------------------------------------------------------------*/

#include "baci_discretization_fem_general_extract_values.hpp"
#include "baci_discretization_geometry_position_array.hpp"
#include "baci_inpar_levelset.hpp"
#include "baci_lib_discret.hpp"
#include "baci_scatra_ele_action.hpp"
#include "baci_scatra_ele_calc_ls.hpp"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ScaTraEleCalcLS<distype>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  //(for now) only first dof set considered
  const std::vector<int>& lm = la[0].lm_;

  // determine and evaluate action
  switch (action)
  {
    case SCATRA::Action::calc_error:
    {
      // extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> phizero = discretization.GetState("phiref");
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phizero == Teuchos::null or phinp == Teuchos::null)
        dserror("Cannot get state vector 'phizero' and/ or 'phinp'!");

      std::vector<CORE::LINALG::Matrix<nen_, 1>> ephizero(my::numscal_);

      CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phinp, my::ephinp_, lm);
      CORE::FE::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phizero, ephizero, lm);

      // check if length suffices
      if (elevec1_epetra.length() < 1) dserror("Result vector too short");

      CalErrorComparedToAnalytSolution(ele, ephizero, params, elevec1_epetra);

      break;
    }
    default:
    {
      my::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
  }  // switch(action)

  return 0;
}


/*---------------------------------------------------------------------*
 |  calculate error compared to analytical solution    rasthofer 04/14 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLS<distype>::CalErrorComparedToAnalytSolution(
    const DRT::Element* ele, const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephizero,
    Teuchos::ParameterList& params, CORE::LINALG::SerialDenseVector& errors)
{
  // get element volume
  const double vol = my::EvalShapeFuncAndDerivsAtEleCenter();

  // get elemet length
  // cast dimension to a double varibale -> pow()
  const double dim = static_cast<double>(nsd_);
  const double h = std::pow(vol, 1.0 / dim);

  // integration points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  CORE::FE::IntPointsAndWeights<nsd_> intpoints(
      SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  const INPAR::SCATRA::CalcErrorLevelSet errortype =
      CORE::UTILS::GetAsEnum<INPAR::SCATRA::CalcErrorLevelSet>(params, "calcerrorflag");
  switch (errortype)
  {
    case INPAR::SCATRA::calcerror_initial_field:
    {
      // start loop over integration points
      for (int iquad = 0; iquad < intpoints.IP().nquad; iquad++)
      {
        const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        const double phizero =
            my::funct_.Dot(ephizero[0]);  // only one scalar available for level set
        double smoothH_exact = 0.0;
        SmoothHeavisideFunction(h, phizero, smoothH_exact);

        const double phinp =
            my::funct_.Dot(my::ephinp_[0]);  // only one scalar available for level set
        double smoothH = 0.0;
        SmoothHeavisideFunction(h, phinp, smoothH);

        // add error
        errors[0] += std::abs(smoothH_exact - smoothH) * fac;
      }  // end of loop over integration points

      errors[1] += vol;
    }
    break;
    default:
      dserror("Unknown analytical solution!");
      break;
  }  // switch(errortype)

  return;
}


/*----------------------------------------------------------------------*
 | smoothed heaviside function                          rasthofer 04/12 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLS<distype>::SmoothHeavisideFunction(
    const double charelelength, const double phi, double& smoothH)
{
  // assume interface thickness
  const double epsilon = charelelength;

  if (phi < -epsilon)
    smoothH = 0.0;
  else if (phi > epsilon)
    smoothH = 1.0;
  else
    smoothH = 0.5 * (1.0 + phi / epsilon + sin(phi * M_PI / epsilon) / M_PI);

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<CORE::FE::CellType::nurbs27>;

BACI_NAMESPACE_CLOSE
