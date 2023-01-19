/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for level sets

\level 2

*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_ls.H"

#include "scatra_ele_action.H"

#include "inpar_levelset.H"

#include "lib_utils.H"

#include "geometry_position_array.H"
#include "lib_discret.H"


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcLS<distype>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
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

      std::vector<LINALG::Matrix<my::nen_, 1>> ephizero(my::numscal_);

      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phinp, my::ephinp_, lm);
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phizero, ephizero, lm);

      // check if length suffices
      if (elevec1_epetra.Length() < 1) dserror("Result vector too short");

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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcLS<distype>::CalErrorComparedToAnalytSolution(
    const DRT::Element* ele, const std::vector<LINALG::Matrix<my::nen_, 1>>& ephizero,
    Teuchos::ParameterList& params, Epetra_SerialDenseVector& errors)
{
  // get element volume
  const double vol = my::EvalShapeFuncAndDerivsAtEleCenter();

  // get elemet length
  // cast dimension to a double varibale -> pow()
  const double dim = double(my::nsd_);
  const double h = std::pow(vol, 1.0 / dim);

  // integration points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
      SCATRA::DisTypeToGaussRuleForExactSol<distype>::rule);

  const INPAR::SCATRA::CalcError errortype =
      DRT::INPUT::get<INPAR::SCATRA::CalcError>(params, "calcerrorflag");
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
template <DRT::Element::DiscretizationType distype>
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
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcLS<DRT::Element::nurbs27>;
