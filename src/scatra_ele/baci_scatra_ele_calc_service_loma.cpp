/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for loma

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "baci_discretization_geometry_position_array.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_utils.hpp"
#include "baci_mat_arrhenius_temp.hpp"
#include "baci_mat_list.hpp"
#include "baci_mat_material.hpp"
#include "baci_mat_sutherland.hpp"
#include "baci_mat_tempdepwater.hpp"
#include "baci_scatra_ele.hpp"
#include "baci_scatra_ele_action.hpp"
#include "baci_scatra_ele_calc_loma.hpp"
#include "baci_scatra_ele_parameter_timint.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
int DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  const std::vector<int>& lm = la[0].lm_;

  // determine and evaluate action
  if (action == SCATRA::Action::calc_dissipation or action == SCATRA::Action::calc_mean_Cai)
  {
    if (my::scatraparatimint_->IsGenAlpha())
    {
      // extract additional local values from global vector
      Teuchos::RCP<const Epetra_Vector> phiam = discretization.GetState("phiam");
      if (phiam == Teuchos::null) dserror("Cannot get state vector 'phiam'");
      DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phiam, ephiam_, lm);
    }
  }

  if (action == SCATRA::Action::calc_dissipation or action == SCATRA::Action::calc_mean_Cai)
  {
    // get thermodynamic pressure (for CalcDissipation)
    thermpressnp_ = params.get<double>("thermodynamic pressure");
    thermpressdt_ = params.get<double>("time derivative of thermodynamic pressure");
    if (my::scatraparatimint_->IsGenAlpha())
      thermpressam_ = params.get<double>("thermodynamic pressure at n+alpha_M");
  }

  switch (action)
  {
    case SCATRA::Action::calc_domain_and_bodyforce:
    {
      // NOTE: add integral values only for elements which are NOT ghosted!
      if (ele->Owner() == discretization.Comm().MyPID())
        // calculate domain and bodyforce integral
        CalculateDomainAndBodyforce(elevec1_epetra, ele);

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


/*----------------------------------------------------------------------*
 |  calculate domain integral                                   vg 01/09|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::CalculateDomainAndBodyforce(
    CORE::LINALG::SerialDenseVector& scalars, const DRT::Element* ele)
{
  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  my::BodyForce(ele);

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // get bodyforce in gausspoint
    const double rhs = my::bodyforce_[0].Dot(my::funct_);

    // calculate integrals of domain and bodyforce
    for (unsigned i = 0; i < nen_; i++)
    {
      scalars[0] += fac * my::funct_(i);
    }
    scalars[1] += fac * rhs;

  }  // loop over integration points

  return;
}  // ScaTraEleCalcLoma::CalculateDomain


/*-----------------------------------------------------------------------------------------*
 | extract element based or nodal values and return extracted values of phinp   fang 02/15 |
 *-----------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::ExtractElementAndNodeValues(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  // add loma-specific values
  if (my::scatraparatimint_->IsGenAlpha())
  {
    // extract local values from global vector
    Teuchos::RCP<const Epetra_Vector> phiam = discretization.GetState("phiam");
    if (phiam == Teuchos::null) dserror("Cannot get state vector 'phiam'");
    DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phiam, ephiam_, la[0].lm_);
  }

  // get thermodynamic pressure
  thermpressnp_ = params.get<double>("thermodynamic pressure");
  thermpressdt_ = params.get<double>("time derivative of thermodynamic pressure");
  if (my::scatraparatimint_->IsGenAlpha())
    thermpressam_ = params.get<double>("thermodynamic pressure at n+alpha_M");

  // call base class routine
  my::ExtractElementAndNodeValues(ele, params, discretization, la);

  return;
}


/*-----------------------------------------------------------------------------*
 | get density at integration point                                 fang 02/15 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
double DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::GetDensity(const DRT::Element* ele,
    Teuchos::RCP<const MAT::Material> material, Teuchos::ParameterList& params, const double tempnp)
{
  // initialization
  double density(0.);

  // check whether temperature is positive
  if (tempnp < 0.0) dserror("Negative temperature in ScaTra low-Mach-number routine 'GetDensity'!");

  if (material->MaterialType() == INPAR::MAT::m_sutherland)
  {
    // get thermodynamic pressure
    const double thermpress = params.get<double>("thermpress");

    density = Teuchos::rcp_dynamic_cast<const MAT::Sutherland>(material)->ComputeDensity(
        tempnp, thermpress);
  }
  else if (material->MaterialType() == INPAR::MAT::m_tempdepwater)
  {
    density = Teuchos::rcp_dynamic_cast<const MAT::TempDepWater>(material)->ComputeDensity(tempnp);
  }
  else if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    // get thermodynamic pressure
    const double thermpress = params.get<double>("thermpress");

    const MAT::MatList* actmat = static_cast<const MAT::MatList*>(material.get());

    const int lastmatid = actmat->NumMat() - 1;

    // compute density based on temperature and thermodynamic pressure
    if (actmat->MaterialById(lastmatid)->MaterialType() == INPAR::MAT::m_arrhenius_temp)
      density = static_cast<const MAT::ArrheniusTemp*>(actmat->MaterialById(lastmatid).get())
                    ->ComputeDensity(tempnp, thermpress);

    else
      dserror(
          "Type of material found in material list not supported, should be Arrhenius-type "
          "temperature!");
  }
  else
    dserror("Invalid material type!");

  return density;
}  // DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::GetDensity


/*-----------------------------------------------------------------------------*
 | calculate viscous part of subgrid-scale velocity                 fang 02/15 |
 *-----------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::CalcSubgrVelocityVisc(
    CORE::LINALG::Matrix<nsd_, 1>& epsilonvel)
{
  double prefac = 1.0 / 3.0;
  my::derxy2_.Scale(prefac);

  if (nsd_ == 3)
  {
    for (unsigned i = 0; i < nen_; ++i)
    {
      double sum = (my::derxy2_(0, i) + my::derxy2_(1, i) + my::derxy2_(2, i)) / prefac;

      epsilonvel(0) +=
          ((sum + my::derxy2_(0, i)) * my::evelnp_(0, i) + my::derxy2_(3, i) * my::evelnp_(1, i) +
              my::derxy2_(4, i) * my::evelnp_(2, i)) /
          2.0;
      epsilonvel(1) +=
          (my::derxy2_(3, i) * my::evelnp_(0, i) + (sum + my::derxy2_(1, i)) * my::evelnp_(1, i) +
              my::derxy2_(5, i) * my::evelnp_(2, i)) /
          2.0;
      epsilonvel(2) +=
          (my::derxy2_(4, i) * my::evelnp_(0, i) + my::derxy2_(5, i) * my::evelnp_(1, i) +
              (sum + my::derxy2_(2, i)) * my::evelnp_(2, i)) /
          2.0;
    }
  }

  else if (nsd_ == 2)
  {
    for (unsigned i = 0; i < nen_; ++i)
    {
      double sum = (my::derxy2_(0, i) + my::derxy2_(1, i)) / prefac;

      epsilonvel(0) +=
          ((sum + my::derxy2_(0, i)) * my::evelnp_(0, i) + my::derxy2_(2, i) * my::evelnp_(1, i)) /
          2.0;
      epsilonvel(1) +=
          (my::derxy2_(2, i) * my::evelnp_(0, i) + (sum + my::derxy2_(1, i)) * my::evelnp_(1, i)) /
          2.0;
    }
  }

  else
    dserror("Epsilon(u) is not implemented for the 1D case!");

  my::derxy2_.Scale(1.0 / prefac);

  return;
}  // DRT::ELEMENTS::ScaTraEleCalcLoma<distype>::CalcSubgrVelocityVisc


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcLoma<CORE::FE::CellType::nurbs27>;

BACI_NAMESPACE_CLOSE
