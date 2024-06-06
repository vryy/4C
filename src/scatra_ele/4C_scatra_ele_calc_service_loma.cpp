/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for loma

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_discretization_fem_general_extract_values.hpp"
#include "4C_discretization_geometry_position_array.hpp"
#include "4C_lib_discret.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_sutherland.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_calc_loma.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::evaluate_action(Core::Elements::Element* ele,
    Teuchos::ParameterList& params, Discret::Discretization& discretization,
    const ScaTra::Action& action, Core::Elements::Element::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  const std::vector<int>& lm = la[0].lm_;

  // determine and evaluate action
  if (action == ScaTra::Action::calc_dissipation or action == ScaTra::Action::calc_mean_Cai)
  {
    if (my::scatraparatimint_->IsGenAlpha())
    {
      // extract additional local values from global vector
      Teuchos::RCP<const Epetra_Vector> phiam = discretization.GetState("phiam");
      if (phiam == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'phiam'");
      Core::FE::ExtractMyValues<Core::LinAlg::Matrix<nen_, 1>>(*phiam, ephiam_, lm);
    }
  }

  if (action == ScaTra::Action::calc_dissipation or action == ScaTra::Action::calc_mean_Cai)
  {
    // get thermodynamic pressure (for calc_dissipation)
    thermpressnp_ = params.get<double>("thermodynamic pressure");
    thermpressdt_ = params.get<double>("time derivative of thermodynamic pressure");
    if (my::scatraparatimint_->IsGenAlpha())
      thermpressam_ = params.get<double>("thermodynamic pressure at n+alpha_M");
  }

  switch (action)
  {
    case ScaTra::Action::calc_domain_and_bodyforce:
    {
      // NOTE: add integral values only for elements which are NOT ghosted!
      if (ele->Owner() == discretization.Comm().MyPID())
        // calculate domain and bodyforce integral
        calculate_domain_and_bodyforce(elevec1_epetra, ele);

      break;
    }
    default:
    {
      my::evaluate_action(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);
      break;
    }
  }  // switch(action)

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate domain integral                                   vg 01/09|
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::calculate_domain_and_bodyforce(
    Core::LinAlg::SerialDenseVector& scalars, const Core::Elements::Element* ele)
{
  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  my::body_force(ele);

  // integration points and weights
  const Core::FE::IntPointsAndWeights<nsd_> intpoints(ScaTra::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

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
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::extract_element_and_node_values(
    Core::Elements::Element* ele, Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Elements::Element::LocationArray& la)
{
  // add loma-specific values
  if (my::scatraparatimint_->IsGenAlpha())
  {
    // extract local values from global vector
    Teuchos::RCP<const Epetra_Vector> phiam = discretization.GetState("phiam");
    if (phiam == Teuchos::null) FOUR_C_THROW("Cannot get state vector 'phiam'");
    Core::FE::ExtractMyValues<Core::LinAlg::Matrix<nen_, 1>>(*phiam, ephiam_, la[0].lm_);
  }

  // get thermodynamic pressure
  thermpressnp_ = params.get<double>("thermodynamic pressure");
  thermpressdt_ = params.get<double>("time derivative of thermodynamic pressure");
  if (my::scatraparatimint_->IsGenAlpha())
    thermpressam_ = params.get<double>("thermodynamic pressure at n+alpha_M");

  // call base class routine
  my::extract_element_and_node_values(ele, params, discretization, la);

  return;
}


/*-----------------------------------------------------------------------------*
 | get density at integration point                                 fang 02/15 |
 *-----------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
double Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::get_density(
    const Core::Elements::Element* ele, Teuchos::RCP<const Core::Mat::Material> material,
    Teuchos::ParameterList& params, const double tempnp)
{
  // initialization
  double density(0.);

  // check whether temperature is positive
  if (tempnp < 0.0)
    FOUR_C_THROW("Negative temperature in ScaTra low-Mach-number routine 'GetDensity'!");

  if (material->MaterialType() == Core::Materials::m_sutherland)
  {
    // get thermodynamic pressure
    const double thermpress = params.get<double>("thermpress");

    density = Teuchos::rcp_dynamic_cast<const Mat::Sutherland>(material)->ComputeDensity(
        tempnp, thermpress);
  }
  else
    FOUR_C_THROW("Invalid material type!");

  return density;
}  // Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::GetDensity


/*-----------------------------------------------------------------------------*
 | calculate viscous part of subgrid-scale velocity                 fang 02/15 |
 *-----------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::calc_subgr_velocity_visc(
    Core::LinAlg::Matrix<nsd_, 1>& epsilonvel)
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
    FOUR_C_THROW("Epsilon(u) is not implemented for the 1D case!");

  my::derxy2_.Scale(1.0 / prefac);

  return;
}  // Discret::ELEMENTS::ScaTraEleCalcLoma<distype>::calc_subgr_velocity_visc


// template classes

// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::line3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::quad4>;
// template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::nurbs9>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::hex8>;
// template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::tet10>;
// template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::pyramid5>;
// template class Discret::ELEMENTS::ScaTraEleCalcLoma<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
