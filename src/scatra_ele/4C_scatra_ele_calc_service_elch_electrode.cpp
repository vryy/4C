// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_discretization.hpp"
#include "4C_fem_general_extract_values.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_ele_calc_elch_electrode.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::Elements::ScaTraEleCalcElchElectrode<distype, probdim>::evaluate_action(
    Core::Elements::Element* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, const ScaTra::Action& action,
    Core::Elements::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case ScaTra::Action::calc_elch_electrode_soc_and_c_rate:
    {
      calculate_electrode_soc_and_c_rate(ele, discretization, la, elevec1_epetra);

      break;
    }
    case ScaTra::Action::calc_elch_elctrode_mean_concentration:
    {
      calculate_mean_electrode_concentration(ele, discretization, la, elevec1_epetra);

      break;
    }

    default:
    {
      myelch::evaluate_action(ele, params, discretization, action, la, elemat1_epetra,
          elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch(action)

  return 0;
}


/*----------------------------------------------------------------------------------------------------------*
 | validity check with respect to input parameters, degrees of freedom, number of scalars etc. fang
 02/15 |
 *----------------------------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchElectrode<distype, probdim>::check_elch_element_parameter(
    Core::Elements::Element* ele  //!< current element
)
{
  // safety checks
  if (ele->material()->material_type() != Core::Materials::m_electrode)
    FOUR_C_THROW("Invalid material type!");

  if (my::numscal_ != 1) FOUR_C_THROW("Invalid number of transported scalars!");
}  // Discret::Elements::ScaTraEleCalcElchElectrode<distype>::check_elch_element_parameter


/*----------------------------------------------------------------------*
 | get conductivity                                          fang 02/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchElectrode<distype, probdim>::get_conductivity(
    const enum Inpar::ElCh::EquPot equpot,  //!< type of closing equation for electric potential
    double& sigma_all,                      //!< conductivity of electrolyte solution
    std::vector<double>& sigma,  //!< conductivity or a single ion + overall electrolyte solution
    bool effCond)
{
  // use precomputed conductivity
  sigma_all = diff_manager()->get_cond();
}  // Discret::Elements::ScaTraEleCalcElchElectrode<distype>::get_conductivity


/*----------------------------------------------------------------------*
 | calculate weighted current density                        fang 07/16 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchElectrode<distype, probdim>::calculate_current(
    Core::LinAlg::Matrix<nsd_, 1>& q,        //!< flux of species k
    const Inpar::ScaTra::FluxType fluxtype,  //!< type fo flux
    const double fac                         //!< integration factor
)
{
  /*
  Actually, we compute here a weighted (and integrated) form of the current density!
  On time integration level, these contributions are then used to calculate
  an L2-projected representation of the current density.
  Thus, this method here DOES NOT YET provide current density values that are ready to use!

  /                           \
  |                    /   \  |
  | w, -sigma * nabla | phi | |
  |                    \   /  |
  \                           /
  */

  switch (fluxtype)
  {
    case Inpar::ScaTra::flux_diffusive:
    case Inpar::ScaTra::flux_total:
    {
      // ohmic contribution to current density
      q.update(-diff_manager()->get_cond(), var_manager()->grad_pot());
      break;
    }

    default:
    {
      FOUR_C_THROW("Invalid flux type!");
      break;
    }
  }
}  // Discret::Elements::ScaTraEleCalcElchElectrode<distype>::calculate_current


/*----------------------------------------------------------------------*
 | calculate electrode state of charge and C rate            fang 01/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchElectrode<distype,
    probdim>::calculate_electrode_soc_and_c_rate(const Core::Elements::Element* const&
                                                     ele,  //!< the element we are dealing with
    const Core::FE::Discretization& discretization,        //!< discretization
    Core::Elements::LocationArray& la,                     //!< location array
    Core::LinAlg::SerialDenseVector& scalars  //!< result vector for scalar integrals to be computed
)
{
  // safety check
  if (my::numscal_ != 1)
    FOUR_C_THROW("Electrode state of charge can only be computed for one transported scalar!");

  // get global state vectors
  const Teuchos::RCP<const Core::LinAlg::Vector<double>> phinp = discretization.get_state("phinp");
  if (phinp == Teuchos::null) FOUR_C_THROW("Cannot get state vector \"phinp\"!");
  const Teuchos::RCP<const Core::LinAlg::Vector<double>> phidtnp =
      discretization.get_state("phidtnp");
  if (phidtnp == Teuchos::null) FOUR_C_THROW("Cannot get state vector \"phidtnp\"!");

  // extract local nodal values from global state vectors
  Core::FE::extract_my_values(*phinp, my::ephinp_, la[0].lm_);
  static std::vector<Core::LinAlg::Matrix<nen_, 1>> ephidtnp(2);
  Core::FE::extract_my_values(*phidtnp, ephidtnp, la[0].lm_);

  // initialize variables for integrals of concentration, its time derivative, and domain
  double intconcentration(0.);
  double intconcentrationtimederiv(0.);
  double intdomain(0.);

  // integration points and weights
  const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      ScaTra::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

    // calculate integrals of concentration and its time derivative
    for (unsigned vi = 0; vi < nen_; ++vi)
    {
      const double vi_fac = my::funct_(vi) * fac;

      // integral of concentration
      intconcentration += my::ephinp_[0](vi) * vi_fac;

      // integral of time derivative of concentration
      intconcentrationtimederiv += ephidtnp[0](vi) * vi_fac;
    }

    // domain integral
    intdomain += fac;
  }  // loop over integration points

  // safety check
  if (not my::scatrapara_->is_ale() and scalars.length() != 3)
    FOUR_C_THROW(
        "Result vector for electrode state of charge and C rate computation has invalid length!");

  // write results for integrals of concentration, its time derivative, and domain into result
  // vector
  scalars(0) = intconcentration;
  scalars(1) = intconcentrationtimederiv;
  scalars(2) = intdomain;

  // additional computations in case of ALE
  if (my::scatrapara_->is_ale())
  {
    const int ndsvel = my::scatrapara_->nds_vel();
    // extract velocities
    const Teuchos::RCP<const Core::LinAlg::Vector<double>> vel =
        discretization.get_state(ndsvel, "velocity field");
    if (vel == Teuchos::null) FOUR_C_THROW("Cannot get state vector \"velocity field\"!");
    Core::FE::extract_my_values(*vel, my::evelnp_, la[ndsvel].lm_);

    // initialize additional variables for integrals related to velocity divergence
    double intdivv(0.);
    double intcdivv(0.);
    double intvgradc(0.);

    // loop over integration points
    for (int iquad = 0; iquad < intpoints.ip().nquad; ++iquad)
    {
      // evaluate values of shape functions and domain integration factor at current integration
      // point
      const double fac = my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

      // compute internal variables at current integration point
      var_manager()->set_internal_variables_elch_electrode_soc_and_c_rate(
          my::funct_, my::derxy_, my::ephinp_, my::ephin_, my::econvelnp_, my::ehist_);

      // compute velocity and its divergence
      static Core::LinAlg::Matrix<nsd_, 1> v;
      v.multiply(my::evelnp_, my::funct_);
      double divv(0.);
      my::get_divergence(divv, my::evelnp_);

      // integral of velocity divergence
      const double divv_fac = divv * fac;
      intdivv += divv_fac;

      // integral of concentration times velocity divergence
      intcdivv += my::scatravarmanager_->phinp(0) * divv_fac;

      // integral of velocity times concentration gradient
      intvgradc += v.dot(my::scatravarmanager_->grad_phi(0)) * fac;
    }  // loop over integration points

    // safety check
    if (scalars.length() != 6)
      FOUR_C_THROW(
          "Result vector for electrode state of charge and C rate computation has invalid length!");

    // write results for integrals related to velocity divergence into result vector
    scalars(3) = intdivv;
    scalars(4) = intcdivv;
    scalars(5) = intvgradc;
  }
}  // Discret::Elements::ScaTraEleCalcElchElectrode<distype>::calculate_electrode_soc_and_c_rate


/*---------------------------------------------------------------------*
 | calculate weighted mass flux (no reactive flux so far)   fang 02/15 |
 *---------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchElectrode<distype, probdim>::calculate_flux(
    Core::LinAlg::Matrix<nsd_, 1>& q,        //!< flux of species k
    const Inpar::ScaTra::FluxType fluxtype,  //!< type fo flux
    const int k                              //!< index of current scalar
)
{
  /*
    Actually, we compute here a weighted (and integrated) form of the fluxes!
    On time integration level, these contributions are then used to calculate
    an L2-projected representation of fluxes.
    Thus, this method here DOES NOT YET provide flux values that are ready to use!
  */

  // add convective flux contribution
  switch (fluxtype)
  {
    case Inpar::ScaTra::flux_diffusive:
    case Inpar::ScaTra::flux_total:
    {
      // diffusive flux contribution
      q.update(-diff_manager()->get_isotropic_diff(k), var_manager()->grad_phi(k));
      break;
    }

    default:
    {
      FOUR_C_THROW("received illegal flag inside flux evaluation for whole domain");
      break;
    }
  }
}  // Discret::Elements::ScaTraEleCalcElchElectrode<distype>::calculate_flux


/*----------------------------------------------------------------------------------------*
 | calculate error of numerical solution with respect to analytical solution   fang 10/16 |
 *----------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchElectrode<distype,
    probdim>::cal_error_compared_to_analyt_solution(const Core::Elements::Element*
                                                        ele,  //!< element
    Teuchos::ParameterList& params,                           //!< parameter list
    Core::LinAlg::SerialDenseVector& errors  //!< vector containing L2 and H1 error norms
)
{
  // call base class routine
  myelch::cal_error_compared_to_analyt_solution(ele, params, errors);
}  // Discret::Elements::ScaTraEleCalcElchElectrode<distype>::cal_error_compared_to_analyt_solution


/*------------------------------------------------------------------------------*
 | set internal variables for electrodes                             fang 02/15 |
 *------------------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchElectrode<distype,
    probdim>::set_internal_variables_for_mat_and_rhs()
{
  // set internal variables
  var_manager()->set_internal_variables_elch_electrode(
      my::funct_, my::derxy_, my::ephinp_, my::ephin_, my::econvelnp_, my::ehist_);
}  // Discret::Elements::ScaTraEleCalcElchElectrode<distype>::set_internal_variables_for_mat_and_rhs()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::Elements::ScaTraEleCalcElchElectrode<distype,
    probdim>::calculate_mean_electrode_concentration(const Core::Elements::Element* const& ele,
    const Core::FE::Discretization& discretization, Core::Elements::LocationArray& la,
    Core::LinAlg::SerialDenseVector& conc)
{
  // for complete 1D simulation of battery:
  // Micro state must exist for electrolyte -> set value to 00
  for (int node = 0; node < static_cast<int>(nen_); ++node) conc(node) = 0.0;
}

// template classes
// 1D elements
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::line2, 1>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::line2, 2>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::line2, 3>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::line3, 1>;

// 2D elements
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::tri3, 2>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::tri3, 3>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::tri6, 2>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::quad4, 2>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::quad4, 3>;
// template class
// Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::quad8>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::quad9, 2>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::nurbs9, 2>;

// 3D elements
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::hex8, 3>;
// template class
// Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::hex20>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::hex27, 3>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::tet4, 3>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::tet10, 3>;
// template class
// Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::wedge6>;
template class Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::pyramid5, 3>;
// template class
// Discret::Elements::ScaTraEleCalcElchElectrode<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
