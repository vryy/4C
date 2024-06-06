/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge
within thermodynamic electrodes

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_scatra_ele_calc_elch_electrode_sti_thermo.hpp"

#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>*
Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcElchElectrodeSTIThermo<distype>>(
            new ScaTraEleCalcElchElectrodeSTIThermo<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | extract quantities for element evaluation                 fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    distype>::extract_element_and_node_values(Core::Elements::Element* ele,
    Teuchos::ParameterList& params, Discret::Discretization& discretization,
    Core::Elements::Element::LocationArray& la)
{
  // call base class routine to extract scatra-related quantities
  myelch::extract_element_and_node_values(ele, params, discretization, la);

  // call base class routine to extract thermo-related quantitites
  mythermo::extract_element_and_node_values(ele, params, discretization, la);
}


/*----------------------------------------------------------------------*
 | get material parameters                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::get_material_params(
    const Core::Elements::Element* ele, std::vector<double>& densn, std::vector<double>& densnp,
    std::vector<double>& densam, double& visc, const int iquad)
{
  // Set GP values to mat_electrode
  myelectrode::utils()->mat_electrode(
      ele->Material(), var_manager()->Phinp(0), var_manager()->Temp(), myelectrode::diff_manager());

  // get parameters of secondary, thermodynamic electrolyte material
  Teuchos::RCP<const Core::Mat::Material> material = ele->Material(1);
  materialtype_ = material->MaterialType();
  if (materialtype_ == Core::Materials::m_soret) mythermo::mat_soret(material);
}  // Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::get_material_params


/*--------------------------------------------------------------------------*
 | calculate element matrix and element right-hand side vector   fang 11/15 |
 *--------------------------------------------------------------------------*/

template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::calc_mat_and_rhs(
    Core::LinAlg::SerialDenseMatrix& emat, Core::LinAlg::SerialDenseVector& erhs, const int k,
    const double fac, const double timefacfac, const double rhsfac, const double taufac,
    const double timetaufac, const double rhstaufac, Core::LinAlg::Matrix<nen_, 1>& tauderpot,
    double& rhsint)
{
  // call base class routine for isothermal problems
  myelectrode::calc_mat_and_rhs(
      emat, erhs, k, fac, timefacfac, rhsfac, taufac, timetaufac, rhstaufac, tauderpot, rhsint);

  if (materialtype_ == Core::Materials::m_soret)
  {
    // matrix and vector contributions arising from additional, thermodynamic term for Soret effect
    mythermo::calc_mat_soret(emat, timefacfac, var_manager()->Phinp(0),
        myelectrode::diff_manager()->GetIsotropicDiff(0),
        myelectrode::diff_manager()->get_conc_deriv_iso_diff_coef(0, 0), var_manager()->Temp(),
        var_manager()->GradTemp(), my::funct_, my::derxy_);
    mythermo::calc_rhs_soret(erhs, var_manager()->Phinp(0),
        myelectrode::diff_manager()->GetIsotropicDiff(0), rhsfac, var_manager()->Temp(),
        var_manager()->GradTemp(), my::derxy_);
  }
}


/*----------------------------------------------------------------------*
 | evaluate action for off-diagonal system matrix block      fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::EvaluateActionOD(
    Core::Elements::Element* ele, Teuchos::ParameterList& params,
    Discret::Discretization& discretization, const ScaTra::Action& action,
    Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case ScaTra::Action::calc_scatra_mono_odblock_scatrathermo:
    {
      sysmat_od_scatra_thermo(ele, elemat1_epetra);

      break;
    }

    default:
    {
      // call base class routine
      my::EvaluateActionOD(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch(action)

  return 0;
}


/*------------------------------------------------------------------------------------------------------*
 | fill element matrix with linearizations of discrete scatra residuals w.r.t. thermo dofs   fang
 11/15 |
 *------------------------------------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<distype>::sysmat_od_scatra_thermo(
    Core::Elements::Element* ele, Core::LinAlg::SerialDenseMatrix& emat)
{
  // integration points and weights
  Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(ScaTra::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions, their derivatives, and domain integration factor at current
    // integration point
    const double fac = my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

    // evaluate overall integration factor
    const double timefacfac = my::scatraparatimint_->TimeFac() * fac;

    // evaluate internal variables at current integration point
    set_internal_variables_for_mat_and_rhs();

    // evaluate material parameters at current integration point
    double dummy(0.);
    std::vector<double> dummyvec(my::numscal_, 0.);
    get_material_params(ele, dummyvec, dummyvec, dummyvec, dummy, iquad);

    // calculating the off diagonal for the temperature derivative of concentration and electric
    // potential
    mythermo::calc_mat_diff_thermo_od(emat, my::numdofpernode_, timefacfac, var_manager()->InvF(),
        var_manager()->GradPhi(0), var_manager()->GradPot(),
        myelectrode::diff_manager()->get_temp_deriv_iso_diff_coef(0, 0),
        myelectrode::diff_manager()->GetTempDerivCond(0), my::funct_, my::derxy_, 1.);

    if (materialtype_ == Core::Materials::m_soret)
    {
      // provide element matrix with linearizations of Soret term in discrete scatra residuals
      // w.r.t. thermo dofs
      mythermo::calc_mat_soret_od(emat, timefacfac, var_manager()->Phinp(0),
          myelectrode::diff_manager()->GetIsotropicDiff(0), var_manager()->Temp(),
          var_manager()->GradTemp(), my::funct_, my::derxy_);
    }
  }
}


/*------------------------------------------------------------------------------*
 | set internal variables for element evaluation                     fang 11/15 |
 *------------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    distype>::set_internal_variables_for_mat_and_rhs()
{
  // set internal variables for element evaluation
  var_manager()->set_internal_variables(my::funct_, my::derxy_, my::ephinp_, my::ephin_,
      mythermo::etempnp_, my::econvelnp_, my::ehist_);
}

/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<
    distype>::ScaTraEleCalcElchElectrodeSTIThermo(const int numdofpernode, const int numscal,
    const std::string& disname)
    :  // constructors of base classes
      ScaTraEleCalcElchElectrode<distype>::ScaTraEleCalcElchElectrode(
          numdofpernode, numscal, disname),
      ScaTraEleSTIThermo<distype>::ScaTraEleSTIThermo(numscal)
{
  // safety check
  if (numscal != 1 or numdofpernode != 2)
    FOUR_C_THROW("Invalid number of transported scalars or degrees of freedom per node!");

  // replace internal variable manager for isothermal electrodes by internal variable manager for
  // thermodynamic electrodes
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerElchElectrodeSTIThermo<nsd_, nen_>(
          my::numscal_, myelch::elchparams_));
}


// template classes
// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::line3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::quad4>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::nurbs9>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::hex8>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::tet10>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::pyramid5>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcElchElectrodeSTIThermo<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
