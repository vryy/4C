/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for thermodynamic diffusion-conduction ion-transport equations

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_scatra_ele_calc_elch_diffcond_sti_thermo.hpp"

#include "4C_scatra_ele_parameter_elch.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>*
Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcElchDiffCondSTIThermo<distype>>(
            new ScaTraEleCalcElchDiffCondSTIThermo<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 | extract quantities for element evaluation                 fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<
    distype>::extract_element_and_node_values(Core::Elements::Element* ele,  //!< current element
    Teuchos::ParameterList& params,                                          //!< parameter list
    Core::FE::Discretization& discretization,                                //!< discretization
    Core::Elements::Element::LocationArray& la                               //!< location array
)
{
  // call base class routine to extract scatra-related quantities
  mydiffcond::extract_element_and_node_values(ele, params, discretization, la);

  // call base class routine to extract thermo-related quantitites
  mythermo::extract_element_and_node_values(ele, params, discretization, la);
}


/*----------------------------------------------------------------------*
 | get material parameters                                   fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::get_material_params(
    const Core::Elements::Element* ele,  //!< current element
    std::vector<double>& densn,          //!< density at t_(n)
    std::vector<double>& densnp,         //!< density at t_(n+1) or t_(n+alpha_F)
    std::vector<double>& densam,         //!< density at t_(n+alpha_M)
    double& visc,                        //!< fluid viscosity
    const int iquad                      //!< ID of current integration point
)
{
  // call base class routine to get parameters of primary, isothermal electrolyte material
  mydiffcond::get_material_params(ele, densn, densnp, densam, visc, iquad);

  // get parameters of secondary, thermodynamic electrolyte material
  Teuchos::RCP<const Core::Mat::Material> material = ele->Material(1);
  materialtype_ = material->MaterialType();
  if (materialtype_ == Core::Materials::m_soret) mythermo::mat_soret(material);
}  // Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::get_material_params


/*--------------------------------------------------------------------------*
 | calculate element matrix and element right-hand side vector   fang 11/15 |
 *--------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::calc_mat_and_rhs(
    Core::LinAlg::SerialDenseMatrix& emat,  //!< element matrix
    Core::LinAlg::SerialDenseVector& erhs,  //!< element right-hand side vector
    const int k,                            //!< index of current scalar
    const double fac,                       //!< domain integration factor
    const double timefacfac,  //!< domain integration factor times time integration factor
    const double rhsfac,      //!< domain integration factor times time integration factor for
                              //!< right-hand side vector
    const double taufac,      //!< domain integration factor times stabilization parameter
    const double timetaufac,  //!< domain integration factor times stabilization parameter times
                              //!< time integration factor
    const double rhstaufac,  //!< domain integration factor times stabilization parameter times time
                             //!< integration factor for right-hand side vector
    Core::LinAlg::Matrix<nen_, 1>&
        tauderpot,  //!< derivatives of stabilization parameter w.r.t. electric potential
    double& rhsint  //!< body force value
)
{
  // call base class routine for isothermal problems
  mydiffcond::calc_mat_and_rhs(
      emat, erhs, k, fac, timefacfac, rhsfac, taufac, timetaufac, rhstaufac, tauderpot, rhsint);

  if (materialtype_ == Core::Materials::m_soret)
  {
    // extract variables and parameters
    const double& concentration = var_manager()->Phinp(0);
    const Core::LinAlg::Matrix<nsd_, 1>& gradtemp = var_manager()->GradTemp();
    const double& kappa = mydiffcond::diff_manager()->GetCond();
    const double& kappaderiv = mydiffcond::diff_manager()->GetConcDerivCond(0);
    const double faraday = Discret::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
    const double invffval = mydiffcond::diff_manager()->InvFVal(0) / faraday;
    const double& invfval = mydiffcond::diff_manager()->InvFVal(0);
    const double& R = Discret::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
    const double& t = mydiffcond::diff_manager()->GetTransNum(0);
    const double& tderiv = mydiffcond::diff_manager()->GetDerivTransNum(0, 0);

    // matrix and vector contributions arising from additional, thermodynamic term in expression for
    // current density
    for (int vi = 0; vi < static_cast<int>(nen_); ++vi)
    {
      // recurring indices
      const int rowconc(vi * 2);
      const int rowpot(vi * 2 + 1);

      // gradient of test function times gradient of temperature
      double laplawfrhs_temp(0.);
      my::get_laplacian_weak_form_rhs(laplawfrhs_temp, gradtemp, vi);

      for (int ui = 0; ui < static_cast<int>(nen_); ++ui)
      {
        // recurring index
        const int colconc(ui * 2);

        // linearizations of contributions for concentration residuals w.r.t. concentration dofs
        emat(rowconc, colconc) -=
            timefacfac * my::funct_(ui) * laplawfrhs_temp * pow(invfval, 2.0) *
            (tderiv * kappa * R * log(concentration) + t * kappaderiv * R * log(concentration) +
                t * kappa * R / concentration);

        // linearizations of contributions for electric potential residuals w.r.t. concentrationdofs
        emat(rowpot, colconc) -= timefacfac * my::funct_(ui) * laplawfrhs_temp * invffval *
                                 (kappaderiv * R * log(concentration) + kappa * R / concentration);

        // linearizations w.r.t. electric potential dofs are zero
      }

      // contribution for concentration residual
      erhs[rowconc] +=
          rhsfac * laplawfrhs_temp * kappa * pow(invfval, 2.0) * t * R * log(concentration);

      // contribution for electric potential residual
      erhs[rowpot] += rhsfac * laplawfrhs_temp * kappa * invffval * R * log(concentration);
    }

    // matrix and vector contributions arising from additional, thermodynamic term for Soret effect
    mythermo::calc_mat_soret(emat, timefacfac, var_manager()->Phinp(0),
        mydiffcond::diff_manager()->GetIsotropicDiff(0),
        mydiffcond::diff_manager()->get_conc_deriv_iso_diff_coef(0, 0), var_manager()->Temp(),
        var_manager()->GradTemp(), my::funct_, my::derxy_);
    mythermo::calc_rhs_soret(erhs, var_manager()->Phinp(0),
        mydiffcond::diff_manager()->GetIsotropicDiff(0), rhsfac, var_manager()->Temp(),
        var_manager()->GradTemp(), my::derxy_);
  }
}


/*----------------------------------------------------------------------*
 | evaluate action for off-diagonal system matrix block      fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::EvaluateActionOD(
    Core::Elements::Element* ele,                     //!< current element
    Teuchos::ParameterList& params,                   //!< parameter list
    Core::FE::Discretization& discretization,         //!< discretization
    const ScaTra::Action& action,                     //!< action parameter
    Core::Elements::Element::LocationArray& la,       //!< location array
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
    Core::LinAlg::SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
    Core::LinAlg::SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
    Core::LinAlg::SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
)
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
void Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::sysmat_od_scatra_thermo(
    Core::Elements::Element* ele,          //!< current element
    Core::LinAlg::SerialDenseMatrix& emat  //!< element matrix
)
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

    if (materialtype_ == Core::Materials::m_soret)
    {
      // extract variables and parameters
      const double& concentration = var_manager()->Phinp(0);
      const Core::LinAlg::Matrix<nsd_, 1>& gradconc = var_manager()->GradPhi(0);
      const double faraday =
          Discret::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->Faraday();
      const double& invffval = mydiffcond::diff_manager()->InvFVal(0) / faraday;
      const double& invfval = mydiffcond::diff_manager()->InvFVal(0);
      const double& kappa = mydiffcond::diff_manager()->GetCond();
      const double& R =
          Discret::ELEMENTS::ScaTraEleParameterElch::Instance("scatra")->GasConstant();
      const double& t = mydiffcond::diff_manager()->GetTransNum(0);

      // matrix contributions arising from additional, thermodynamic term in expression for current
      // density
      for (int vi = 0; vi < static_cast<int>(nen_); ++vi)
      {
        // recurring indices
        const int rowconc(vi * 2);
        const int rowpot(vi * 2 + 1);

        // gradient of test function times gradient of concentration
        double laplawfrhs_conc = 0.0;
        my::get_laplacian_weak_form_rhs(laplawfrhs_conc, gradconc, vi);

        for (int ui = 0; ui < static_cast<int>(nen_); ++ui)
        {
          // gradient of test function times gradient of shape function
          double laplawf(0.);
          my::get_laplacian_weak_form(laplawf, vi, ui);

          // gradient of test function times derivative of current density w.r.t. temperature
          const double di_dT =
              2.0 * kappa * (1.0 - t) * invfval * R / concentration * laplawfrhs_conc;

          // formal, symbolic derivative of current density w.r.t. temperature gradient
          const double di_dgradT = kappa * invfval * R * log(concentration);

          // linearizations of contributions for concentration residuals w.r.t. thermo dofs
          emat(rowconc, ui) -= timefacfac * (t * invfval * di_dT * my::funct_(ui) +
                                                t * invfval * di_dgradT * laplawf);

          // linearizations of contributions for electric potential residuals w.r.t. thermo dofs
          emat(rowpot, ui) -= timefacfac * (laplawf * kappa * invffval * R * log(concentration) +
                                               my::funct_(ui) * kappa * (1.0 - t) * invffval * R /
                                                   concentration * 2.0 * laplawfrhs_conc);
        }
      }

      // provide element matrix with linearizations of Soret term in discrete scatra residuals
      // w.r.t. thermo dofs
      mythermo::calc_mat_soret_od(emat, timefacfac, var_manager()->Phinp(0),
          mydiffcond::diff_manager()->GetIsotropicDiff(0), var_manager()->Temp(),
          var_manager()->GradTemp(), my::funct_, my::derxy_);
    }

    // calculating the off diagonal for the temperature derivative of concentration and electric
    // potential
    mythermo::calc_mat_diff_thermo_od(emat, my::numdofpernode_, timefacfac, var_manager()->InvF(),
        var_manager()->GradPhi(0), var_manager()->GradPot(),
        myelectrode::diff_manager()->get_temp_deriv_iso_diff_coef(0, 0),
        myelectrode::diff_manager()->GetTempDerivCond(0), my::funct_, my::derxy_, 1.0);
  }
}


/*------------------------------------------------------------------------------*
 | set internal variables for element evaluation                     fang 11/15 |
 *------------------------------------------------------------------------------*/
template <Core::FE::CellType distype>
void Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<
    distype>::set_internal_variables_for_mat_and_rhs()
{
  // set internal variables for element evaluation
  var_manager()->set_internal_variables(my::funct_, my::derxy_, mythermo::etempnp_, my::ephinp_,
      my::ephin_, my::econvelnp_, my::ehist_);
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 11/15 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<distype>::ScaTraEleCalcElchDiffCondSTIThermo(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructors of base classes
      mydiffcond::ScaTraEleCalcElchDiffCond(numdofpernode, numscal, disname),
      mythermo::ScaTraEleSTIThermo(numscal)
{
  // safety check
  if (numscal != 1 or numdofpernode != 2)
    FOUR_C_THROW("Invalid number of transported scalars or degrees of freedom per node!");

  // replace internal variable manager for isothermal diffusion-conduction formulation by internal
  // variable manager for thermodynamic diffusion-conduction formulation
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerElchDiffCondSTIThermo<nsd_, nen_>(
          my::numscal_, myelch::elchparams_, mydiffcond::diffcondparams_));
}


// template classes
// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::line2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::line3>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::tri3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::tri6>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::quad4>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::quad9>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::nurbs9>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::hex8>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::hex27>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::tet4>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::tet10>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::pyramid5>;
// template class
// Discret::ELEMENTS::ScaTraEleCalcElchDiffCondSTIThermo<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
