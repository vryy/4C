/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements for diffusion-conduction ion-transport equations within a
multi-scale framework

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_scatra_ele_calc_elch_diffcond_multiscale.hpp"

#include "4C_mat_elchmat.hpp"
#include "4C_mat_elchphase.hpp"
#include "4C_mat_newman_multiscale.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>*
Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>>(
            new ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>(
                numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype,
    probdim>::ScaTraEleCalcElchDiffCondMultiScale(const int numdofpernode, const int numscal,
    const std::string& disname)
    : Discret::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::ScaTraEleCalcElchDiffCond(
          numdofpernode, numscal, disname)
{
  // replace diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElchDiffCondMultiScale(my::numscal_));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype,
    probdim>::calc_mat_and_rhs_multi_scale(const Core::Elements::Element* const ele,
    Core::LinAlg::SerialDenseMatrix& emat, Core::LinAlg::SerialDenseVector& erhs, const int k,
    const int iquad, const double timefacfac, const double rhsfac)
{
  // extract multi-scale Newman material
  const auto elchmat = Teuchos::rcp_dynamic_cast<const Mat::ElchMat>(ele->Material());
  const auto elchphase =
      Teuchos::rcp_dynamic_cast<const Mat::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(0)));
  const auto newmanmultiscale =
      Teuchos::rcp_dynamic_cast<Mat::NewmanMultiScale>(elchphase->MatById(elchphase->MatID(0)));

  // initialize variables for micro-scale coupling flux and derivatives of micro-scale coupling flux
  // w.r.t. macro-scale state variables
  double q_micro(0.);
  std::vector<double> dq_dphi_micro(3, 0.);

  // initialize vector with macro-scale state variables
  std::vector<double> phinp(3, 0.0);
  phinp[0] = my::scatravarmanager_->Phinp(0);
  phinp[1] = my::funct_.Dot(my::ephinp_[1]);
  phinp[2] = my::funct_.Dot(my::ephinp_[2]);

  const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      ScaTra::DisTypeToOptGaussRule<distype>::rule);

  const double detF = my::eval_det_f_at_int_point(ele, intpoints, iquad);

  // evaluate multi-scale Newman material
  newmanmultiscale->evaluate(iquad, phinp, q_micro, dq_dphi_micro, detF,
      not Discret::ELEMENTS::ScaTraEleParameterStd::Instance("scatra")->partitioned_multi_scale());

  // calculate gradient of electric potential inside electrode
  Core::LinAlg::Matrix<nsd_, 1> gradpot_ed(true);
  gradpot_ed.Multiply(my::derxy_, my::ephinp_[2]);

  // evaluate and assemble macro-scale matrix and vector contributions:
  //
  // 1.) additional source term on left-hand side of mass conservation equation for lithium
  // concentration inside electrolyte:
  //   [...] - A_s * q = 0
  //
  // 2.) additional source term on left-hand side of charge conservation equation for electric
  // potential inside electrolyte:
  //   [...] - A_s * q = 0
  //
  // 3.) additional charge conservation equation for electric potential inside electrode:
  //   1/F * nabla cdot (-sigma nabla phi_ed) + A_s * q = 0
  //                    |___________________|
  //                           = i_ed
  //
  const double specific_micro_scale_surface_area =
      newmanmultiscale->specific_micro_scale_surface_area(detF);
  const double dq_dc_el = timefacfac * dq_dphi_micro[0] * specific_micro_scale_surface_area;
  const double dq_dpot_el = timefacfac * dq_dphi_micro[1] * specific_micro_scale_surface_area;
  const double dq_dpot_ed = timefacfac * dq_dphi_micro[2] * specific_micro_scale_surface_area;
  const double q = rhsfac * q_micro * specific_micro_scale_surface_area;
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    // matrix contributions
    const double vi_dq_dc_el = my::funct_(vi) * dq_dc_el;
    const double vi_dq_dpot_el = my::funct_(vi) * dq_dpot_el;
    const double vi_dq_dpot_ed = my::funct_(vi) * dq_dpot_ed;
    const int fvi = vi * my::numdofpernode_;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const double vi_dq_dc_el_ui = vi_dq_dc_el * my::funct_(ui);
      const double vi_dq_dpot_el_ui = vi_dq_dpot_el * my::funct_(ui);
      const double vi_dq_dpot_ed_ui = vi_dq_dpot_ed * my::funct_(ui);
      const int fui = ui * my::numdofpernode_;

      emat(fvi, fui) += vi_dq_dc_el_ui;
      emat(fvi, fui + 1) += vi_dq_dpot_el_ui;
      emat(fvi, fui + 2) += vi_dq_dpot_ed_ui;

      emat(fvi + 1, fui) += vi_dq_dc_el_ui;
      emat(fvi + 1, fui + 1) += vi_dq_dpot_el_ui;
      emat(fvi + 1, fui + 2) += vi_dq_dpot_ed_ui;

      double laplawf(0.);
      my::get_laplacian_weak_form(laplawf, ui, vi);
      emat(fvi + 2, fui) += -vi_dq_dc_el_ui;
      emat(fvi + 2, fui + 1) += -vi_dq_dpot_el_ui;
      emat(fvi + 2, fui + 2) += timefacfac * mydiffcond::var_manager()->InvF() *
                                    diff_manager()->GetPhasePoroTort(0) *
                                    newmanmultiscale->electronic_cond(iquad) * laplawf -
                                vi_dq_dpot_ed_ui;
    }

    // vector contributions
    const double vi_rhsterm = my::funct_(vi) * q;

    erhs[fvi] -= vi_rhsterm;

    erhs[fvi + 1] -= vi_rhsterm;

    double laplawfrhs_gradpot(0.0);
    my::get_laplacian_weak_form_rhs(laplawfrhs_gradpot, gradpot_ed, vi);
    erhs[fvi + 2] -= rhsfac * mydiffcond::var_manager()->InvF() *
                         diff_manager()->GetPhasePoroTort(0) *
                         newmanmultiscale->electronic_cond(iquad) * laplawfrhs_gradpot -
                     vi_rhsterm;
  }
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>::sysmat(
    Core::Elements::Element* ele, Core::LinAlg::SerialDenseMatrix& emat,
    Core::LinAlg::SerialDenseVector& erhs, Core::LinAlg::SerialDenseVector& subgrdiff)
{
  // call base class routine
  mydiffcond::sysmat(ele, emat, erhs, subgrdiff);

  // integration points and weights
  const Core::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      ScaTra::DisTypeToOptGaussRule<distype>::rule);

  // loop over all integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // preparations
    const double fac = my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);
    this->set_internal_variables_for_mat_and_rhs();
    std::vector<double> dummy(1, 0.);
    this->get_material_params(ele, dummy, dummy, dummy, dummy[0], iquad);

    calc_mat_and_rhs_multi_scale(ele, emat, erhs, 0, iquad, my::scatraparatimint_->TimeFac() * fac,
        my::scatraparatimint_->TimeFacRhs() * fac);
  }
}

// template classes
// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::line2, 1>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::line2, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::line2, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::line3, 1>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::tri3, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::tri3, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::tri6, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::quad4, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::quad4, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::quad9, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::nurbs9,
    2>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::hex8, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::hex27, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::tet4, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::tet10, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<Core::FE::CellType::pyramid5,
    3>;

FOUR_C_NAMESPACE_CLOSE
