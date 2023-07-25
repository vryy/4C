/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements for diffusion-conduction ion-transport equations within a
multi-scale framework

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "baci_scatra_ele_calc_elch_diffcond_multiscale.H"

#include "baci_scatra_ele_parameter_std.H"
#include "baci_scatra_ele_parameter_timint.H"

#include "baci_mat_elchmat.H"
#include "baci_mat_elchphase.H"
#include "baci_mat_newman_multiscale.H"
#include "baci_utils_singleton_owner.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>>(
            new ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>(
                numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype,
    probdim>::ScaTraEleCalcElchDiffCondMultiScale(const int numdofpernode, const int numscal,
    const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype, probdim>::ScaTraEleCalcElchDiffCond(
          numdofpernode, numscal, disname)
{
  // replace diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElchDiffCondMultiScale(my::numscal_));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>::CalcMatAndRhsMultiScale(
    const DRT::Element* const ele, Epetra_SerialDenseMatrix& emat, Epetra_SerialDenseVector& erhs,
    const int k, const int iquad, const double timefacfac, const double rhsfac)
{
  // extract multi-scale Newman material
  const Teuchos::RCP<const MAT::ElchMat> elchmat =
      Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(ele->Material());
  const Teuchos::RCP<const MAT::ElchPhase> elchphase =
      Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(0)));
  const Teuchos::RCP<MAT::NewmanMultiScale> newmanmultiscale =
      Teuchos::rcp_dynamic_cast<MAT::NewmanMultiScale>(elchphase->MatById(elchphase->MatID(0)));

  // initialize variables for micro-scale coupling flux and derivatives of micro-scale coupling flux
  // w.r.t. macro-scale state variables
  double q_micro(0.);
  std::vector<double> dq_dphi_micro(3, 0.);

  // initialize vector with macro-scale state variables
  std::vector<double> phinp(3, 0.0);
  phinp[0] = my::scatravarmanager_->Phinp(0);
  phinp[1] = my::funct_.Dot(my::ephinp_[1]);
  phinp[2] = my::funct_.Dot(my::ephinp_[2]);

  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  const double detF = my::EvalDetFAtIntPoint(ele, intpoints, iquad);

  // evaluate multi-scale Newman material
  newmanmultiscale->Evaluate(iquad, phinp, q_micro, dq_dphi_micro, detF,
      not DRT::ELEMENTS::ScaTraEleParameterStd::Instance("scatra")->PartitionedMultiScale());

  // calculate gradient of electric potential inside electrode
  CORE::LINALG::Matrix<nsd_, 1> gradpot_ed(true);
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
      newmanmultiscale->SpecificMicroScaleSurfaceArea(detF);
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
      my::GetLaplacianWeakForm(laplawf, ui, vi);
      emat(fvi + 2, fui) += -vi_dq_dc_el_ui;
      emat(fvi + 2, fui + 1) += -vi_dq_dpot_el_ui;
      emat(fvi + 2, fui + 2) += timefacfac * mydiffcond::VarManager()->InvF() *
                                    DiffManager()->GetPhasePoroTort(0) * DiffManager()->GetSigma() *
                                    laplawf -
                                vi_dq_dpot_ed_ui;
    }

    // vector contributions
    const double vi_rhsterm = my::funct_(vi) * q;

    erhs[fvi] -= vi_rhsterm;

    erhs[fvi + 1] -= vi_rhsterm;

    double laplawfrhs_gradpot(0.0);
    my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot, gradpot_ed, vi);
    erhs[fvi + 2] -= rhsfac * mydiffcond::VarManager()->InvF() *
                         DiffManager()->GetPhasePoroTort(0) * DiffManager()->GetSigma() *
                         laplawfrhs_gradpot -
                     vi_rhsterm;
  }
}

/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>::Sysmat(DRT::Element* ele,
    Epetra_SerialDenseMatrix& emat, Epetra_SerialDenseVector& erhs,
    Epetra_SerialDenseVector& subgrdiff)
{
  // call base class routine
  mydiffcond::Sysmat(ele, emat, erhs, subgrdiff);

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over all integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // preparations
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);
    this->SetInternalVariablesForMatAndRHS();
    std::vector<double> dummy(1, 0.);
    this->GetMaterialParams(ele, dummy, dummy, dummy, dummy[0], iquad);

    CalcMatAndRhsMultiScale(ele, emat, erhs, 0, iquad, my::scatraparatimint_->TimeFac() * fac,
        my::scatraparatimint_->TimeFacRhs() * fac);
  }
}

// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tet10, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::pyramid5, 3>;
