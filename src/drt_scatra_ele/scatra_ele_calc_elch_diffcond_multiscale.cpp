/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements for diffusion-conduction ion-transport equations within a
multi-scale framework

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_elch_diffcond_multiscale.H"

#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_mat/elchmat.H"
#include "../drt_mat/elchphase.H"
#include "../drt_mat/newman_multiscale.H"

/*----------------------------------------------------------------------*
 | singleton access method                                   fang 07/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype>*
DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname,
    const ScaTraEleCalcElchDiffCondMultiScale* delete_me)
{
  static std::map<std::string, ScaTraEleCalcElchDiffCondMultiScale<distype>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleCalcElchDiffCondMultiScale<distype>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleCalcElchDiffCondMultiScale<distype>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 | singleton destruction                                     fang 07/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype>::Done()
{
  // delete singleton
  Instance(0, 0, "", this);

  return;
}


/*----------------------------------------------------------------------*
 | private constructor for singletons                        fang 07/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype>::ScaTraEleCalcElchDiffCondMultiScale(
    const int numdofpernode, const int numscal,
    const std::string& disname)
    :  // constructor of base class
      DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::ScaTraEleCalcElchDiffCond(
          numdofpernode, numscal, disname)
{
  // replace diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElchDiffCondMultiScale(my::numscal_));

  return;
}


/*-----------------------------------------------------------------------------------------------------------------------*
 | macro-scale matrix and vector contributions arising from macro-micro coupling in multi-scale
 simulations   fang 07/17 |
 *-----------------------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype>::CalcMatAndRhsMultiScale(
    const DRT::Element* const ele,   //!< element
    Epetra_SerialDenseMatrix& emat,  //!< element matrix
    Epetra_SerialDenseVector& erhs,  //!< element right-hand side vector
    const int k,                     //!< species index
    const int iquad,                 //!< Gauss point index
    const double timefacfac,         //!< domain integration factor times time integration factor
    const double rhsfac  //!< domain integration factor times time integration factor for right-hand
                         //!< side vector
)
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
  std::vector<double> phinp(3, 0.);
  phinp[0] = my::scatravarmanager_->Phinp(0);
  phinp[1] = my::funct_.Dot(my::ephinp_[1]);
  phinp[2] = my::funct_.Dot(my::ephinp_[2]);

  // evaluate multi-scale Newman material
  newmanmultiscale->Evaluate(iquad, phinp, q_micro, dq_dphi_micro,
      not DRT::ELEMENTS::ScaTraEleParameterStd::Instance("scatra")->PartitionedMultiScale());

  // calculate gradient of electric potential inside electrode
  LINALG::Matrix<my::nsd_, 1> gradpot_ed(true);
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
  const double dq_dc_el = timefacfac * dq_dphi_micro[0] * newmanmultiscale->A_s();
  const double dq_dpot_el = timefacfac * dq_dphi_micro[1] * newmanmultiscale->A_s();
  const double dq_dpot_ed = timefacfac * dq_dphi_micro[2] * newmanmultiscale->A_s();
  const double q = rhsfac * q_micro * newmanmultiscale->A_s();
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    // matrix contributions
    const double vi_dq_dc_el = my::funct_(vi) * dq_dc_el;
    const double vi_dq_dpot_el = my::funct_(vi) * dq_dpot_el;
    const double vi_dq_dpot_ed = my::funct_(vi) * dq_dpot_ed;
    const int fvi = vi * my::numdofpernode_;

    for (unsigned ui = 0; ui < my::nen_; ++ui)
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

    double laplawfrhs_gradpot(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot, gradpot_ed, vi);
    erhs[fvi + 2] -= rhsfac * mydiffcond::VarManager()->InvF() *
                         DiffManager()->GetPhasePoroTort(0) * DiffManager()->GetSigma() *
                         laplawfrhs_gradpot -
                     vi_rhsterm;
  }

  return;
}


/*------------------------------------------------------------------------*
 | compute element matrix and element right-hand side vector   fang 07/17 |
 *------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype>::Sysmat(
    DRT::Element* ele,                   //!< element
    Epetra_SerialDenseMatrix& emat,      //!< element matrix
    Epetra_SerialDenseVector& erhs,      //!< element right-hand side vector
    Epetra_SerialDenseVector& subgrdiff  //!< subgrid diffusivity vector
)
{
  // call base class routine
  mydiffcond::Sysmat(ele, emat, erhs, subgrdiff);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
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

  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::nurbs27>;
