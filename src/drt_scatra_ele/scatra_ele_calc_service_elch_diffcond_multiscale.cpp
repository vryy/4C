/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch_diffcond_multiscale.cpp

\brief evaluation of ScaTra elements for diffusion-conduction ion-transport equations within a multi-scale framework

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_elch_diffcond_multiscale.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/elchmat.H"
#include "../drt_mat/elchphase.H"
#include "../drt_mat/newman_multiscale.H"

/*----------------------------------------------------------------------*
 | calculate electrode state of charge                       fang 08/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype>::CalculateElectrodeSOC(
    const DRT::Element*         ele,              //!< element
    Teuchos::ParameterList&     params,           //!< parameter list
    DRT::Discretization&        discretization,   //!< discretization
    const std::vector<int>&     lm,               //!< location vector
    Epetra_SerialDenseVector&   scalars           //!< result vector for scalar integrals to be computed
    )
{
  // safety check
  if(my::numscal_ != 1)
    dserror("Electrode state of charge can only be computed for one transported scalar!");

  // extract multi-scale material
  const Teuchos::RCP<const MAT::ElchMat> elchmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(ele->Material());
  const Teuchos::RCP<const MAT::ElchPhase> elchphase = Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(0)));
  const Teuchos::RCP<MAT::NewmanMultiScale> newmanmultiscale = Teuchos::rcp_dynamic_cast<MAT::NewmanMultiScale>(elchphase->MatById(elchphase->MatID(0)));

  // initialize variables for concentration and domain integrals
  double intconcentration(0.);
  double intdomain(0.);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for(int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    // calculate concentration integral
    intconcentration += newmanmultiscale->EvaluateMeanConcentration(iquad)*fac;

    // calculate domain integral
    intdomain += fac;
  } // loop over integration points

  // safety check
  if(scalars.Length() != 2)
    dserror("Result vector for electrode state of charge computation has invalid length!");

  // write results for concentration and domain integrals into result vector
  scalars(0) = intconcentration;
  scalars(1) = intdomain;

  return;
} // DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype>::CalculateElectrodeSOC


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 07/17 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype>::EvaluateAction(
    DRT::Element*                  ele,
    Teuchos::ParameterList&        params,
    DRT::Discretization&           discretization,
    const SCATRA::Action&          action,
    DRT::Element::LocationArray&   la,
    Epetra_SerialDenseMatrix&      elemat1_epetra,
    Epetra_SerialDenseMatrix&      elemat2_epetra,
    Epetra_SerialDenseVector&      elevec1_epetra,
    Epetra_SerialDenseVector&      elevec2_epetra,
    Epetra_SerialDenseVector&      elevec3_epetra
    )
{
  // extract multi-scale material
  const Teuchos::RCP<const MAT::ElchMat> elchmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(ele->Material());
  const Teuchos::RCP<const MAT::ElchPhase> elchphase = Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(0)));
  const Teuchos::RCP<MAT::NewmanMultiScale> newmanmultiscale = Teuchos::rcp_dynamic_cast<MAT::NewmanMultiScale>(elchphase->MatById(elchphase->MatID(0)));

  // determine and evaluate action
  switch(action)
  {
    case SCATRA::micro_scale_initialize:
    {
      const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for(int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
        // initialize micro scale in multi-scale simulations
        newmanmultiscale->Initialize(ele->Id(),iquad);

      break;
    }

    case SCATRA::micro_scale_prepare_time_step:
    {
      // extract state variables at element nodes
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_,1> >(*discretization.GetState("phinp"),my::ephinp_,la[0].lm_);

      const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for(int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
      {
        // evaluate shape functions at Gauss point
        this->EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

        // evaluate state variables at Gauss point
        this->SetInternalVariablesForMatAndRHS();

        // initialize vector with macro-scale state variables
        std::vector<double> phinp(3,0.);
        phinp[0] = my::scatravarmanager_->Phinp(0);
        phinp[1] = my::funct_.Dot(my::ephinp_[1]);
        phinp[2] = my::funct_.Dot(my::ephinp_[2]);

        // prepare time step on micro scale
        newmanmultiscale->PrepareTimeStep(iquad,phinp);
      }

      break;
    }

    case SCATRA::micro_scale_update:
    {
      const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for(int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
        // update multi-scale scalar transport material
        newmanmultiscale->Update(iquad);

      break;
    }

    case SCATRA::micro_scale_output:
    {
      const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for(int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
        // create output on micro scale
        newmanmultiscale->Output(iquad);

      break;
    }

    case SCATRA::micro_scale_read_restart:
    {
      const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for(int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
        // read restart on micro scale
        newmanmultiscale->ReadRestart(iquad);

      break;
    }

    default:
    {
      mydiffcond::EvaluateAction(
          ele,
          params,
          discretization,
          action,
          la,
          elemat1_epetra,
          elemat2_epetra,
          elevec1_epetra,
          elevec2_epetra,
          elevec3_epetra
          );

      break;
    }
  } // switch(action)

  return -1;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<DRT::Element::nurbs27>;
