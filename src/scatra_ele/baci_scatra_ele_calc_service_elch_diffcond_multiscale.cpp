/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements for diffusion-conduction ion-transport equations within a
multi-scale framework

\level 2

*/
/*--------------------------------------------------------------------------*/
#include <Teuchos_SerialDenseSolver.hpp>
#include "baci_scatra_ele_calc_elch_diffcond_multiscale.H"
#include "baci_scatra_ele_parameter_std.H"

#include "baci_lib_discret.H"
#include "baci_lib_utils.H"

#include "baci_mat_elchmat.H"
#include "baci_mat_elchphase.H"
#include "baci_mat_newman_multiscale.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype,
    probdim>::CalculateElectrodeSOCAndCRate(const DRT::Element* const& ele,
    const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseVector& scalars)
{
  // safety check
  if (my::numscal_ != 1)
    dserror("Electrode state of charge can only be computed for one transported scalar!");

  // extract multi-scale material
  auto elchmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(ele->Material());
  auto elchphase =
      Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(0)));
  auto newmanmultiscale =
      Teuchos::rcp_dynamic_cast<MAT::NewmanMultiScale>(elchphase->MatById(elchphase->MatID(0)));

  // initialize variables for integrals of concentration, its time derivative, and domain
  double intconcentration(0.0);
  double intconcentrationtimederiv(0.0);
  double intdomain(0.0);

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // calculate integral of concentration
    intconcentration += newmanmultiscale->EvaluateMeanConcentration(iquad) * fac;

    // calculate integral of time derivative of concentration
    intconcentrationtimederiv +=
        newmanmultiscale->EvaluateMeanConcentrationTimeDerivative(iquad) * fac;

    // calculate integral of domain
    intdomain += fac;
  }  // loop over integration points

  // safety check
  if (scalars.length() != 3 and scalars.length() != 6)
    dserror("Result vector for electrode state of charge computation has invalid length!");

  // write results for concentration and domain integrals into result vector
  scalars(0) = intconcentration;
  scalars(1) = intconcentrationtimederiv;
  scalars(2) = intdomain;

  // set ale quantities to zero
  if (scalars.length() == 6) scalars(3) = scalars(4) = scalars(5) = 0.0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype,
    probdim>::CalculateMeanElectrodeConcentration(const DRT::Element* const& ele,
    const DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseVector& conc)
{
  // safety check
  if (my::numscal_ != 1)
    dserror("Electrode state of charge can only be computed for one transported scalar!");

  // extract multi-scale material
  auto elchmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(ele->Material());
  auto elchphase =
      Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(0)));
  auto newmanmultiscale =
      Teuchos::rcp_dynamic_cast<MAT::NewmanMultiScale>(elchphase->MatById(elchphase->MatID(0)));

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  if (intpoints.IP().nquad != nen_)
    dserror("number of element nodes must equal number of Gauss points for reasonable projection");

  // matrix of shape functions evaluated at Gauss points
  CORE::LINALG::SerialDenseMatrix N(nen_, nen_);

  // Gauss point concentration of electrode
  CORE::LINALG::SerialDenseMatrix conc_gp(nen_, 1);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // calculate mean concentration at Gauss point and store in vector
    const double concentration_gp = newmanmultiscale->EvaluateMeanConcentration(iquad);
    conc_gp(iquad, 0) = concentration_gp;

    // build matrix of shape functions
    for (int node = 0; node < static_cast<int>(nen_); ++node) N(iquad, node) = my::funct_(node, 0);
  }

  // conc_gp = N * conc --> conc = N^-1 * conc_gp
  using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
  using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> invert;
  invert.setMatrix(Teuchos::rcpFromRef(N));
  invert.invert();
  conc.multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, N, conc_gp, 0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>::CalculateScalars(
    const DRT::Element* ele, CORE::LINALG::SerialDenseVector& scalars, const bool inverting,
    const bool calc_grad_phi)
{
  my::CalculateScalars(ele, scalars, inverting, calc_grad_phi);

  // extract multi-scale material
  auto elchmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(ele->Material());
  auto elchphase =
      Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(0)));
  auto newmanmultiscale =
      Teuchos::rcp_dynamic_cast<MAT::NewmanMultiScale>(elchphase->MatById(elchphase->MatID(0)));

  // initialize variables for integrals of concentration, its time derivative, and domain
  double intconcentration(0.0);

  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // calculate integral of concentration
    intconcentration += newmanmultiscale->EvaluateMeanConcentration(iquad) * fac;
  }

  scalars(scalars.length() - 1) = intconcentration;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcElchDiffCondMultiScale<distype, probdim>::EvaluateAction(
    DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // extract multi-scale material
  auto elchmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(ele->Material());
  auto elchphase =
      Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(0)));
  auto newmanmultiscale =
      Teuchos::rcp_dynamic_cast<MAT::NewmanMultiScale>(elchphase->MatById(elchphase->MatID(0)));

  // determine and evaluate action
  switch (action)
  {
    case SCATRA::Action::micro_scale_initialize:
    {
      const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
        // initialize micro scale in multi-scale simulations
        newmanmultiscale->Initialize(ele->Id(), iquad, my::scatrapara_->IsAle());

      break;
    }

    case SCATRA::Action::micro_scale_prepare_time_step:
    case SCATRA::Action::micro_scale_solve:
    {
      // extract state variables at element nodes
      DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(
          *discretization.GetState("phinp"), my::ephinp_, la[0].lm_);

      const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
      {
        // evaluate shape functions at Gauss point
        this->EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

        // evaluate state variables at Gauss point
        this->SetInternalVariablesForMatAndRHS();

        // initialize vector with macro-scale state variables
        std::vector<double> phinp(3, 0.);
        phinp[0] = my::scatravarmanager_->Phinp(0);
        phinp[1] = my::funct_.Dot(my::ephinp_[1]);
        phinp[2] = my::funct_.Dot(my::ephinp_[2]);

        if (action == SCATRA::Action::micro_scale_prepare_time_step)
        {
          // prepare time step on micro scale
          newmanmultiscale->PrepareTimeStep(iquad, phinp);
        }
        else
        {
          // solve micro scale
          std::vector<double> dummy(3, 0.0);
          const double detF = my::EvalDetFAtIntPoint(ele, intpoints, iquad);
          newmanmultiscale->Evaluate(iquad, phinp, dummy[0], dummy, detF);
        }
      }

      break;
    }

    case SCATRA::Action::micro_scale_update:
    {
      const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
        // update multi-scale scalar transport material
        newmanmultiscale->Update(iquad);

      break;
    }

    case SCATRA::Action::micro_scale_output:
    {
      const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
        // create output on micro scale
        newmanmultiscale->Output(iquad);

      break;
    }

    case SCATRA::Action::micro_scale_read_restart:
    {
      const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);

      // loop over all Gauss points
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
        // read restart on micro scale
        newmanmultiscale->ReadRestart(iquad);

      break;
    }
    case SCATRA::Action::micro_scale_set_time:
    {
      const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);
      for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
        newmanmultiscale->SetTimeStepping(
            iquad, params.get<double>("dt"), params.get<double>("time"), params.get<int>("step"));
    }

    default:
    {
      mydiffcond::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra,
          elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch(action)

  return -1;
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
