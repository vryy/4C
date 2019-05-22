/*--------------------------------------------------------------------------*/
/*!

\brief evaluation of scatra elements for reinitialization equation

\level 2

\maintainer Christoph Ager
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele_calc_lsreinit.H"

#include "scatra_ele.H"
#include "scatra_ele_action.H"

#include "scatra_ele_parameter_lsreinit.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_utils.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, unsigned probDim>
int DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype, probDim>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  const std::vector<int>& lm = la[0].lm_;

  // determine and evaluate action
  switch (action)
  {
    case SCATRA::calc_mat_and_rhs_lsreinit_correction_step:
    {
      // extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> phizero = discretization.GetState("phizero");
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phizero == Teuchos::null or phinp == Teuchos::null)
        dserror("Cannot get state vector 'phizero' and/ or 'phinp'!");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phinp, my::ephinp_, lm);
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phizero, ephizero_, lm);

      //------------------------------------------------------
      // Step 1: precompute element penalty parameter
      //------------------------------------------------------

      // for correction a penalty parameter
      //                /          (phinp - phi_0)
      //               | H'(phi_0) ---------------- dOmega
      //              /                  time
      //  lambda = - ------------------------------------------
      //                /           2
      //               | (H'(phi_0))  dOmega
      //              /
      // is defined for each element

      // penalty parameter to be computed
      double penalty = 0.0;
      // calculate penalty parameter for element
      CalcElePenaltyParameter(penalty);

      //------------------------------------------------------
      // Step 2: calculate matrix and rhs
      //------------------------------------------------------

      SysmatCorrection(penalty, elemat1_epetra, elevec1_epetra);

      break;
    }
    case SCATRA::calc_node_based_reinit_velocity:
    {
      // extract local values from the global vectors
      Teuchos::RCP<const Epetra_Vector> phizero = discretization.GetState("phizero");
      Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
      if (phizero == Teuchos::null or phinp == Teuchos::null)
        dserror("Cannot get state vector 'phizero' and/ or 'phinp'!");
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phinp, my::ephinp_, lm);
      DRT::UTILS::ExtractMyValues<LINALG::Matrix<my::nen_, 1>>(*phizero, ephizero_, lm);

      // get current direction
      const int dir = params.get<int>("direction");

      SysmatNodalVel(dir, elemat1_epetra, elevec1_epetra);

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
 | setup element evaluation                                  fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, unsigned probDim>
int DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype, probDim>::SetupCalc(
    DRT::Element* ele, DRT::Discretization& discretization)
{
  // reset all managers to their default values (I feel better this way)
  DiffManager()->Reset();
  VarManager()->Reset();

  // clear all unused variables
  my::edispnp_.Clear();
  my::weights_.Clear();
  my::evelnp_.Clear();
  my::eaccnp_.Clear();
  my::eprenp_.Clear();

  // call base class routine
  return my::SetupCalc(ele, discretization);
}


/*----------------------------------------------------------------------*
 | calculate system matrix and rhs for correction step  rasthofer 12/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, unsigned probDim>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype, probDim>::SysmatCorrection(
    const double penalty,            ///< element penalty parameter
    Epetra_SerialDenseMatrix& emat,  ///< element matrix to calculate
    Epetra_SerialDenseVector& erhs   ///< element rhs to calculate
)
{
  //----------------------------------------------------------------------
  // calculation of element volume for characteristic element length
  //----------------------------------------------------------------------
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints_tau(
      SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau, 0);

  //----------------------------------------------------------------------
  // calculation of characteristic element length
  //----------------------------------------------------------------------

  // get gradient of initial phi at element center
  LINALG::Matrix<my::nsd_, 1> gradphizero(true);
  gradphizero.Multiply(my::derxy_, ephizero_[0]);

  // get characteristic element length
  const double charelelength = CalcCharEleLengthReinit(vol, gradphizero);

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // initial phi at Gauss point
    double phizero = 0.0;
    phizero = my::funct_.Dot(ephizero_[0]);
    // and corresponding gradient
    gradphizero.Clear();
    gradphizero.Multiply(my::derxy_, ephizero_[0]);
    double norm_gradphizero = gradphizero.Norm2();

    // derivative of sign function at Gauss point
    double deriv_sign = 0.0;
    DerivSignFunction(deriv_sign, charelelength, phizero);

    // scalar at integration point at time step n+1
    const double phinp = my::funct_.Dot(my::ephinp_[0]);
    Teuchos::rcp_dynamic_cast<
        DRT::ELEMENTS::ScaTraEleInternalVariableManagerLsReinit<my::nsd_, my::nen_>>(
        my::scatravarmanager_)
        ->SetPhinp(0, phinp);
    Teuchos::rcp_dynamic_cast<
        DRT::ELEMENTS::ScaTraEleInternalVariableManagerLsReinit<my::nsd_, my::nen_>>(
        my::scatravarmanager_)
        ->SetHist(0, 0.0);

    //------------------------------------------------
    // element matrix
    //------------------------------------------------

    my::CalcMatMass(emat, 0, fac, 1.0);

    //------------------------------------------------
    // element rhs
    //------------------------------------------------

    // predictor for phinp
    // caution: this function can be used here, since gen-alpha is excluded
    if (my::scatraparatimint_->IsGenAlpha()) dserror("Not supported by this implementation!");
    // note: this function computes phinp at integration point
    my::CalcRHSLinMass(erhs, 0, 0.0, -fac, 1.0, 1.0);  // sign has to be changed!!!!

    // penalty term
    CalcRHSPenalty(erhs, fac, penalty, deriv_sign, norm_gradphizero);

  }  // end: loop all Gauss points

  return;
}


/*-------------------------------------------------------------------------------*
 | calculation of element-wise denominator of penalty parameter  rasthofer 12/13 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, unsigned probDim>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype, probDim>::CalcElePenaltyParameter(
    double& penalty)
{
  // safety check
  if (lsreinitparams_->SignType() != INPAR::SCATRA::signtype_SussmanFatemi1999)
    dserror("Penalty method only for smoothed sign function: SussmanFatemi1999!");

  // denominator
  double ele_dom = 0.0;
  // nominator
  double ele_nom = 0.0;

  //----------------------------------------------------------------------
  // calculation of element volume for characteristic element length
  //----------------------------------------------------------------------
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints_tau(
      SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints_tau, 0);

  //----------------------------------------------------------------------
  // calculation of characteristic element length
  //----------------------------------------------------------------------

  // get gradient of initial phi at element center
  LINALG::Matrix<my::nsd_, 1> gradphizero(true);
  gradphizero.Multiply(my::derxy_, ephizero_[0]);

  // get characteristic element length
  const double charelelength = CalcCharEleLengthReinit(vol, gradphizero);

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // initial phi at Gauss point
    double phizero = 0.0;
    phizero = my::funct_.Dot(ephizero_[0]);
    // and corresponding gradient
    gradphizero.Clear();
    gradphizero.Multiply(my::derxy_, ephizero_[0]);
    double norm_gradphizero = gradphizero.Norm2();

    // derivative of sign function at Gauss point
    double deriv_sign = 0.0;
    DerivSignFunction(deriv_sign, charelelength, phizero);

    // current phi at Gauss point
    double phinp = 0.0;
    phinp = my::funct_.Dot(my::ephinp_[0]);

    // get sign function
    double signphi = 0.0;
    // gradient of current scalar
    LINALG::Matrix<my::nsd_, 1> gradphi(true);
    gradphi.Multiply(my::derxy_, my::ephinp_[0]);
    // get norm
    const double gradphi_norm = gradphi.Norm2();
    SignFunction(signphi, charelelength, phizero, gradphizero, phinp, gradphi);

    // get velocity at element center
    LINALG::Matrix<my::nsd_, 1> convelint(true);
    if (gradphi_norm > 1e-8) convelint.Update(signphi / gradphi_norm, gradphi);
    // convective term
    //    double conv_phi = convelint.Dot(gradphi);

    // add Gauss point contribution to denominator
    // TODO: mit norm_gradphizero Sussman-Style
    ele_dom += (fac * deriv_sign * deriv_sign * norm_gradphizero);
    ele_nom -= (fac * deriv_sign * (phinp - phizero) / my::scatraparatimint_->Time());
    //    ele_nom -= (fac * deriv_sign * (-conv_phi + signphi)); // gecheckt
  }  // end: loop all Gauss points

  // compute penalty parameter
  if (std::abs(ele_dom) > 1.0e-9) penalty = ele_nom / ele_dom;

  return;
}


/*------------------------------------------------------------------- *
 |  calculation of penalty term on rhs                rasthofer 12/13 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, unsigned probDim>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype, probDim>::CalcRHSPenalty(
    Epetra_SerialDenseVector& erhs, const double fac, const double penalty, const double deriv_sign,
    const double norm_gradphizero)
{
  double vpenalty = fac * my::scatraparatimint_->Dt() * penalty * deriv_sign * norm_gradphizero;

  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_;

    erhs[fvi] += vpenalty * my::funct_(vi);
  }

  return;
}


/*-------------------------------------------------------------------------*
 | calculate system matrix and rhs for velocity projection rasthofer 12/13 |
 *-------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, unsigned probDim>
void DRT::ELEMENTS::ScaTraEleCalcLsReinit<distype, probDim>::SysmatNodalVel(
    const int dir,                   ///< current spatial direction
    Epetra_SerialDenseMatrix& emat,  ///< element matrix to calculate
    Epetra_SerialDenseVector& erhs   ///< element rhs to calculate
)
{
  //----------------------------------------------------------------------
  // calculation of element volume for characteristic element length
  //----------------------------------------------------------------------
  // use one-point Gauss rule to do calculations at the element center
  DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints_center(
      SCATRA::DisTypeToStabGaussRule<distype>::rule);

  // volume of the element (2D: element surface area; 1D: element length)
  // (Integration of f(x) = 1 gives exactly the volume/surface/length of element)
  const double vol = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints_center, 0);

  //----------------------------------------------------------------------
  // calculation of characteristic element length
  //----------------------------------------------------------------------

  // get gradient of initial phi at element center
  LINALG::Matrix<my::nsd_, 1> gradphizero(true);
  gradphizero.Multiply(my::derxy_, ephizero_[0]);

  // get characteristic element length
  const double charelelength = CalcCharEleLengthReinit(vol, gradphizero);

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // initial phi at Gauss point
    double phizero = 0.0;
    phizero = my::funct_.Dot(ephizero_[0]);
    // and corresponding gradient
    gradphizero.Clear();
    gradphizero.Multiply(my::derxy_, ephizero_[0]);

    // current phi at Gauss point
    double phinp = 0.0;
    phinp = my::funct_.Dot(my::ephinp_[0]);
    // gradient of current scalar
    LINALG::Matrix<my::nsd_, 1> gradphi(true);
    gradphi.Multiply(my::derxy_, my::ephinp_[0]);
    // get norm
    const double gradphi_norm = gradphi.Norm2();

    // TODO: remove
    //    if (std::abs(my::ephinp_[0](0,0)-my::ephinp_[0](1,0))>1.0e-10 or
    //        std::abs(my::ephinp_[0](2,0)-my::ephinp_[0](3,0))>1.0e-10 or
    //        std::abs(my::ephinp_[0](4,0)-my::ephinp_[0](5,0))>1.0e-10 or
    //        std::abs(my::ephinp_[0](6,0)-my::ephinp_[0](7,0))>1.0e-10)
    //    {
    //        std::cout << my::ephinp_[0] << std::endl;
    //        dserror("ENDE");
    //    }

    // get velocity at element center
    LINALG::Matrix<my::nsd_, 1> convelint(true);
    if (lsreinitparams_->ReinitType() == INPAR::SCATRA::reinitaction_sussman)
    {
      // get sign function
      double signphi = 0.0;
      SignFunction(signphi, charelelength, phizero, gradphizero, phinp, gradphi);

      if (gradphi_norm > 1e-8) convelint.Update(signphi / gradphi_norm, gradphi);
    }

    //------------------------------------------------
    // element matrix
    //------------------------------------------------
    my::CalcMatMass(emat, 0, fac, 1.0);

    //------------------------------------------------
    // add dissipation for smooth fields
    //------------------------------------------------
    // should not be used together with lumping
    // prevented by dserror in lsreinit parameters
    if (lsreinitparams_->ProjectDiff() > 0.0)
    {
      const double diff = (lsreinitparams_->ProjectDiff()) * charelelength * charelelength;
      my::diffmanager_->SetIsotropicDiff(diff, 0);
      my::CalcMatDiff(emat, 0, fac);
    }

    //------------------------------------------------
    // element rhs
    //------------------------------------------------
    // distinguish reinitialization
    switch (lsreinitparams_->ReinitType())
    {
      case INPAR::SCATRA::reinitaction_sussman:
      {
        my::CalcRHSHistAndSource(erhs, 0, fac, convelint(dir, 0));
        break;
      }
      case INPAR::SCATRA::reinitaction_ellipticeq:
      {
        my::CalcRHSHistAndSource(erhs, 0, fac, gradphi(dir, 0));
        break;
      }
      default:
        break;
    }
  }  // loop Gauss points

  // do lumping: row sum
  if (lsreinitparams_->Lumping())
  {
    for (unsigned vi = 0; vi < my::nen_; ++vi)
    {
      const int fvi = vi * my::numdofpernode_;

      double sum = 0.0;
      // loop all columns
      for (unsigned ui = 0; ui < my::nen_; ++ui)
      {
        const int fui = ui * my::numdofpernode_;
        sum += emat(fvi, fui);
        // reset
        emat(fvi, fui) = 0.0;
      }

      emat(fvi, fvi) = sum;
    }
  }

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::quad8,2>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::hex20,3>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::wedge6,3>;
template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcLsReinit<DRT::Element::nurbs27,3>;
