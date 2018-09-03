/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_elch.cpp

\brief evaluation of ScaTra elements for ion-transport equation

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_elch.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"
#include "scatra_ele_utils_elch.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcElch<distype>::ScaTraEleCalcElch(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname),
      elchparams_(DRT::ELEMENTS::ScaTraEleParameterElch::Instance(
          disname)),  // parameter class for electrochemistry problems
      utils_(DRT::ELEMENTS::ScaTraEleUtilsElch<distype>::Instance(numdofpernode, numscal, disname))
{
  // replace standard scatra diffusion manager by elch diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElch(my::numscal_));

  // safety check
  if (not my::scatraparatimint_->IsIncremental())
    dserror(
        "Since the ion-transport equations are non-linear, it can be solved only incrementally!!");

  return;
}


/*----------------------------------------------------------------------*
 | Action type: Evaluate                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Evaluate(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra)
{
  // call base class routine
  my::Evaluate(ele, params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
      elevec2_epetra, elevec3_epetra);

  // for certain ELCH problem formulations we have to provide
  // additional flux terms / currents across Dirichlet boundaries
  if (elchparams_->BoundaryFluxCoupling())
    CorrectionForFluxAcrossDC(discretization, la[0].lm_, elemat1_epetra, elevec1_epetra);

  return 0;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 ehrl  08/08|
*----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::Sysmat(
    DRT::Element* ele,                   ///< the element whose matrix is calculated
    Epetra_SerialDenseMatrix& emat,      ///< element matrix to calculate
    Epetra_SerialDenseVector& erhs,      ///< element rhs to calculate
    Epetra_SerialDenseVector& subgrdiff  ///< subgrid-diff.-scaling vector
)
{
  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  const double vol = my::EvalShapeFuncAndDerivsAtEleCenter();

  //-----------------------------------------------------------------------------------------------
  // calculate material and stabilization parameters (one per transported scalar) at element center
  //-----------------------------------------------------------------------------------------------
  // density at t_(n) (one per transported scalar)
  std::vector<double> densn(my::numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F) (one per transported scalar)
  std::vector<double> densnp(my::numscal_, 1.0);
  // density at t_(n+alpha_M) (one per transported scalar)
  std::vector<double> densam(my::numscal_, 1.0);

  // fluid viscosity
  double visc(0.0);

  // stabilization variables
  std::vector<double> tau(my::numscal_, 0.);
  std::vector<LINALG::Matrix<my::nen_, 1>> tauderpot(
      my::numscal_, LINALG::Matrix<my::nen_, 1>(true));

  if (not my::scatrapara_->MatGP() or not my::scatrapara_->TauGP())
  {
    // set internal variables at element center
    SetInternalVariablesForMatAndRHS();

    // material parameters at element center
    GetMaterialParams(ele, densn, densnp, densam, visc);

    if (not my::scatrapara_->TauGP()) PrepareStabilization(tau, tauderpot, densnp, vol);
  }

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    SetInternalVariablesForMatAndRHS();

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (my::scatrapara_->MatGP()) GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

    //-------------------------------------------------------------------------------------
    // calculate stabilization parameters (one per transported scalar) at integration point
    //-------------------------------------------------------------------------------------
    if (my::scatrapara_->TauGP()) PrepareStabilization(tau, tauderpot, densnp, vol);

    //-----------------------------------------------------------------------------------
    // calculate contributions to element matrix and right-hand side at integration point
    //-----------------------------------------------------------------------------------
    const double timefacfac = my::scatraparatimint_->TimeFac() * fac;
    const double rhsfac = my::scatraparatimint_->TimeFacRhs() * fac;

    // loop all scalars
    // deal with a system of transported scalars
    for (int k = 0; k < my::numscal_; ++k)
    {
      const double taufac = tau[k] * fac;
      const double timetaufac = my::scatraparatimint_->TimeFac() * taufac;
      const double rhstaufac = my::scatraparatimint_->TimeFacRhsTau() * taufac;

      // compute rhs containing bodyforce (divided by specific heat capacity) and,
      // for temperature equation, the time derivative of thermodynamic pressure,
      // if not constant, and for temperature equation of a reactive
      // equation system, the reaction-rate term
      double rhsint(0.0);
      my::GetRhsInt(rhsint, densnp[k], k);

      // Compute element matrix and rhs
      CalcMatAndRhs(emat, erhs, k, fac, timefacfac, rhsfac, taufac, timetaufac, rhstaufac,
          tauderpot[k], rhsint);
    }  // end loop over scalar

    // Compute element matrix and rhs
    CalcMatAndRhsOutsideScalarLoop(emat, erhs, fac, timefacfac, rhsfac);
  }

  return;
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Potential equation ENC                                       ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalcMatPotEquENC(
    Epetra_SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int k,                     //!< index of current scalar
    const double fac,                //!< domain-integration factor
    const double alphaf              //!< time factor for ENC
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    for (unsigned ui = 0; ui < my::nen_; ++ui)
    {
      // linearization of the transference number in the conduction term (transport equation)
      //
      // (w, sum(z_k c_k))
      //
      // electroneutrality condition (only derivative w.r.t. concentration c_k)
      emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
          alphaf * DiffManager()->GetValence(k) * fac * my::funct_(vi) * my::funct_(ui);
    }
  }

  return;
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Potential equation ENC                                         ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalcRhsPotEquENC(
    Epetra_SerialDenseVector& erhs,  //!< element vector to be filled
    const int k,                     //!< index of current scalar
    const double fac,                //!< domain-integration factor
    const double conint              //!< concentration at GP
)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
    // electroneutrality condition
    // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
    erhs[vi * my::numdofpernode_ + my::numscal_] -=
        DiffManager()->GetValence(k) * fac * my::funct_(vi) * conint;

  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::pyramid5>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs27>;
