/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of ScaTra elements for ion-transport equation

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "4C_scatra_ele_calc_elch.hpp"

#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_scatra_ele_utils_elch.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcElch<distype, probdim>::ScaTraEleCalcElch(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(
          numdofpernode, numscal, disname),
      elchparams_(DRT::ELEMENTS::ScaTraEleParameterElch::Instance(
          disname)),  // parameter class for electrochemistry problems
      utils_(DRT::ELEMENTS::ScaTraEleUtilsElch<distype>::Instance(numdofpernode, numscal, disname))
{
  // replace standard scatra diffusion manager by elch diffusion manager
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElch(my::numscal_));

  // safety check
  if (not my::scatraparatimint_->IsIncremental())
    FOUR_C_THROW(
        "Since the ion-transport equations are non-linear, it can be solved only incrementally!!");
}


/*----------------------------------------------------------------------*
 | Action type: Evaluate                                     ehrl 01/14 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcElch<distype, probdim>::Evaluate(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // call base class routine
  my::Evaluate(ele, params, discretization, la, elemat1_epetra, elemat2_epetra, elevec1_epetra,
      elevec2_epetra, elevec3_epetra);

  // for certain ELCH problem formulations we have to provide
  // additional flux terms / currents across Dirichlet boundaries
  if (elchparams_->boundary_flux_coupling())
    correction_for_flux_across_dc(discretization, la[0].lm_, elemat1_epetra, elevec1_epetra);

  return 0;
}


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs (public)                 ehrl  08/08|
*----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype, probdim>::sysmat(
    DRT::Element* ele,                          ///< the element whose matrix is calculated
    CORE::LINALG::SerialDenseMatrix& emat,      ///< element matrix to calculate
    CORE::LINALG::SerialDenseVector& erhs,      ///< element rhs to calculate
    CORE::LINALG::SerialDenseVector& subgrdiff  ///< subgrid-diff.-scaling vector
)
{
  //----------------------------------------------------------------------
  // calculation of element volume both for tau at ele. cent. and int. pt.
  //----------------------------------------------------------------------
  const double vol = my::eval_shape_func_and_derivs_at_ele_center();

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
  std::vector<CORE::LINALG::Matrix<nen_, 1>> tauderpot(
      my::numscal_, CORE::LINALG::Matrix<nen_, 1>(true));

  if (not my::scatrapara_->MatGP() or not my::scatrapara_->TauGP())
  {
    // set internal variables at element center
    set_internal_variables_for_mat_and_rhs();

    // material parameters at element center
    get_material_params(ele, densn, densnp, densam, visc);

    if (not my::scatrapara_->TauGP()) prepare_stabilization(tau, tauderpot, densnp, vol);
  }

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::eval_shape_func_and_derivs_at_int_point(intpoints, iquad);

    set_internal_variables_for_mat_and_rhs();

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (my::scatrapara_->MatGP()) get_material_params(ele, densn, densnp, densam, visc, iquad);

    //-------------------------------------------------------------------------------------
    // calculate stabilization parameters (one per transported scalar) at integration point
    //-------------------------------------------------------------------------------------
    if (my::scatrapara_->TauGP()) prepare_stabilization(tau, tauderpot, densnp, vol);

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
      my::get_rhs_int(rhsint, densnp[k], k);

      // Compute element matrix and rhs
      calc_mat_and_rhs(emat, erhs, k, fac, timefacfac, rhsfac, taufac, timetaufac, rhstaufac,
          tauderpot[k], rhsint);
    }  // end loop over scalar

    // Compute element matrix and rhs
    calc_mat_and_rhs_outside_scalar_loop(emat, erhs, fac, timefacfac, rhsfac);
  }
}


/*----------------------------------------------------------------------------------*
|  CalcMat: Potential equation ENC                                       ehrl  02/14|
*-----------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype, probdim>::calc_mat_pot_equ_enc(
    CORE::LINALG::SerialDenseMatrix& emat,  //!< element matrix to be filled
    const int k,                            //!< index of current scalar
    const double fac,                       //!< domain-integration factor
    const double alphaf                     //!< time factor for ENC
)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      // linearization of the transference number in the conduction term (transport equation)
      //
      // (w, sum(z_k c_k))
      //
      // electroneutrality condition (only derivative w.r.t. concentration c_k)
      emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
          alphaf * diff_manager()->GetValence(k) * fac * my::funct_(vi) * my::funct_(ui);
    }
  }
}


/*-------------------------------------------------------------------------------------*
 |  CalcRhs: Potential equation ENC                                         ehrl 11/13 |
 *-------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype, probdim>::calc_rhs_pot_equ_enc(
    CORE::LINALG::SerialDenseVector& erhs,  //!< element vector to be filled
    const int k,                            //!< index of current scalar
    const double fac,                       //!< domain-integration factor
    const double conint                     //!< concentration at GP
)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    // electroneutrality condition
    // for incremental formulation, there is the residuum on the rhs! : 0-sum(z_k c_k)
    erhs[vi * my::numdofpernode_ + my::numscal_] -=
        diff_manager()->GetValence(k) * fac * my::funct_(vi) * conint;
  }
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcElch<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
