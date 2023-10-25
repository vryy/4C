/*----------------------------------------------------------------------*/
/*! \file

\brief scatra_ele_calc_cardiac_monodomain.cpp

\level 2

 *----------------------------------------------------------------------*/


#include "baci_scatra_ele_calc_cardiac_monodomain.H"

#include "baci_inpar_cardiac_monodomain.H"
#include "baci_lib_discret.H"
#include "baci_lib_element.H"
#include "baci_lib_globalproblem.H"
#include "baci_mat_list.H"
#include "baci_mat_myocard.H"
#include "baci_scatra_ele_parameter_std.H"
#include "baci_scatra_ele_parameter_timint.H"
#include "baci_utils_singleton_owner.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::ScaTraEleCalcCardiacMonodomain(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(
          numdofpernode, numscal, disname),
      DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::ScaTraEleCalcAniso(
          numdofpernode, numscal, disname),
      DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype, probdim>::ScaTraEleCalcAdvReac(
          numdofpernode, numscal, disname)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcCardiacMonodomain<distype, probdim>>(
            new ScaTraEleCalcCardiacMonodomain<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    ljag 06/14 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::Materials(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    double& densn,                                     //!< density at t_(n)
    double& densnp,                                    //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                                    //!< density at t_(n+alpha_M)
    double& visc,                                      //!< fluid viscosity
    const int iquad                                    //!< id of current gauss point

)
{
  // safety check
  if (material->MaterialType() != INPAR::MAT::m_myocard) dserror("Material type is not supported");

  // safety check
  Teuchos::RCP<MAT::Myocard> actmat =
      Teuchos::rcp_dynamic_cast<MAT::Myocard>(Teuchos::rcp_const_cast<MAT::Material>(material));
  if (actmat->GetNumberOfGP() != 1 and not my::scatrapara_->MatGP())
  {
    actmat->SetGP(1);
    actmat->ResizeInternalStateVariables();
  }
  MatMyocard(material, k, densn, densnp, densam, visc, iquad);

  return;
}


/*----------------------------------------------------------------------*
 |  Material ScaTra                                          ljag 06/14 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::MatMyocard(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    double& densn,                                     //!< density at t_(n)
    double& densnp,                                    //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                                    //!< density at t_(n+alpha_M)
    double& visc,                                      //!< fluid viscosity
    const int iquad                                    //!< id of current gauss point (default = -1)
)
{
  const Teuchos::RCP<const MAT::Myocard>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::Myocard>(material);

  // dynamic cast to Advanced_Reaction-specific reaction manager
  Teuchos::RCP<ScaTraEleReaManagerAdvReac> advreamanager =
      Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerAdvReac>(my::reamanager_);

  // dynamic cast to anisotropic diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerAniso<nsd_>> diffmanageraniso =
      Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerAniso<nsd_>>(my::diffmanager_);

  // get constant diffusivity
  CORE::LINALG::Matrix<nsd_, nsd_> difftensor(true);
  actmat->Diffusivity(difftensor);

  diffmanageraniso->SetAnisotropicDiff(difftensor, k);

  // clear
  advreamanager->Clear(my::numscal_);

  if (my::scatrapara_->SemiImplicit())
  {
    // get membrane potential at n at integration point
    const double phin = my::scatravarmanager_->Phin(k);
    const double phinp = my::scatravarmanager_->Phinp(k);
    // get reaction coefficient
    double react = -actmat->ReaCoeffN(phin, my::scatraparatimint_->Dt(), iquad);
    if (my::scatraparatimint_->IsGenAlpha())
      react *= my::scatraparatimint_->Dt() / my::scatraparatimint_->TimeFac();
    advreamanager->AddToReaBodyForce(react, k);
    advreamanager->AddToReaBodyForceDerivMatrix(0.0, k, k);
    actmat->ReaCoeff(phinp, my::scatraparatimint_->Dt(), iquad);
  }
  else
  {
    // get membrane potential at n+1 or n+alpha_F at integration point
    const double phinp = my::scatravarmanager_->Phinp(k);
    // get reaction coefficient
    advreamanager->AddToReaBodyForce(
        -actmat->ReaCoeff(phinp, my::scatraparatimint_->Dt(), iquad), k);
    advreamanager->AddToReaBodyForceDerivMatrix(
        -actmat->ReaCoeffDeriv(phinp, my::scatraparatimint_->Dt(), iquad), k, k);
  }

  return;
}  // ScaTraEleCalcCardiacMonodomain<distype>::MatMyocard


/*----------------------------------------------------------------------*
|  calculate system matrix and rhs for ep                 hoermann 06/16|
*----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::Sysmat(
    DRT::Element* ele,                          ///< the element whose matrix is calculated
    CORE::LINALG::SerialDenseMatrix& emat,      ///< element matrix to calculate
    CORE::LINALG::SerialDenseVector& erhs,      ///< element rhs to calculate
    CORE::LINALG::SerialDenseVector& subgrdiff  ///< subgrid-diff.-scaling vector
)
{
  // density at t_(n) (one per transported scalar)
  std::vector<double> densn(my::numscal_, 1.0);
  // density at t_(n+1) or t_(n+alpha_F) (one per transported scalar)
  std::vector<double> densnp(my::numscal_, 1.0);
  // density at t_(n+alpha_M) (one per transported scalar)
  std::vector<double> densam(my::numscal_, 1.0);

  // fluid viscosity
  double visc(0.0);

  // calculation of material parameter at element center
  if (not my::scatrapara_->MatGP())
  {
    advreac::EvalShapeFuncAndDerivsAtEleCenter();
    // set Gauss point variables needed for evaluation of mat and rhs
    my::SetInternalVariablesForMatAndRHS();
    advreac::GetMaterialParams(ele, densn, densnp, densam, visc);
  }
  // calculation of material at integration points (different number of integration points possible)
  else
  {
    int deg = 0;
    if (ele->Degree() == 1)
      deg = 4 * ele->Degree();
    else
      deg = 3 * ele->Degree();
    const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
        SCATRA::DisTypeToMatGaussRule<distype>::GetGaussRule(deg));

    // loop over integration points
    for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

      // set gauss point variables needed for evaluation of mat and rhs
      my::SetInternalVariablesForMatAndRHS();

      // get material parameters (evaluation at integration point)
      advreac::GetMaterialParams(ele, densn, densnp, densam, visc, iquad);

      // loop all scalars
      for (int k = 0; k < my::numscal_; ++k)  // deal with a system of transported scalars
      {
        double rhsint(0.0);
        advreac::GetRhsInt(rhsint, densnp[k], k);

        CORE::LINALG::Matrix<nen_, 1> dummy(true);
        const double timefacfac = my::scatraparatimint_->TimeFac() * fac;

        // reactive terms on integration point on rhs
        my::ComputeRhsInt(rhsint, densam[k], densnp[k], my::scatravarmanager_->Hist(k));

        // standard Galerkin transient, old part of rhs and bodyforce term
        my::CalcRHSHistAndSource(erhs, k, fac, rhsint);

        // element matrix: reactive term
        advreac::CalcMatReact(emat, k, timefacfac, 0., 0., densnp[k], dummy, dummy);
      }
    }
  }

  //----------------------------------------------------------------------
  // integration loop for one element
  //----------------------------------------------------------------------
  // integration points and weights
  const CORE::DRT::UTILS::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // set gauss point variables needed for evaluation of mat and rhs
    my::SetInternalVariablesForMatAndRHS();

    // loop all scalars
    for (int k = 0; k < my::numscal_; ++k)  // deal with a system of transported scalars
    {
      // compute rhs containing bodyforce
      double rhsint(0.0);
      advreac::GetRhsInt(rhsint, densnp[k], k);

      // integration factors
      const double timefacfac = my::scatraparatimint_->TimeFac() * fac;

      //----------------------------------------------------------------
      // 1) element matrix: stationary terms
      //----------------------------------------------------------------

      // calculation of diffusive element matrix
      aniso::CalcMatDiff(emat, k, timefacfac);

      //----------------------------------------------------------------
      // 2) element matrix: instationary terms
      //----------------------------------------------------------------

      if (not my::scatraparatimint_->IsStationary()) my::CalcMatMass(emat, k, fac, densam[k]);

      //----------------------------------------------------------------
      // 3) element matrix: reactive term
      //----------------------------------------------------------------

      CORE::LINALG::Matrix<nen_, 1> dummy(true);
      if (not my::scatrapara_->MatGP())
        advreac::CalcMatReact(emat, k, timefacfac, 0., 0., densnp[k], dummy, dummy);

      //----------------------------------------------------------------
      // 5) element right hand side
      //----------------------------------------------------------------
      //----------------------------------------------------------------
      // computation of bodyforce (and potentially history) term,
      // residual, integration factors and standard Galerkin transient
      // term (if required) on right hand side depending on respective
      // (non-)incremental stationary or time-integration scheme
      //----------------------------------------------------------------
      double rhsfac = my::scatraparatimint_->TimeFacRhs() * fac;

      if (my::scatraparatimint_->IsIncremental() and not my::scatraparatimint_->IsStationary())
        my::CalcRHSLinMass(erhs, k, rhsfac, fac, densam[k], densnp[k]);


      if (not my::scatrapara_->MatGP())
      {
        my::ComputeRhsInt(rhsint, densam[k], densnp[k], my::scatravarmanager_->Hist(k));
        // standard Galerkin transient, old part of rhs and bodyforce term
        my::CalcRHSHistAndSource(erhs, k, fac, rhsint);
      }

      //----------------------------------------------------------------
      // standard Galerkin terms on right hand side
      //----------------------------------------------------------------

      // diffusive term
      aniso::CalcRHSDiff(erhs, k, rhsfac);

    }  // end loop all scalars
  }    // end loop Gauss points

  return;
}


/*----------------------------------------------------------------------*
 | extract element based or nodal values                 hoermann 06/16 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<distype, probdim>::ExtractElementAndNodeValues(
    DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la)
{
  my::ExtractElementAndNodeValues(ele, params, discretization, la);

  // extract additional local values from global vector
  Teuchos::RCP<const Epetra_Vector> phin = discretization.GetState("phin");
  if (phin == Teuchos::null) dserror("Cannot get state vector 'phin'");
  DRT::UTILS::ExtractMyValues<CORE::LINALG::Matrix<nen_, 1>>(*phin, my::ephin_, la[0].lm_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::quad4, 3>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::hex8, 3>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::tet10, 3>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::pyramid5, 3>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcCardiacMonodomain<CORE::FE::CellType::nurbs27>;
