/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge
within isothermal electrodes

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "baci_scatra_ele_calc_elch_electrode.H"

#include "baci_mat_material.H"
#include "baci_scatra_ele_parameter_std.H"
#include "baci_scatra_ele_parameter_timint.H"
#include "baci_scatra_ele_utils_elch_electrode.H"
#include "baci_utils_singleton_owner.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcElchElectrode<distype, probdim>>(
            new ScaTraEleCalcElchElectrode<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::ScaTraEleCalcElchElectrode(
    const int numdofpernode, const int numscal, const std::string& disname)
    : myelch::ScaTraEleCalcElch(numdofpernode, numscal, disname)
{
  // replace elch diffusion manager by diffusion manager for electrodes
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerElchElectrode(my::numscal_));

  // replace elch internal variable manager by internal variable manager for electrodes
  my::scatravarmanager_ =
      Teuchos::rcp(new ScaTraEleInternalVariableManagerElchElectrode<nsd_, nen_>(
          my::numscal_, myelch::elchparams_));

  // replace elch utility class by utility class for electrodes
  myelch::utils_ = DRT::ELEMENTS::ScaTraEleUtilsElchElectrode<distype>::Instance(
      numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalcMatAndRhs(
    CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs, const int k,
    const double fac, const double timefacfac, const double rhsfac, const double taufac,
    const double timetaufac, const double rhstaufac, CORE::LINALG::Matrix<nen_, 1>& tauderpot,
    double& rhsint)
{
  //----------------------------------------------------------------------
  // 1) element matrix: instationary terms arising from transport equation
  //----------------------------------------------------------------------

  if (not my::scatraparatimint_->IsStationary())
    // 1a) element matrix: standard Galerkin mass term
    my::CalcMatMass(emat, k, fac, 1.);

  //--------------------------------------------------------------------
  // 2) element matrix: stationary terms arising from transport equation
  //--------------------------------------------------------------------

  // 2a) element matrix: standard Galerkin diffusive term
  my::CalcMatDiff(emat, k, timefacfac);

  // 2b) element matrix: additional term arising from concentration dependency of diffusion
  // coefficient
  CalcMatDiffCoeffLin(emat, k, timefacfac, VarManager()->GradPhi(k), 1.);

  // 2c) element matrix: conservative part of convective term, needed for deforming electrodes,
  //                     i.e., for scalar-structure interaction
  double vdiv(0.);
  if (my::scatrapara_->IsConservative())
  {
    my::GetDivergence(vdiv, my::evelnp_);
    my::CalcMatConvAddCons(emat, k, timefacfac, vdiv, 1.);
  }

  //----------------------------------------------------------------------------
  // 3) element right hand side vector (negative residual of nonlinear problem):
  //    terms arising from transport equation
  //----------------------------------------------------------------------------

  // 3a) element rhs: standard Galerkin contributions from non-history part of instationary term if
  // needed
  if (not my::scatraparatimint_->IsStationary()) my::CalcRHSLinMass(erhs, k, rhsfac, fac, 1., 1.);

  // 3b) element rhs: standard Galerkin contributions from rhsint vector (contains body force vector
  // and history vector) need to adapt rhsint vector to time integration scheme first
  my::ComputeRhsInt(rhsint, 1., 1., VarManager()->Hist(k));
  my::CalcRHSHistAndSource(erhs, k, fac, rhsint);

  // 3c) element rhs: standard Galerkin diffusion term
  my::CalcRHSDiff(erhs, k, rhsfac);

  // 3d) element rhs: conservative part of convective term, needed for deforming electrodes,
  //                  i.e., for scalar-structure interaction
  if (my::scatrapara_->IsConservative())
    CalcRhsConservativePartOfConvectiveTerm(erhs, k, rhsfac, vdiv);

  //----------------------------------------------------------------------------
  // 4) element matrix: stationary terms arising from potential equation
  // 5) element right hand side vector (negative residual of nonlinear problem):
  //    terms arising from potential equation
  //----------------------------------------------------------------------------
  // see function CalcMatAndRhsOutsideScalarLoop()
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalcMatAndRhsOutsideScalarLoop(
    CORE::LINALG::SerialDenseMatrix& emat, CORE::LINALG::SerialDenseVector& erhs, const double fac,
    const double timefacfac, const double rhsfac)
{
  //--------------------------------------------------------------------
  // 4) element matrix: stationary terms arising from potential equation
  //--------------------------------------------------------------------

  // element matrix: standard Galerkin terms from potential equation
  CalcMatPotEquDiviOhm(emat, timefacfac, VarManager()->InvF(), VarManager()->GradPot(), 1.);

  //----------------------------------------------------------------------------
  // 5) element right hand side vector (negative residual of nonlinear problem):
  //    terms arising from potential equation
  //----------------------------------------------------------------------------

  // element rhs: standard Galerkin terms from potential equation
  CalcRhsPotEquDiviOhm(erhs, rhsfac, VarManager()->InvF(), VarManager()->GradPot(), 1.);

  // safety check
  if (my::bodyforce_[my::numscal_].Dot(my::funct_) != 0.0)
    dserror("body force not implemented for potential equation");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalcDiffODMesh(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const int ndofpernodemesh,
    const double diffcoeff, const double fac, const double rhsfac, const double J,
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi, const CORE::LINALG::Matrix<nsd_, 1>& convelint,
    const CORE::LINALG::Matrix<1, nsd_ * nen_>& dJ_dmesh)
{
  // safety check
  if (k != 0) dserror("Invalid species index!");

  // call base class routine to compute linearizations of diffusion term w.r.t. structural
  // displacements
  my::CalcDiffODMesh(
      emat, 0, ndofpernodemesh, diffcoeff, fac, rhsfac, J, gradphi, convelint, dJ_dmesh);

  // call base class routine again to compute linearizations of Ohmic overpotential w.r.t.
  // structural displacements
  my::CalcDiffODMesh(emat, 1, ndofpernodemesh, VarManager()->InvF() * DiffManager()->GetCond(), fac,
      rhsfac, J, VarManager()->GradPot(), convelint, dJ_dmesh);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalcMatDiffCoeffLin(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac,
    const CORE::LINALG::Matrix<nsd_, 1>& gradphi, const double scalar)
{
  // linearization of diffusion coefficient in ionic diffusion term (transport equation):
  //
  // (nabla w, D(D(c)) nabla c)
  //
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      double laplawfrhs_gradphi(0.);
      my::GetLaplacianWeakFormRHS(laplawfrhs_gradphi, gradphi, vi);

      emat(vi * my::numdofpernode_ + k, ui * my::numdofpernode_ + k) +=
          scalar * timefacfac * DiffManager()->GetConcDerivIsoDiffCoef(k, k) * laplawfrhs_gradphi *
          my::funct_(ui);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalcMatPotEquDiviOhm(
    CORE::LINALG::SerialDenseMatrix& emat, const double timefacfac, const double invf,
    const CORE::LINALG::Matrix<nsd_, 1>& gradpot, const double scalar)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      double laplawf(0.);
      my::GetLaplacianWeakForm(laplawf, ui, vi);

      // linearization of the ohmic term
      //
      // (grad w, 1/F kappa D(grad pot))
      //
      emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + my::numscal_) +=
          scalar * timefacfac * invf * DiffManager()->GetCond() * laplawf;

      for (int iscal = 0; iscal < my::numscal_; ++iscal)
      {
        double laplawfrhs_gradpot(0.);
        my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot, gradpot, vi);

        // linearization of the ohmic term with respect to conductivity
        //
        // (grad w, 1/F kappa D(grad pot))
        //
        emat(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + iscal) +=
            scalar * timefacfac * invf * DiffManager()->GetConcDerivCond(iscal) * my::funct_(ui) *
            laplawfrhs_gradpot;
      }
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype,
    probdim>::CalcRhsConservativePartOfConvectiveTerm(CORE::LINALG::SerialDenseVector& erhs,
    const int k, const double rhsfac, const double vdiv)
{
  double vrhs = rhsfac * my::scatravarmanager_->Phinp(k) * vdiv;
  for (unsigned vi = 0; vi < nen_; ++vi) erhs[vi * my::numdofpernode_ + k] -= vrhs * my::funct_(vi);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalcRhsPotEquDiviOhm(
    CORE::LINALG::SerialDenseVector& erhs, const double rhsfac, const double invf,
    const CORE::LINALG::Matrix<nsd_, 1>& gradpot, const double scalar)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    double laplawfrhs_gradpot(0.);
    my::GetLaplacianWeakFormRHS(laplawfrhs_gradpot, gradpot, vi);

    erhs[vi * my::numdofpernode_ + my::numscal_] -=
        scalar * rhsfac * invf * DiffManager()->GetCond() * laplawfrhs_gradpot;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::GetMaterialParams(
    const DRT::Element* ele, std::vector<double>& densn, std::vector<double>& densnp,
    std::vector<double>& densam, double& visc, const int iquad)
{
  // get material
  Teuchos::RCP<const MAT::Material> material = ele->Material();

  // evaluate electrode material
  if (material->MaterialType() == INPAR::MAT::m_electrode)
  {
    Utils()->MatElectrode(
        material, VarManager()->Phinp(0), VarManager()->Temperature(), DiffManager());
  }
  else
    dserror("Material type not supported!");
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::hex8, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tet10, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::pyramid5, 3>;
