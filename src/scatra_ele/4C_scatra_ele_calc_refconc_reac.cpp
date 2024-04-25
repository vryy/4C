/*----------------------------------------------------------------------*/
/*! \file
\brief main file containing routines for calculation of scatra element formulated in reference
concentrations and with advanced reaction terms

\level 3

 *----------------------------------------------------------------------*/


#include "4C_scatra_ele_calc_refconc_reac.hpp"

#include "4C_mat_list_reactions.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::ScaTraEleCalcRefConcReac(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode, numscal, disname),
      DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(
          numdofpernode, numscal, disname),
      j_(1.0),
      c_inv_(true),
      d_jd_x_(true)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>*
DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcRefConcReac<distype>>(
            new ScaTraEleCalcRefConcReac<distype>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}

//!
/*----------------------------------------------------------------------------*
 |  Set reac. body force, reaction coefficient and derivatives     thon 02/16 |
 *---------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::SetAdvancedReactionTerms(
    const int k,                                            //!< index of current scalar
    const Teuchos::RCP<MAT::MatListReactions> matreaclist,  //!< index of current scalar
    const double* gpcoord                                   //!< current Gauss-point coordinates
)
{
  const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = advreac::ReaManager();

  remanager->AddToReaBodyForce(
      matreaclist->CalcReaBodyForceTerm(k, my::scatravarmanager_->Phinp(), gpcoord, 1.0 / j_) * j_,
      k);

  matreaclist->CalcReaBodyForceDerivMatrix(
      k, remanager->GetReaBodyForceDerivVector(k), my::scatravarmanager_->Phinp(), gpcoord);
}

/*------------------------------------------------------------------------------------------*
 |  calculation of convective element matrix: add conservative contributions     thon 02/16 |
 *------------------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::CalcMatConvAddCons(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac, const double vdiv,
    const double densnp)
{
  FOUR_C_THROW(
      "If you want to calculate the reference concentrations the CONVFORM must be 'convective'!");
}

/*------------------------------------------------------------------------------*
 | set internal variables                                           thon 02/16  |
 *------------------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::SetInternalVariablesForMatAndRHS()
{
  // do the usual and...
  advreac::SetInternalVariablesForMatAndRHS();

  /////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////
  // spatial node coordinates
  CORE::LINALG::Matrix<nsd_, nen_> xyze(my::xyze_);
  xyze += my::edispnp_;

  //! transposed jacobian "dx/ds"
  CORE::LINALG::Matrix<nsd_, nsd_> dxds(true);
  dxds.MultiplyNT(my::deriv_, xyze);

  // deformation gradtient dx/dX = dx/ds * ds/dX = dx/ds * (dX/ds)^(-1)
  CORE::LINALG::Matrix<nsd_, nsd_> F(true);
  F.MultiplyTT(dxds, my::xij_);

  // inverse of jacobian "dx/dX"
  CORE::LINALG::Matrix<nsd_, nsd_> F_inv(true);
  j_ = F_inv.Invert(F);

  // calculate inverse of cauchy-green stress tensor
  c_inv_.MultiplyNT(F_inv, F_inv);

  ////////////////////////////////////////////////////////////////////////////////////////////////
  // calculate derivative dJ/dX by finite differences
  ////////////////////////////////////////////////////////////////////////////////////////////////
  const double epsilon = 1.0e-8;

  for (unsigned i = 0; i < 3; i++)
  {
    CORE::LINALG::Matrix<nsd_, nen_> xyze_epsilon(my::xyze_);
    for (unsigned j = 0; j < nen_; ++j) xyze_epsilon(i, j) = xyze_epsilon(i, j) + epsilon;

    CORE::LINALG::Matrix<nsd_, nsd_> xjm_epsilon(true);
    xjm_epsilon.MultiplyNT(my::deriv_, xyze_epsilon);

    CORE::LINALG::Matrix<nsd_, nsd_> xij_epsilon(true);
    xij_epsilon.Invert(xjm_epsilon);

    // dx/dX = dx/ds * ds/dX = dx/ds * (dX/ds)^(-1)
    CORE::LINALG::Matrix<nsd_, nsd_> F_epsilon(true);
    F_epsilon.MultiplyTT(dxds, xij_epsilon);

    // inverse of transposed jacobian "ds/dX"
    const double J_epsilon = F_epsilon.Determinant();
    const double dJdX_i = (J_epsilon - j_) / epsilon;

    d_jd_x_(i, 0) = dJdX_i;
  }

  return;
}

/*------------------------------------------------------------------- *
 |  calculation of diffusive element matrix                thon 02/16 |
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::CalcMatDiff(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac)
{
  CORE::LINALG::Matrix<nsd_, nsd_> Diff_tens(c_inv_);
  Diff_tens.Scale(my::diffmanager_->GetIsotropicDiff(k));

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * my::numdofpernode_ + k;

      double laplawf = 0.0;
      //      GetLaplacianWeakForm(laplawf,Diff_tens,vi,ui);
      for (unsigned j = 0; j < nsd_; j++)
      {
        for (unsigned i = 0; i < nsd_; i++)
        {
          laplawf += my::derxy_(j, vi) * Diff_tens(j, i) * my::derxy_(i, ui);
        }
      }

      emat(fvi, fui) += timefacfac * laplawf;
    }
  }



  CORE::LINALG::Matrix<nsd_, nsd_> Diff_tens2(c_inv_);
  Diff_tens2.Scale(my::diffmanager_->GetIsotropicDiff(k) / j_);

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    double laplawf2 = 0.0;
    for (unsigned j = 0; j < nsd_; j++)
    {
      for (unsigned i = 0; i < nsd_; i++)
      {
        laplawf2 += my::derxy_(j, vi) * Diff_tens2(j, i) * d_jd_x_(i);
      }
    }

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * my::numdofpernode_ + k;

      emat(fvi, fui) -= timefacfac * laplawf2 * my::funct_(ui);
    }
  }

  return;
}

/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term on right hand side     ehrl 11/13 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
void DRT::ELEMENTS::ScaTraEleCalcRefConcReac<distype>::CalcRHSDiff(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac)
{
  /////////////////////////////////////////////////////////////////////
  // \D* \grad c_0 \times \grad \phi ...
  /////////////////////////////////////////////////////////////////////
  CORE::LINALG::Matrix<nsd_, nsd_> Diff_tens(c_inv_);
  Diff_tens.Scale(my::diffmanager_->GetIsotropicDiff(k));

  const CORE::LINALG::Matrix<nsd_, 1>& gradphi = my::scatravarmanager_->GradPhi(k);

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    double laplawf(0.0);
    //    GetLaplacianWeakFormRHS(laplawf,Diff_tens,gradphi,vi);
    for (unsigned j = 0; j < nsd_; j++)
    {
      for (unsigned i = 0; i < nsd_; i++)
      {
        laplawf += my::derxy_(j, vi) * Diff_tens(j, i) * gradphi(i);
      }
    }

    erhs[fvi] -= rhsfac * laplawf;
  }

  /////////////////////////////////////////////////////////////////////
  // ... + \D* c_0/J * \grad J \times \grad \phi
  /////////////////////////////////////////////////////////////////////
  CORE::LINALG::Matrix<nsd_, nsd_> Diff_tens2(c_inv_);
  Diff_tens2.Scale(my::diffmanager_->GetIsotropicDiff(k) / j_ * my::scatravarmanager_->Phinp(k));

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    double laplawf2(0.0);
    //    GetLaplacianWeakFormRHS(laplawf2,Diff_tens2,dJdX_,vi);
    for (unsigned j = 0; j < nsd_; j++)
    {
      for (unsigned i = 0; i < nsd_; i++)
      {
        laplawf2 += my::derxy_(j, vi) * Diff_tens2(j, i) * d_jd_x_(i);
      }
    }

    erhs[fvi] += rhsfac * laplawf2;
  }

  return;
}



// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::quad4>;
// template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::hex8>;
// template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::tet10>;
// template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::nurbs9>;
// template class
// DRT::ELEMENTS::ScaTraEleCalcRefConcReac<CORE::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
