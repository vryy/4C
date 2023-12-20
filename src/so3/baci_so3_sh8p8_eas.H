/*----------------------------------------------------------------------*/
/*! \file
\brief 8-node solid shell element
\level 2
*/

/*----------------------------------------------------------------------*/
#ifndef BACI_SO3_SH8P8_EAS_H
#define BACI_SO3_SH8P8_EAS_H

/*----------------------------------------------------------------------*/
#include "baci_config.H"

#include "baci_linalg_fixedsizematrix_voigt_notation.H"
#include "baci_so3_sh8p8.H"

#include <Teuchos_SerialDenseSolver.hpp>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int NUMEAS_T>
void DRT::ELEMENTS::So_sh8p8::EasUpdateIncrementally(CORE::LINALG::SerialDenseMatrix*& oldfeas,
    CORE::LINALG::SerialDenseMatrix*& oldKaainv, CORE::LINALG::SerialDenseMatrix*& oldKad,
    CORE::LINALG::SerialDenseMatrix*& oldKap, Teuchos::RCP<CORE::LINALG::SerialDenseVector>& feas,
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Kaa,
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Kad,
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Kap, CORE::LINALG::SerialDenseMatrix*& alpha,
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& M, DRT::Container& data,
    const CORE::LINALG::Matrix<NUMDISP_, 1>& dispi, const CORE::LINALG::Matrix<NUMPRES_, 1>& presi)
{
  // retrieve history
  alpha = data.Get<CORE::LINALG::SerialDenseMatrix>("alpha");  // get old alpha
  // evaluate current (updated) EAS alphas (from history variables)
  // get stored EAS history
  oldfeas = data.Get<CORE::LINALG::SerialDenseMatrix>("feas");
  oldKaainv = data.Get<CORE::LINALG::SerialDenseMatrix>("invKaa");
  oldKad = data.Get<CORE::LINALG::SerialDenseMatrix>("Kda");  // actually k_ad
  oldKap = data.Get<CORE::LINALG::SerialDenseMatrix>("Kap");
  if (!alpha || !oldKaainv || !oldKad || !oldKap || !oldfeas) dserror("Missing EAS history-data");

  // feas^{k+1} := feas^k + k_ad^k . Ddisp^k + k_ap^k . Dpres^k
  CORE::LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMEAS_T, NUMDISP_, 1>(
      1.0, oldfeas->values(), 1.0, oldKad->values(), dispi.values());
  CORE::LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMEAS_T, NUMPRES_, 1>(
      1.0, oldfeas->values(), 1.0, oldKap->values(), presi.values());

  // alpha^{k+1} := alpha^k + Dalpha^k
  //                alpha^k - k_aa^{k;-1} . feas^{k+1}
  CORE::LINALG::DENSEFUNCTIONS::multiply<double, NUMEAS_T, NUMEAS_T, 1>(
      1.0, *alpha, -1.0, *oldKaainv, *oldfeas);

  // EAS portion of internal forces, also called enhacement vector s or Rtilde
  feas = Teuchos::rcp(new CORE::LINALG::SerialDenseVector(NUMEAS_T));
  // EAS matrix K_{alpha alpha}, also called Dtilde
  Kaa = Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(NUMEAS_T, NUMEAS_T));
  // EAS matrix K_{alpha disp}
  Kad = Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(NUMEAS_T, NUMDISP_));
  // EAS matrix K_{alpha pres}
  Kap = Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(NUMEAS_T, NUMPRES_));

  // M-operator, ie EAS shape functions
  M = Teuchos::rcp(new CORE::LINALG::SerialDenseMatrix(MAT::NUM_STRESS_3D, NUMEAS_T));

  // see you
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int NUMEAS_T>
void DRT::ELEMENTS::So_sh8p8::EasMaterialiseShapeFcts(
    const Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& M, const double& detJ0, const double& detJ,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& T0invT,
    const CORE::LINALG::SerialDenseMatrix& Mloc)
{
  // map local M to global, also enhancement is refered to element origin
  // M = detJ0/detJ T0^{-T} . M
  CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D, NUMEAS_T>(
      M->values(), detJ0 / detJ, T0invT.values(), Mloc.values());

  // watch out
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int NUMEAS_T>
void DRT::ELEMENTS::So_sh8p8::EasAddStrain(CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain,
    const Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& M,
    const CORE::LINALG::SerialDenseMatrix* alpha)
{
  // add enhanced strains = M . alpha to GL strains to "unlock" element
  CORE::LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, NUMEAS_T, 1>(
      1.0, glstrain.values(), 1.0, M->values(), alpha->values());
  // ready
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int NUMEAS_T>
void DRT::ELEMENTS::So_sh8p8::EasConstraintAndTangent(
    Teuchos::RCP<CORE::LINALG::SerialDenseVector>& feas,
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Kaa,
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Kad,
    Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Kap,
    const CORE::LINALG::Matrix<NUMDIM_, NUMDIM_>& defgradD,
    const CORE::LINALG::Matrix<NUMDIM_, NUMDIM_>& invrgtstrD,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& rcgbyrgtstr,
    const double& detdefgrad, const CORE::LINALG::Matrix<NUMDFGR_, 1>& tinvdefgrad,
    const CORE::LINALG::Matrix<NUMDFGR_, NUMDFGR_>& WmT,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& cmat,
    const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress, const double& effpressure,
    const double& detJ_w, const CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISP_>& cb,
    const CORE::LINALG::Matrix<NUMDFGR_, NUMDISP_>& defgradbydisp,
    const CORE::LINALG::Matrix<NUMPRES_, 1>& prshfct,
    const Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& M)
{
  // derivative of assumed right stretch tensor w.r.t. EAS parameters
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMEAS_T> rgtstrbyalpha;
  CORE::LINALG::DENSEFUNCTIONS::multiplyNN<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
      NUMEAS_T>(rgtstrbyalpha.values(), 2.0, rcgbyrgtstr.values(), M->values());

  // derivative of pseudo identity with respect to EAS parameters
  // I^{assd}_{CB,e} = U^{d;-1}_{CD} . U^{ass}_{DB,e}
  // WARNING: I^{assd}_{CB} and I^{assd}_{CB,e} might be non-symmetric in CB
  CORE::LINALG::Matrix<NUMDFGR_, NUMEAS_T> pseudoidentity;
  for (int e = 0; e < NUMEAS_T; ++e)
  {
    for (int CB = 0; CB < NUMDFGR_; ++CB)
    {
      const int C = VOIGT9ROW_INCONSISTENT_[CB];
      const int B = VOIGT9COL_INCONSISTENT_[CB];
      double pseudoidentity_CBe = 0.0;
      for (int D = 0; D < NUMDIM_; ++D)
      {
        const int DB = CORE::LINALG::VOIGT::IndexMappings::SymToVoigt6(D, B);
        const double DBfact = (D == B) ? 1.0 : 0.5;
        pseudoidentity_CBe += invrgtstrD(C, D) * DBfact * rgtstrbyalpha(DB, e);
      }
      pseudoidentity(CB, e) = pseudoidentity_CBe;
    }
  }

  // derivative of def.grad. with respect to e EAS parameters alpha^e
  // F_{aB,e} = F^d_{aC} . U^{d;-1}_{CD} . U^{ass}_{DB,e}
  //            = F^d_{aC} . I^{assd}_{CB,e}
  CORE::LINALG::Matrix<NUMDFGR_, NUMEAS_T> defgradbyalpha;
  for (int e = 0; e < NUMEAS_T; ++e)
  {
    for (int aB = 0; aB < NUMDFGR_; ++aB)
    {
      const int a = VOIGT9ROW_INCONSISTENT_[aB];
      const int B = VOIGT9COL_INCONSISTENT_[aB];
      double defgradbyalpha_aBe = 0.0;
      for (int C = 0; C < NUMDIM_; ++C)
      {
        const int CB = VOIGT3X3NONSYM_INCONSISTENT_[NUMDIM_ * C + B];
        defgradbyalpha_aBe += defgradD(a, C) * pseudoidentity(CB, e);
      }
      defgradbyalpha(aB, e) = defgradbyalpha_aBe;
    }
  }

  // M^T = M
  // cmat = cmat
  // detJ * w = detJ_w
  // fv = tinvdefgrad
  // B_F = defgradbydisp
  // M_F = defgradbyalpha
  // H = prshfct
  // H.p = effpressure
  // ( fv . fv^T + Wm ) = WmT
  // B = bop
  // cmat . B = cb
  // sigma = stress

  // temporary c . M
  CORE::LINALG::Matrix<MAT::NUM_STRESS_3D, NUMEAS_T> cM;
  CORE::LINALG::DENSEFUNCTIONS::multiplyNN<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D,
      NUMEAS_T>(cM.values(), cmat.values(), M->values());
  // temporary ( fv . fv^T + Wm ) . M_F
  CORE::LINALG::Matrix<NUMDFGR_, NUMEAS_T> ffwmf;
  ffwmf.MultiplyNN(WmT, defgradbyalpha);
  // temporary M_F^T . fv
  CORE::LINALG::Matrix<NUMEAS_T, 1> mff;
  mff.MultiplyTN(defgradbyalpha, tinvdefgrad);
  // temporary integration factor
  const double fac = effpressure * detdefgrad * detJ_w;

  // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
  //                     - (M_F^T . ( fv . fv^T + Wm ) . M_F) * (H . p) * detF * detJ * w
  CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMEAS_T, MAT::NUM_STRESS_3D, NUMEAS_T>(
      1.0, Kaa->values(), detJ_w, M->values(), cM.values());
  CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMEAS_T, NUMDFGR_, NUMEAS_T>(
      1.0, Kaa->values(), -fac, defgradbyalpha.values(), ffwmf.values());
  // integrate Kad: Kad += (M^T . cmat . B) * detJ * w(gp)
  //                     - (M_F^T . ( fv . fv^T + Wm ) . B_F) * (H . p) * detF * detJ * w
  CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMEAS_T, MAT::NUM_STRESS_3D, NUMDISP_>(
      1.0, Kad->values(), detJ_w, M->values(), cb.values());
  CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMEAS_T, NUMDFGR_, NUMDISP_>(
      1.0, Kad->values(), -fac, ffwmf.values(), defgradbydisp.values());
  // integrate Kap: Kap += - (M_F^T . fv . H)  * detF * detJ * w
  CORE::LINALG::DENSEFUNCTIONS::multiplyNT<double, NUMEAS_T, 1, NUMPRES_>(
      1.0, Kap->values(), -detdefgrad * detJ_w, mff.values(), prshfct.values());
  // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
  //                       - (M_F^T . fv) * (H . p)  * detF * detJ * w
  CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMEAS_T, MAT::NUM_STRESS_3D, 1>(
      1.0, feas->values(), detJ_w, M->values(), stress.values());
  CORE::LINALG::DENSEFUNCTIONS::update<double, NUMEAS_T, 1>(
      1.0, feas->values(), -fac, mff.values());

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int NUMEAS_T>
void DRT::ELEMENTS::So_sh8p8::EasCondensation(
    CORE::LINALG::Matrix<NUMDISP_, 1>* force,               ///< element internal force vector
    CORE::LINALG::Matrix<NUMDISP_, NUMDISP_>* stiffmatrix,  // element stiffness matrix
    CORE::LINALG::Matrix<NUMDISP_, NUMPRES_>* gradmatrix,   // element gradient matrix
    CORE::LINALG::Matrix<NUMPRES_, 1>* incomp,              ///< incompressibility residual
    CORE::LINALG::Matrix<NUMPRES_, NUMDISP_>* dargmatrix,   // 'transposed' element gradient matrix
    CORE::LINALG::Matrix<NUMPRES_, NUMPRES_>* stabmatrix,   // element stabilisation matrix
    CORE::LINALG::SerialDenseMatrix*& oldfeas, CORE::LINALG::SerialDenseMatrix*& oldKaainv,
    CORE::LINALG::SerialDenseMatrix*& oldKad, CORE::LINALG::SerialDenseMatrix*& oldKap,
    const Teuchos::RCP<CORE::LINALG::SerialDenseVector>& feas,
    const Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Kaa,
    const Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Kad,
    const Teuchos::RCP<CORE::LINALG::SerialDenseMatrix>& Kap)
{
  // we need the inverse of Kaa
  using ordinalType = CORE::LINALG::SerialDenseMatrix::ordinalType;
  using scalarType = CORE::LINALG::SerialDenseMatrix::scalarType;
  Teuchos::SerialDenseSolver<ordinalType, scalarType> solve_for_inverseKaa;
  solve_for_inverseKaa.setMatrix(Kaa);
  solve_for_inverseKaa.invert();

  // temporary Kda.Kaa^{-1}
  CORE::LINALG::Matrix<NUMDISP_, NUMEAS_T> KdaKaa(false);
  CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMDISP_, NUMEAS_T, NUMEAS_T>(
      KdaKaa.values(), Kad->values(), Kaa->values());
  // temporary Kpa.Kaa^{-1}
  CORE::LINALG::Matrix<NUMPRES_, NUMEAS_T> KpaKaa(false);
  CORE::LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMPRES_, NUMEAS_T, NUMEAS_T>(
      KpaKaa.values(), Kap->values(), Kaa->values());

  // EAS stiffness matrix is: Kdd - Kda . Kaa^-1 . Kad
  if (stiffmatrix != nullptr)
    CORE::LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMDISP_, NUMEAS_T, NUMDISP_>(
        1.0, stiffmatrix->values(), -1.0, KdaKaa.values(), Kad->values());
  // EAS stiffness matrix is: Kdp - Kda . Kaa^-1 . Kap
  if (gradmatrix != nullptr)
    CORE::LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMDISP_, NUMEAS_T, NUMPRES_>(
        1.0, gradmatrix->values(), -1.0, KdaKaa.values(), Kap->values());
  // EAS stiffness matrix is: Kpd - Kpa . Kaa^-1 . Kad
  if (dargmatrix != nullptr)
    CORE::LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMPRES_, NUMEAS_T, NUMDISP_>(
        1.0, dargmatrix->values(), -1.0, KpaKaa.values(), Kad->values());
  // EAS stiffness matrix is: Kpp - Kpa . Kaa^-1 . Kap
  if (stabmatrix != nullptr)
    CORE::LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMPRES_, NUMEAS_T, NUMPRES_>(
        1.0, stabmatrix->values(), -1.0, KpaKaa.values(), Kap->values());

  // EAS internal force is: fint - Kda^T . Kaa^-1 . feas
  if (force != nullptr)
    CORE::LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMDISP_, NUMEAS_T, 1>(
        1.0, force->values(), -1.0, KdaKaa.values(), feas->values());
  // EAS incompressibility is: fint - Kpa^T . Kaa^-1 . feas
  if (incomp != nullptr)
    CORE::LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMPRES_, NUMEAS_T, 1>(
        1.0, incomp->values(), -1.0, KpaKaa.values(), feas->values());

  // store current EAS data in history
  CORE::LINALG::DENSEFUNCTIONS::update<double, NUMEAS_T, NUMEAS_T>(*oldKaainv, *Kaa);
  CORE::LINALG::DENSEFUNCTIONS::update<double, NUMEAS_T, NUMDISP_>(*oldKad, *Kad);
  CORE::LINALG::DENSEFUNCTIONS::update<double, NUMEAS_T, NUMPRES_>(*oldKap, *Kap);
  CORE::LINALG::DENSEFUNCTIONS::update<double, NUMEAS_T, 1>(*oldfeas, *feas);

  // done
  return;
}


/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // SO3_SH8P8_EAS_H
