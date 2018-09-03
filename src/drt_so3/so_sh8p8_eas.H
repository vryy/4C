/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8_eas.H
\brief
\level 2
\maintainer Alexander Seitz
*/

/*----------------------------------------------------------------------*/
#ifndef SO_SH8P8_EAS_H
#define SO_SH8P8_EAS_H

/*----------------------------------------------------------------------*/
#include "so_sh8p8.H"
#include "Epetra_SerialDenseSolver.h"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int NUMEAS_T>
void DRT::ELEMENTS::So_sh8p8::EasUpdateIncrementally(Epetra_SerialDenseMatrix*& oldfeas,
    Epetra_SerialDenseMatrix*& oldKaainv, Epetra_SerialDenseMatrix*& oldKad,
    Epetra_SerialDenseMatrix*& oldKap, Teuchos::RCP<Epetra_SerialDenseVector>& feas,
    Teuchos::RCP<Epetra_SerialDenseMatrix>& Kaa, Teuchos::RCP<Epetra_SerialDenseMatrix>& Kad,
    Teuchos::RCP<Epetra_SerialDenseMatrix>& Kap, Epetra_SerialDenseMatrix*& alpha,
    Teuchos::RCP<Epetra_SerialDenseMatrix>& M, DRT::Container& data,
    const LINALG::Matrix<NUMDISP_, 1>& dispi, const LINALG::Matrix<NUMPRES_, 1>& presi)
{
  // retrieve history
  alpha = data.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // get old alpha
  // evaluate current (updated) EAS alphas (from history variables)
  // get stored EAS history
  oldfeas = data.GetMutable<Epetra_SerialDenseMatrix>("feas");
  oldKaainv = data.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
  oldKad = data.GetMutable<Epetra_SerialDenseMatrix>("Kda");  // actually k_ad
  oldKap = data.GetMutable<Epetra_SerialDenseMatrix>("Kap");
  if (!alpha || !oldKaainv || !oldKad || !oldKap || !oldfeas) dserror("Missing EAS history-data");

  // feas^{k+1} := feas^k + k_ad^k . Ddisp^k + k_ap^k . Dpres^k
  LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMEAS_T, NUMDISP_, 1>(
      1.0, oldfeas->A(), 1.0, oldKad->A(), dispi.A());
  LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMEAS_T, NUMPRES_, 1>(
      1.0, oldfeas->A(), 1.0, oldKap->A(), presi.A());

  // alpha^{k+1} := alpha^k + Dalpha^k
  //                alpha^k - k_aa^{k;-1} . feas^{k+1}
  LINALG::DENSEFUNCTIONS::multiply<double, NUMEAS_T, NUMEAS_T, 1>(
      1.0, *alpha, -1.0, *oldKaainv, *oldfeas);

  // EAS portion of internal forces, also called enhacement vector s or Rtilde
  feas = Teuchos::rcp(new Epetra_SerialDenseVector(NUMEAS_T));
  // EAS matrix K_{alpha alpha}, also called Dtilde
  Kaa = Teuchos::rcp(new Epetra_SerialDenseMatrix(NUMEAS_T, NUMEAS_T));
  // EAS matrix K_{alpha disp}
  Kad = Teuchos::rcp(new Epetra_SerialDenseMatrix(NUMEAS_T, NUMDISP_));
  // EAS matrix K_{alpha pres}
  Kap = Teuchos::rcp(new Epetra_SerialDenseMatrix(NUMEAS_T, NUMPRES_));

  // M-operator, ie EAS shape functions
  M = Teuchos::rcp(new Epetra_SerialDenseMatrix(MAT::NUM_STRESS_3D, NUMEAS_T));

  // see you
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int NUMEAS_T>
void DRT::ELEMENTS::So_sh8p8::EasMaterialiseShapeFcts(Teuchos::RCP<Epetra_SerialDenseMatrix> M,
    const double& detJ0, const double& detJ,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& T0invT,
    const Epetra_SerialDenseMatrix& Mloc)
{
  // map local M to global, also enhancement is refered to element origin
  // M = detJ0/detJ T0^{-T} . M
  LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D, NUMEAS_T>(
      M->A(), detJ0 / detJ, T0invT.A(), Mloc.A());

  // watch out
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int NUMEAS_T>
void DRT::ELEMENTS::So_sh8p8::EasAddStrain(LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& glstrain,
    const Teuchos::RCP<Epetra_SerialDenseMatrix>& M, const Epetra_SerialDenseMatrix* alpha)
{
  // add enhanced strains = M . alpha to GL strains to "unlock" element
  LINALG::DENSEFUNCTIONS::multiply<double, MAT::NUM_STRESS_3D, NUMEAS_T, 1>(
      1.0, glstrain.A(), 1.0, M->A(), alpha->A());
  // ready
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int NUMEAS_T>
void DRT::ELEMENTS::So_sh8p8::EasConstraintAndTangent(Teuchos::RCP<Epetra_SerialDenseVector>& feas,
    Teuchos::RCP<Epetra_SerialDenseMatrix>& Kaa, Teuchos::RCP<Epetra_SerialDenseMatrix>& Kad,
    Teuchos::RCP<Epetra_SerialDenseMatrix>& Kap, const LINALG::Matrix<NUMDIM_, NUMDIM_>& defgradD,
    const LINALG::Matrix<NUMDIM_, NUMDIM_>& invrgtstrD,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& rcgbyrgtstr,
    const double& detdefgrad, const LINALG::Matrix<NUMDFGR_, 1>& tinvdefgrad,
    const LINALG::Matrix<NUMDFGR_, NUMDFGR_>& WmT,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D>& cmat,
    const LINALG::Matrix<MAT::NUM_STRESS_3D, 1>& stress, const double& effpressure,
    const double& detJ_w, const LINALG::Matrix<MAT::NUM_STRESS_3D, NUMDISP_>& cb,
    const LINALG::Matrix<NUMDFGR_, NUMDISP_>& defgradbydisp,
    const LINALG::Matrix<NUMPRES_, 1>& prshfct, const Teuchos::RCP<Epetra_SerialDenseMatrix>& M)
{
  // derivative of assumed right stretch tensor w.r.t. EAS parameters
  LINALG::Matrix<MAT::NUM_STRESS_3D, NUMEAS_T> rgtstrbyalpha;
  LINALG::DENSEFUNCTIONS::multiplyNN<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D, NUMEAS_T>(
      rgtstrbyalpha.A(), 2.0, rcgbyrgtstr.A(), M->A());

  // derivative of pseudo identity with respect to EAS parameters
  // I^{assd}_{CB,e} = U^{d;-1}_{CD} . U^{ass}_{DB,e}
  // WARNING: I^{assd}_{CB} and I^{assd}_{CB,e} might be non-symmetric in CB
  LINALG::Matrix<NUMDFGR_, NUMEAS_T> pseudoidentity;
  for (int e = 0; e < NUMEAS_T; ++e)
  {
    for (int CB = 0; CB < NUMDFGR_; ++CB)
    {
      const int C = VOIGT9ROW_[CB];
      const int B = VOIGT9COL_[CB];
      double pseudoidentity_CBe = 0.0;
      for (int D = 0; D < NUMDIM_; ++D)
      {
        const int DB = VOIGT3X3SYM_[NUMDIM_ * D + B];
        const double DBfact = (D == B) ? 1.0 : 0.5;
        pseudoidentity_CBe += invrgtstrD(C, D) * DBfact * rgtstrbyalpha(DB, e);
      }
      pseudoidentity(CB, e) = pseudoidentity_CBe;
    }
  }

  // derivative of def.grad. with respect to e EAS parameters alpha^e
  // F_{aB,e} = F^d_{aC} . U^{d;-1}_{CD} . U^{ass}_{DB,e}
  //            = F^d_{aC} . I^{assd}_{CB,e}
  LINALG::Matrix<NUMDFGR_, NUMEAS_T> defgradbyalpha;
  for (int e = 0; e < NUMEAS_T; ++e)
  {
    for (int aB = 0; aB < NUMDFGR_; ++aB)
    {
      const int a = VOIGT9ROW_[aB];
      const int B = VOIGT9COL_[aB];
      double defgradbyalpha_aBe = 0.0;
      for (int C = 0; C < NUMDIM_; ++C)
      {
        const int CB = VOIGT3X3_[NUMDIM_ * C + B];
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
  LINALG::Matrix<MAT::NUM_STRESS_3D, NUMEAS_T> cM;
  LINALG::DENSEFUNCTIONS::multiplyNN<double, MAT::NUM_STRESS_3D, MAT::NUM_STRESS_3D, NUMEAS_T>(
      cM.A(), cmat.A(), M->A());
  // temporary ( fv . fv^T + Wm ) . M_F
  LINALG::Matrix<NUMDFGR_, NUMEAS_T> ffwmf;
  ffwmf.MultiplyNN(WmT, defgradbyalpha);
  // temporary M_F^T . fv
  LINALG::Matrix<NUMEAS_T, 1> mff;
  mff.MultiplyTN(defgradbyalpha, tinvdefgrad);
  // temporary integration factor
  const double fac = effpressure * detdefgrad * detJ_w;
#if 0
  std::cout << "effpressure=" << effpressure << std::endl;
  std::cout << "detJ_w=" << detJ_w << std::endl;
#endif


#if 0
  {
  LINALG::Matrix<NUMEAS_T,NUMEAS_T> Kaa1(true);
  LINALG::DENSEFUNCTIONS::multiplyTN<NUMEAS_T,MAT::NUM_STRESS_3D,NUMEAS_T>(1.0,Kaa1.A(), detJ_w,M->A(),cM.A());
  std::cout << "Kaa1=" << Kaa1 << std::endl;
  LINALG::Matrix<NUMEAS_T,NUMEAS_T> Kaa2(true);
  LINALG::DENSEFUNCTIONS::multiplyTN<NUMEAS_T,NUMDFGR_,NUMEAS_T>(1.0,Kaa2.A(), -fac,defgradbyalpha.A(),ffwmf.A());
  std::cout << "Kaa2=" << Kaa2 << std::endl;
  LINALG::Matrix<NUMEAS_T,NUMDISP_> Kad1(true);
  LINALG::DENSEFUNCTIONS::multiplyTN<NUMEAS_T,MAT::NUM_STRESS_3D,NUMDISP_>(1.0,Kad1.A(), detJ_w,M->A(),cb.A());
  std::cout << "Kad1=" << Kad1 << std::endl;
  std::cout << "M=" << *M << std::endl;
  std::cout << "C.B=" << cb << std::endl;
  }
#endif

  // integrate Kaa: Kaa += (M^T . cmat . M) * detJ * w(gp)
  //                     - (M_F^T . ( fv . fv^T + Wm ) . M_F) * (H . p) * detF * detJ * w
  LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMEAS_T, MAT::NUM_STRESS_3D, NUMEAS_T>(
      1.0, Kaa->A(), detJ_w, M->A(), cM.A());
  LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMEAS_T, NUMDFGR_, NUMEAS_T>(
      1.0, Kaa->A(), -fac, defgradbyalpha.A(), ffwmf.A());
  // integrate Kad: Kad += (M^T . cmat . B) * detJ * w(gp)
  //                     - (M_F^T . ( fv . fv^T + Wm ) . B_F) * (H . p) * detF * detJ * w
  LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMEAS_T, MAT::NUM_STRESS_3D, NUMDISP_>(
      1.0, Kad->A(), detJ_w, M->A(), cb.A());
  LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMEAS_T, NUMDFGR_, NUMDISP_>(
      1.0, Kad->A(), -fac, ffwmf.A(), defgradbydisp.A());
  // integrate Kap: Kap += - (M_F^T . fv . H)  * detF * detJ * w
  LINALG::DENSEFUNCTIONS::multiplyNT<double, NUMEAS_T, 1, NUMPRES_>(
      1.0, Kap->A(), -detdefgrad * detJ_w, mff.A(), prshfct.A());
  // integrate feas: feas += (M^T . sigma) * detJ *wp(gp)
  //                       - (M_F^T . fv) * (H . p)  * detF * detJ * w
  LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMEAS_T, MAT::NUM_STRESS_3D, 1>(
      1.0, feas->A(), detJ_w, M->A(), stress.A());
  LINALG::DENSEFUNCTIONS::update<double, NUMEAS_T, 1>(1.0, feas->A(), -fac, mff.A());

  // bye
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
template <int NUMEAS_T>
void DRT::ELEMENTS::So_sh8p8::EasCondensation(
    LINALG::Matrix<NUMDISP_, 1>* force,               ///< element internal force vector
    LINALG::Matrix<NUMDISP_, NUMDISP_>* stiffmatrix,  // element stiffness matrix
    LINALG::Matrix<NUMDISP_, NUMPRES_>* gradmatrix,   // element gradient matrix
    LINALG::Matrix<NUMPRES_, 1>* incomp,              ///< incompressibility residual
    LINALG::Matrix<NUMPRES_, NUMDISP_>* dargmatrix,   // 'transposed' element gradient matrix
    LINALG::Matrix<NUMPRES_, NUMPRES_>* stabmatrix,   // element stabilisation matrix
    Epetra_SerialDenseMatrix*& oldfeas, Epetra_SerialDenseMatrix*& oldKaainv,
    Epetra_SerialDenseMatrix*& oldKad, Epetra_SerialDenseMatrix*& oldKap,
    const Teuchos::RCP<Epetra_SerialDenseVector>& feas,
    const Teuchos::RCP<Epetra_SerialDenseMatrix>& Kaa,
    const Teuchos::RCP<Epetra_SerialDenseMatrix>& Kad,
    const Teuchos::RCP<Epetra_SerialDenseMatrix>& Kap)
{
#if 0
  std::cout << " -  - - - - -  " << std::endl;
  std::cout << "Kaa=" << *Kaa << std::endl;
#endif

  // we need the inverse of Kaa
  Epetra_SerialDenseSolver solve_for_inverseKaa;
  solve_for_inverseKaa.SetMatrix(*Kaa);
  solve_for_inverseKaa.Invert();

#if 0
  std::cout << "Kaa=" << *Kaa << std::endl;
  std::cout << "Kap=" << *Kap << std::endl;
  std::cout << "Kad=" << std::scientific << *Kad << std::endl;
#endif

  // temporary Kda.Kaa^{-1}
  LINALG::Matrix<NUMDISP_, NUMEAS_T> KdaKaa(false);
  LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMDISP_, NUMEAS_T, NUMEAS_T>(
      KdaKaa.A(), Kad->A(), Kaa->A());
  // temporary Kpa.Kaa^{-1}
  LINALG::Matrix<NUMPRES_, NUMEAS_T> KpaKaa(false);
  LINALG::DENSEFUNCTIONS::multiplyTN<double, NUMPRES_, NUMEAS_T, NUMEAS_T>(
      KpaKaa.A(), Kap->A(), Kaa->A());

#if 0
  LINALG::Matrix<NUMDISP_,NUMDISP_> KdaKaaKad(true);
  LINALG::DENSEFUNCTIONS::multiplyNN<double,NUMDISP_,NUMEAS_T,NUMDISP_>(1.0,KdaKaaKad.A(), -1.0,KdaKaa.A(),Kad->A());
  std::cout << "KdaKaaKad=" << KdaKaaKad << std::endl;
  LINALG::Matrix<NUMDISP_,NUMPRES_> KdaKaaKap(true);
  LINALG::DENSEFUNCTIONS::multiplyNN<double,NUMDISP_,NUMEAS_T,NUMPRES_>(1.0,KdaKaaKap.A(), -1.0,KdaKaa.A(),Kap->A());
  std::cout << "KdaKaaKap=" << KdaKaaKap << std::endl;
  LINALG::Matrix<NUMPRES_,NUMDISP_> KpaKaaKad(true);
  LINALG::DENSEFUNCTIONS::multiplyNN<double,NUMPRES_,NUMEAS_T,NUMDISP_>(1.0,KpaKaaKad.A(), -1.0,KpaKaa.A(),Kad->A());
  std::cout << "KpaKaaKad=" << KpaKaaKad << std::endl;
  LINALG::Matrix<NUMPRES_,NUMPRES_> KpaKaaKap(true);
  LINALG::DENSEFUNCTIONS::multiplyNN<double,NUMPRES_,NUMEAS_T,NUMPRES_>(1.0,KpaKaaKap.A(), -1.0,KpaKaa.A(),Kap->A());
  std::cout << "KpaKaaKap=" << KpaKaaKap << std::endl;
  std::cout << "stiffmatrix=" << *stiffmatrix << std::endl;
  std::cout << "gradmatrix=" << *gradmatrix << std::endl;
  std::cout << "stabmatrix=" << *stabmatrix << std::endl;
#endif

  // EAS stiffness matrix is: Kdd - Kda . Kaa^-1 . Kad
  if (stiffmatrix != NULL)
    LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMDISP_, NUMEAS_T, NUMDISP_>(
        1.0, stiffmatrix->A(), -1.0, KdaKaa.A(), Kad->A());
  // EAS stiffness matrix is: Kdp - Kda . Kaa^-1 . Kap
  if (gradmatrix != NULL)
    LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMDISP_, NUMEAS_T, NUMPRES_>(
        1.0, gradmatrix->A(), -1.0, KdaKaa.A(), Kap->A());
  // EAS stiffness matrix is: Kpd - Kpa . Kaa^-1 . Kad
  if (dargmatrix != NULL)
    LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMPRES_, NUMEAS_T, NUMDISP_>(
        1.0, dargmatrix->A(), -1.0, KpaKaa.A(), Kad->A());
  // EAS stiffness matrix is: Kpp - Kpa . Kaa^-1 . Kap
  if (stabmatrix != NULL)
    LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMPRES_, NUMEAS_T, NUMPRES_>(
        1.0, stabmatrix->A(), -1.0, KpaKaa.A(), Kap->A());

  // EAS internal force is: fint - Kda^T . Kaa^-1 . feas
  if (force != NULL)
    LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMDISP_, NUMEAS_T, 1>(
        1.0, force->A(), -1.0, KdaKaa.A(), feas->A());
  // EAS incompressibility is: fint - Kpa^T . Kaa^-1 . feas
  if (incomp != NULL)
    LINALG::DENSEFUNCTIONS::multiplyNN<double, NUMPRES_, NUMEAS_T, 1>(
        1.0, incomp->A(), -1.0, KpaKaa.A(), feas->A());

  // store current EAS data in history
  LINALG::DENSEFUNCTIONS::update<double, NUMEAS_T, NUMEAS_T>(*oldKaainv, *Kaa);
  LINALG::DENSEFUNCTIONS::update<double, NUMEAS_T, NUMDISP_>(*oldKad, *Kad);
  LINALG::DENSEFUNCTIONS::update<double, NUMEAS_T, NUMPRES_>(*oldKap, *Kap);
  LINALG::DENSEFUNCTIONS::update<double, NUMEAS_T, 1>(*oldfeas, *feas);

  // done
  return;
}


/*----------------------------------------------------------------------*/
#endif  // #ifndef SO_SH8P8_EAS_H
