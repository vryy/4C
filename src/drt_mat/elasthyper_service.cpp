/*----------------------------------------------------------------------*/
/*!

\brief Contains the declaration of service functions for hyperelastic materials

\level 1

\maintainer Fabian Braeu

*/
/*----------------------------------------------------------------------*/

#include "elasthyper_service.H"
#include "material_service.H"
#include "../headers/definitions.h"

void MAT::ElastHyperEvaluate(const LINALG::Matrix<3, 3>* defgrd,
    const LINALG::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    LINALG::Matrix<6, 1>* stress, LINALG::Matrix<6, 6>* cmat, int eleGID,
    const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum,
    MAT::SummandProperties& properties, bool checkpolyconvexity)
{
  static LINALG::Matrix<6, 1> id2(false);
  id2.Clear();
  static LINALG::Matrix<6, 1> C_strain(false);
  C_strain.Clear();
  static LINALG::Matrix<6, 1> iC_strain(false);
  iC_strain.Clear();
  static LINALG::Matrix<3, 1> prinv(false);
  prinv.Clear();
  static LINALG::Matrix<3, 1> dPI(false);
  dPI.Clear();
  static LINALG::Matrix<6, 1> ddPII(false);
  ddPII.Clear();

  // Evaluate identity tensor
  IdentityMatrix(id2);

  // Evalutate Right Cauchy-Green strain tensor in strain-like Voigt notation
  EvaluateRightCauchyGreenStrainLikeVoigt(*glstrain, C_strain);

  // Invert Right Cauchy Green Strain tensor
  VStrainUtils::InverseTensor(C_strain, iC_strain);

  // Evaluate principle invariants
  InvariantsPrincipal<MAT::VoigtNotation::strain>(prinv, C_strain);

  // Evaluate derivatives of potsum w.r.t the principal invariants
  ElastHyperEvaluateInvariantDerivatives(prinv, dPI, ddPII, potsum, properties, eleGID);

  // check if system is polyconvex (set "POLYCONVEX 1" in material input-line)
  if (checkpolyconvexity)
    ElastHyperCheckPolyconvexity(*defgrd, prinv, dPI, ddPII, params, eleGID, properties);


  // clear stress and cmat (for safety reasons)
  stress->Clear();
  cmat->Clear();

  // Evaluate isotropic stress response
  ElastHyperAddIsotropicStressCmat(*stress, *cmat, C_strain, iC_strain, prinv, dPI, ddPII);

  if (properties.coeffStretchesPrinc || properties.coeffStretchesMod)
  {
    ElastHyperAddResponseStretches(*cmat, *stress, C_strain, potsum, properties, eleGID);
  }

  // Evaluate anisotropic stress response from summands with principle invariants formulation
  if (properties.anisoprinc)
    ElastHyperAddAnisotropicPrinc(*stress, *cmat, C_strain, params, eleGID, potsum);

  // Evaluate anisotropic stress response from summands with modified invariants formulation
  if (properties.anisomod)
    ElastHyperAddAnisotropicMod(*stress, *cmat, C_strain, iC_strain, prinv, eleGID, potsum);
}

void MAT::EvaluateRightCauchyGreenStrainLikeVoigt(
    const LINALG::Matrix<6, 1>& E_strain, LINALG::Matrix<6, 1>& C_strain)
{
  // C = 2*E+I
  C_strain.Update(2.0, E_strain, 0.0);

  // Add Identity
  for (unsigned i = 0; i < 3; ++i) C_strain(i) += 1.0;
}

void MAT::ElastHyperEvaluateInvariantDerivatives(const LINALG::Matrix<3, 1>& prinv,
    LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII,
    const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum,
    const MAT::SummandProperties& properties, int eleGID)
{
  // derivatives of principla materials
  if (properties.isoprinc)
  {
    // loop map of associated potential summands
    for (auto& p : potsum)
    {
      p->AddDerivativesPrincipal(dPI, ddPII, prinv, eleGID);
    }
  }

  // derivatives of decoupled (volumetric or isochoric) materials
  if (properties.isomod)
  {
    static LINALG::Matrix<3, 1> modinv(true);
    modinv.Clear();
    static LINALG::Matrix<3, 1> dPmodI(true);
    dPmodI.Clear();
    static LINALG::Matrix<6, 1> ddPmodII(true);
    ddPmodII.Clear();

    // Evaluate modified invariants
    MAT::InvariantsModified(modinv, prinv);

    for (auto& p : potsum)
    {
      p->AddDerivativesModified(dPmodI, ddPmodII, modinv, eleGID);
    }

    // convert decoupled derivatives to principal derivatives
    MAT::ConvertModToPrinc(prinv, dPmodI, ddPmodII, dPI, ddPII);
  }
}

void MAT::ConvertModToPrinc(const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPmodI,
    const LINALG::Matrix<6, 1>& ddPmodII, LINALG::Matrix<3, 1>& dPI, LINALG::Matrix<6, 1>& ddPII)
{
  // Conversions to dPI
  dPI(0) += std::pow(prinv(2), -1. / 3.) * dPmodI(0);
  dPI(1) += std::pow(prinv(2), -2. / 3.) * dPmodI(1);
  dPI(2) += 0.5 * std::pow(prinv(2), -0.5) * dPmodI(2) -
            1. / 3. * prinv(0) * std::pow(prinv(2), -4. / 3.) * dPmodI(0) -
            2. / 3. * prinv(1) * std::pow(prinv(2), -5. / 3.) * dPmodI(1);

  // Conversions to ddPII
  ddPII(0) += std::pow(prinv(2), -2. / 3.) * ddPmodII(0);
  ddPII(1) += std::pow(prinv(2), -4. / 3.) * ddPmodII(1);
  ddPII(2) += (1. / 9.) * std::pow(prinv(2), -8. / 3.) * prinv(0) * prinv(0) * ddPmodII(0) +
              (4. / 9.) * prinv(0) * prinv(1) * std::pow(prinv(2), -3.) * ddPmodII(5) -
              (1. / 3.) * std::pow(prinv(2), -11. / 6.) * prinv(0) * ddPmodII(4) +
              (4. / 9.) * std::pow(prinv(2), -7. / 3.) * prinv(0) * dPmodI(0) +
              (4. / 9.) * std::pow(prinv(2), -10. / 3.) * prinv(1) * prinv(1) * ddPmodII(1) -
              (2. / 3.) * std::pow(prinv(2), -13. / 6.) * prinv(1) * ddPmodII(3) +
              (10. / 9.) * std::pow(prinv(2), -8. / 3.) * prinv(1) * dPmodI(1) +
              0.25 * std::pow(prinv(2), -1.) * ddPmodII(2) -
              0.25 * std::pow(prinv(2), -1.5) * dPmodI(2);
  ddPII(3) += -(1. / 3.) * std::pow(prinv(2), -2.) * prinv(0) * ddPmodII(5) -
              (2. / 3.) * std::pow(prinv(2), -7. / 3.) * prinv(1) * ddPmodII(1) +
              0.5 * std::pow(prinv(2), -7. / 6.) * ddPmodII(3) -
              (2. / 3.) * std::pow(prinv(2), -5. / 3.) * dPmodI(1);
  ddPII(4) += -(1. / 3.) * std::pow(prinv(2), -5. / 3.) * prinv(0) * ddPmodII(0) -
              (2. / 3.) * std::pow(prinv(2), -2.) * prinv(1) * ddPmodII(5) +
              0.5 * std::pow(prinv(2), -5. / 6.) * ddPmodII(4) -
              (1. / 3.) * std::pow(prinv(2), -4. / 3.) * dPmodI(0);
  ddPII(5) += std::pow(prinv(2), -1.) * ddPmodII(5);
}

void MAT::ElastHyperAddIsotropicStressCmat(LINALG::Matrix<6, 1>& S_stress,
    LINALG::Matrix<6, 6>& cmat, const LINALG::Matrix<6, 1>& C_strain,
    const LINALG::Matrix<6, 1>& iC_strain, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& dPI, const LINALG::Matrix<6, 1>& ddPII)
{
  // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  static LINALG::Matrix<3, 1> gamma(true);
  // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  static LINALG::Matrix<8, 1> delta(true);
  // 2nd order identity tensor
  static LINALG::Matrix<6, 1> id2(false);
  // Right Cauchy-Green tensor in stress-like Voigt notation
  static LINALG::Matrix<6, 1> C_stress(false);
  // Inverse Right Cauchy-Green tensor in stress-like Voigt notation
  static LINALG::Matrix<6, 1> iC_stress(false);
  // 4th order identity tensor (rows and colums are stress-like)
  static LINALG::Matrix<6, 6> id4sharp(false);
  MAT::FourthOrderIdentityMatrix<MAT::VoigtNotation::stress, MAT::VoigtNotation::stress>(id4sharp);

  // initialize matrices
  MAT::IdentityMatrix(id2);
  MAT::VStrainUtils::ToStressLike(C_strain, C_stress);
  MAT::VStrainUtils::ToStressLike(iC_strain, iC_stress);

  // compose coefficients
  CalculateGammaDelta(gamma, delta, prinv, dPI, ddPII);

  // 2nd Piola Kirchhoff stress
  S_stress.Update(gamma(0), id2, 1.0);
  S_stress.Update(gamma(1), C_stress, 1.0);
  S_stress.Update(gamma(2), iC_stress, 1.0);

  // constitutive tensor
  // contribution: Id \otimes Id
  cmat.MultiplyNT(delta(0), id2, id2, 1.0);
  // contribution: Id \otimes C + C \otimes Id
  cmat.MultiplyNT(delta(1), id2, C_stress, 1.0);
  cmat.MultiplyNT(delta(1), C_stress, id2, 1.0);
  // contribution: Id \otimes Cinv + Cinv \otimes Id
  cmat.MultiplyNT(delta(2), id2, iC_stress, 1.0);
  cmat.MultiplyNT(delta(2), iC_stress, id2, 1.0);
  // contribution: C \otimes C
  cmat.MultiplyNT(delta(3), C_stress, C_stress, 1.0);
  // contribution: C \otimes Cinv + Cinv \otimes C
  cmat.MultiplyNT(delta(4), C_stress, iC_stress, 1.0);
  cmat.MultiplyNT(delta(4), iC_stress, C_stress, 1.0);
  // contribution: Cinv \otimes Cinv
  cmat.MultiplyNT(delta(5), iC_stress, iC_stress, 1.0);
  // contribution: Cinv \odot Cinv
  AddtoCmatHolzapfelProduct(cmat, iC_stress, delta(6));
  // contribution: Id4^#
  cmat.Update(delta(7), id4sharp, 1.0);
}

void MAT::ElastHyperAddResponseStretches(LINALG::Matrix<6, 6>& cmat, LINALG::Matrix<6, 1>& S_stress,
    const LINALG::Matrix<6, 1>& C_strain,
    const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum,
    const SummandProperties& properties, int eleGID)
{
  // get principal stretches and directions
  LINALG::Matrix<3, 1> prstr;
  LINALG::Matrix<3, 3> prdir;
  StretchesPrincipal(prstr, prdir, C_strain);
  // modified stretches
  LINALG::Matrix<3, 1> modstr;
  StretchesModified(modstr, prstr);
  // determinant of deformation gradient
  const double detdefgrad = prstr(0) * prstr(1) * prstr(2);

  // get coefficients
  LINALG::Matrix<3, 1> gamma_(true);
  LINALG::Matrix<6, 1> delta_(true);
  if (properties.coeffStretchesPrinc)
  {
    // loop map of associated potential summands
    for (const auto& p : potsum)
    {
      p->AddCoefficientsStretchesPrincipal(gamma_, delta_, prstr);
    }
  }
  if (properties.coeffStretchesMod)
  {
    // reciprocal of cubic root of determinant of deformation gradient (convenience)
    const double detdefgrad13 = std::pow(detdefgrad, -1.0 / 3.0);
    // retrieve coefficients with respect to modified principal stretches
    static LINALG::Matrix<3, 1> modgamma(true);
    modgamma.Clear();
    static LINALG::Matrix<6, 1> moddelta(true);
    moddelta.Clear();
    {
      // loop map of associated potential summands
      for (const auto& p : potsum)
      {
        p->AddCoefficientsStretchesModified(modgamma, moddelta, modstr);
      }
    }
    // convert modified coefficients to oridinary counterparts
    //
    // derivatives of modified pr. stretches WRT pr. stretches
    static LINALG::Matrix<3, 3> modbypr(false);
    for (int al = 0; al < 3; ++al)
    {
      for (int be = 0; be < 3; ++be)
      {
        modbypr(al, be) = -modstr(al) / modstr(be);
      }
      modbypr(al, al) += 3.0;
    }
    modbypr.Scale(detdefgrad13 / 3.0);
    // determine unmodified coefficients gamma and add them
    gamma_.MultiplyTN(1.0, modbypr, modgamma, 1.0);
    // determine unmodified coefficients delta and add them
    //
    // rewrite mod.coeff. as 2-tensor
    static LINALG::Matrix<3, 3> moddeltat(false);
    moddeltat(0, 0) = moddelta(0);
    moddeltat(1, 1) = moddelta(1);
    moddeltat(2, 2) = moddelta(2);
    moddeltat(0, 1) = moddeltat(1, 0) = moddelta(3);
    moddeltat(1, 2) = moddeltat(2, 1) = moddelta(4);
    moddeltat(2, 0) = moddeltat(0, 2) = moddelta(5);
    // Psi_{,barlam barlam} barlam_{,lam} barlam_{,lam}
    static LINALG::Matrix<3, 3> aux(false);
    aux.MultiplyTN(modbypr, moddeltat);
    static LINALG::Matrix<3, 3> deltat(false);
    deltat.MultiplyNN(aux, modbypr);
    // Psi_{,barlam} barlam_{,lam lam}
    for (int be = 0; be < 3; ++be)
    {
      for (int ga = 0; ga < 3; ++ga)
      {
        double deltat_bega = 0.0;
        for (int al = 0; al < 3; ++al)
        {
          deltat_bega += -modgamma(al) * modbypr(al, be) / (3.0 * prstr(ga));
          if (ga == al) deltat_bega += -modgamma(al) * detdefgrad13 / (3.0 * prstr(be));
          if (be == ga)
            deltat_bega += modgamma(al) * detdefgrad13 * prstr(al) / (3.0 * prstr(be) * prstr(be));
        }
        deltat(be, ga) += deltat_bega;
      }
    }
    // add to delta
    // Psi_{lam lam} = Psi_{,barlam barlam} barlam_{,lam} barlam_{,lam}
    //               + Psi_{,barlam} barlam_{,lam lam}
    delta_(0) += deltat(0, 0);
    delta_(1) += deltat(1, 1);
    delta_(2) += deltat(2, 2);
    delta_(3) += deltat(0, 1);
    delta_(4) += deltat(1, 2);
    delta_(5) += deltat(2, 0);
  }

  // principal 2nd Piola--Kirchhoff stress tensor, cf [1] Eq (6.47)
  static LINALG::Matrix<3, 1> prsts(true);
  prsts.Clear();
  for (int al = 0; al < 3; ++al)
  {
    // PK2 principal stresses
    prsts(al) = gamma_(al) / prstr(al);
    // PK2 tensor in Voigt notation
    S_stress(0) += prsts(al) * prdir(0, al) * prdir(0, al);  // S^11
    S_stress(1) += prsts(al) * prdir(1, al) * prdir(1, al);  // S^22
    S_stress(2) += prsts(al) * prdir(2, al) * prdir(2, al);  // S^33
    S_stress(3) += prsts(al) * prdir(0, al) * prdir(1, al);  // S^12
    S_stress(4) += prsts(al) * prdir(1, al) * prdir(2, al);  // S^23
    S_stress(5) += prsts(al) * prdir(2, al) * prdir(0, al);  // S^31
  }

  // integration factor prfact_{al be}
  static LINALG::Matrix<6, 1> prfact1(true);
  prfact1.Clear();
  static LINALG::Matrix<6, 1> prfact2(true);
  prfact2.Clear();
  for (int albe = 0; albe < 6; ++albe)
  {
    const int al = VOIGT6ROW_[albe];
    const int be = VOIGT6COL_[albe];
    double prfact1_albe = delta_(albe) / (prstr(al) * prstr(be));
    if (albe < 3) prfact1_albe -= gamma_(al) / (prstr(be) * prstr(al) * prstr(al));
    prfact1(albe) = prfact1_albe;
    if (al != be)
    {
      if (fabs(prstr(al) - prstr(be)) < EPS6)
        prfact2(albe) = (prfact1(be) - prfact1(albe)) / 2.0;
      else
        prfact2(albe) = (prsts(be) - prsts(al)) / (prstr(be) * prstr(be) - prstr(al) * prstr(al));
    }
  }

  // add elasticity 4-tensor, cf Holzapfel [1] Eq (6.180),(6.196)
  for (int kl = 0; kl < 6; ++kl)
  {
    const int k = VOIGT6ROW_[kl];
    const int l = VOIGT6COL_[kl];
    for (int ij = 0; ij < 6; ++ij)
    {
      const int i = VOIGT6ROW_[ij];
      const int j = VOIGT6COL_[ij];
      double c_ijkl = 0.0;
      for (int albe = 0; albe < 6; ++albe)
      {
        const int al = VOIGT6ROW_[albe];
        const int be = VOIGT6COL_[albe];
        const double fact1 = prfact1(albe);
        c_ijkl += fact1 * prdir(i, al) * prdir(j, al) * prdir(k, be) * prdir(l, be);
        if (albe >= 3)
        {  // al!=be
          c_ijkl += fact1 * prdir(i, be) * prdir(j, be) * prdir(k, al) * prdir(l, al);
          const double fact2 = prfact2(albe);
          c_ijkl += fact2 * prdir(i, al) * prdir(j, be) * prdir(k, al) * prdir(l, be) +
                    fact2 * prdir(i, al) * prdir(j, be) * prdir(k, be) * prdir(l, al) +
                    fact2 * prdir(i, be) * prdir(j, al) * prdir(k, be) * prdir(l, al) +
                    fact2 * prdir(i, be) * prdir(j, al) * prdir(k, al) * prdir(l, be);
        }
      }
      cmat(ij, kl) += c_ijkl;
    }
  }
}

void MAT::ElastHyperAddAnisotropicPrinc(LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
    const LINALG::Matrix<6, 1>& C_strain, Teuchos::ParameterList& params, int eleGID,
    const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum)
{
  // Loop over all summands and add aniso stress
  // ToDo: This should be solved in analogy to the solution in elast_remodelfiber.cpp
  // ToDo: i.e. by evaluating the derivatives of the potsum w.r.t. the anisotropic invariants
  for (auto& p : potsum) p->AddStressAnisoPrincipal(C_strain, cmat, S_stress, params, eleGID);
}

void MAT::ElastHyperAddAnisotropicMod(LINALG::Matrix<6, 1>& S_stress, LINALG::Matrix<6, 6>& cmat,
    const LINALG::Matrix<6, 1>& C_strain, const LINALG::Matrix<6, 1>& iC_strain,
    const LINALG::Matrix<3, 1>& prinv, int eleGID,
    const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum)
{
  static LINALG::Matrix<6, 1> iC_stress(false);
  VStrainUtils::ToStressLike(iC_strain, iC_stress);
  // Loop over all summands and add aniso stress
  // ToDo: This should be solved in analogy to the solution in elast_remodelfiber.cpp
  // ToDo: i.e. by evaluating the derivatives of the potsum w.r.t. the anisotropic invariants
  for (auto& p : potsum)
    p->AddStressAnisoModified(C_strain, iC_stress, cmat, S_stress, prinv(2), eleGID);
}

void MAT::CalculateGammaDelta(LINALG::Matrix<3, 1>& gamma, LINALG::Matrix<8, 1>& delta,
    const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPI,
    const LINALG::Matrix<6, 1>& ddPII)
{
  // according to Holzapfel-Nonlinear Solid Mechanics p. 216
  gamma(0) = 2. * (dPI(0) + prinv(0) * dPI(1));
  gamma(1) = -2. * dPI(1);
  gamma(2) = 2. * prinv(2) * dPI(2);

  // according to Holzapfel-Nonlinear Solid Mechanics p. 261
  delta(0) = 4. * (ddPII(0) + 2. * prinv(0) * ddPII(5) + dPI(1) + prinv(0) * prinv(0) * ddPII(1));
  delta(1) = -4. * (ddPII(5) + prinv(0) * ddPII(1));
  delta(2) = 4. * (prinv(2) * ddPII(4) + prinv(0) * prinv(2) * ddPII(3));
  delta(3) = 4. * ddPII(1);
  delta(4) = -4. * prinv(2) * ddPII(3);
  delta(5) = 4. * (prinv(2) * dPI(2) + prinv(2) * prinv(2) * ddPII(2));
  delta(6) = -4. * prinv(2) * dPI(2);
  delta(7) = -4. * dPI(1);
}

void MAT::ElastHyperProperties(
    const std::vector<Teuchos::RCP<MAT::ELASTIC::Summand>>& potsum, SummandProperties& properties)
{
  for (auto& p : potsum)
  {
    p->SpecifyFormulation(properties.isoprinc, properties.isomod, properties.anisoprinc,
        properties.anisomod, properties.viscoGeneral);

    properties.coeffStretchesPrinc |= p->HaveCoefficientsStretchesPrincipal();
    properties.coeffStretchesMod |= p->HaveCoefficientsStretchesModified();
  }
}

void MAT::ElastHyperCheckPolyconvexity(const LINALG::Matrix<3, 3>& defgrd,
    const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& dPI,
    const LINALG::Matrix<6, 1>& ddPII, Teuchos::ParameterList& params, const int eleGID,
    const SummandProperties& properties)
{
  // This polyconvexity-test is just implemented for isotropic hyperelastic-materials
  // --> error if anisotropic material is tested (plastic and viscoelastic materials should not get
  // in here)
  if (properties.anisoprinc || properties.anisomod)
    dserror(
        "This polyconvexity-check is just implemented for isotropic "
        "hyperelastic-materials (do not use for anistropic materials).");

  // principal invariants (i)
  // first strain energy derivative dPI (i)
  // second strain energy derivative ddPII (i)

  // J = sqrt(I_3) = modinv(2)
  double J = std::pow(prinv(2), 1. / 2.);

  // defgrd = F (i)
  // dfgrd = F in Voigt - Notation
  static LINALG::Matrix<9, 1> dfgrd(true);
  Matrix3x3to9x1(defgrd, dfgrd);

  // Cof(F) = J*F^(-T)
  static LINALG::Matrix<3, 3> CoFacF(true);  // Cof(F) in Matrix-Notation
  static LINALG::Matrix<9, 1> CofF(true);    // Cof(F) in Voigt-Notation
  CoFacF.Invert(defgrd);
  CoFacF.Scale(J);
  // sort in Voigt-Notation and invert!
  Matrix3x3to9x1(CoFacF, CofF);

  // id4 (9x9)
  static LINALG::Matrix<9, 9> ID4(true);
  for (int i = 0; i < 9; i++)
    for (int j = 0; j < 9; j++)
      if (i == j) ID4(i, j) = 1.0;

  // Frechet Derivative according to Ebbing, PhD-thesis page 79, Eq: (5.31)
  static LINALG::Matrix<19, 19> FreD(true);
  FreD.Clear();

  // single matrices of Frechet Derivative:

  // d²P/dFdF
  // = 4 d^2\Psi/dI_1dI_1 F \otimes F + 2 \d\Psi/dI_1 *II
  static LINALG::Matrix<9, 9> FreDFF(true);
  FreDFF.Clear();
  FreDFF.MultiplyNT(4 * ddPII(0), dfgrd, dfgrd, 1.0);
  FreDFF.Update(2 * dPI(0), ID4, 1.0);

  // d²P/d(cofF)d(cofF)
  // = = 4 d^2\Psi/dI_2dI_2 cof(F) \otimes cof(F) + 2 \d\Psi/dI_2 *II
  static LINALG::Matrix<9, 9> FreDcFcF(true);
  FreDcFcF.Clear();
  FreDcFcF.MultiplyNT(4 * ddPII(1), CofF, CofF, 1.0);
  FreDcFcF.Update(2 * dPI(1), ID4, 1.0);

  // d²P/d(detF)d(detF)
  // = 2*d \Psi/dI_3 + 4*I_3*d²\Psi/dI_3dI_3
  double FreDJJ(true);
  FreDJJ += 2 * dPI(2) + 4 * prinv(2) * ddPII(2);

  // d²P/d(cofF)dF
  // = 4*d\Psi/dI_1dI_2 F /otimes CofF
  static LINALG::Matrix<9, 9> FreDcFF(true);
  FreDcFF.Clear();
  FreDcFF.MultiplyNT(4 * ddPII(5), dfgrd, CofF, 1.0);

  // d²P/d(detF)d(cofF)
  // = 4*J*d^2 \Psi /dI_2 dI_3 \mat{CofF}
  static LINALG::Matrix<9, 1> FreDcFJ(true);
  FreDcFF.Clear();
  FreDcFJ.Update(4 * J * ddPII(3), CofF, 1.0);

  // d²P/d(detF) dF = d²P/dF d(detF)
  // = 4*J*d^2 \Psi /dI_1 dI_3 \mat{F}
  static LINALG::Matrix<9, 1> FreDFJ(true);
  FreDcFF.Clear();
  FreDFJ.Update(4 * J * ddPII(4), dfgrd, 1.0);

  // Sort values in Frechet Derivative

  // FreD = [FreDFF   FreDcFF    FreDFJ
  //         FreDcFF  FreDcFcF   FreDcFJ
  //         FreDFJ   FreDcFJ    FreDJJ]
  for (int i = 0; i < 9; i++)
    for (int j = 0; j < 9; j++)
    {
      FreD(i, j) = FreDFF(i, j);
      FreD(i, j + 9) = FreDcFF(i, j);
      FreD(i + 9, j) = FreDcFF(i, j);
      FreD(i + 9, j + 9) = FreDcFcF(i, j);
    }

  for (int i = 0; i < 9; i++)
  {
    FreD(i + 9, 18) = FreDcFJ(i);
    FreD(18, i + 9) = FreDcFJ(i);
    FreD(i, 18) = FreDFJ(i);
    FreD(18, i) = FreDFJ(i);
  }

  FreD(18, 18) = FreDJJ;

  // EigenValues of Frechet Derivative
  static LINALG::Matrix<19, 19> EWFreD(true);  // EW on diagonal
  static LINALG::Matrix<19, 19> EVFreD(true);
  LINALG::SYEV(FreD, EWFreD, EVFreD);

  // Just positive EigenValues --> System is polyconvex
  for (int i = 0; i < 19; i++)
    for (int j = 0; j < 19; j++)
      if (i == j)  // values on diagonal = EigenValues
        if (EWFreD(i, i) <
            (-1.0e-10 * EWFreD.NormInf()))  // do not test < 0, but reasonable small value
        {
          std::cout << "\nWARNING: Your system is not polyconvex!" << std::endl;
          std::cout << "Polyconvexity fails at: Element-Id: " << eleGID
                    << " and Gauß-Point: " << params.get<int>("gp") << std::endl;
          std::cout << "Eigenvalues of the Frechet Derivative are: " << EWFreD << std::endl;
        }
}