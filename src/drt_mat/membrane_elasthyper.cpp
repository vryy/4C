/*!----------------------------------------------------------------------
\brief hyperelastic toolbox for membranes assuming incompressibility and plane stress

The input line should read
MAT 0 MAT_Membrane_ElastHyper NUMMAT 2 MATIDS 1 2 DENS 0

\level 3

\maintainer  Sebastian Fuchs

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
#include "membrane_elasthyper.H"

#include "../drt_matelast/elast_summand.H"

/*----------------------------------------------------------------------*
 | constructor                                           sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
MAT::PAR::Membrane_ElastHyper::Membrane_ElastHyper(Teuchos::RCP<MAT::PAR::Material> matdata)
    : MAT::PAR::ElastHyper(matdata)
{
  return;
}  // MAT::PAR::Membrane_ElastHyper::Membrane_ElastHyper

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Membrane_ElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Membrane_ElastHyper(this));
}  // MAT::PAR::Membrane_ElastHyper::CreateMaterial


MAT::Membrane_ElastHyperType MAT::Membrane_ElastHyperType::instance_;


DRT::ParObject* MAT::Membrane_ElastHyperType::Create(const std::vector<char>& data)
{
  MAT::Membrane_ElastHyper* memelhy = new MAT::Membrane_ElastHyper();
  memelhy->Unpack(data);

  return memelhy;
}  // MAT::Membrane_ElastHyperType::Create

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
MAT::Membrane_ElastHyper::Membrane_ElastHyper() : MAT::ElastHyper(), fibervecs_(true)
{
  return;
}  // MAT::Membrane_ElastHyper::Membrane_ElastHyper()

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
MAT::Membrane_ElastHyper::Membrane_ElastHyper(MAT::PAR::Membrane_ElastHyper* params)
    : MAT::ElastHyper(params), fibervecs_(true)
{
  return;
}  // MAT::Membrane_ElastHyper::Membrane_ElastHyper()

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  MAT::ElastHyper::Pack(data);

  AddtoPack(data, fibervecs_);

  return;
}  // MAT::Membrane_ElastHyper::Pack()

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  MAT::ElastHyper::Unpack(basedata);

  ExtractfromPack(position, data, fibervecs_);

  return;
}  // MAT::Membrane_ElastHyper::Unpack()

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  // call setup of base class
  MAT::ElastHyper::Setup(numgp, linedef);

  GetFiberVecs(fibervecs_);

  return;
}  // MAT::Membrane_ElastHyper::Setup()

/*----------------------------------------------------------------------*
 | hyperelastic stress response plus elasticity tensor   sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::Evaluate(
    LINALG::Matrix<3, 3>& cauchygreen,  ///< right Cauchy-Green tensor
    Teuchos::ParameterList& params,     ///< Container for additional information
    LINALG::Matrix<3, 3>&
        Q_trafo,  ///< Trafo from local membrane orthonormal coordinates to global coordinates
    LINALG::Matrix<3, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses in stress-like voigt notation
    LINALG::Matrix<3, 3>* cmat,    ///< Constitutive matrix
    const int eleGID               ///< Element GID
)
{
  // blank resulting quantities
  stress->Clear();
  cmat->Clear();

  // kinematic quantities and identity tensors
  LINALG::Matrix<3, 1> id2(true);
  LINALG::Matrix<3, 3> id4sharp(true);
  LINALG::Matrix<3, 1> rcg(true);
  double rcg33;
  LINALG::Matrix<3, 1> icg(true);
  EvaluateKinQuant(cauchygreen, id2, id4sharp, rcg, rcg33, icg);

  // evaluate isotropic 2nd Piola-Kirchhoff stress and constitutive tensor
  LINALG::Matrix<3, 1> stress_iso(true);
  LINALG::Matrix<3, 3> cmat_iso(true);
  EvaluateIsotropicStressCmat(stress_iso, cmat_iso, id2, id4sharp, rcg, rcg33, icg, eleGID);

  // update 2nd Piola-Kirchhoff stress and constitutive tensor
  stress->Update(1.0, stress_iso, 1.0);
  cmat->Update(1.0, cmat_iso, 1.0);

  // evaluate anisotropic 2nd Piola-Kirchhoff stress and constitutive tensor
  if (summandProperties_.anisoprinc)
  {
    LINALG::Matrix<3, 1> stress_aniso(true);
    LINALG::Matrix<3, 3> cmat_aniso(true);
    EvaluateAnisotropicStressCmat(stress_aniso, cmat_aniso, Q_trafo, rcg, rcg33, params, eleGID);

    // update 2nd Piola-Kirchhoff stress and constitutive tensor
    stress->Update(1.0, stress_aniso, 1.0);
    cmat->Update(1.0, cmat_aniso, 1.0);
  }
  if (summandProperties_.anisomod)
  {
    dserror("anisomod_ not implemented for membrane elasthyper materials!");
  }

  return;
}  // MAT::Membrane_ElastHyper::Evaluate

/*----------------------------------------------------------------------*
 | evaluate strain energy function                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::StrainEnergy(
    LINALG::Matrix<3, 3>& cauchygreen,  ///< right Cauchy-Green tensor
    double& psi,                        ///< Strain energy function
    const int eleGID                    ///< Element GID
)
{
  // kinematic quantities and identity tensors
  LINALG::Matrix<3, 1> id2(true);
  LINALG::Matrix<3, 3> id4sharp(true);
  LINALG::Matrix<3, 1> rcg(true);
  double rcg33;
  LINALG::Matrix<3, 1> icg(true);
  EvaluateKinQuant(cauchygreen, id2, id4sharp, rcg, rcg33, icg);

  // Green-Lagrange strains matrix E = 0.5 * (Cauchy-Green - Identity)
  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  LINALG::Matrix<6, 1> glstrain(true);
  glstrain(0) = 0.5 * (rcg(0) - 1.0);
  glstrain(1) = 0.5 * (rcg(1) - 1.0);
  glstrain(2) = 0.5 * (rcg33 - 1.0);
  glstrain(3) = rcg(2);

  // principal isotropic invariants
  LINALG::Matrix<3, 1> prinv_iso(true);
  InvariantsPrincipal(prinv_iso, rcg, rcg33);

  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    // note that modified invariants equal the principal invariants as detF=J=1 (incompressibility)
    potsum_[p]->AddStrainEnergy(psi, prinv_iso, prinv_iso, glstrain, eleGID);
  }

  return;
}  // MAT::Membrane_ElastHyper::StrainEnergy

/*----------------------------------------------------------------------*
 | calculate kinematic quantities and identity tensors   sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::EvaluateKinQuant(
    const LINALG::Matrix<3, 3>& cauchygreen,  ///< right Cauchy-Green tensor
    LINALG::Matrix<3, 1>& id2,                ///< cartesian identity 2-tensor I_{AB}
    LINALG::Matrix<3, 3>& id4sharp,           ///< cartesian identity 4-tensor
    LINALG::Matrix<3, 1>& rcg,  ///< right Cauchy-Green in stress-like 3-Voigt notation
    double& rcg33,              ///< principal stretch in thickness direction
    LINALG::Matrix<3, 1>& icg   ///< inverse right Cauchy-Green in stress-like 3-Voigt notation
)
{
  // build Cartesian identity 2-tensor I_{AB}
  for (int i = 0; i < 2; i++) id2(i) = 1.0;

  // set Cartesian identity 4-tensor in 3-Voigt matrix notation
  // this is a fully 'contra-variant' identity tensor, ie I^{ABCD}
  // REMARK: rows are stress-like 3-Voigt
  //         columns are stress-like 3-Voigt
  for (int i = 0; i < 2; i++) id4sharp(i, i) = 1.0;
  for (int i = 2; i < 3; i++) id4sharp(i, i) = 0.5;

  // right Cauchy-Green
  // REMARK: stress-like 3-Voigt vector
  rcg(0) = cauchygreen(0, 0);
  rcg(1) = cauchygreen(1, 1);
  rcg(2) = cauchygreen(0, 1);

  // component in thickness direction of membrane
  // assuming incompressibility (J=detF=1)
  rcg33 = 1.0 / (rcg(0) * rcg(1) - std::pow(rcg(2), 2.0));

  // inverse right Cauchy-Green
  // REMARK: stress-like 3-Voigt vector
  icg(0) = rcg(1) * rcg33;
  icg(1) = rcg(0) * rcg33;
  icg(2) = -rcg(2) * rcg33;

  return;
}  // MAT::Membrane_ElastHyper::EvaluateKinQuant

/*----------------------------------------------------------------------*
 | calculate isotropic stress and elasticity tensor      sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::EvaluateIsotropicStressCmat(LINALG::Matrix<3, 1>& stress_iso,
    LINALG::Matrix<3, 3>& cmat_iso, const LINALG::Matrix<3, 1>& id2,
    const LINALG::Matrix<3, 3>& id4sharp, const LINALG::Matrix<3, 1>& rcg, const double& rcg33,
    const LINALG::Matrix<3, 1>& icg, int eleGID)
{
  // principal isotropic invariants
  LINALG::Matrix<3, 1> prinv_iso(true);
  InvariantsPrincipal(prinv_iso, rcg, rcg33);

  // 1st and 2nd derivative of the isotropic strain energy function
  LINALG::Matrix<2, 1> dPI_iso(true);
  LINALG::Matrix<3, 1> ddPII_iso(true);
  EvaluateInvariantDerivatives(prinv_iso, dPI_iso, ddPII_iso, eleGID);

  // stress and constitutive tensor factors according to Fakhreddine2011 equation (11,15)
  LINALG::Matrix<3, 1> gamma_iso(true);
  LINALG::Matrix<8, 1> delta_iso(true);
  CalculateGammaDelta(gamma_iso, delta_iso, prinv_iso, dPI_iso, ddPII_iso, rcg33);

  // isotropic 2nd Piola Kirchhoff stress
  stress_iso.Update(gamma_iso(0), id2, 1.0);
  stress_iso.Update(gamma_iso(1), rcg, 1.0);
  stress_iso.Update(gamma_iso(2), icg, 1.0);

  // isotropic constitutive tensor
  // contribution: Id \otimes Id
  cmat_iso.MultiplyNT(delta_iso(0), id2, id2, 0.0);
  // contribution: Id \otimes C + C \otimes Id
  cmat_iso.MultiplyNT(delta_iso(1), id2, rcg, 1.0);
  cmat_iso.MultiplyNT(delta_iso(1), rcg, id2, 1.0);
  // contribution: Id \otimes Cinv + Cinv \otimes Id
  cmat_iso.MultiplyNT(delta_iso(2), id2, icg, 1.0);
  cmat_iso.MultiplyNT(delta_iso(2), icg, id2, 1.0);
  // contribution: C \otimes C
  cmat_iso.MultiplyNT(delta_iso(3), rcg, rcg, 1.0);
  // contribution: C \otimes Cinv + Cinv \otimes C
  cmat_iso.MultiplyNT(delta_iso(4), rcg, icg, 1.0);
  cmat_iso.MultiplyNT(delta_iso(4), icg, rcg, 1.0);
  // contribution: Cinv \otimes Cinv
  cmat_iso.MultiplyNT(delta_iso(5), icg, icg, 1.0);
  // contribution: Cinv \odot Cinv
  cmat_iso(0, 0) += delta_iso(6) * icg(0) * icg(0);
  cmat_iso(0, 1) += delta_iso(6) * icg(2) * icg(2);
  cmat_iso(0, 2) += delta_iso(6) * icg(0) * icg(2);
  cmat_iso(1, 0) += delta_iso(6) * icg(2) * icg(2);
  cmat_iso(1, 1) += delta_iso(6) * icg(1) * icg(1);
  cmat_iso(1, 2) += delta_iso(6) * icg(1) * icg(2);
  cmat_iso(2, 0) += delta_iso(6) * icg(0) * icg(2);
  cmat_iso(2, 1) += delta_iso(6) * icg(1) * icg(2);
  cmat_iso(2, 2) += delta_iso(6) * 0.5 * (icg(0) * icg(1) + icg(2) * icg(2));
  // contribution: Id4^#
  cmat_iso.Update(delta_iso(7), id4sharp, 1.0);

  return;
}  // MAT::Membrane_ElastHyper::EvaluateIsotropicStressCmat

/*----------------------------------------------------------------------*
 | calculate principal invariants                        sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::InvariantsPrincipal(
    LINALG::Matrix<3, 1>& prinv,      ///< principal invariants
    const LINALG::Matrix<3, 1>& rcg,  ///< right cauchy-green in stress-like 3-Voigt notation
    const double& rcg33               ///< principal stretch in thickness direction
)
{
  prinv(0) = rcg(0) + rcg(1) + rcg33;
  prinv(1) =
      0.5 * (std::pow(prinv(0), 2.0) - (std::pow(rcg(0), 2.0) + std::pow(rcg(1), 2.0) +
                                           std::pow(rcg33, 2.0) + 2.0 * std::pow(rcg(2), 2.0)));
  prinv(2) = 1.0;  // incompressibility condition

  return;
}  // MAT::Membrane_ElastHyper::InvariantsPrincipal

/*----------------------------------------------------------------------*
 | calculate derivatives with respect to invariants                     |
 | from all materials of the elasthyper-toolbox          sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::EvaluateInvariantDerivatives(const LINALG::Matrix<3, 1>& prinv,
    LINALG::Matrix<2, 1>& dPI, LINALG::Matrix<3, 1>& ddPII, int eleGID)
{
  LINALG::Matrix<3, 1> dPI_full(true);
  LINALG::Matrix<6, 1> ddPII_full(true);

  // REMARK: since incompressibility (J=1) is assumed principal and modified invariants are equal
  // no transformation between principal and modified invariants needed below
  if (prinv(2) != 1.0)
    dserror("Incompressibility assumption not fulfilled in membrane hyperelastic material!");

  // derivatives of principal materials
  if (summandProperties_.isoprinc)
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->AddDerivativesPrincipal(dPI_full, ddPII_full, prinv, eleGID);
    }
  }

  // derivatives of decoupled (volumetric or isochoric) materials
  if (summandProperties_.isomod)
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->AddDerivativesModified(dPI_full, ddPII_full, prinv, eleGID);
    }
  }

  // the derivatives dPI_full(2), ddPII_full(2), ddPII_full(3) and ddPII_full(4) are w.r.t. prinv(2)
  // = 1.0 and thus not needed in this formulation
  dPI(0) = dPI_full(0);
  dPI(1) = dPI_full(1);

  ddPII(0) = ddPII_full(0);
  ddPII(1) = ddPII_full(1);
  ddPII(2) = ddPII_full(5);

  return;
}  // MAT::Membrane_ElastHyper::EvaluateInvariantDerivatives

/*----------------------------------------------------------------------*
 | calculate stress and constitutive tensor factors      sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::CalculateGammaDelta(LINALG::Matrix<3, 1>& gamma,
    LINALG::Matrix<8, 1>& delta, const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<2, 1>& dPI,
    const LINALG::Matrix<3, 1>& ddPII, const double rcg33)
{
  // according to Fakhreddine2011 equation (11)
  gamma(0) = 2.0 * (dPI(0) + prinv(0) * dPI(1));
  gamma(1) = -2.0 * dPI(1);
  gamma(2) = -rcg33 * gamma(0) - rcg33 * rcg33 * gamma(1);

  // according to Fakhreddine2011 equation (15)
  delta(0) = 4.0 * (ddPII(0) + 2.0 * prinv(0) * ddPII(2) + dPI(1) + prinv(0) * prinv(0) * ddPII(1));
  delta(1) = -4.0 * (ddPII(2) + prinv(0) * ddPII(1));
  delta(2) = -4.0 * rcg33 *
             (ddPII(0) + prinv(0) * ddPII(2) + dPI(1) +
                 (prinv(0) - rcg33) * (ddPII(2) + prinv(0) * ddPII(1)));
  delta(3) = 4.0 * ddPII(1);
  delta(4) = 4.0 * rcg33 * (ddPII(2) + (prinv(0) - rcg33) * ddPII(1));
  delta(5) = -2.0 * gamma(2) + 4.0 * rcg33 * rcg33 *
                                   (ddPII(0) + 2.0 * (prinv(0) - rcg33) * ddPII(2) +
                                       std::pow((prinv(0) - rcg33), 2.0) * ddPII(1));
  delta(6) = -2.0 * gamma(2);
  delta(7) = 2.0 * gamma(1);

  return;
}  // MAT::Membrane_ElastHyper::CalculateGammaDelta

/*----------------------------------------------------------------------*
 | calculate anisotropic stress and elasticity tensor    sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::Membrane_ElastHyper::EvaluateAnisotropicStressCmat(LINALG::Matrix<3, 1>& stress_aniso,
    LINALG::Matrix<3, 3>& cmat_aniso, LINALG::Matrix<3, 3>& Q_trafo,
    const LINALG::Matrix<3, 1>& rcg, const double& rcg33, Teuchos::ParameterList& params,
    int eleGID)
{
  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    // skip for materials without fiber
    if (fibervecs_[p].Norm2() == 0) continue;

    // fibervector in orthonormal frame on membrane surface
    LINALG::Matrix<3, 1> fibervector(true);
    fibervector.MultiplyTN(1.0, Q_trafo, fibervecs_[p], 0.0);

    // set new fibervector in anisotropic material
    potsum_[p]->SetFiberVecs(fibervector);

    // three dimensional right Cauchy-Green
    // REMARK: strain-like 6-Voigt vector
    // NOTE: rcg is a stress-like 3-Voigt vector
    LINALG::Matrix<6, 1> rcg_full(true);
    rcg_full(0) = rcg(0);
    rcg_full(1) = rcg(1);
    rcg_full(2) = rcg33;
    rcg_full(3) = 2.0 * rcg(2);

    // three dimensional anisotropic stress and constitutive tensor
    LINALG::Matrix<6, 1> stress_aniso_full(true);
    LINALG::Matrix<6, 6> cmat_aniso_full(true);

    potsum_[p]->AddStressAnisoPrincipal(
        rcg_full, cmat_aniso_full, stress_aniso_full, params, eleGID);

    // reduced anisotropic stress and constitutive tensor
    LINALG::Matrix<3, 1> stress_aniso_red(true);
    stress_aniso_red(0) = stress_aniso_full(0);
    stress_aniso_red(1) = stress_aniso_full(1);
    stress_aniso_red(2) = stress_aniso_full(3);

    LINALG::Matrix<3, 3> cmat_aniso_red(true);
    cmat_aniso_red(0, 0) = cmat_aniso_full(0, 0);
    cmat_aniso_red(0, 1) = cmat_aniso_full(0, 1);
    cmat_aniso_red(0, 2) = cmat_aniso_full(0, 3);
    cmat_aniso_red(1, 0) = cmat_aniso_full(1, 0);
    cmat_aniso_red(1, 1) = cmat_aniso_full(1, 1);
    cmat_aniso_red(1, 2) = cmat_aniso_full(1, 3);
    cmat_aniso_red(2, 0) = cmat_aniso_full(3, 0);
    cmat_aniso_red(2, 1) = cmat_aniso_full(3, 1);
    cmat_aniso_red(2, 2) = cmat_aniso_full(3, 3);

    // anisotropic 2nd Piola Kirchhoff stress
    stress_aniso.Update(1.0, stress_aniso_red, 1.0);

    // anisotropic constitutive tensor
    cmat_aniso.Update(1.0, cmat_aniso_red, 1.0);
  }

  return;
}  // MAT::Membrane_ElastHyper::EvaluateAnisotropicStressCmat
