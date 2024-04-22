/*----------------------------------------------------------------------*/
/*! \file
\brief hyperelastic toolbox for membranes assuming incompressibility and plane stress

The input line should read
MAT 0 MAT_Membrane_ElastHyper NUMMAT 2 MATIDS 1 2 DENS 0

\level 3


*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                               sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
#include "baci_mat_membrane_elasthyper.hpp"

#include "baci_mat_membrane_elasthyper_service.hpp"
#include "baci_matelast_summand.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | constructor                                           sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
MAT::PAR::MembraneElastHyper::MembraneElastHyper(Teuchos::RCP<MAT::PAR::Material> matdata)
    : MAT::PAR::ElastHyper(matdata)
{
  return;
}  // MAT::PAR::MembraneElastHyper::MembraneElastHyper

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::MembraneElastHyper::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MembraneElastHyper(this));
}  // MAT::PAR::MembraneElastHyper::CreateMaterial


MAT::MembraneElastHyperType MAT::MembraneElastHyperType::instance_;


CORE::COMM::ParObject* MAT::MembraneElastHyperType::Create(const std::vector<char>& data)
{
  MAT::MembraneElastHyper* memelhy = new MAT::MembraneElastHyper();
  memelhy->Unpack(data);

  return memelhy;
}  // MAT::Membrane_ElastHyperType::Create

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
MAT::MembraneElastHyper::MembraneElastHyper() : MAT::ElastHyper(), fibervecs_(true)
{
  return;
}  // MAT::MembraneElastHyper::MembraneElastHyper()

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
MAT::MembraneElastHyper::MembraneElastHyper(MAT::PAR::MembraneElastHyper* params)
    : MAT::ElastHyper(params), fibervecs_(true)
{
  return;
}  // MAT::MembraneElastHyper::MembraneElastHyper()

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::MembraneElastHyper::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // add base class Element
  MAT::ElastHyper::Pack(data);

  AddtoPack(data, fibervecs_);

  return;
}  // MAT::MembraneElastHyper::Pack()

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::MembraneElastHyper::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  MAT::ElastHyper::Unpack(basedata);

  ExtractfromPack(position, data, fibervecs_);

  return;
}  // MAT::MembraneElastHyper::Unpack()

/*----------------------------------------------------------------------*
 |                                                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::MembraneElastHyper::Setup(int numgp, INPUT::LineDefinition* linedef)
{
  // call setup of base class
  MAT::ElastHyper::Setup(numgp, linedef);

  GetFiberVecs(fibervecs_);

  return;
}  // MAT::MembraneElastHyper::Setup()

/*----------------------------------------------------------------------*
 | hyperelastic stress response plus elasticity tensor   sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::MembraneElastHyper::EvaluateMembrane(const CORE::LINALG::Matrix<3, 3>& defgrd,
    const CORE::LINALG::Matrix<3, 3>& cauchygreen, Teuchos::ParameterList& params,
    const CORE::LINALG::Matrix<3, 3>& Q_trafo, CORE::LINALG::Matrix<3, 1>& stress,
    CORE::LINALG::Matrix<3, 3>& cmat, const int gp, const int eleGID)
{
  // blank resulting quantities
  stress.Clear();
  cmat.Clear();

  // kinematic quantities and identity tensors
  CORE::LINALG::Matrix<3, 1> id2(true);
  CORE::LINALG::Matrix<3, 3> id4sharp(true);
  CORE::LINALG::Matrix<3, 1> rcg(true);
  double rcg33;
  CORE::LINALG::Matrix<3, 1> icg(true);
  MembraneElastHyperEvaluateKinQuant(cauchygreen, id2, id4sharp, rcg, rcg33, icg);

  // evaluate isotropic 2nd Piola-Kirchhoff stress and constitutive tensor
  CORE::LINALG::Matrix<3, 1> stress_iso(true);
  CORE::LINALG::Matrix<3, 3> cmat_iso(true);
  MembraneElastHyperEvaluateIsotropicStressCmat(stress_iso, cmat_iso, id2, id4sharp, rcg, rcg33,
      icg, gp, eleGID, potsum_, summandProperties_);

  // update 2nd Piola-Kirchhoff stress and constitutive tensor
  stress.Update(1.0, stress_iso, 1.0);
  cmat.Update(1.0, cmat_iso, 1.0);

  // evaluate anisotropic 2nd Piola-Kirchhoff stress and constitutive tensor
  if (summandProperties_.anisoprinc)
  {
    CORE::LINALG::Matrix<3, 1> stress_aniso(true);
    CORE::LINALG::Matrix<3, 3> cmat_aniso(true);
    EvaluateAnisotropicStressCmat(
        stress_aniso, cmat_aniso, Q_trafo, rcg, rcg33, params, gp, eleGID);

    // update 2nd Piola-Kirchhoff stress and constitutive tensor
    stress.Update(1.0, stress_aniso, 1.0);
    cmat.Update(1.0, cmat_aniso, 1.0);
  }
  if (summandProperties_.anisomod)
  {
    FOUR_C_THROW("anisomod_ not implemented for membrane elasthyper materials!");
  }

  return;
}  // MAT::MembraneElastHyper::Evaluate

/*----------------------------------------------------------------------*
 | evaluate strain energy function                       sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::MembraneElastHyper::StrainEnergy(
    CORE::LINALG::Matrix<3, 3>& cauchygreen, double& psi, const int gp, const int eleGID)
{
  // kinematic quantities and identity tensors
  CORE::LINALG::Matrix<3, 1> id2(true);
  CORE::LINALG::Matrix<3, 3> id4sharp(true);
  CORE::LINALG::Matrix<3, 1> rcg(true);
  double rcg33;
  CORE::LINALG::Matrix<3, 1> icg(true);
  MembraneElastHyperEvaluateKinQuant(cauchygreen, id2, id4sharp, rcg, rcg33, icg);

  // Green-Lagrange strains matrix E = 0.5 * (Cauchy-Green - Identity)
  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  CORE::LINALG::Matrix<6, 1> glstrain(true);
  glstrain(0) = 0.5 * (rcg(0) - 1.0);
  glstrain(1) = 0.5 * (rcg(1) - 1.0);
  glstrain(2) = 0.5 * (rcg33 - 1.0);
  glstrain(3) = rcg(2);

  // principal isotropic invariants
  CORE::LINALG::Matrix<3, 1> prinv_iso(true);
  MembraneElastHyperInvariantsPrincipal(prinv_iso, rcg, rcg33);

  // loop map of associated potential summands
  for (auto& p : potsum_)
  {
    // note that modified invariants equal the principal invariants as detF=J=1 (incompressibility)
    p->AddStrainEnergy(psi, prinv_iso, prinv_iso, glstrain, gp, eleGID);
  }
}  // MAT::MembraneElastHyper::StrainEnergy

/*----------------------------------------------------------------------*
 | calculate anisotropic stress and elasticity tensor    sfuchs 08/2017 |
 *----------------------------------------------------------------------*/
void MAT::MembraneElastHyper::EvaluateAnisotropicStressCmat(
    CORE::LINALG::Matrix<3, 1>& stress_aniso, CORE::LINALG::Matrix<3, 3>& cmat_aniso,
    const CORE::LINALG::Matrix<3, 3>& Q_trafo, const CORE::LINALG::Matrix<3, 1>& rcg,
    const double& rcg33, Teuchos::ParameterList& params, const int gp, int eleGID)
{
  // loop map of associated potential summands
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    // skip for materials without fiber
    if (fibervecs_[p].Norm2() == 0) continue;

    // fibervector in orthonormal frame on membrane surface
    CORE::LINALG::Matrix<3, 1> fibervector(true);
    fibervector.MultiplyTN(1.0, Q_trafo, fibervecs_[p], 0.0);

    // set new fibervector in anisotropic material
    potsum_[p]->SetFiberVecs(fibervector);

    // three dimensional right Cauchy-Green
    // REMARK: strain-like 6-Voigt vector
    // NOTE: rcg is a stress-like 3-Voigt vector
    CORE::LINALG::Matrix<6, 1> rcg_full(true);
    rcg_full(0) = rcg(0);
    rcg_full(1) = rcg(1);
    rcg_full(2) = rcg33;
    rcg_full(3) = 2.0 * rcg(2);

    // three dimensional anisotropic stress and constitutive tensor
    CORE::LINALG::Matrix<6, 1> stress_aniso_full(true);
    CORE::LINALG::Matrix<6, 6> cmat_aniso_full(true);

    potsum_[p]->AddStressAnisoPrincipal(
        rcg_full, cmat_aniso_full, stress_aniso_full, params, gp, eleGID);

    // reduced anisotropic stress and constitutive tensor
    CORE::LINALG::Matrix<3, 1> stress_aniso_red(true);
    stress_aniso_red(0) = stress_aniso_full(0);
    stress_aniso_red(1) = stress_aniso_full(1);
    stress_aniso_red(2) = stress_aniso_full(3);

    CORE::LINALG::Matrix<3, 3> cmat_aniso_red(true);
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
}  // MAT::MembraneElastHyper::EvaluateAnisotropicStressCmat

FOUR_C_NAMESPACE_CLOSE
