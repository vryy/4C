/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an interface class for materials of the (visco)elasthyper toolbox.

\level 1
*/
/*----------------------------------------------------------------------*/

#include "matelast_aniso_structuraltensor_strategy.H"
#include "matelast_anisoactivestress_evolution.H"
#include "matelast_coup13apow.H"
#include "matelast_coup1pow.H"
#include "matelast_coup2pow.H"
#include "matelast_coup3pow.H"
#include "matelast_coupSaintVenantKirchhoff.H"
#include "matelast_coupanisoexpo.H"
#include "matelast_coupanisoexposhear.H"
#include "matelast_coupanisoexpotwocoup.H"
#include "matelast_coupanisoneohooke.H"
#include "matelast_coupanisoneohooke_VarProp.H"
#include "matelast_coupanisopow.H"
#include "matelast_coupblatzko.H"
#include "matelast_coupexppol.H"
#include "matelast_couplogmixneohooke.H"
#include "matelast_couplogneohooke.H"
#include "matelast_coupmooneyrivlin.H"
#include "matelast_coupneohooke.H"
#include "matelast_coupsimopister.H"
#include "matelast_couptransverselyisotropic.H"
#include "matelast_coupvarga.H"
#include "matelast_iso1pow.H"
#include "matelast_iso2pow.H"
#include "matelast_isoanisoexpo.H"
#include "matelast_isoexpopow.H"
#include "matelast_isomooneyrivlin.H"
#include "matelast_isomuscle_blemker.H"
#include "matelast_isoneohooke.H"
#include "matelast_isoogden.H"
#include "matelast_isotestmaterial.H"
#include "matelast_isovarga.H"
#include "matelast_isovolaaagasser.H"
#include "matelast_isoyeoh.H"
#include "matelast_remodelfiber.H"
#include "matelast_vologden.H"
#include "matelast_volpenalty.H"
#include "matelast_volpow.H"
#include "matelast_volsussmanbathe.H"
#include "matelast_visco_coupmyocard.H"
#include "matelast_visco_fract.H"
#include "matelast_visco_generalizedgenmax.H"
#include "matelast_visco_genmax.H"
#include "matelast_visco_isoratedep.H"

#include "lib_globalproblem.H"
#include "lib_linedefinition.H"
#include "mat_par_bundle.H"
#include "matelast_summand.H"

Teuchos::RCP<MAT::ELASTIC::Summand> MAT::ELASTIC::Summand::Factory(int matnum)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");

  // yet another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      DRT::Problem::Instance(probinst)->Materials()->ById(matnum);

  // construct structural tensor strategy for anisotropic materials
  switch (curmat->Type())
  {
    case INPAR::MAT::mes_isoanisoexpo:
    case INPAR::MAT::mes_isomuscleblemker:
    case INPAR::MAT::mes_coupanisoexpo:
    case INPAR::MAT::mes_coupanisoexpoactive:
    case INPAR::MAT::mes_coupanisoexpotwocoup:
    case INPAR::MAT::mes_coupanisoneohooke:
    case INPAR::MAT::mes_coupanisopow:
    case INPAR::MAT::mes_coupanisoneohooke_varprop:
    {
      break;
    }
    default:
      break;
  }

  switch (curmat->Type())
  {
    case INPAR::MAT::mes_anisoactivestress_evolution:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::AnisoActiveStress_Evolution(curmat));
      auto* params =
          dynamic_cast<MAT::ELASTIC::PAR::AnisoActiveStress_Evolution*>(curmat->Parameter());
      return Teuchos::rcp(new AnisoActiveStress_Evolution(params));
    }
    case INPAR::MAT::mes_coupanisoexpoactive:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoExpoActive(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupAnisoExpoActive*>(curmat->Parameter());
      return Teuchos::rcp(new CoupAnisoExpoActive(params));
    }
    case INPAR::MAT::mes_coupanisoexpo:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoExpo(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupAnisoExpo*>(curmat->Parameter());
      return Teuchos::rcp(new CoupAnisoExpo(params));
    }
    case INPAR::MAT::mes_coupanisoexposhear:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoExpoShear(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupAnisoExpoShear*>(curmat->Parameter());
      return Teuchos::rcp(new CoupAnisoExpoShear(params));
    }
    case INPAR::MAT::mes_coupanisoexpotwocoup:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup*>(curmat->Parameter());
      return Teuchos::rcp(new CoupAnisoExpoTwoCoup(params));
    }
    case INPAR::MAT::mes_coupanisoneohooke:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoNeoHooke(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupAnisoNeoHooke*>(curmat->Parameter());
      return Teuchos::rcp(new CoupAnisoNeoHooke(params));
    }
    case INPAR::MAT::mes_coupanisoneohooke_varprop:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp(curmat));
      auto* params =
          dynamic_cast<MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp*>(curmat->Parameter());
      return Teuchos::rcp(new CoupAnisoNeoHooke_VarProp(params));
    }
    case INPAR::MAT::mes_coupanisopow:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoPow(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupAnisoPow*>(curmat->Parameter());
      return Teuchos::rcp(new CoupAnisoPow(params));
    }
    case INPAR::MAT::mes_coupblatzko:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupBlatzKo(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupBlatzKo*>(curmat->Parameter());
      return Teuchos::rcp(new CoupBlatzKo(params));
    }
    case INPAR::MAT::mes_coupexppol:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupExpPol(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupExpPol*>(curmat->Parameter());
      return Teuchos::rcp(new CoupExpPol(params));
    }
    case INPAR::MAT::mes_couplogneohooke:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupLogNeoHooke(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupLogNeoHooke*>(curmat->Parameter());
      return Teuchos::rcp(new CoupLogNeoHooke(params));
    }
    case INPAR::MAT::mes_couplogmixneohooke:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupLogMixNeoHooke(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupLogMixNeoHooke*>(curmat->Parameter());
      return Teuchos::rcp(new CoupLogMixNeoHooke(params));
    }
    case INPAR::MAT::mes_coupmooneyrivlin:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupMooneyRivlin(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupMooneyRivlin*>(curmat->Parameter());
      return Teuchos::rcp(new CoupMooneyRivlin(params));
    }
    case INPAR::MAT::mes_coupmyocard:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupMyocard(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupMyocard*>(curmat->Parameter());
      return Teuchos::rcp(new CoupMyocard(params));
    }
    case INPAR::MAT::mes_coupneohooke:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupNeoHooke(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupNeoHooke*>(curmat->Parameter());
      return Teuchos::rcp(new CoupNeoHooke(params));
    }
    case INPAR::MAT::mes_coup1pow:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::Coup1Pow(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::Coup1Pow*>(curmat->Parameter());
      return Teuchos::rcp(new Coup1Pow(params));
    }
    case INPAR::MAT::mes_coup2pow:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::Coup2Pow(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::Coup2Pow*>(curmat->Parameter());
      return Teuchos::rcp(new Coup2Pow(params));
    }
    case INPAR::MAT::mes_coup3pow:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::Coup3Pow(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::Coup3Pow*>(curmat->Parameter());
      return Teuchos::rcp(new Coup3Pow(params));
    }
    case INPAR::MAT::mes_coup13apow:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::Coup13aPow(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::Coup13aPow*>(curmat->Parameter());
      return Teuchos::rcp(new Coup13aPow(params));
    }
    case INPAR::MAT::mes_coupsimopister:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupSimoPister(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupSimoPister*>(curmat->Parameter());
      return Teuchos::rcp(new CoupSimoPister(params));
    }
    case INPAR::MAT::mes_coupSVK:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupSVK(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupSVK*>(curmat->Parameter());
      return Teuchos::rcp(new CoupSVK(params));
    }
    case INPAR::MAT::mes_couptransverselyisotropic:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupTransverselyIsotropic(curmat));
      auto* params =
          dynamic_cast<MAT::ELASTIC::PAR::CoupTransverselyIsotropic*>(curmat->Parameter());
      return Teuchos::rcp(new CoupTransverselyIsotropic(params));
    }
    case INPAR::MAT::mes_coupvarga:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupVarga(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::CoupVarga*>(curmat->Parameter());
      return Teuchos::rcp(new CoupVarga(params));
    }
    case INPAR::MAT::mes_fract:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::Fract(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::Fract*>(curmat->Parameter());
      return Teuchos::rcp(new Fract(params));
    }
    case INPAR::MAT::mes_genmax:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::GenMax(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::GenMax*>(curmat->Parameter());
      return Teuchos::rcp(new GenMax(params));
    }
    case INPAR::MAT::mes_generalizedgenmax:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::GeneralizedGenMax(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::GeneralizedGenMax*>(curmat->Parameter());
      return Teuchos::rcp(new GeneralizedGenMax(params));
    }
    case INPAR::MAT::mes_isoanisoexpo:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoAnisoExpo(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::IsoAnisoExpo*>(curmat->Parameter());
      return Teuchos::rcp(new IsoAnisoExpo(params));
    }
    case INPAR::MAT::mes_isoexpopow:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoExpoPow(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::IsoExpoPow*>(curmat->Parameter());
      return Teuchos::rcp(new IsoExpoPow(params));
    }
    case INPAR::MAT::mes_isomooneyrivlin:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoMooneyRivlin(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::IsoMooneyRivlin*>(curmat->Parameter());
      return Teuchos::rcp(new IsoMooneyRivlin(params));
    }
    case INPAR::MAT::mes_isomuscleblemker:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoMuscleBlemker(curmat));
      MAT::ELASTIC::PAR::IsoMuscleBlemker* params =
          static_cast<MAT::ELASTIC::PAR::IsoMuscleBlemker*>(curmat->Parameter());
      return Teuchos::rcp(new IsoMuscleBlemker(params));
    }
    case INPAR::MAT::mes_isoneohooke:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoNeoHooke(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::IsoNeoHooke*>(curmat->Parameter());
      return Teuchos::rcp(new IsoNeoHooke(params));
    }
    case INPAR::MAT::mes_isoogden:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoOgden(curmat));
      auto* params = static_cast<MAT::ELASTIC::PAR::IsoOgden*>(curmat->Parameter());
      return Teuchos::rcp(new IsoOgden(params));
    }
    case INPAR::MAT::mes_iso1pow:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::Iso1Pow(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::Iso1Pow*>(curmat->Parameter());
      return Teuchos::rcp(new Iso1Pow(params));
    }
    case INPAR::MAT::mes_iso2pow:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::Iso2Pow(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::Iso2Pow*>(curmat->Parameter());
      return Teuchos::rcp(new Iso2Pow(params));
    }
    case INPAR::MAT::mes_isoratedep:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoRateDep(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::IsoRateDep*>(curmat->Parameter());
      return Teuchos::rcp(new IsoRateDep(params));
    }
    case INPAR::MAT::mes_isotestmaterial:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoTestMaterial(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::IsoTestMaterial*>(curmat->Parameter());
      return Teuchos::rcp(new IsoTestMaterial(params));
    }
    case INPAR::MAT::mes_isovarga:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoVarga(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::IsoVarga*>(curmat->Parameter());
      return Teuchos::rcp(new IsoVarga(params));
    }
    case INPAR::MAT::mes_isovolaaagasser:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoVolAAAGasser(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::IsoVolAAAGasser*>(curmat->Parameter());
      return Teuchos::rcp(new IsoVolAAAGasser(params));
    }
    case INPAR::MAT::mes_isoyeoh:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::IsoYeoh(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::IsoYeoh*>(curmat->Parameter());
      return Teuchos::rcp(new IsoYeoh(params));
    }
    case INPAR::MAT::mes_remodelfiber:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::RemodelFiber(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::RemodelFiber*>(curmat->Parameter());
      return Teuchos::rcp(new RemodelFiber(params));
    }
    case INPAR::MAT::mes_vologden:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::VolOgden(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::VolOgden*>(curmat->Parameter());
      return Teuchos::rcp(new VolOgden(params));
    }
    case INPAR::MAT::mes_volpenalty:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::VolPenalty(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::VolPenalty*>(curmat->Parameter());
      return Teuchos::rcp(new VolPenalty(params));
    }
    case INPAR::MAT::mes_volpow:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::VolPow(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::VolPow*>(curmat->Parameter());
      return Teuchos::rcp(new VolPow(params));
    }
    case INPAR::MAT::mes_volsussmanbathe:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::VolSussmanBathe(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::VolSussmanBathe*>(curmat->Parameter());
      return Teuchos::rcp(new VolSussmanBathe(params));
    }
    case INPAR::MAT::mes_viscobranch:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::ViscoBranch(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::ViscoBranch*>(curmat->Parameter());
      return Teuchos::rcp(new ViscoBranch(params));
    }
    case INPAR::MAT::mes_viscopart:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::ELASTIC::PAR::ViscoPart(curmat));
      auto* params = dynamic_cast<MAT::ELASTIC::PAR::ViscoPart*>(curmat->Parameter());
      return Teuchos::rcp(new ViscoPart(params));
    }
    default:
      dserror("cannot deal with type %d", curmat->Type());
  }
  return Teuchos::null;
}

void MAT::ELASTIC::Summand::AddShearMod(bool& haveshearmod, double& shearmod) const
{
  dserror("MAT::ELASTIC::Summand::AddShearMod: Add Shear Modulus not implemented - do so!");
}

int MAT::ELASTIC::Summand::UniqueParObjectId() const { return -1; }

void MAT::ELASTIC::Summand::Pack(DRT::PackBuffer& data) const { return; }

void MAT::ELASTIC::Summand::Unpack(const std::vector<char>& data) { return; };


// Function which reads in the given fiber value due to the FIBER1 nomenclature
void MAT::ELASTIC::Summand::ReadFiber(DRT::INPUT::LineDefinition* linedef,
    const std::string& specifier, CORE::LINALG::Matrix<3, 1>& fiber_vector)
{
  std::vector<double> fiber1;
  linedef->ExtractDoubleVector(specifier, fiber1);
  double f1norm = 0.;
  // normalization
  for (int i = 0; i < 3; ++i)
  {
    f1norm += fiber1[i] * fiber1[i];
  }
  f1norm = sqrt(f1norm);

  // fill final fiber vector
  for (int i = 0; i < 3; ++i) fiber_vector(i) = fiber1[i] / f1norm;
}

// Function which reads in the given fiber value due to the CIR-AXI-RAD nomenclature
void MAT::ELASTIC::Summand::ReadRadAxiCir(
    DRT::INPUT::LineDefinition* linedef, CORE::LINALG::Matrix<3, 3>& locsys)
{
  // read local (cylindrical) cosy-directions at current element
  // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
  CORE::LINALG::Matrix<3, 1> fiber_rad;
  CORE::LINALG::Matrix<3, 1> fiber_axi;
  CORE::LINALG::Matrix<3, 1> fiber_cir;

  ReadFiber(linedef, "RAD", fiber_rad);
  ReadFiber(linedef, "AXI", fiber_axi);
  ReadFiber(linedef, "CIR", fiber_cir);

  for (int i = 0; i < 3; ++i)
  {
    locsys(i, 0) = fiber_rad(i);
    locsys(i, 1) = fiber_axi(i);
    locsys(i, 2) = fiber_cir(i);
  }
}

void MAT::ELASTIC::Summand::EvaluateFirstDerivativesAniso(CORE::LINALG::Matrix<2, 1>& dPI_aniso,
    CORE::LINALG::Matrix<3, 3> const& rcg, int gp, int eleGID)
{
  bool isoprinc, isomod, anisoprinc, anisomod, viscogeneral;
  SpecifyFormulation(isoprinc, isomod, anisoprinc, anisomod, viscogeneral);
  if (anisoprinc or anisomod)
  {
    dserror(
        "This anisotropic material does not support the first derivative of the free-energy "
        "function with respect to the anisotropic invariants. You need to implement them.");
  }
}

void MAT::ELASTIC::Summand::EvaluateSecondDerivativesAniso(CORE::LINALG::Matrix<3, 1>& ddPII_aniso,
    CORE::LINALG::Matrix<3, 3> const& rcg, int gp, int eleGID)
{
  bool isoprinc, isomod, anisoprinc, anisomod, viscogeneral;
  SpecifyFormulation(isoprinc, isomod, anisoprinc, anisomod, viscogeneral);
  if (anisoprinc or anisomod)
  {
    dserror(
        "This anisotropic material does not support the second derivative of the free-energy "
        "function with respect to the anisotropic invariants. You need to implement them.");
  }
}