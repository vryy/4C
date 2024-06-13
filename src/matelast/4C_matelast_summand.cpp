/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an interface class for materials of the (visco)elasthyper toolbox.

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_summand.hpp"

#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_matelast_anisoactivestress_evolution.hpp"
#include "4C_matelast_coup13apow.hpp"
#include "4C_matelast_coup1pow.hpp"
#include "4C_matelast_coup2pow.hpp"
#include "4C_matelast_coup3pow.hpp"
#include "4C_matelast_coupanisoexpo.hpp"
#include "4C_matelast_coupanisoexposhear.hpp"
#include "4C_matelast_coupanisoexpotwocoup.hpp"
#include "4C_matelast_coupanisoneohooke.hpp"
#include "4C_matelast_coupanisoneohooke_VarProp.hpp"
#include "4C_matelast_coupanisopow.hpp"
#include "4C_matelast_coupblatzko.hpp"
#include "4C_matelast_coupexppol.hpp"
#include "4C_matelast_couplogmixneohooke.hpp"
#include "4C_matelast_couplogneohooke.hpp"
#include "4C_matelast_coupmooneyrivlin.hpp"
#include "4C_matelast_coupneohooke.hpp"
#include "4C_matelast_coupSaintVenantKirchhoff.hpp"
#include "4C_matelast_coupsimopister.hpp"
#include "4C_matelast_couptransverselyisotropic.hpp"
#include "4C_matelast_coupvarga.hpp"
#include "4C_matelast_iso1pow.hpp"
#include "4C_matelast_iso2pow.hpp"
#include "4C_matelast_isoanisoexpo.hpp"
#include "4C_matelast_isoexpopow.hpp"
#include "4C_matelast_isomooneyrivlin.hpp"
#include "4C_matelast_isomuscle_blemker.hpp"
#include "4C_matelast_isoneohooke.hpp"
#include "4C_matelast_isoogden.hpp"
#include "4C_matelast_isotestmaterial.hpp"
#include "4C_matelast_isovarga.hpp"
#include "4C_matelast_isoyeoh.hpp"
#include "4C_matelast_remodelfiber.hpp"
#include "4C_matelast_visco_coupmyocard.hpp"
#include "4C_matelast_visco_fract.hpp"
#include "4C_matelast_visco_generalizedgenmax.hpp"
#include "4C_matelast_visco_genmax.hpp"
#include "4C_matelast_visco_isoratedep.hpp"
#include "4C_matelast_vologden.hpp"
#include "4C_matelast_volpenalty.hpp"
#include "4C_matelast_volpow.hpp"
#include "4C_matelast_volsussmanbathe.hpp"

FOUR_C_NAMESPACE_OPEN

Teuchos::RCP<Mat::Elastic::Summand> Mat::Elastic::Summand::Factory(int matnum)
{
  // for the sake of safety
  if (Global::Problem::Instance()->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");

  // yet another safety check
  if (Global::Problem::Instance()->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  // retrieve problem instance to read from
  const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(matnum);

  // construct structural tensor strategy for anisotropic materials
  switch (curmat->Type())
  {
    case Core::Materials::mes_isoanisoexpo:
    case Core::Materials::mes_isomuscleblemker:
    case Core::Materials::mes_coupanisoexpo:
    case Core::Materials::mes_coupanisoexpoactive:
    case Core::Materials::mes_coupanisoexpotwocoup:
    case Core::Materials::mes_coupanisoneohooke:
    case Core::Materials::mes_coupanisopow:
    case Core::Materials::mes_coupanisoneohooke_varprop:
    {
      break;
    }
    default:
      break;
  }

  switch (curmat->Type())
  {
    case Core::Materials::mes_anisoactivestress_evolution:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::AnisoActiveStressEvolution*>(curmat);
      return Teuchos::rcp(new AnisoActiveStressEvolution(params));
    }
    case Core::Materials::mes_coupanisoexpoactive:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupAnisoExpoActive*>(curmat);
      return Teuchos::rcp(new CoupAnisoExpoActive(params));
    }
    case Core::Materials::mes_coupanisoexpo:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupAnisoExpo*>(curmat);
      return Teuchos::rcp(new CoupAnisoExpo(params));
    }
    case Core::Materials::mes_coupanisoexposhear:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupAnisoExpoShear*>(curmat);
      return Teuchos::rcp(new CoupAnisoExpoShear(params));
    }
    case Core::Materials::mes_coupanisoexpotwocoup:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupAnisoExpoTwoCoup*>(curmat);
      return Teuchos::rcp(new CoupAnisoExpoTwoCoup(params));
    }
    case Core::Materials::mes_coupanisoneohooke:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupAnisoNeoHooke*>(curmat);
      return Teuchos::rcp(new CoupAnisoNeoHooke(params));
    }
    case Core::Materials::mes_coupanisoneohooke_varprop:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupAnisoNeoHookeVarProp*>(curmat);
      return Teuchos::rcp(new CoupAnisoNeoHookeVarProp(params));
    }
    case Core::Materials::mes_coupanisopow:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupAnisoPow*>(curmat);
      return Teuchos::rcp(new CoupAnisoPow(params));
    }
    case Core::Materials::mes_coupblatzko:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupBlatzKo*>(curmat);
      return Teuchos::rcp(new CoupBlatzKo(params));
    }
    case Core::Materials::mes_coupexppol:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupExpPol*>(curmat);
      return Teuchos::rcp(new CoupExpPol(params));
    }
    case Core::Materials::mes_couplogneohooke:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupLogNeoHooke*>(curmat);
      return Teuchos::rcp(new CoupLogNeoHooke(params));
    }
    case Core::Materials::mes_couplogmixneohooke:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupLogMixNeoHooke*>(curmat);
      return Teuchos::rcp(new CoupLogMixNeoHooke(params));
    }
    case Core::Materials::mes_coupmooneyrivlin:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupMooneyRivlin*>(curmat);
      return Teuchos::rcp(new CoupMooneyRivlin(params));
    }
    case Core::Materials::mes_coupmyocard:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupMyocard*>(curmat);
      return Teuchos::rcp(new CoupMyocard(params));
    }
    case Core::Materials::mes_coupneohooke:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupNeoHooke*>(curmat);
      return Teuchos::rcp(new CoupNeoHooke(params));
    }
    case Core::Materials::mes_coup1pow:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::Coup1Pow*>(curmat);
      return Teuchos::rcp(new Coup1Pow(params));
    }
    case Core::Materials::mes_coup2pow:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::Coup2Pow*>(curmat);
      return Teuchos::rcp(new Coup2Pow(params));
    }
    case Core::Materials::mes_coup3pow:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::Coup3Pow*>(curmat);
      return Teuchos::rcp(new Coup3Pow(params));
    }
    case Core::Materials::mes_coup13apow:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::Coup13aPow*>(curmat);
      return Teuchos::rcp(new Coup13aPow(params));
    }
    case Core::Materials::mes_coupsimopister:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupSimoPister*>(curmat);
      return Teuchos::rcp(new CoupSimoPister(params));
    }
    case Core::Materials::mes_coupSVK:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupSVK*>(curmat);
      return Teuchos::rcp(new CoupSVK(params));
    }
    case Core::Materials::mes_couptransverselyisotropic:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupTransverselyIsotropic*>(curmat);
      return Teuchos::rcp(new CoupTransverselyIsotropic(params));
    }
    case Core::Materials::mes_coupvarga:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::CoupVarga*>(curmat);
      return Teuchos::rcp(new CoupVarga(params));
    }
    case Core::Materials::mes_fract:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::Fract*>(curmat);
      return Teuchos::rcp(new Fract(params));
    }
    case Core::Materials::mes_genmax:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::GenMax*>(curmat);
      return Teuchos::rcp(new GenMax(params));
    }
    case Core::Materials::mes_generalizedgenmax:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::GeneralizedGenMax*>(curmat);
      return Teuchos::rcp(new GeneralizedGenMax(params));
    }
    case Core::Materials::mes_isoanisoexpo:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::IsoAnisoExpo*>(curmat);
      return Teuchos::rcp(new IsoAnisoExpo(params));
    }
    case Core::Materials::mes_isoexpopow:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::IsoExpoPow*>(curmat);
      return Teuchos::rcp(new IsoExpoPow(params));
    }
    case Core::Materials::mes_isomooneyrivlin:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::IsoMooneyRivlin*>(curmat);
      return Teuchos::rcp(new IsoMooneyRivlin(params));
    }
    case Core::Materials::mes_isomuscleblemker:
    {
      Mat::Elastic::PAR::IsoMuscleBlemker* params =
          static_cast<Mat::Elastic::PAR::IsoMuscleBlemker*>(curmat);
      return Teuchos::rcp(new IsoMuscleBlemker(params));
    }
    case Core::Materials::mes_isoneohooke:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::IsoNeoHooke*>(curmat);
      return Teuchos::rcp(new IsoNeoHooke(params));
    }
    case Core::Materials::mes_isoogden:
    {
      auto* params = static_cast<Mat::Elastic::PAR::IsoOgden*>(curmat);
      return Teuchos::rcp(new IsoOgden(params));
    }
    case Core::Materials::mes_iso1pow:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::Iso1Pow*>(curmat);
      return Teuchos::rcp(new Iso1Pow(params));
    }
    case Core::Materials::mes_iso2pow:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::Iso2Pow*>(curmat);
      return Teuchos::rcp(new Iso2Pow(params));
    }
    case Core::Materials::mes_isoratedep:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::IsoRateDep*>(curmat);
      return Teuchos::rcp(new IsoRateDep(params));
    }
    case Core::Materials::mes_isotestmaterial:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::IsoTestMaterial*>(curmat);
      return Teuchos::rcp(new IsoTestMaterial(params));
    }
    case Core::Materials::mes_isovarga:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::IsoVarga*>(curmat);
      return Teuchos::rcp(new IsoVarga(params));
    }
    case Core::Materials::mes_isoyeoh:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::IsoYeoh*>(curmat);
      return Teuchos::rcp(new IsoYeoh(params));
    }
    case Core::Materials::mes_remodelfiber:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::RemodelFiber*>(curmat);
      return Teuchos::rcp(new RemodelFiber(params));
    }
    case Core::Materials::mes_vologden:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::VolOgden*>(curmat);
      return Teuchos::rcp(new VolOgden(params));
    }
    case Core::Materials::mes_volpenalty:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::VolPenalty*>(curmat);
      return Teuchos::rcp(new VolPenalty(params));
    }
    case Core::Materials::mes_volpow:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::VolPow*>(curmat);
      return Teuchos::rcp(new VolPow(params));
    }
    case Core::Materials::mes_volsussmanbathe:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::VolSussmanBathe*>(curmat);
      return Teuchos::rcp(new VolSussmanBathe(params));
    }
    case Core::Materials::mes_viscobranch:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::ViscoBranch*>(curmat);
      return Teuchos::rcp(new ViscoBranch(params));
    }
    case Core::Materials::mes_viscopart:
    {
      auto* params = dynamic_cast<Mat::Elastic::PAR::ViscoPart*>(curmat);
      return Teuchos::rcp(new ViscoPart(params));
    }
    default:
      FOUR_C_THROW("cannot deal with type %d", curmat->Type());
  }
  return Teuchos::null;
}

void Mat::Elastic::Summand::AddShearMod(bool& haveshearmod, double& shearmod) const
{
  FOUR_C_THROW("Mat::Elastic::Summand::AddShearMod: Add Shear Modulus not implemented - do so!");
}

int Mat::Elastic::Summand::UniqueParObjectId() const { return -1; }

void Mat::Elastic::Summand::Pack(Core::Communication::PackBuffer& data) const { return; }

void Mat::Elastic::Summand::Unpack(const std::vector<char>& data) { return; };


// Function which reads in the given fiber value due to the FIBER1 nomenclature
void Mat::Elastic::Summand::ReadFiber(Input::LineDefinition* linedef, const std::string& specifier,
    Core::LinAlg::Matrix<3, 1>& fiber_vector)
{
  std::vector<double> fiber1;
  linedef->extract_double_vector(specifier, fiber1);
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
void Mat::Elastic::Summand::ReadRadAxiCir(
    Input::LineDefinition* linedef, Core::LinAlg::Matrix<3, 3>& locsys)
{
  // read local (cylindrical) cosy-directions at current element
  // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
  Core::LinAlg::Matrix<3, 1> fiber_rad;
  Core::LinAlg::Matrix<3, 1> fiber_axi;
  Core::LinAlg::Matrix<3, 1> fiber_cir;

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

void Mat::Elastic::Summand::evaluate_first_derivatives_aniso(Core::LinAlg::Matrix<2, 1>& dPI_aniso,
    Core::LinAlg::Matrix<3, 3> const& rcg, int gp, int eleGID)
{
  bool isoprinc, isomod, anisoprinc, anisomod, viscogeneral;
  SpecifyFormulation(isoprinc, isomod, anisoprinc, anisomod, viscogeneral);
  if (anisoprinc or anisomod)
  {
    FOUR_C_THROW(
        "This anisotropic material does not support the first derivative of the free-energy "
        "function with respect to the anisotropic invariants. You need to implement them.");
  }
}

void Mat::Elastic::Summand::evaluate_second_derivatives_aniso(
    Core::LinAlg::Matrix<3, 1>& ddPII_aniso, Core::LinAlg::Matrix<3, 3> const& rcg, int gp,
    int eleGID)
{
  bool isoprinc, isomod, anisoprinc, anisomod, viscogeneral;
  SpecifyFormulation(isoprinc, isomod, anisoprinc, anisomod, viscogeneral);
  if (anisoprinc or anisomod)
  {
    FOUR_C_THROW(
        "This anisotropic material does not support the second derivative of the free-energy "
        "function with respect to the anisotropic invariants. You need to implement them.");
  }
}
FOUR_C_NAMESPACE_CLOSE
