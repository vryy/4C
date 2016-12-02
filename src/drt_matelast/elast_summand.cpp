/*----------------------------------------------------------------------*/
/*!
\file elast_summand.cpp

\brief
Interface class for materials of (visco)elasthyper toolbox.

\level 1

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/


#include "elast_couplogneohooke.H"
#include "elast_couplogmixneohooke.H"
#include "elast_coupexppol.H"
#include "elast_coupneohooke.H"
#include "elast_coupblatzko.H"
#include "elast_coupmooneyrivlin.H"
#include "elast_isoneohooke.H"
#include "elast_isoyeoh.H"
#include "elast_iso1pow.H"
#include "elast_iso2pow.H"
#include "elast_coup1pow.H"
#include "elast_coup2pow.H"
#include "elast_coup13apow.H"
#include "elast_isoexpopow.H"
#include "elast_isomooneyrivlin.H"
#include "elast_isotestmaterial.H"
#include "elast_coupSaintVenantKirchhoff.H"
#include "elast_coupsimopister.H"
#include "elast_remodelfiber.H"
#include "elast_volsussmanbathe.H"
#include "elast_volpenalty.H"
#include "elast_vologden.H"
#include "elast_volpow.H"
#include "elast_coupanisoexpo.H"
#include "elast_coupanisoexpotwocoup.H"
#include "elast_coupanisoneohooke.H"
#include "elast_coupanisoneohooke_VarProp.H"
#include "elast_coupanisopow.H"
#include "elast_isoanisoexpo.H"
#include "elast_coupvarga.H"
#include "elast_isovarga.H"
#include "elast_isovolHUdependentneohooke.H"
#include "elast_isovolaaagasser.H"
#include "visco_coupmyocard.H"
#include "visco_isoratedep.H"
#include "visco_genmax.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_linedefinition.H"
#include "elast_anisoactivestress_evolution.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<MAT::ELASTIC::Summand> MAT::ELASTIC::Summand::Factory(int matnum)
{
  // for the sake of safety
  if (DRT::Problem::Instance()->Materials() == Teuchos::null)
    dserror("Sorry dude, cannot work out problem instance.");

  // yet another safety check
  if (DRT::Problem::Instance()->Materials()->Num() == 0)
    dserror("Sorry dude, no materials defined.");

  // retrieve problem instance to read from
  const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat = DRT::Problem::Instance(probinst)->Materials()->ById(matnum);

  // check what was read
  //cout << *curmat << endl;

  switch (curmat->Type())
  {
  case INPAR::MAT::mes_couplogneohooke:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupLogNeoHooke(curmat));
    MAT::ELASTIC::PAR::CoupLogNeoHooke* params = static_cast<MAT::ELASTIC::PAR::CoupLogNeoHooke*>(curmat->Parameter());
    return Teuchos::rcp(new CoupLogNeoHooke(params));
  }
  case INPAR::MAT::mes_coupSVK:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupSVK(curmat));
    MAT::ELASTIC::PAR::CoupSVK* params = static_cast<MAT::ELASTIC::PAR::CoupSVK*>(curmat->Parameter());
    return Teuchos::rcp(new CoupSVK(params));
  }
  case INPAR::MAT::mes_coupsimopister:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupSimoPister(curmat));
    MAT::ELASTIC::PAR::CoupSimoPister* params = static_cast<MAT::ELASTIC::PAR::CoupSimoPister*>(curmat->Parameter());
    return Teuchos::rcp(new CoupSimoPister(params));
  }
  case INPAR::MAT::mes_couplogmixneohooke:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupLogMixNeoHooke(curmat));
    MAT::ELASTIC::PAR::CoupLogMixNeoHooke* params = static_cast<MAT::ELASTIC::PAR::CoupLogMixNeoHooke*>(curmat->Parameter());
    return Teuchos::rcp(new CoupLogMixNeoHooke(params));
  }
  case INPAR::MAT::mes_coupexppol:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupExpPol(curmat));
    MAT::ELASTIC::PAR::CoupExpPol* params = static_cast<MAT::ELASTIC::PAR::CoupExpPol*>(curmat->Parameter());
    return Teuchos::rcp(new CoupExpPol(params));
  }
  case INPAR::MAT::mes_coupneohooke:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupNeoHooke(curmat));
    MAT::ELASTIC::PAR::CoupNeoHooke* params = static_cast<MAT::ELASTIC::PAR::CoupNeoHooke*>(curmat->Parameter());
    return Teuchos::rcp(new CoupNeoHooke(params));
  }
  case INPAR::MAT::mes_coupblatzko:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupBlatzKo(curmat));
    MAT::ELASTIC::PAR::CoupBlatzKo* params = static_cast<MAT::ELASTIC::PAR::CoupBlatzKo*>(curmat->Parameter());
    return Teuchos::rcp(new CoupBlatzKo(params));
  }
  case INPAR::MAT::mes_coupmooneyrivlin:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupMooneyRivlin(curmat));
    MAT::ELASTIC::PAR::CoupMooneyRivlin* params = static_cast<MAT::ELASTIC::PAR::CoupMooneyRivlin*>(curmat->Parameter());
    return Teuchos::rcp(new CoupMooneyRivlin(params));
  }
  case INPAR::MAT::mes_isoneohooke:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoNeoHooke(curmat));
    MAT::ELASTIC::PAR::IsoNeoHooke* params = static_cast<MAT::ELASTIC::PAR::IsoNeoHooke*>(curmat->Parameter());
    return Teuchos::rcp(new IsoNeoHooke(params));
  }
  case INPAR::MAT::mes_isoyeoh:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoYeoh(curmat));
    MAT::ELASTIC::PAR::IsoYeoh* params = static_cast<MAT::ELASTIC::PAR::IsoYeoh*>(curmat->Parameter());
    return Teuchos::rcp(new IsoYeoh(params));
  }
  case INPAR::MAT::mes_iso1pow:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::Iso1Pow(curmat));
    MAT::ELASTIC::PAR::Iso1Pow* params = static_cast<MAT::ELASTIC::PAR::Iso1Pow*>(curmat->Parameter());
    return Teuchos::rcp(new Iso1Pow(params));
  }
  case INPAR::MAT::mes_iso2pow:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::Iso2Pow(curmat));
    MAT::ELASTIC::PAR::Iso2Pow* params = static_cast<MAT::ELASTIC::PAR::Iso2Pow*>(curmat->Parameter());
    return Teuchos::rcp(new Iso2Pow(params));
  }
  case INPAR::MAT::mes_coup1pow:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::Coup1Pow(curmat));
    MAT::ELASTIC::PAR::Coup1Pow* params = static_cast<MAT::ELASTIC::PAR::Coup1Pow*>(curmat->Parameter());
    return Teuchos::rcp(new Coup1Pow(params));
  }
  case INPAR::MAT::mes_coup2pow:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::Coup2Pow(curmat));
    MAT::ELASTIC::PAR::Coup2Pow* params = static_cast<MAT::ELASTIC::PAR::Coup2Pow*>(curmat->Parameter());
    return Teuchos::rcp(new Coup2Pow(params));
  }
  case INPAR::MAT::mes_coup13apow:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::Coup13aPow(curmat));
    MAT::ELASTIC::PAR::Coup13aPow* params = static_cast<MAT::ELASTIC::PAR::Coup13aPow*>(curmat->Parameter());
    return Teuchos::rcp(new Coup13aPow(params));
  }
  case INPAR::MAT::mes_isoexpopow:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoExpoPow(curmat));
    MAT::ELASTIC::PAR::IsoExpoPow* params = static_cast<MAT::ELASTIC::PAR::IsoExpoPow*>(curmat->Parameter());
    return Teuchos::rcp(new IsoExpoPow(params));
  }
  case INPAR::MAT::mes_isomooneyrivlin:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoMooneyRivlin(curmat));
    MAT::ELASTIC::PAR::IsoMooneyRivlin* params = static_cast<MAT::ELASTIC::PAR::IsoMooneyRivlin*>(curmat->Parameter());
    return Teuchos::rcp(new IsoMooneyRivlin(params));
  }
  case INPAR::MAT::mes_isotestmaterial:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoTestMaterial(curmat));
    MAT::ELASTIC::PAR::IsoTestMaterial* params = static_cast<MAT::ELASTIC::PAR::IsoTestMaterial*>(curmat->Parameter());
    return Teuchos::rcp(new IsoTestMaterial(params));
  }
  case INPAR::MAT::mes_isovolHUdependentneohooke:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke(curmat));
    MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke* params = static_cast<MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke*>(curmat->Parameter());
    return Teuchos::rcp(new IsoVolHUDependentNeoHooke(params));
  }
  case INPAR::MAT::mes_isovolaaagasser:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoVolAAAGasser(curmat));
    MAT::ELASTIC::PAR::IsoVolAAAGasser* params = static_cast<MAT::ELASTIC::PAR::IsoVolAAAGasser*>(curmat->Parameter());
    return Teuchos::rcp(new IsoVolAAAGasser(params));
  }
  case INPAR::MAT::mes_volsussmanbathe:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::VolSussmanBathe(curmat));
    MAT::ELASTIC::PAR::VolSussmanBathe* params = static_cast<MAT::ELASTIC::PAR::VolSussmanBathe*>(curmat->Parameter());
    return Teuchos::rcp(new VolSussmanBathe(params));
  }
  case INPAR::MAT::mes_remodelfiber:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::RemodelFiber(curmat));
    MAT::ELASTIC::PAR::RemodelFiber* params = static_cast<MAT::ELASTIC::PAR::RemodelFiber*>(curmat->Parameter());
    return Teuchos::rcp(new RemodelFiber(params));
  }
  case INPAR::MAT::mes_volpenalty:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::VolPenalty(curmat));
    MAT::ELASTIC::PAR::VolPenalty* params = static_cast<MAT::ELASTIC::PAR::VolPenalty*>(curmat->Parameter());
    return Teuchos::rcp(new VolPenalty(params));
  }
  case INPAR::MAT::mes_vologden:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::VolOgden(curmat));
    MAT::ELASTIC::PAR::VolOgden* params = static_cast<MAT::ELASTIC::PAR::VolOgden*>(curmat->Parameter());
    return Teuchos::rcp(new VolOgden(params));
  }
  case INPAR::MAT::mes_volpow:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::VolPow(curmat));
    MAT::ELASTIC::PAR::VolPow* params = static_cast<MAT::ELASTIC::PAR::VolPow*>(curmat->Parameter());
    return Teuchos::rcp(new VolPow(params));
  }
  case INPAR::MAT::mes_anisoactivestress_evolution:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::AnisoActiveStress_Evolution(curmat));
    MAT::ELASTIC::PAR::AnisoActiveStress_Evolution* params = static_cast<MAT::ELASTIC::PAR::AnisoActiveStress_Evolution*>(curmat->Parameter());
    return Teuchos::rcp(new AnisoActiveStress_Evolution(params));
  }
  case INPAR::MAT::mes_coupanisoexpoactive:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoExpoActive(curmat));
    MAT::ELASTIC::PAR::CoupAnisoExpoActive* params = static_cast<MAT::ELASTIC::PAR::CoupAnisoExpoActive*>(curmat->Parameter());
    return Teuchos::rcp(new CoupAnisoExpoActive(params));
  }
  case INPAR::MAT::mes_coupanisoexpo:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoExpo(curmat));
    MAT::ELASTIC::PAR::CoupAnisoExpo* params = static_cast<MAT::ELASTIC::PAR::CoupAnisoExpo*>(curmat->Parameter());
    return Teuchos::rcp(new CoupAnisoExpo(params));
  }
  case INPAR::MAT::mes_coupanisopow:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoPow(curmat));
    MAT::ELASTIC::PAR::CoupAnisoPow* params = static_cast<MAT::ELASTIC::PAR::CoupAnisoPow*>(curmat->Parameter());
    return Teuchos::rcp(new CoupAnisoPow(params));
  }
  case INPAR::MAT::mes_coupanisoexpotwocoup:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup(curmat));
    MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup* params = static_cast<MAT::ELASTIC::PAR::CoupAnisoExpoTwoCoup*>(curmat->Parameter());
    return Teuchos::rcp(new CoupAnisoExpoTwoCoup(params));
  }
  case INPAR::MAT::mes_coupanisoneohooke:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoNeoHooke(curmat));
    MAT::ELASTIC::PAR::CoupAnisoNeoHooke* params = static_cast<MAT::ELASTIC::PAR::CoupAnisoNeoHooke*>(curmat->Parameter());
    return Teuchos::rcp(new CoupAnisoNeoHooke(params));
  }
  case INPAR::MAT::mes_coupanisoneohooke_varprop:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp(curmat));
    MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp* params = static_cast<MAT::ELASTIC::PAR::CoupAnisoNeoHooke_VarProp*>(curmat->Parameter());
    return Teuchos::rcp(new CoupAnisoNeoHooke_VarProp(params));
  }
  case INPAR::MAT::mes_isoanisoexpo:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoAnisoExpo(curmat));
    MAT::ELASTIC::PAR::IsoAnisoExpo* params = static_cast<MAT::ELASTIC::PAR::IsoAnisoExpo*>(curmat->Parameter());
    return Teuchos::rcp(new IsoAnisoExpo(params));
  }
  case INPAR::MAT::mes_coupvarga:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupVarga(curmat));
    MAT::ELASTIC::PAR::CoupVarga* params = static_cast<MAT::ELASTIC::PAR::CoupVarga*>(curmat->Parameter());
    return Teuchos::rcp(new CoupVarga(params));
  }
  case INPAR::MAT::mes_isovarga:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoVarga(curmat));
    MAT::ELASTIC::PAR::IsoVarga* params = static_cast<MAT::ELASTIC::PAR::IsoVarga*>(curmat->Parameter());
    return Teuchos::rcp(new IsoVarga(params));
  }
  case INPAR::MAT::mes_coupmyocard:
    {
      if (curmat->Parameter() == NULL)
        curmat->SetParameter(new MAT::ELASTIC::PAR::CoupMyocard(curmat));
      MAT::ELASTIC::PAR::CoupMyocard* params = static_cast<MAT::ELASTIC::PAR::CoupMyocard*>(curmat->Parameter());
      return Teuchos::rcp(new CoupMyocard(params));
    }
  case INPAR::MAT::mes_isoratedep:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoRateDep(curmat));
    MAT::ELASTIC::PAR::IsoRateDep* params = static_cast<MAT::ELASTIC::PAR::IsoRateDep*>(curmat->Parameter());
    return Teuchos::rcp(new IsoRateDep(params));
  }
  case INPAR::MAT::mes_genmax:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::GenMax(curmat));
    MAT::ELASTIC::PAR::GenMax* params = static_cast<MAT::ELASTIC::PAR::GenMax*>(curmat->Parameter());
    return Teuchos::rcp(new GenMax(params));
  }
  default:
    dserror("cannot deal with type %d", curmat->Type());
  }
  return Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Summand::AddShearMod(
  bool& haveshearmod,
  double& shearmod
  ) const
{
  dserror("MAT::ELASTIC::Summand::AddShearMod: Add Shear Modulus not implemented - do so!");
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int MAT::ELASTIC::Summand::UniqueParObjectId() const
{
  return -1;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Summand::Pack(DRT::PackBuffer& data) const
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Summand::Unpack(const std::vector<char>& data)
{
  return;
};

/*----------------------------------------------------------------------*
 *
 * Function for computing the structural tensor for anisotropic materials
 * via a dyadic product of the current fiber direction
 *
 *----------------------------------------------------------------------*/
void MAT::ELASTIC::Summand::SetupStructuralTensor(
    LINALG::Matrix<3,1>  &fiber_vector,
    LINALG::Matrix<6,1>  &structural_tensor
)
{
  for (int i = 0; i < 3; ++i)
    structural_tensor(i) = fiber_vector(i)*fiber_vector(i);
  structural_tensor(3) = fiber_vector(0)*fiber_vector(1);
  structural_tensor(4) = fiber_vector(1)*fiber_vector(2);
  structural_tensor(5) = fiber_vector(0)*fiber_vector(2);
}

/*----------------------------------------------------------------------*
 * Function which reads in the given fiber value due to the
 * FIBER1 nomenclature
 *----------------------------------------------------------------------*/
void MAT::ELASTIC::Summand::ReadFiber(
    DRT::INPUT::LineDefinition* linedef,
    std::string specifier,
    LINALG::Matrix<3,1> &fiber_vector
)
{
  std::vector<double> fiber1;
  linedef->ExtractDoubleVector(specifier,fiber1);
  double f1norm=0.;
  //normalization
  for (int i = 0; i < 3; ++i)
  {
    f1norm += fiber1[i]*fiber1[i];
  }
  f1norm = sqrt(f1norm);

  // fill final fiber vector
  for (int i = 0; i < 3; ++i)
    fiber_vector(i) = fiber1[i]/f1norm;
}

/*----------------------------------------------------------------------*
 * Function which reads in the given fiber value due to the
 * CIR-AXI-RAD nomenclature
 *----------------------------------------------------------------------*/
void MAT::ELASTIC::Summand::ReadRadAxiCir(
    DRT::INPUT::LineDefinition* linedef,
    LINALG::Matrix<3,3>& locsys
)
{
    // read local (cylindrical) cosy-directions at current element
    // basis is local cosy with third vec e3 = circumferential dir and e2 = axial dir
    LINALG::Matrix<3,1> fiber_rad;
    LINALG::Matrix<3,1> fiber_axi;
    LINALG::Matrix<3,1> fiber_cir;

    ReadFiber(linedef,"RAD",fiber_rad);
    ReadFiber(linedef,"AXI",fiber_axi);
    ReadFiber(linedef,"CIR",fiber_cir);

    for (int i=0; i<3; ++i)
    {
      locsys(i,0) = fiber_rad(i);
      locsys(i,1) = fiber_axi(i);
      locsys(i,2) = fiber_cir(i);
    }
}
