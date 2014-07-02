/*----------------------------------------------------------------------*/
/*!
\file elast_summand.cpp

\brief Interface class for complex materials at Gauss points

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/


#include "elast_couplogneohooke.H"
#include "elast_coupexppol.H"
#include "elast_coupneohooke.H"
#include "elast_coupblatzko.H"
#include "elast_coupmooneyrivlin.H"
#include "elast_isoneohooke.H"
#include "elast_isoyeoh.H"
#include "elast_isoquad.H"
#include "elast_isocub.H"
#include "elast_iso1pow.H"
#include "elast_iso2pow.H"
#include "elast_coup1pow.H"
#include "elast_coup2pow.H"
#include "elast_isoexpopow.H"
#include "elast_isomooneyrivlin.H"
#include "elast_volsussmanbathe.H"
#include "elast_volpenalty.H"
#include "elast_vologden.H"
#include "elast_coupanisoexpo.H"
#include "elast_coupanisoexpotwocoup.H"
#include "elast_coupanisoneohooke.H"
#include "elast_coupanisoneohooke_ActiveStress.H"
#include "elast_coupanisoneohooke_VarProp.H"
#include "elast_coupanisopow.H"
#include "elast_isoanisoexpo.H"
#include "elast_isoanisoexpodispersion.H"
#include "elast_coupvarga.H"
#include "elast_isovarga.H"
#include "elast_isovolHUdependentneohooke.H"
#include "elast_isovolaaagasser.H"
#include "visco_isoratedep.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_linedefinition.H"

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
  case INPAR::MAT::mes_isoquad:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoQuad(curmat));
    MAT::ELASTIC::PAR::IsoQuad* params = static_cast<MAT::ELASTIC::PAR::IsoQuad*>(curmat->Parameter());
    return Teuchos::rcp(new IsoQuad(params));
  }
  case INPAR::MAT::mes_isocub:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoCub(curmat));
    MAT::ELASTIC::PAR::IsoCub* params = static_cast<MAT::ELASTIC::PAR::IsoCub*>(curmat->Parameter());
    return Teuchos::rcp(new IsoCub(params));
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
  case INPAR::MAT::mes_coupanisoneohooke_activestress:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoNeoHooke_ActiveStress(curmat));
    MAT::ELASTIC::PAR::CoupAnisoNeoHooke_ActiveStress* params = static_cast<MAT::ELASTIC::PAR::CoupAnisoNeoHooke_ActiveStress*>(curmat->Parameter());
    return Teuchos::rcp(new CoupAnisoNeoHooke_ActiveStress(params));
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
  case INPAR::MAT::mes_isoanisoexpodispersion:
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
  case INPAR::MAT::mes_isoratedep:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoRateDep(curmat));
    MAT::ELASTIC::PAR::IsoRateDep* params = static_cast<MAT::ELASTIC::PAR::IsoRateDep*>(curmat->Parameter());
    return Teuchos::rcp(new IsoRateDep(params));
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
