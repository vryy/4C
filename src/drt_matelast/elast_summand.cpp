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
#ifdef CCADISCRET

#include "elast_couplogneohooke.H"
#include "elast_coupneohooke.H"
#include "elast_coupblatzko.H"
#include "elast_coupmooneyrivlin.H"
#include "elast_holzapfel_cardiac.H"
#include "elast_isoneohooke.H"
#include "elast_isoyeoh.H"
#include "elast_isoquad.H"
#include "elast_isocub.H"
#include "elast_iso1pow.H"
#include "elast_iso2pow.H"
#include "elast_coup1pow.H"
#include "elast_coup2pow.H"
#include "elast_isoexpo.H"
#include "elast_isomooneyrivlin.H"
#include "elast_volsussmanbathe.H"
#include "elast_volpenalty.H"
#include "elast_vologden.H"
#include "elast_coupanisoexpotwo.H"
#include "elast_coupanisoneohooketwo.H"
#include "elast_coupvarga.H"
#include "elast_isovarga.H"
#include "elast_isovolHUdependentneohooke.H"
#include "elast_isovolaaagasser.H"
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
  case INPAR::MAT::mes_holzapfel_cardiac:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::Holzapfel_Cardiac(curmat));
    MAT::ELASTIC::PAR::Holzapfel_Cardiac* params = static_cast<MAT::ELASTIC::PAR::Holzapfel_Cardiac*>(curmat->Parameter());
    return Teuchos::rcp(new Holzapfel_Cardiac(params));
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
  case INPAR::MAT::mes_isoexpo:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::IsoExpo(curmat));
    MAT::ELASTIC::PAR::IsoExpo* params = static_cast<MAT::ELASTIC::PAR::IsoExpo*>(curmat->Parameter());
    return Teuchos::rcp(new IsoExpo(params));
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
  case INPAR::MAT::mes_coupanisoexpotwo:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoExpoTwo(curmat));
    MAT::ELASTIC::PAR::CoupAnisoExpoTwo* params = static_cast<MAT::ELASTIC::PAR::CoupAnisoExpoTwo*>(curmat->Parameter());
    return Teuchos::rcp(new CoupAnisoExpoTwo(params));
  }
  case INPAR::MAT::mes_coupanisoneohooketwo:
  {
    if (curmat->Parameter() == NULL)
      curmat->SetParameter(new MAT::ELASTIC::PAR::CoupAnisoNeoHookeTwo(curmat));
    MAT::ELASTIC::PAR::CoupAnisoNeoHookeTwo* params = static_cast<MAT::ELASTIC::PAR::CoupAnisoNeoHookeTwo*>(curmat->Parameter());
    return Teuchos::rcp(new CoupAnisoNeoHookeTwo(params));
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
void MAT::ELASTIC::Summand::Unpack(const vector<char>& data)
{
  return;
};

#endif
