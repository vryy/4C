/*----------------------------------------------------------------------*/
/*!
 \brief scatra material for transport within multiphase porous medium

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/



#include <vector>
#include "scatra_mat_multiporo.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraMatMultiPoroFluid::ScatraMatMultiPoroFluid(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ScatraMat(matdata),
      phaseID_(matdata->GetInt("PHASEID")),
      delta_(matdata->GetDouble("DELTA")),
      min_sat_(matdata->GetDouble("MIN_SAT"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMatMultiPoroFluid::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMatMultiPoroFluid(this));
}


MAT::ScatraMatMultiPoroFluidType MAT::ScatraMatMultiPoroFluidType::instance_;

DRT::ParObject* MAT::ScatraMatMultiPoroFluidType::Create(const std::vector<char>& data)
{
  MAT::ScatraMatMultiPoroFluid* scatra_mat = new MAT::ScatraMatMultiPoroFluid();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatMultiPoroFluid::ScatraMatMultiPoroFluid() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatMultiPoroFluid::ScatraMatMultiPoroFluid(MAT::PAR::ScatraMatMultiPoroFluid* params)
    : ScatraMat(params), params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatMultiPoroFluid::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // add base class material
  ScatraMat::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatMultiPoroFluid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ScatraMatMultiPoroFluid*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraMat::Unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraMatMultiPoroVolFrac::ScatraMatMultiPoroVolFrac(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : ScatraMat(matdata), phaseID_(matdata->GetInt("PHASEID")), delta_(matdata->GetDouble("DELTA"))

{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMatMultiPoroVolFrac::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMatMultiPoroVolFrac(this));
}


MAT::ScatraMatMultiPoroVolFracType MAT::ScatraMatMultiPoroVolFracType::instance_;

DRT::ParObject* MAT::ScatraMatMultiPoroVolFracType::Create(const std::vector<char>& data)
{
  MAT::ScatraMatMultiPoroVolFrac* scatra_mat = new MAT::ScatraMatMultiPoroVolFrac();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatMultiPoroVolFrac::ScatraMatMultiPoroVolFrac() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatMultiPoroVolFrac::ScatraMatMultiPoroVolFrac(
    MAT::PAR::ScatraMatMultiPoroVolFrac* params)
    : ScatraMat(params), params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatMultiPoroVolFrac::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // add base class material
  ScatraMat::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatMultiPoroVolFrac::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ScatraMatMultiPoroVolFrac*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraMat::Unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

MAT::PAR::ScatraMatMultiPoroSolid::ScatraMatMultiPoroSolid(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ScatraMat(matdata), delta_(matdata->GetDouble("DELTA"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMatMultiPoroSolid::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMatMultiPoroSolid(this));
}

MAT::ScatraMatMultiPoroSolidType MAT::ScatraMatMultiPoroSolidType::instance_;

DRT::ParObject* MAT::ScatraMatMultiPoroSolidType::Create(const std::vector<char>& data)
{
  MAT::ScatraMatMultiPoroSolid* scatra_mat = new MAT::ScatraMatMultiPoroSolid();
  scatra_mat->Unpack(data);
  return scatra_mat;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatMultiPoroSolid::ScatraMatMultiPoroSolid() : params_(NULL) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatMultiPoroSolid::ScatraMatMultiPoroSolid(MAT::PAR::ScatraMatMultiPoroSolid* params)
    : ScatraMat(params), params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatMultiPoroSolid::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // add base class material
  ScatraMat::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatMultiPoroSolid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ScatraMatMultiPoroSolid*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraMat::Unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

MAT::PAR::ScatraMatMultiPoroTemperature::ScatraMatMultiPoroTemperature(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : ScatraMat(matdata),
      numfluidphases_(matdata->GetInt("NUMFLUIDPHASES")),
      numvolfrac_(matdata->GetInt("NUMVOLFRAC")),
      cp_fluid_(*(matdata->Get<std::vector<double>>("CP_FLUID"))),
      cp_volfrac_(*(matdata->Get<std::vector<double>>("CP_VOLFRAC"))),
      cp_solid_(matdata->GetDouble("CP_SOLID")),
      kappa_fluid_(*(matdata->Get<std::vector<double>>("KAPPA_FLUID"))),
      kappa_volfrac_(*(matdata->Get<std::vector<double>>("KAPPA_VOLFRAC"))),
      kappa_solid_(matdata->GetDouble("KAPPA_SOLID"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMatMultiPoroTemperature::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMatMultiPoroTemperature(this));
}

MAT::ScatraMatMultiPoroTemperatureType MAT::ScatraMatMultiPoroTemperatureType::instance_;

DRT::ParObject* MAT::ScatraMatMultiPoroTemperatureType::Create(const std::vector<char>& data)
{
  MAT::ScatraMatMultiPoroTemperature* scatra_mat = new MAT::ScatraMatMultiPoroTemperature();
  scatra_mat->Unpack(data);
  return scatra_mat;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatMultiPoroTemperature::ScatraMatMultiPoroTemperature() : params_(NULL) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatMultiPoroTemperature::ScatraMatMultiPoroTemperature(
    MAT::PAR::ScatraMatMultiPoroTemperature* params)
    : ScatraMat(params), params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatMultiPoroTemperature::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // add base class material
  ScatraMat::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatMultiPoroTemperature::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ScatraMatMultiPoroTemperature*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraMat::Unpack(basedata);
}
