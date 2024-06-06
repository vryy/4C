/*----------------------------------------------------------------------*/
/*! \file
 \brief scatra material for transport within multiphase porous medium

   \level 3

 *----------------------------------------------------------------------*/



#include "4C_mat_scatra_multiporo.hpp"

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ScatraMatMultiPoroFluid::ScatraMatMultiPoroFluid(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : ScatraMat(matdata),
      phaseID_(matdata->Get<int>("PHASEID")),
      delta_(matdata->Get<double>("DELTA")),
      min_sat_(matdata->Get<double>("MIN_SAT")),
      relative_mobility_funct_id_(matdata->Get<int>("RELATIVE_MOBILITY_FUNCTION_ID"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMatMultiPoroFluid::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMatMultiPoroFluid(this));
}


Mat::ScatraMatMultiPoroFluidType Mat::ScatraMatMultiPoroFluidType::instance_;

Core::Communication::ParObject* Mat::ScatraMatMultiPoroFluidType::Create(
    const std::vector<char>& data)
{
  Mat::ScatraMatMultiPoroFluid* scatra_mat = new Mat::ScatraMatMultiPoroFluid();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMatMultiPoroFluid::ScatraMatMultiPoroFluid() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMatMultiPoroFluid::ScatraMatMultiPoroFluid(Mat::PAR::ScatraMatMultiPoroFluid* params)
    : ScatraMat(params), params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroFluid::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // add base class material
  ScatraMat::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroFluid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::ScatraMatMultiPoroFluid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraMat::Unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ScatraMatMultiPoroVolFrac::ScatraMatMultiPoroVolFrac(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : ScatraMat(matdata),
      phaseID_(matdata->Get<int>("PHASEID")),
      delta_(matdata->Get<double>("DELTA")),
      relative_mobility_funct_id_(matdata->Get<int>("RELATIVE_MOBILITY_FUNCTION_ID"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMatMultiPoroVolFrac::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMatMultiPoroVolFrac(this));
}


Mat::ScatraMatMultiPoroVolFracType Mat::ScatraMatMultiPoroVolFracType::instance_;

Core::Communication::ParObject* Mat::ScatraMatMultiPoroVolFracType::Create(
    const std::vector<char>& data)
{
  Mat::ScatraMatMultiPoroVolFrac* scatra_mat = new Mat::ScatraMatMultiPoroVolFrac();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMatMultiPoroVolFrac::ScatraMatMultiPoroVolFrac() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMatMultiPoroVolFrac::ScatraMatMultiPoroVolFrac(
    Mat::PAR::ScatraMatMultiPoroVolFrac* params)
    : ScatraMat(params), params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroVolFrac::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // add base class material
  ScatraMat::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroVolFrac::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::ScatraMatMultiPoroVolFrac*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraMat::Unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Mat::PAR::ScatraMatMultiPoroSolid::ScatraMatMultiPoroSolid(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : ScatraMat(matdata), delta_(matdata->Get<double>("DELTA"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMatMultiPoroSolid::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMatMultiPoroSolid(this));
}

Mat::ScatraMatMultiPoroSolidType Mat::ScatraMatMultiPoroSolidType::instance_;

Core::Communication::ParObject* Mat::ScatraMatMultiPoroSolidType::Create(
    const std::vector<char>& data)
{
  Mat::ScatraMatMultiPoroSolid* scatra_mat = new Mat::ScatraMatMultiPoroSolid();
  scatra_mat->Unpack(data);
  return scatra_mat;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMatMultiPoroSolid::ScatraMatMultiPoroSolid() : params_(nullptr) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMatMultiPoroSolid::ScatraMatMultiPoroSolid(Mat::PAR::ScatraMatMultiPoroSolid* params)
    : ScatraMat(params), params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroSolid::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // add base class material
  ScatraMat::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroSolid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::ScatraMatMultiPoroSolid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraMat::Unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Mat::PAR::ScatraMatMultiPoroTemperature::ScatraMatMultiPoroTemperature(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : ScatraMat(matdata),
      numfluidphases_(matdata->Get<int>("NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE")),
      numvolfrac_(matdata->Get<int>("NUMVOLFRAC")),
      cp_fluid_((matdata->Get<std::vector<double>>("CP_FLUID"))),
      cp_volfrac_((matdata->Get<std::vector<double>>("CP_VOLFRAC"))),
      cp_solid_(matdata->Get<double>("CP_SOLID")),
      kappa_fluid_((matdata->Get<std::vector<double>>("KAPPA_FLUID"))),
      kappa_volfrac_((matdata->Get<std::vector<double>>("KAPPA_VOLFRAC"))),
      kappa_solid_(matdata->Get<double>("KAPPA_SOLID"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMatMultiPoroTemperature::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMatMultiPoroTemperature(this));
}

Mat::ScatraMatMultiPoroTemperatureType Mat::ScatraMatMultiPoroTemperatureType::instance_;

Core::Communication::ParObject* Mat::ScatraMatMultiPoroTemperatureType::Create(
    const std::vector<char>& data)
{
  Mat::ScatraMatMultiPoroTemperature* scatra_mat = new Mat::ScatraMatMultiPoroTemperature();
  scatra_mat->Unpack(data);
  return scatra_mat;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMatMultiPoroTemperature::ScatraMatMultiPoroTemperature() : params_(nullptr) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMatMultiPoroTemperature::ScatraMatMultiPoroTemperature(
    Mat::PAR::ScatraMatMultiPoroTemperature* params)
    : ScatraMat(params), params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroTemperature::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // add base class material
  ScatraMat::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroTemperature::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::ScatraMatMultiPoroTemperature*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraMat::Unpack(basedata);
}

FOUR_C_NAMESPACE_CLOSE
