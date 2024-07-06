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
    const Core::Mat::PAR::Parameter::Data& matdata)
    : ScatraMat(matdata),
      phaseID_(matdata.parameters.get<int>("PHASEID")),
      delta_(matdata.parameters.get<double>("DELTA")),
      min_sat_(matdata.parameters.get<double>("MIN_SAT")),
      relative_mobility_funct_id_(matdata.parameters.get<int>("RELATIVE_MOBILITY_FUNCTION_ID"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMatMultiPoroFluid::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMatMultiPoroFluid(this));
}


Mat::ScatraMatMultiPoroFluidType Mat::ScatraMatMultiPoroFluidType::instance_;

Core::Communication::ParObject* Mat::ScatraMatMultiPoroFluidType::create(
    const std::vector<char>& data)
{
  Mat::ScatraMatMultiPoroFluid* scatra_mat = new Mat::ScatraMatMultiPoroFluid();
  scatra_mat->unpack(data);
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
void Mat::ScatraMatMultiPoroFluid::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // add base class material
  ScatraMat::pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroFluid::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::ScatraMatMultiPoroFluid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  // extract base class material
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  ScatraMat::unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ScatraMatMultiPoroVolFrac::ScatraMatMultiPoroVolFrac(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : ScatraMat(matdata),
      phaseID_(matdata.parameters.get<int>("PHASEID")),
      delta_(matdata.parameters.get<double>("DELTA")),
      relative_mobility_funct_id_(matdata.parameters.get<int>("RELATIVE_MOBILITY_FUNCTION_ID"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMatMultiPoroVolFrac::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMatMultiPoroVolFrac(this));
}


Mat::ScatraMatMultiPoroVolFracType Mat::ScatraMatMultiPoroVolFracType::instance_;

Core::Communication::ParObject* Mat::ScatraMatMultiPoroVolFracType::create(
    const std::vector<char>& data)
{
  Mat::ScatraMatMultiPoroVolFrac* scatra_mat = new Mat::ScatraMatMultiPoroVolFrac();
  scatra_mat->unpack(data);
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
void Mat::ScatraMatMultiPoroVolFrac::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // add base class material
  ScatraMat::pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroVolFrac::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::ScatraMatMultiPoroVolFrac*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  // extract base class material
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  ScatraMat::unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Mat::PAR::ScatraMatMultiPoroSolid::ScatraMatMultiPoroSolid(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : ScatraMat(matdata), delta_(matdata.parameters.get<double>("DELTA"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMatMultiPoroSolid::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMatMultiPoroSolid(this));
}

Mat::ScatraMatMultiPoroSolidType Mat::ScatraMatMultiPoroSolidType::instance_;

Core::Communication::ParObject* Mat::ScatraMatMultiPoroSolidType::create(
    const std::vector<char>& data)
{
  Mat::ScatraMatMultiPoroSolid* scatra_mat = new Mat::ScatraMatMultiPoroSolid();
  scatra_mat->unpack(data);
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
void Mat::ScatraMatMultiPoroSolid::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // add base class material
  ScatraMat::pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroSolid::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::ScatraMatMultiPoroSolid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  // extract base class material
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  ScatraMat::unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

Mat::PAR::ScatraMatMultiPoroTemperature::ScatraMatMultiPoroTemperature(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : ScatraMat(matdata),
      numfluidphases_(matdata.parameters.get<int>("NUMFLUIDPHASES_IN_MULTIPHASEPORESPACE")),
      numvolfrac_(matdata.parameters.get<int>("NUMVOLFRAC")),
      cp_fluid_((matdata.parameters.get<std::vector<double>>("CP_FLUID"))),
      cp_volfrac_((matdata.parameters.get<std::vector<double>>("CP_VOLFRAC"))),
      cp_solid_(matdata.parameters.get<double>("CP_SOLID")),
      kappa_fluid_((matdata.parameters.get<std::vector<double>>("KAPPA_FLUID"))),
      kappa_volfrac_((matdata.parameters.get<std::vector<double>>("KAPPA_VOLFRAC"))),
      kappa_solid_(matdata.parameters.get<double>("KAPPA_SOLID"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMatMultiPoroTemperature::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMatMultiPoroTemperature(this));
}

Mat::ScatraMatMultiPoroTemperatureType Mat::ScatraMatMultiPoroTemperatureType::instance_;

Core::Communication::ParObject* Mat::ScatraMatMultiPoroTemperatureType::create(
    const std::vector<char>& data)
{
  Mat::ScatraMatMultiPoroTemperature* scatra_mat = new Mat::ScatraMatMultiPoroTemperature();
  scatra_mat->unpack(data);
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
void Mat::ScatraMatMultiPoroTemperature::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // add base class material
  ScatraMat::pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatMultiPoroTemperature::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::ScatraMatMultiPoroTemperature*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  // extract base class material
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  ScatraMat::unpack(basedata);
}

FOUR_C_NAMESPACE_CLOSE
