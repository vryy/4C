/*----------------------------------------------------------------------*/
/*! \file
\brief Linear law (pressure-dependent) for the density and the viscosity

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_fluid_linear_density_viscosity.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::LinearDensityViscosity::LinearDensityViscosity(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      refdensity_(matdata.parameters.get<double>("REFDENSITY")),
      refviscosity_(matdata.parameters.get<double>("REFVISCOSITY")),
      refpressure_(matdata.parameters.get<double>("REFPRESSURE")),
      coeffdensity_(matdata.parameters.get<double>("COEFFDENSITY")),
      coeffviscosity_(matdata.parameters.get<double>("COEFFVISCOSITY")),
      gamma_(matdata.parameters.get<double>("GAMMA"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::LinearDensityViscosity::create_material()
{
  return Teuchos::rcp(new Mat::LinearDensityViscosity(this));
}


Mat::LinearDensityViscosityType Mat::LinearDensityViscosityType::instance_;


Core::Communication::ParObject* Mat::LinearDensityViscosityType::Create(
    const std::vector<char>& data)
{
  Mat::LinearDensityViscosity* fluid = new Mat::LinearDensityViscosity();
  fluid->Unpack(data);
  return fluid;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::LinearDensityViscosity::LinearDensityViscosity() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::LinearDensityViscosity::LinearDensityViscosity(Mat::PAR::LinearDensityViscosity* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::LinearDensityViscosity::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::LinearDensityViscosity::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::LinearDensityViscosity*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::LinearDensityViscosity::ComputeDensity(const double press) const
{
  const double density = RefDensity() * (1.0 + CoeffDensity() * (press - RefPressure()));

  return density;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::LinearDensityViscosity::ComputeViscosity(const double press) const
{
  const double viscosity = RefViscosity() * (1.0 + CoeffViscosity() * (press - RefPressure()));

  return viscosity;
}

FOUR_C_NAMESPACE_CLOSE
