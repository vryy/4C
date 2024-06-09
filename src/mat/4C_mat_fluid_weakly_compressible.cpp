/*----------------------------------------------------------------------------*/
/*! \file
\brief Weakly compressible fluid

\level 1

*/
/*----------------------------------------------------------------------------*/

#include "4C_mat_fluid_weakly_compressible.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::WeaklyCompressibleFluid::WeaklyCompressibleFluid(
    Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(matdata->Get<double>("VISCOSITY")),
      refdensity_(matdata->Get<double>("REFDENSITY")),
      refpressure_(matdata->Get<double>("REFPRESSURE")),
      comprcoeff_(matdata->Get<double>("COMPRCOEFF"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::WeaklyCompressibleFluid::create_material()
{
  return Teuchos::rcp(new Mat::WeaklyCompressibleFluid(this));
}


Mat::WeaklyCompressibleFluidType Mat::WeaklyCompressibleFluidType::instance_;


Core::Communication::ParObject* Mat::WeaklyCompressibleFluidType::Create(
    const std::vector<char>& data)
{
  Mat::WeaklyCompressibleFluid* fluid = new Mat::WeaklyCompressibleFluid();
  fluid->Unpack(data);
  return fluid;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::WeaklyCompressibleFluid::WeaklyCompressibleFluid() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::WeaklyCompressibleFluid::WeaklyCompressibleFluid(Mat::PAR::WeaklyCompressibleFluid* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::WeaklyCompressibleFluid::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

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
void Mat::WeaklyCompressibleFluid::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::WeaklyCompressibleFluid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::WeaklyCompressibleFluid::ComputeDensity(const double press) const
{
  const double density = RefDensity() + ComprCoeff() * (press - RefPressure());

  return density;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::WeaklyCompressibleFluid::ComputePressure(const double dens) const
{
  const double pressure = RefPressure() + 1.0 / ComprCoeff() * (dens - RefDensity());

  return pressure;
}

FOUR_C_NAMESPACE_CLOSE
