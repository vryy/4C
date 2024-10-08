/*----------------------------------------------------------------------*/
/*! \file
\brief Linear law (pressure-dependent) for the density and the viscosity

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_fluid_linear_density_viscosity.hpp"

#include "4C_comm_pack_helpers.hpp"
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
  return Teuchos::make_rcp<Mat::LinearDensityViscosity>(this);
}


Mat::LinearDensityViscosityType Mat::LinearDensityViscosityType::instance_;


Core::Communication::ParObject* Mat::LinearDensityViscosityType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::LinearDensityViscosity* fluid = new Mat::LinearDensityViscosity();
  fluid->unpack(buffer);
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
void Mat::LinearDensityViscosity::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::LinearDensityViscosity::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::LinearDensityViscosity*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::LinearDensityViscosity::compute_density(const double press) const
{
  const double density = ref_density() * (1.0 + coeff_density() * (press - ref_pressure()));

  return density;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::LinearDensityViscosity::compute_viscosity(const double press) const
{
  const double viscosity = ref_viscosity() * (1.0 + coeff_viscosity() * (press - ref_pressure()));

  return viscosity;
}

FOUR_C_NAMESPACE_CLOSE
