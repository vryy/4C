/*----------------------------------------------------------------------*/
/*! \file
\brief temperature-dependent gas according to Sutherland law

\level 2

*/
/*----------------------------------------------------------------------*/


#include "4C_mat_sutherland.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::Sutherland::Sutherland(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      refvisc_(matdata.parameters.get<double>("REFVISC")),
      reftemp_(matdata.parameters.get<double>("REFTEMP")),
      suthtemp_(matdata.parameters.get<double>("SUTHTEMP")),
      shc_(matdata.parameters.get<double>("SHC")),
      pranum_(matdata.parameters.get<double>("PRANUM")),
      thermpress_(matdata.parameters.get<double>("THERMPRESS")),
      gasconst_(matdata.parameters.get<double>("GASCON"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::Sutherland::create_material()
{
  return Teuchos::rcp(new Mat::Sutherland(this));
}


Mat::SutherlandType Mat::SutherlandType::instance_;


Core::Communication::ParObject* Mat::SutherlandType::create(const std::vector<char>& data)
{
  Mat::Sutherland* sutherland = new Mat::Sutherland();
  sutherland->unpack(data);
  return sutherland;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::Sutherland::Sutherland() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::Sutherland::Sutherland(Mat::PAR::Sutherland* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::Sutherland::pack(Core::Communication::PackBuffer& data) const
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
void Mat::Sutherland::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid and recover params_
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
        params_ = static_cast<Mat::PAR::Sutherland*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::Sutherland::compute_viscosity(const double temp) const
{
  // previous implementation using "pow"-function appears to be extremely
  // time-consuming sometimes, at least on the computing cluster
  // const double visc =
  // std::pow((temp/RefTemp()),1.5)*((RefTemp()+SuthTemp())/(temp+SuthTemp()))*RefVisc();
  const double visc = sqrt((temp / ref_temp()) * (temp / ref_temp()) * (temp / ref_temp())) *
                      ((ref_temp() + suth_temp()) / (temp + suth_temp())) * ref_visc();

  return visc;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::Sutherland::compute_diffusivity(const double temp) const
{
  // previous implementation using "pow"-function appears to be extremely
  // time-consuming sometimes, at least on the computing cluster
  // const double diffus =
  // std::pow((temp/RefTemp()),1.5)*((RefTemp()+SuthTemp())/(temp+SuthTemp()))*RefVisc()/PraNum();
  const double diffus = sqrt((temp / ref_temp()) * (temp / ref_temp()) * (temp / ref_temp())) *
                        ((ref_temp() + suth_temp()) / (temp + suth_temp())) * ref_visc() /
                        pra_num();

  return diffus;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::Sutherland::compute_density(const double temp, const double thermpress) const
{
  const double density = thermpress / (gas_const() * temp);

  return density;
}

FOUR_C_NAMESPACE_CLOSE
