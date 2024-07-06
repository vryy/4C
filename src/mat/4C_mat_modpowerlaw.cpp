/*----------------------------------------------------------------------*/
/*! \file
\brief
Nonlinear viscosity according to a modified power law


\level 3
*/
/*----------------------------------------------------------------------*/


#include "4C_mat_modpowerlaw.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ModPowerLaw::ModPowerLaw(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      m_cons_(matdata.parameters.get<double>("MCONS")),
      delta_(matdata.parameters.get<double>("DELTA")),
      a_exp_(matdata.parameters.get<double>("AEXP")),
      density_(matdata.parameters.get<double>("DENSITY"))
{
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::ModPowerLaw::create_material()
{
  return Teuchos::rcp(new Mat::ModPowerLaw(this));
}

Mat::ModPowerLawType Mat::ModPowerLawType::instance_;


Core::Communication::ParObject* Mat::ModPowerLawType::create(const std::vector<char>& data)
{
  Mat::ModPowerLaw* powLaw = new Mat::ModPowerLaw();
  powLaw->unpack(data);
  return powLaw;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ModPowerLaw::ModPowerLaw() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ModPowerLaw::ModPowerLaw(Mat::PAR::ModPowerLaw* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ModPowerLaw::pack(Core::Communication::PackBuffer& data) const
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
void Mat::ModPowerLaw::unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::ModPowerLaw*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
