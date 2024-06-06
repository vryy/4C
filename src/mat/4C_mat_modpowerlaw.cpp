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
Mat::PAR::ModPowerLaw::ModPowerLaw(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata),
      m_cons_(matdata->Get<double>("MCONS")),
      delta_(matdata->Get<double>("DELTA")),
      a_exp_(matdata->Get<double>("AEXP")),
      density_(matdata->Get<double>("DENSITY"))
{
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::ModPowerLaw::create_material()
{
  return Teuchos::rcp(new Mat::ModPowerLaw(this));
}

Mat::ModPowerLawType Mat::ModPowerLawType::instance_;


Core::Communication::ParObject* Mat::ModPowerLawType::Create(const std::vector<char>& data)
{
  Mat::ModPowerLaw* powLaw = new Mat::ModPowerLaw();
  powLaw->Unpack(data);
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
void Mat::ModPowerLaw::Pack(Core::Communication::PackBuffer& data) const
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ModPowerLaw::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
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
        params_ = static_cast<Mat::PAR::ModPowerLaw*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
