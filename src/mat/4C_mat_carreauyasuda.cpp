/*----------------------------------------------------------------------*/
/*! \file
\brief
Former file of Ursula Mayer

\level 3


*/
/*----------------------------------------------------------------------*/


#include "4C_mat_carreauyasuda.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::CarreauYasuda::CarreauYasuda(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      nu_0_(matdata.parameters.get<double>("NU_0")),
      nu_inf_(matdata.parameters.get<double>("NU_INF")),
      lambda_(matdata.parameters.get<double>("LAMBDA")),
      a_param_(matdata.parameters.get<double>("APARAM")),
      b_param_(matdata.parameters.get<double>("BPARAM")),
      density_(matdata.parameters.get<double>("DENSITY"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::CarreauYasuda::create_material()
{
  return Teuchos::rcp(new Mat::CarreauYasuda(this));
}


Mat::CarreauYasudaType Mat::CarreauYasudaType::instance_;


Core::Communication::ParObject* Mat::CarreauYasudaType::Create(const std::vector<char>& data)
{
  Mat::CarreauYasuda* carYas = new Mat::CarreauYasuda();
  carYas->unpack(data);
  return carYas;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::CarreauYasuda::CarreauYasuda() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::CarreauYasuda::CarreauYasuda(Mat::PAR::CarreauYasuda* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::CarreauYasuda::pack(Core::Communication::PackBuffer& data) const
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
void Mat::CarreauYasuda::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
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
        params_ = static_cast<Mat::PAR::CarreauYasuda*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
