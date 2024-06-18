/*----------------------------------------------------------------------*/
/*! \file
\brief non-Newtonian fluid of Herschel-Bulkley type

\level 2

*/
/*----------------------------------------------------------------------*/


#include "4C_mat_herschelbulkley.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::HerschelBulkley::HerschelBulkley(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      tau0_(matdata.parameters.get<double>("TAU_0")),
      kfac_(matdata.parameters.get<double>("KFAC")),
      nexp_(matdata.parameters.get<double>("NEXP")),
      mexp_(matdata.parameters.get<double>("MEXP")),
      lolimshearrate_(matdata.parameters.get<double>("LOLIMSHEARRATE")),
      uplimshearrate_(matdata.parameters.get<double>("UPLIMSHEARRATE")),
      density_(matdata.parameters.get<double>("DENSITY"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::HerschelBulkley::create_material()
{
  return Teuchos::rcp(new Mat::HerschelBulkley(this));
}


Mat::HerschelBulkleyType Mat::HerschelBulkleyType::instance_;


Core::Communication::ParObject* Mat::HerschelBulkleyType::Create(const std::vector<char>& data)
{
  Mat::HerschelBulkley* herbul = new Mat::HerschelBulkley();
  herbul->unpack(data);
  return herbul;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::HerschelBulkley::HerschelBulkley() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::HerschelBulkley::HerschelBulkley(Mat::PAR::HerschelBulkley* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::HerschelBulkley::pack(Core::Communication::PackBuffer& data) const
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
void Mat::HerschelBulkley::unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::HerschelBulkley*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
