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


Core::Communication::ParObject* Mat::HerschelBulkleyType::create(const std::vector<char>& data)
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
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::HerschelBulkley::unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::HerschelBulkley*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
