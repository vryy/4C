/*----------------------------------------------------------------------*/
/*! \file
\brief class for handling of micro-macro transitions

\level 3

*/
/*----------------------------------------------------------------------*/


#include "4C_mat_micromaterial.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


// Be careful when adding new member functions of MicroMaterial that
// relate to MicroMaterialGP (which is NOT in the filter). See also
// comments in micromaterial_evaluate.cpp which is separated from
// this file for the very reason.


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::MicroMaterial::MicroMaterial(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata),
      microfile_((matdata->Get<std::string>("MICROFILE"))),
      microdisnum_(matdata->Get<int>("MICRODIS_NUM")),
      initvol_(matdata->Get<double>("INITVOL"))
{
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::MicroMaterial::create_material()
{
  return Teuchos::rcp(new Mat::MicroMaterial(this));
}


Mat::MicroMaterialType Mat::MicroMaterialType::instance_;


Core::Communication::ParObject* Mat::MicroMaterialType::Create(const std::vector<char>& data)
{
  Mat::MicroMaterial* micro = new Mat::MicroMaterial();
  micro->Unpack(data);
  return micro;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::MicroMaterial::MicroMaterial() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::MicroMaterial::MicroMaterial(Mat::PAR::MicroMaterial* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MicroMaterial::Pack(Core::Communication::PackBuffer& data) const
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
void Mat::MicroMaterial::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::MicroMaterial*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
