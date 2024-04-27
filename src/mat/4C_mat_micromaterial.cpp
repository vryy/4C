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
MAT::PAR::MicroMaterial::MicroMaterial(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      microfile_((matdata->Get<std::string>("MICROFILE"))),
      microdisnum_(matdata->Get<int>("MICRODIS_NUM")),
      initvol_(matdata->Get<double>("INITVOL"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::MicroMaterial::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MicroMaterial(this));
}


MAT::MicroMaterialType MAT::MicroMaterialType::instance_;


CORE::COMM::ParObject* MAT::MicroMaterialType::Create(const std::vector<char>& data)
{
  MAT::MicroMaterial* micro = new MAT::MicroMaterial();
  micro->Unpack(data);
  return micro;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MicroMaterial::MicroMaterial() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MicroMaterial::MicroMaterial(MAT::PAR::MicroMaterial* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MicroMaterial::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
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
void MAT::MicroMaterial::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::MicroMaterial*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
