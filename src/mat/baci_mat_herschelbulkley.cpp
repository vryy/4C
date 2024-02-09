/*----------------------------------------------------------------------*/
/*! \file
\brief non-Newtonian fluid of Herschel-Bulkley type

\level 2

*/
/*----------------------------------------------------------------------*/


#include "baci_mat_herschelbulkley.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::HerschelBulkley::HerschelBulkley(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      tau0_(matdata->GetDouble("TAU_0")),
      kfac_(matdata->GetDouble("KFAC")),
      nexp_(matdata->GetDouble("NEXP")),
      mexp_(matdata->GetDouble("MEXP")),
      lolimshearrate_(matdata->GetDouble("LOLIMSHEARRATE")),
      uplimshearrate_(matdata->GetDouble("UPLIMSHEARRATE")),
      density_(matdata->GetDouble("DENSITY"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::HerschelBulkley::CreateMaterial()
{
  return Teuchos::rcp(new MAT::HerschelBulkley(this));
}


MAT::HerschelBulkleyType MAT::HerschelBulkleyType::instance_;


CORE::COMM::ParObject* MAT::HerschelBulkleyType::Create(const std::vector<char>& data)
{
  MAT::HerschelBulkley* herbul = new MAT::HerschelBulkley();
  herbul->Unpack(data);
  return herbul;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::HerschelBulkley::HerschelBulkley() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::HerschelBulkley::HerschelBulkley(MAT::PAR::HerschelBulkley* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::HerschelBulkley::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::HerschelBulkley::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
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
        params_ = static_cast<MAT::PAR::HerschelBulkley*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

BACI_NAMESPACE_CLOSE
