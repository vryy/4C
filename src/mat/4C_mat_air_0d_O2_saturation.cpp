/*----------------------------------------------------------------------*/
/*! \file

\brief Gives relevant quantities of O2 saturation of air, used for scatra in reduced dimensional
airway elements framework (transport in elements and between air and blood)


\level 3
*/
/*----------------------------------------------------------------------*/


#include "4C_mat_air_0d_O2_saturation.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Air0dO2Saturation::Air0dO2Saturation(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata),
      atmospheric_p_(matdata->Get<double>("AtmosphericPressure")),
      nO2_per_VO2_(matdata->Get<double>("NumberOfO2PerVO2"))
{
}

Teuchos::RCP<CORE::MAT::Material> MAT::PAR::Air0dO2Saturation::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Air0dO2Saturation(this));
}


MAT::Air0dO2SaturationType MAT::Air0dO2SaturationType::instance_;


CORE::COMM::ParObject* MAT::Air0dO2SaturationType::Create(const std::vector<char>& data)
{
  MAT::Air0dO2Saturation* air_0d_O2_sat = new MAT::Air0dO2Saturation();
  air_0d_O2_sat->Unpack(data);
  return air_0d_O2_sat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Air0dO2Saturation::Air0dO2Saturation() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Air0dO2Saturation::Air0dO2Saturation(MAT::PAR::Air0dO2Saturation* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Air0dO2Saturation::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::Air0dO2Saturation::Unpack(const std::vector<char>& data)
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
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Air0dO2Saturation*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
