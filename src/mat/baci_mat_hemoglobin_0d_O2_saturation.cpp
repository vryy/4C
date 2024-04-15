/*----------------------------------------------------------------------*/
/*! \file
\brief Gives relevant quantities of O2 saturation of blood (hemoglobin), used for scatra in reduced
dimensional airway elements framework (transport in elements and between air and blood)


\level 3
*/
/*----------------------------------------------------------------------*/


#include "baci_mat_hemoglobin_0d_O2_saturation.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Hemoglobin0dO2Saturation::Hemoglobin0dO2Saturation(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      per_volume_blood_(*matdata->Get<double>("PerVolumeBlood")),
      o2_sat_per_vol_blood_(*matdata->Get<double>("O2SaturationPerVolBlood")),
      p_half_(*matdata->Get<double>("PressureHalf")),
      power_(*matdata->Get<double>("Power")),
      nO2_per_VO2_(*matdata->Get<double>("NumberOfO2PerVO2"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Hemoglobin0dO2Saturation::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Hemoglobin0dO2Saturation(this));
}


MAT::Hemoglobin0dO2SaturationType MAT::Hemoglobin0dO2SaturationType::instance_;


CORE::COMM::ParObject* MAT::Hemoglobin0dO2SaturationType::Create(const std::vector<char>& data)
{
  MAT::Hemoglobin0dO2Saturation* hem_0d_O2_sat = new MAT::Hemoglobin0dO2Saturation();
  hem_0d_O2_sat->Unpack(data);
  return hem_0d_O2_sat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Hemoglobin0dO2Saturation::Hemoglobin0dO2Saturation() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Hemoglobin0dO2Saturation::Hemoglobin0dO2Saturation(MAT::PAR::Hemoglobin0dO2Saturation* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Hemoglobin0dO2Saturation::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::Hemoglobin0dO2Saturation::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::Hemoglobin0dO2Saturation*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
