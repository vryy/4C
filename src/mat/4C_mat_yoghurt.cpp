/*----------------------------------------------------------------------*/
/*! \file
\brief yoghurt-type fluid

\level 2

*/
/*----------------------------------------------------------------------*/


#include "4C_mat_yoghurt.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Yoghurt::Yoghurt(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      shc_(matdata->Get<double>("SHC")),
      density_(matdata->Get<double>("DENSITY")),
      thermcond_(matdata->Get<double>("THERMCOND")),
      strrateexp_(matdata->Get<double>("STRAINRATEEXP")),
      preexcon_(matdata->Get<double>("PREEXCON")),
      actenergy_(matdata->Get<double>("ACTENERGY")),
      gasconst_(matdata->Get<double>("GASCON")),
      delta_(matdata->Get<double>("DELTA"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Yoghurt::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Yoghurt(this));
}


MAT::YoghurtType MAT::YoghurtType::instance_;


CORE::COMM::ParObject* MAT::YoghurtType::Create(const std::vector<char>& data)
{
  MAT::Yoghurt* yoghurt = new MAT::Yoghurt();
  yoghurt->Unpack(data);
  return yoghurt;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Yoghurt::Yoghurt() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Yoghurt::Yoghurt(MAT::PAR::Yoghurt* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Yoghurt::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::Yoghurt::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::Yoghurt*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Yoghurt::ComputeViscosity(const double rateofstrain, const double temp) const
{
  const double visc = PreExCon() * pow(abs(rateofstrain) + Delta(), (StrRateExp() - 1)) *
                      exp(ActEnergy() / (temp * GasConst()));

  return visc;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Yoghurt::ComputeDiffusivity() const
{
  const double diffus = ThermCond() / Shc();

  return diffus;
}

FOUR_C_NAMESPACE_CLOSE
