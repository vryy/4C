/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport material according to Sutherland law with
       Arrhenius-type chemical kinetics (progress variable)

\level 2

*/
/*----------------------------------------------------------------------*/


#include "4C_mat_arrhenius_pv.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ArrheniusPV::ArrheniusPV(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata),
      refvisc_(matdata->Get<double>("REFVISC")),
      reftemp_(matdata->Get<double>("REFTEMP")),
      suthtemp_(matdata->Get<double>("SUTHTEMP")),
      pranum_(matdata->Get<double>("PRANUM")),
      preexcon_(matdata->Get<double>("PREEXCON")),
      tempexp_(matdata->Get<double>("TEMPEXP")),
      actemp_(matdata->Get<double>("ACTEMP")),
      unbshc_(matdata->Get<double>("UNBSHC")),
      burshc_(matdata->Get<double>("BURSHC")),
      unbtemp_(matdata->Get<double>("UNBTEMP")),
      burtemp_(matdata->Get<double>("BURTEMP")),
      unbdens_(matdata->Get<double>("UNBDENS")),
      burdens_(matdata->Get<double>("BURDENS"))
{
}

Teuchos::RCP<CORE::MAT::Material> MAT::PAR::ArrheniusPV::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ArrheniusPV(this));
}


MAT::ArrheniusPVType MAT::ArrheniusPVType::instance_;


CORE::COMM::ParObject* MAT::ArrheniusPVType::Create(const std::vector<char>& data)
{
  MAT::ArrheniusPV* arrhenius_pv = new MAT::ArrheniusPV();
  arrhenius_pv->Unpack(data);
  return arrhenius_pv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ArrheniusPV::ArrheniusPV() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ArrheniusPV::ArrheniusPV(MAT::PAR::ArrheniusPV* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ArrheniusPV::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::ArrheniusPV::Unpack(const std::vector<char>& data)
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
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ArrheniusPV*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeTemperature(const double provar) const
{
  const double temperature = UnbTemp() + provar * (BurTemp() - UnbTemp());

  return temperature;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeDensity(const double provar) const
{
  // BML hypothesis
  const double density = UnbDens() + provar * (BurDens() - UnbDens());

  // equation of state
  // const double density = UnbDens()*BurDens()/(BurDens() + provar * (UnbDens() - BurDens()));

  return density;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeFactor(const double provar) const
{
  // BML hypothesis
  const double factor = (UnbDens() - BurDens()) / (UnbDens() + provar * (BurDens() - UnbDens()));

  // equation of state
  // const double factor = (UnbDens() - BurDens())/(BurDens() + provar * (UnbDens() - BurDens()));

  return factor;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeShc(const double provar) const
{
  const double shc = UnbShc() + provar * (BurShc() - UnbShc());

  return shc;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeViscosity(const double temp) const
{
  // previous implementation using "pow"-function appears to be extremely
  // time-consuming sometimes, at least on the computing cluster
  // const double visc =
  // std::pow((temp/RefTemp()),1.5)*((RefTemp()+SuthTemp())/(temp+SuthTemp()))*RefVisc();
  const double visc = sqrt((temp / RefTemp()) * (temp / RefTemp()) * (temp / RefTemp())) *
                      ((RefTemp() + SuthTemp()) / (temp + SuthTemp())) * RefVisc();

  return visc;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeDiffusivity(const double temp) const
{
  // previous implementation using "pow"-function appears to be extremely
  // time-consuming sometimes, at least on the computing cluster
  // const double diffus =
  // std::pow((temp/RefTemp()),1.5)*((RefTemp()+SuthTemp())/(temp+SuthTemp()))*RefVisc()/PraNum();
  const double diffus = sqrt((temp / RefTemp()) * (temp / RefTemp()) * (temp / RefTemp())) *
                        ((RefTemp() + SuthTemp()) / (temp + SuthTemp())) * RefVisc() / PraNum();

  return diffus;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeReactionCoeff(const double temp) const
{
  const double reacoeff = -PreExCon() * pow(temp, TempExp()) * exp(-AcTemp() / temp);

  return reacoeff;
}

FOUR_C_NAMESPACE_CLOSE
