/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport material according to Sutherland law with
       Arrhenius-type chemical kinetics (species)

\level 2

*/
/*----------------------------------------------------------------------*/


#include "baci_mat_arrhenius_spec.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ArrheniusSpec::ArrheniusSpec(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      refvisc_(*matdata->Get<double>("REFVISC")),
      reftemp_(*matdata->Get<double>("REFTEMP")),
      suthtemp_(*matdata->Get<double>("SUTHTEMP")),
      schnum_(*matdata->Get<double>("SCHNUM")),
      preexcon_(*matdata->Get<double>("PREEXCON")),
      tempexp_(*matdata->Get<double>("TEMPEXP")),
      actemp_(*matdata->Get<double>("ACTEMP")),
      gasconst_(*matdata->Get<double>("GASCON"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ArrheniusSpec::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ArrheniusSpec(this));
}


MAT::ArrheniusSpecType MAT::ArrheniusSpecType::instance_;


CORE::COMM::ParObject* MAT::ArrheniusSpecType::Create(const std::vector<char>& data)
{
  MAT::ArrheniusSpec* arrhenius_spec = new MAT::ArrheniusSpec();
  arrhenius_spec->Unpack(data);
  return arrhenius_spec;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ArrheniusSpec::ArrheniusSpec() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ArrheniusSpec::ArrheniusSpec(MAT::PAR::ArrheniusSpec* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ArrheniusSpec::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::ArrheniusSpec::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ArrheniusSpec*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusSpec::ComputeViscosity(const double temp) const
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
double MAT::ArrheniusSpec::ComputeDiffusivity(const double temp) const
{
  // previous implementation using "pow"-function appears to be extremely
  // time-consuming sometimes, at least on the computing cluster
  // const double diffus =
  // std::pow((temp/RefTemp()),1.5)*((RefTemp()+SuthTemp())/(temp+SuthTemp()))*RefVisc()/SchNum();
  const double diffus = sqrt((temp / RefTemp()) * (temp / RefTemp()) * (temp / RefTemp())) *
                        ((RefTemp() + SuthTemp()) / (temp + SuthTemp())) * RefVisc() / SchNum();

  return diffus;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusSpec::ComputeDensity(const double temp, const double thermpress) const
{
  const double density = thermpress / (GasConst() * temp);

  return density;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusSpec::ComputeReactionCoeff(const double temp) const
{
  const double reacoeff = -PreExCon() * pow(temp, TempExp()) * exp(-AcTemp() / temp);

  return reacoeff;
}

FOUR_C_NAMESPACE_CLOSE
