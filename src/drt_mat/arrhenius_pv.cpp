/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport material according to Sutherland law with
       Arrhenius-type chemical kinetics (progress variable)

\level 2

*/
/*----------------------------------------------------------------------*/


#include <vector>

#include "arrhenius_pv.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ArrheniusPV::ArrheniusPV(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      refvisc_(matdata->GetDouble("REFVISC")),
      reftemp_(matdata->GetDouble("REFTEMP")),
      suthtemp_(matdata->GetDouble("SUTHTEMP")),
      pranum_(matdata->GetDouble("PRANUM")),
      preexcon_(matdata->GetDouble("PREEXCON")),
      tempexp_(matdata->GetDouble("TEMPEXP")),
      actemp_(matdata->GetDouble("ACTEMP")),
      unbshc_(matdata->GetDouble("UNBSHC")),
      burshc_(matdata->GetDouble("BURSHC")),
      unbtemp_(matdata->GetDouble("UNBTEMP")),
      burtemp_(matdata->GetDouble("BURTEMP")),
      unbdens_(matdata->GetDouble("UNBDENS")),
      burdens_(matdata->GetDouble("BURDENS"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ArrheniusPV::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ArrheniusPV(this));
}


MAT::ArrheniusPVType MAT::ArrheniusPVType::instance_;


DRT::ParObject* MAT::ArrheniusPVType::Create(const std::vector<char>& data)
{
  MAT::ArrheniusPV* arrhenius_pv = new MAT::ArrheniusPV();
  arrhenius_pv->Unpack(data);
  return arrhenius_pv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ArrheniusPV::ArrheniusPV() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ArrheniusPV::ArrheniusPV(MAT::PAR::ArrheniusPV* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ArrheniusPV::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ArrheniusPV::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ArrheniusPV*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
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
