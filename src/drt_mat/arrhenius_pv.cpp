/*----------------------------------------------------------------------*/
/*!
\file arrhenius_pv.cpp

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>

#include "arrhenius_pv.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ArrheniusPV::ArrheniusPV(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  refvisc_(matdata->GetDouble("REFVISC")),
  reftemp_(matdata->GetDouble("REFTEMP")),
  suthtemp_(matdata->GetDouble("SUTHTEMP")),
  pranum_(matdata->GetDouble("PRANUM")),
  preexcon_(matdata->GetDouble("PREEXCON")),
  tempexp_(matdata->GetDouble("TEMPEXP")),
  actemp_(matdata->GetDouble("ACTEMP")),
  flamespeed_(matdata->GetDouble("FLAMESPEED")),
  unbshc_(matdata->GetDouble("UNBSHC")),
  burshc_(matdata->GetDouble("BURSHC")),
  unbtemp_(matdata->GetDouble("UNBTEMP")),
  burtemp_(matdata->GetDouble("BURTEMP")),
  unbdens_(matdata->GetDouble("UNBDENS")),
  burdens_(matdata->GetDouble("BURDENS"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ArrheniusPV::ArrheniusPV()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ArrheniusPV::ArrheniusPV(MAT::PAR::ArrheniusPV* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ArrheniusPV::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ArrheniusPV::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
  MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
  if (mat->Type() == MaterialType())
    params_ = static_cast<MAT::PAR::ArrheniusPV*>(mat);
  else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeAlpha() const
{
  const double alpha = (BurTemp()-UnbTemp())/BurTemp();

  return alpha;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeBeta() const
{
  const double beta = (BurTemp()-UnbTemp())*AcTemp()/(BurTemp()*BurTemp());

  return beta;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeDiffFlameThickness() const
{
  const double delta = RefVisc()/(PraNum()*UnbDens()*FlameSpeed());

  return delta;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeBlintFlameThickness() const
{
  const double delta_B = 2.0*pow((BurTemp()/RefTemp()),1.5)*((RefTemp()+SuthTemp())/(BurTemp()+SuthTemp()))*RefVisc()/(PraNum()*UnbDens()*FlameSpeed());

  return delta_B;
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
double MAT::ArrheniusPV::ComputeTemperature(const double provar) const
{
  const double temperature = UnbTemp() + provar * (BurTemp() - UnbTemp());

  return temperature;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeDensity(const double provar) const
{
  const double density = UnbDens() + provar * (BurDens() - UnbDens());

  return density;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeDiffusivity(const double temp) const
{
  const double diffus = pow((temp/RefTemp()),1.5)*((RefTemp()+SuthTemp())/(temp+SuthTemp()))*RefVisc()/PraNum();

  return diffus;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPV::ComputeReactionCoeff(const double temp,
                                              const double dens) const
{
  const double reacoeff = -PreExCon()*pow(temp,TempExp())*dens*exp(-AcTemp()/temp);

  return reacoeff;
}

#endif
