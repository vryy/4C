/*----------------------------------------------------------------------*/
/*!
\file ferech_pv.cpp

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

#include "ferech_pv.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::FerEchPV::FerEchPV(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  refvisc_(matdata->GetDouble("REFVISC")),
  reftemp_(matdata->GetDouble("REFTEMP")),
  suthtemp_(matdata->GetDouble("SUTHTEMP")),
  pranum_(matdata->GetDouble("PRANUM")),
  reacratecon_(matdata->GetDouble("REACRATECON")),
  pvcrit_(matdata->GetDouble("PVCRIT")),
  unbshc_(matdata->GetDouble("UNBSHC")),
  burshc_(matdata->GetDouble("BURSHC")),
  unbtemp_(matdata->GetDouble("UNBTEMP")),
  burtemp_(matdata->GetDouble("BURTEMP")),
  unbdens_(matdata->GetDouble("UNBDENS")),
  burdens_(matdata->GetDouble("BURDENS")),
  mod_(matdata->GetDouble("MOD"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::FerEchPV::FerEchPV()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::FerEchPV::FerEchPV(MAT::PAR::FerEchPV* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::FerEchPV::Pack(vector<char>& data) const
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
void MAT::FerEchPV::Unpack(const vector<char>& data)
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
    params_ = static_cast<MAT::PAR::FerEchPV*>(mat);
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
double MAT::FerEchPV::ComputeTemperature(const double provar) const
{
  const double temperature = UnbTemp() + provar * (BurTemp() - UnbTemp());

  return temperature;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeDensity(const double provar) const
{
  const double density = UnbDens() + provar * (BurDens() - UnbDens());

  return density;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeShc(const double provar) const
{
  const double shc = UnbShc() + provar * (BurShc() - UnbShc());

  return shc;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeViscosity(const double temp) const
{
  // Sutherland law
  //const double visc = pow((temp/RefTemp()),1.5)*((RefTemp()+SuthTemp())/(temp+SuthTemp()))*RefVisc();

  // modified version by Poinsot and Veynante (1993)
  double visc = RefVisc();

  // original version by Ferziger and Echekki (1993)
  if (Mod() < (1.0-EPS15)) visc *= temp/RefTemp();

  return visc;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeDiffusivity(const double temp) const
{
  // Sutherland law
  //const double diffus = pow((temp/RefTemp()),1.5)*((RefTemp()+SuthTemp())/(temp+SuthTemp()))*RefVisc()/PraNum();

  // modified version by Poinsot and Veynante (1993)
  double diffus = RefVisc()/PraNum();

  // original version by Ferziger and Echekki (1993)
  if (Mod() < (1.0-EPS15)) diffus *= temp/RefTemp();

  return diffus;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::FerEchPV::ComputeReactionCoeff(const double provar) const
{
  double reacoeff = 0.0;

  // Heaviside function loop
  if (provar > PvCrit())
  {
    // original version by Ferziger and Echekki (1993)
    reacoeff = ReacRateCon();

    // modified version by Poinsot and Veynante (1993)
    if (Mod() > EPS15) reacoeff *= UnbDens()/(UnbDens() + provar * (BurDens() - UnbDens()));
  }

  return reacoeff;
}

#endif
