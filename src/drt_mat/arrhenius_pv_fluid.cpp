/*----------------------------------------------------------------------*/
/*!
\file arrhenius_pv_fluid.cpp

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

#include "arrhenius_pv_fluid.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ArrheniusPVFluid::ArrheniusPVFluid(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  refvisc_(matdata->GetDouble("REFVISC")),
  reftemp_(matdata->GetDouble("REFTEMP")),
  suthtemp_(matdata->GetDouble("SUTHTEMP")),
  unbtemp_(matdata->GetDouble("UNBTEMP")),
  burtemp_(matdata->GetDouble("BURTEMP")),
  unbdens_(matdata->GetDouble("UNBDENS")),
  burdens_(matdata->GetDouble("BURDENS"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ArrheniusPVFluid::ArrheniusPVFluid()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ArrheniusPVFluid::ArrheniusPVFluid(MAT::PAR::ArrheniusPVFluid* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ArrheniusPVFluid::Pack(vector<char>& data) const
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
void MAT::ArrheniusPVFluid::Unpack(const vector<char>& data)
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
    params_ = static_cast<MAT::PAR::ArrheniusPVFluid*>(mat);
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
double MAT::ArrheniusPVFluid::ComputeTemperature(const double provar) const
{
  const double temperature = UnbTemp() + provar * (BurTemp() - UnbTemp());

  return temperature;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPVFluid::ComputeDensity(const double provar) const
{
  const double density = UnbDens() + provar * (BurDens() - UnbDens());

  return density;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ArrheniusPVFluid::ComputeViscosity(const double temp) const
{
  const double visc = pow((temp/RefTemp()),1.5)*((RefTemp()+SuthTemp())/(temp+SuthTemp()))*RefVisc();

  return visc;
}

#endif
