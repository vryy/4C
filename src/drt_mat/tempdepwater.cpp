/*----------------------------------------------------------------------*/
/*!
\file tempdepwater.cpp

\brief material for temperature-dependent water according to "VDI Waermeatlas"

\level 3

<pre>
\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/


#include <vector>

#include "tempdepwater.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::TempDepWater::TempDepWater(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  critdens_(matdata->GetDouble("CRITDENS")),
  crittemp_(matdata->GetDouble("CRITTEMP")),
  shc_(matdata->GetDouble("SHC"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::TempDepWater::CreateMaterial()
{
  return Teuchos::rcp(new MAT::TempDepWater(this));
}


MAT::TempDepWaterType MAT::TempDepWaterType::instance_;


DRT::ParObject* MAT::TempDepWaterType::Create( const std::vector<char> & data )
{
  MAT::TempDepWater* tempdepwater = new MAT::TempDepWater();
  tempdepwater->Unpack(data);
  return tempdepwater;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::TempDepWater::TempDepWater()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::TempDepWater::TempDepWater(MAT::PAR::TempDepWater* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::TempDepWater::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

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
void MAT::TempDepWater::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::TempDepWater*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::TempDepWater::ComputeViscosity(const double temp) const
{
  const double A_mu = -3.7188;
  const double B_mu =  578.919;
  const double C_mu = -137.546;

  const double visc = exp(A_mu+(B_mu/(C_mu+temp)))*1e-3;

  return visc;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::TempDepWater::ComputeDiffusivity(const double temp) const
{
  const double A_lambda = -2.4149;
  const double B_lambda =  2.45165e-02;
  const double C_lambda = -7.3121e-05;
  const double D_lambda =  9.9492e-08;
  const double E_lambda = -5.3730e-11;

  // compute diffusivity divided by specific heat capacity at constant pressure
  const double diffus = (A_lambda + B_lambda*temp + C_lambda*temp*temp + D_lambda*temp*temp*temp + E_lambda*temp*temp*temp*temp)/Shc();

  return diffus;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::TempDepWater::ComputeDensity(const double temp) const
{
  const double A_rho =  1094.0233;
  const double B_rho = -1813.2295;
  const double C_rho =  3863.9557;
  const double D_rho = -2479.8130;

  const double oneminustbytc = 1.0 - (temp/CritTemp());

  const double density = CritDens() + A_rho*std::pow(oneminustbytc,0.35) + B_rho*std::pow(oneminustbytc,(2.0/3.0)) + + C_rho*oneminustbytc + D_rho*std::pow(oneminustbytc,(4.0/3.0));

  return density;
}


