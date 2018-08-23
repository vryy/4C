/*----------------------------------------------------------------------*/
/*!
\file fluid_linear_density_viscosity.cpp

\brief Linear law (pressure-dependent) for the density and the viscosity

<pre>
\level 1

\maintainer Andrea La Spina
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#include <vector>
#include "fluid_linear_density_viscosity.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::LinearDensityViscosity::LinearDensityViscosity(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  refdensity_(matdata->GetDouble("REFDENSITY")),
  refviscosity_(matdata->GetDouble("REFVISCOSITY")),
  refpressure_(matdata->GetDouble("REFPRESSURE")),
  coeffdensity_(matdata->GetDouble("COEFFDENSITY")),
  coeffviscosity_(matdata->GetDouble("COEFFVISCOSITY")),
  gamma_(matdata->GetDouble("GAMMA"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::LinearDensityViscosity::CreateMaterial()
{
  return Teuchos::rcp(new MAT::LinearDensityViscosity(this));
}


MAT::LinearDensityViscosityType MAT::LinearDensityViscosityType::instance_;


DRT::ParObject* MAT::LinearDensityViscosityType::Create( const std::vector<char> & data )
{
  MAT::LinearDensityViscosity* fluid = new MAT::LinearDensityViscosity();
  fluid->Unpack(data);
  return fluid;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::LinearDensityViscosity::LinearDensityViscosity()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::LinearDensityViscosity::LinearDensityViscosity(MAT::PAR::LinearDensityViscosity* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::LinearDensityViscosity::Pack(DRT::PackBuffer& data) const
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
void MAT::LinearDensityViscosity::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::LinearDensityViscosity*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::LinearDensityViscosity::ComputeDensity(const double press) const
{
  const double density = RefDensity() * (1.0 + CoeffDensity() * (press - RefPressure()));

  return density;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::LinearDensityViscosity::ComputeViscosity(const double press) const
{
  const double viscosity = RefViscosity() * (1.0 + CoeffViscosity() * (press - RefPressure()));

  return viscosity;
}
