/*----------------------------------------------------------------------------*/
/*! \file
\brief Weakly compressible fluid

\level 1

*/
/*----------------------------------------------------------------------------*/

#include <vector>
#include "fluid_weakly_compressible.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::WeaklyCompressibleFluid::WeaklyCompressibleFluid(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(matdata->GetDouble("VISCOSITY")),
      refdensity_(matdata->GetDouble("REFDENSITY")),
      refpressure_(matdata->GetDouble("REFPRESSURE")),
      comprcoeff_(matdata->GetDouble("COMPRCOEFF"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::WeaklyCompressibleFluid::CreateMaterial()
{
  return Teuchos::rcp(new MAT::WeaklyCompressibleFluid(this));
}


MAT::WeaklyCompressibleFluidType MAT::WeaklyCompressibleFluidType::instance_;


DRT::ParObject* MAT::WeaklyCompressibleFluidType::Create(const std::vector<char>& data)
{
  MAT::WeaklyCompressibleFluid* fluid = new MAT::WeaklyCompressibleFluid();
  fluid->Unpack(data);
  return fluid;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::WeaklyCompressibleFluid::WeaklyCompressibleFluid() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::WeaklyCompressibleFluid::WeaklyCompressibleFluid(MAT::PAR::WeaklyCompressibleFluid* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::WeaklyCompressibleFluid::Pack(DRT::PackBuffer& data) const
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
void MAT::WeaklyCompressibleFluid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
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
        params_ = static_cast<MAT::PAR::WeaklyCompressibleFluid*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::WeaklyCompressibleFluid::ComputeDensity(const double press) const
{
  const double density = RefDensity() + ComprCoeff() * (press - RefPressure());

  return density;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::WeaklyCompressibleFluid::ComputePressure(const double dens) const
{
  const double pressure = RefPressure() + 1.0 / ComprCoeff() * (dens - RefDensity());

  return pressure;
}
