/*----------------------------------------------------------------------------*/
/*! \file
\brief Weakly compressible fluid

\level 1

*/
/*----------------------------------------------------------------------------*/

#include "baci_mat_fluid_weakly_compressible.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::WeaklyCompressibleFluid::WeaklyCompressibleFluid(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(*matdata->Get<double>("VISCOSITY")),
      refdensity_(*matdata->Get<double>("REFDENSITY")),
      refpressure_(*matdata->Get<double>("REFPRESSURE")),
      comprcoeff_(*matdata->Get<double>("COMPRCOEFF"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::WeaklyCompressibleFluid::CreateMaterial()
{
  return Teuchos::rcp(new MAT::WeaklyCompressibleFluid(this));
}


MAT::WeaklyCompressibleFluidType MAT::WeaklyCompressibleFluidType::instance_;


CORE::COMM::ParObject* MAT::WeaklyCompressibleFluidType::Create(const std::vector<char>& data)
{
  MAT::WeaklyCompressibleFluid* fluid = new MAT::WeaklyCompressibleFluid();
  fluid->Unpack(data);
  return fluid;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::WeaklyCompressibleFluid::WeaklyCompressibleFluid() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::WeaklyCompressibleFluid::WeaklyCompressibleFluid(MAT::PAR::WeaklyCompressibleFluid* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::WeaklyCompressibleFluid::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::WeaklyCompressibleFluid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
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

FOUR_C_NAMESPACE_CLOSE
