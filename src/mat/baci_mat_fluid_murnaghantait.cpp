/*----------------------------------------------------------------------*/
/*! \file
\brief Weakly compressible fluid according to Murnaghan-Tait

\level 1

*/
/*----------------------------------------------------------------------*/

#include "baci_mat_fluid_murnaghantait.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::MurnaghanTaitFluid::MurnaghanTaitFluid(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(*matdata->Get<double>("DYNVISCOSITY")),
      refdensity_(*matdata->Get<double>("REFDENSITY")),
      refpressure_(*matdata->Get<double>("REFPRESSURE")),
      refbulkmodulus_(*matdata->Get<double>("REFBULKMODULUS")),
      matparameter_(*matdata->Get<double>("MATPARAMETER")),
      gamma_(*matdata->Get<double>("GAMMA"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::MurnaghanTaitFluid::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MurnaghanTaitFluid(this));
}


MAT::MurnaghanTaitFluidType MAT::MurnaghanTaitFluidType::instance_;


CORE::COMM::ParObject* MAT::MurnaghanTaitFluidType::Create(const std::vector<char>& data)
{
  MAT::MurnaghanTaitFluid* fluid = new MAT::MurnaghanTaitFluid();
  fluid->Unpack(data);
  return fluid;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MurnaghanTaitFluid::MurnaghanTaitFluid() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MurnaghanTaitFluid::MurnaghanTaitFluid(MAT::PAR::MurnaghanTaitFluid* params) : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MurnaghanTaitFluid::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::MurnaghanTaitFluid::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::MurnaghanTaitFluid*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::MurnaghanTaitFluid::ComputeDensity(const double press) const
{
  // $ \rho=\rho_0[\dfrac{n}{K_0}\left(P-P_0)+1\right]^{\dfrac{1}{n}} $
  const double density =
      RefDensity() * std::pow((MatParameter() / RefBulkModulus() * (press - RefPressure()) + 1.0),
                         (1.0 / MatParameter()));

  return density;
}

BACI_NAMESPACE_CLOSE
