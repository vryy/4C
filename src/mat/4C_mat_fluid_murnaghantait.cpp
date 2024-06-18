/*----------------------------------------------------------------------*/
/*! \file
\brief Weakly compressible fluid according to Murnaghan-Tait

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_mat_fluid_murnaghantait.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::MurnaghanTaitFluid::MurnaghanTaitFluid(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      viscosity_(matdata.parameters.get<double>("DYNVISCOSITY")),
      refdensity_(matdata.parameters.get<double>("REFDENSITY")),
      refpressure_(matdata.parameters.get<double>("REFPRESSURE")),
      refbulkmodulus_(matdata.parameters.get<double>("REFBULKMODULUS")),
      matparameter_(matdata.parameters.get<double>("MATPARAMETER")),
      gamma_(matdata.parameters.get<double>("GAMMA"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::MurnaghanTaitFluid::create_material()
{
  return Teuchos::rcp(new Mat::MurnaghanTaitFluid(this));
}


Mat::MurnaghanTaitFluidType Mat::MurnaghanTaitFluidType::instance_;


Core::Communication::ParObject* Mat::MurnaghanTaitFluidType::Create(const std::vector<char>& data)
{
  Mat::MurnaghanTaitFluid* fluid = new Mat::MurnaghanTaitFluid();
  fluid->Unpack(data);
  return fluid;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::MurnaghanTaitFluid::MurnaghanTaitFluid() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::MurnaghanTaitFluid::MurnaghanTaitFluid(Mat::PAR::MurnaghanTaitFluid* params) : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MurnaghanTaitFluid::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MurnaghanTaitFluid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::MurnaghanTaitFluid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double Mat::MurnaghanTaitFluid::ComputeDensity(const double press) const
{
  // $ \rho=\rho_0[\dfrac{n}{K_0}\left(P-P_0)+1\right]^{\dfrac{1}{n}} $
  const double density =
      RefDensity() * std::pow((MatParameter() / RefBulkModulus() * (press - RefPressure()) + 1.0),
                         (1.0 / MatParameter()));

  return density;
}

FOUR_C_NAMESPACE_CLOSE
