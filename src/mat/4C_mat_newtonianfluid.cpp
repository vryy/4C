/*----------------------------------------------------------------------*/
/*! \file
\brief Newtonian fluid material

\level 1

*/
/*----------------------------------------------------------------------*/


#include "4C_mat_newtonianfluid.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::NewtonianFluid::NewtonianFluid(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(matdata->Get<double>("DYNVISCOSITY")),
      density_(matdata->Get<double>("DENSITY")),
      gamma_(matdata->Get<double>("GAMMA"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::NewtonianFluid::create_material()
{
  return Teuchos::rcp(new Mat::NewtonianFluid(this));
}


Mat::NewtonianFluidType Mat::NewtonianFluidType::instance_;


Core::Communication::ParObject* Mat::NewtonianFluidType::Create(const std::vector<char>& data)
{
  Mat::NewtonianFluid* fluid = new Mat::NewtonianFluid();
  fluid->Unpack(data);
  return fluid;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::NewtonianFluid::NewtonianFluid() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::NewtonianFluid::NewtonianFluid(Mat::PAR::NewtonianFluid* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::NewtonianFluid::Pack(Core::Communication::PackBuffer& data) const
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
void Mat::NewtonianFluid::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::NewtonianFluid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
