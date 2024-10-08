/*----------------------------------------------------------------------*/
/*! \file
\brief Newtonian fluid material

\level 1

*/
/*----------------------------------------------------------------------*/


#include "4C_mat_newtonianfluid.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::NewtonianFluid::NewtonianFluid(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      viscosity_(matdata.parameters.get<double>("DYNVISCOSITY")),
      density_(matdata.parameters.get<double>("DENSITY")),
      gamma_(matdata.parameters.get<double>("GAMMA"))
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::NewtonianFluid::create_material()
{
  return Teuchos::RCP(new Mat::NewtonianFluid(this));
}


Mat::NewtonianFluidType Mat::NewtonianFluidType::instance_;


Core::Communication::ParObject* Mat::NewtonianFluidType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::NewtonianFluid* fluid = new Mat::NewtonianFluid();
  fluid->unpack(buffer);
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
void Mat::NewtonianFluid::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::NewtonianFluid::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::NewtonianFluid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

FOUR_C_NAMESPACE_CLOSE
