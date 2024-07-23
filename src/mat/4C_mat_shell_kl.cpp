/*----------------------------------------------------------------------*/
/*! \file
\brief Material for an elastic Kirchhoff-Love shell

\level 3

*----------------------------------------------------------------------*/


#include "4C_mat_shell_kl.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_utils_exceptions.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
Mat::PAR::KirchhoffLoveShell::KirchhoffLoveShell(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      young_modulus_(matdata.parameters.get<double>("YOUNG_MODULUS")),
      poisson_ratio_(matdata.parameters.get<double>("POISSON_RATIO")),
      thickness_(matdata.parameters.get<double>("THICKNESS"))
{
  if (young_modulus_ <= 0.0)
    FOUR_C_THROW("Young's modulus has to be positive. Got %f", young_modulus_);
  if (poisson_ratio_ <= -1.0 or poisson_ratio_ > 0.5)
    FOUR_C_THROW("Poisson's ration has to be in the range (-1,1/2]. Got %f", poisson_ratio_);
  if (thickness_ < 0.0) FOUR_C_THROW("Thickness has to be positive. Got %f", thickness_);
}

/**
 *
 */
Teuchos::RCP<Core::Mat::Material> Mat::PAR::KirchhoffLoveShell::create_material()
{
  return Teuchos::rcp(new Mat::KirchhoffLoveShell(this));
}

/**
 *
 */
Mat::KirchhoffLoveShellType Mat::KirchhoffLoveShellType::instance_;

/**
 *
 */
Core::Communication::ParObject* Mat::KirchhoffLoveShellType::create(const std::vector<char>& data)
{
  auto* shell = new Mat::KirchhoffLoveShell();
  shell->unpack(data);
  return shell;
}

/**
 *
 */
Mat::KirchhoffLoveShell::KirchhoffLoveShell(Mat::PAR::KirchhoffLoveShell* params) : params_(params)
{
}

/**
 *
 */
void Mat::KirchhoffLoveShell::pack(Core::Communication::PackBuffer& data) const
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

/**
 *
 */
void Mat::KirchhoffLoveShell::unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::KirchhoffLoveShell*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
