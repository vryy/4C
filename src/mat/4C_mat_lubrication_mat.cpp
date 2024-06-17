/*--------------------------------------------------------------------------*/
/*! \file
\brief Material model for the lubrication film

\level 3

*/
/*--------------------------------------------------------------------------*/


#include "4C_mat_lubrication_mat.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_lubrication_law.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::LubricationMat::LubricationMat(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      density_(matdata.parameters.Get<double>("DENSITY")),
      lubricationlawID_(matdata.parameters.Get<int>("LUBRICATIONLAWID")),
      lubricationlaw_(nullptr)
{
  // retrieve problem instance to read from
  const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (Global::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (Global::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  auto* curmat = Global::Problem::Instance(probinst)->Materials()->ParameterById(lubricationlawID_);

  switch (curmat->Type())
  {
    case Core::Materials::m_lubrication_law_constant:
    {
      lubricationlaw_ = static_cast<Mat::PAR::LubricationLaw*>(curmat);
      break;
    }
    case Core::Materials::m_lubrication_law_barus:
    {
      lubricationlaw_ = static_cast<Mat::PAR::LubricationLaw*>(curmat);
      break;
    }
    case Core::Materials::m_lubrication_law_roeland:
    {
      lubricationlaw_ = static_cast<Mat::PAR::LubricationLaw*>(curmat);
      break;
    }
    default:
      FOUR_C_THROW("invalid material for lubrication law %d", curmat->Type());
      break;
  }
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::LubricationMat::create_material()
{
  return Teuchos::rcp(new Mat::LubricationMat(this));
}

Mat::LubricationMatType Mat::LubricationMatType::instance_;

Core::Communication::ParObject* Mat::LubricationMatType::Create(const std::vector<char>& data)
{
  Mat::LubricationMat* lubrication_mat = new Mat::LubricationMat();
  lubrication_mat->Unpack(data);
  return lubrication_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::LubricationMat::LubricationMat() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::LubricationMat::LubricationMat(Mat::PAR::LubricationMat* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::LubricationMat::Pack(Core::Communication::PackBuffer& data) const
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
void Mat::LubricationMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
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
        params_ = static_cast<Mat::PAR::LubricationMat*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*
                                                              wirtz 09/16|
*----------------------------------------------------------------------*/
double Mat::LubricationMat::ComputeViscosity(const double press)
{
  double visc = -1.;
  params_->lubricationlaw_->ComputeViscosity(press, visc);

  return visc;
}


/*----------------------------------------------------------------------*
                                                              wirtz 09/16|
*----------------------------------------------------------------------*/
double Mat::LubricationMat::compute_viscosity_deriv(const double press, const double visc)
{
  double visc_dp = -1;
  params_->lubricationlaw_->constitutive_derivatives(press, visc, visc_dp);

  return visc_dp;
}

FOUR_C_NAMESPACE_CLOSE
