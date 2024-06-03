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
MAT::PAR::LubricationMat::LubricationMat(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata),
      density_(matdata->Get<double>("DENSITY")),
      lubricationlawID_(matdata->Get<int>("LUBRICATIONLAWID")),
      lubricationlaw_(nullptr)
{
  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    FOUR_C_THROW("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    FOUR_C_THROW("List of materials in the global problem instance is empty.");

  auto* curmat = GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(lubricationlawID_);

  switch (curmat->Type())
  {
    case CORE::Materials::m_lubrication_law_constant:
    {
      lubricationlaw_ = static_cast<MAT::PAR::LubricationLaw*>(curmat);
      break;
    }
    case CORE::Materials::m_lubrication_law_barus:
    {
      lubricationlaw_ = static_cast<MAT::PAR::LubricationLaw*>(curmat);
      break;
    }
    case CORE::Materials::m_lubrication_law_roeland:
    {
      lubricationlaw_ = static_cast<MAT::PAR::LubricationLaw*>(curmat);
      break;
    }
    default:
      FOUR_C_THROW("invalid material for lubrication law %d", curmat->Type());
      break;
  }
}


Teuchos::RCP<CORE::MAT::Material> MAT::PAR::LubricationMat::create_material()
{
  return Teuchos::rcp(new MAT::LubricationMat(this));
}

MAT::LubricationMatType MAT::LubricationMatType::instance_;

CORE::COMM::ParObject* MAT::LubricationMatType::Create(const std::vector<char>& data)
{
  MAT::LubricationMat* lubrication_mat = new MAT::LubricationMat();
  lubrication_mat->Unpack(data);
  return lubrication_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::LubricationMat::LubricationMat() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::LubricationMat::LubricationMat(MAT::PAR::LubricationMat* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::LubricationMat::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::LubricationMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::LubricationMat*>(mat);
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
double MAT::LubricationMat::ComputeViscosity(const double press)
{
  double visc = -1.;
  params_->lubricationlaw_->ComputeViscosity(press, visc);

  return visc;
}


/*----------------------------------------------------------------------*
                                                              wirtz 09/16|
*----------------------------------------------------------------------*/
double MAT::LubricationMat::compute_viscosity_deriv(const double press, const double visc)
{
  double visc_dp = -1;
  params_->lubricationlaw_->constitutive_derivatives(press, visc, visc_dp);

  return visc_dp;
}

FOUR_C_NAMESPACE_CLOSE
