/*--------------------------------------------------------------------------*/
/*! \file
\brief Material model for the lubrication film

\level 3

*/
/*--------------------------------------------------------------------------*/


#include "baci_mat_lubrication_mat.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_lubrication_law.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::LubricationMat::LubricationMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      density_(*matdata->Get<double>("DENSITY")),
      lubricationlawID_(*matdata->Get<int>("LUBRICATIONLAWID")),
      lubricationlaw_(nullptr)
{
  // retrieve problem instance to read from
  const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();

  // for the sake of safety
  if (GLOBAL::Problem::Instance(probinst)->Materials() == Teuchos::null)
    dserror("List of materials cannot be accessed in the global problem instance.");
  // yet another safety check
  if (GLOBAL::Problem::Instance(probinst)->Materials()->Num() == 0)
    dserror("List of materials in the global problem instance is empty.");

  // retrieve validated input line of material ID in question
  Teuchos::RCP<MAT::PAR::Material> curmat =
      GLOBAL::Problem::Instance(probinst)->Materials()->ById(lubricationlawID_);

  switch (curmat->Type())
  {
    case INPAR::MAT::m_lubrication_law_constant:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationLawConstant(curmat));
      lubricationlaw_ = static_cast<MAT::PAR::LubricationLaw*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_lubrication_law_barus:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationLawBarus(curmat));
      lubricationlaw_ = static_cast<MAT::PAR::LubricationLaw*>(curmat->Parameter());
      break;
    }
    case INPAR::MAT::m_lubrication_law_roeland:
    {
      if (curmat->Parameter() == nullptr)
        curmat->SetParameter(new MAT::PAR::LubricationLawRoeland(curmat));
      lubricationlaw_ = static_cast<MAT::PAR::LubricationLaw*>(curmat->Parameter());
      break;
    }
    default:
      dserror("invalid material for lubrication law %d", curmat->Type());
      break;
  }
}


Teuchos::RCP<MAT::Material> MAT::PAR::LubricationMat::CreateMaterial()
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
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::LubricationMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
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
double MAT::LubricationMat::ComputeViscosityDeriv(const double press, const double visc)
{
  double visc_dp = -1;
  params_->lubricationlaw_->ConstitutiveDerivatives(press, visc, visc_dp);

  return visc_dp;
}

FOUR_C_NAMESPACE_CLOSE
