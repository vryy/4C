/*----------------------------------------------------------------------*/
/*! \file
 \brief scatra material for transport within porous model with special implementations
        for ECM model


\level 3
 *----------------------------------------------------------------------*/


#include "4C_mat_scatra_poro_ecm.hpp"

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ScatraMatPoroECM::ScatraMatPoroECM(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : ScatraReactionMat(matdata), reacscale_(matdata->Get<double>("REACSCALE"))
{
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMatPoroECM::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMatPoroECM(this));
}


Mat::ScatraMatPoroECMType Mat::ScatraMatPoroECMType::instance_;

Core::Communication::ParObject* Mat::ScatraMatPoroECMType::Create(const std::vector<char>& data)
{
  Mat::ScatraMatPoroECM* scatra_mat = new Mat::ScatraMatPoroECM();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMatPoroECM::ScatraMatPoroECM() : params_(nullptr), reaccoeff_(0.0) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMatPoroECM::ScatraMatPoroECM(Mat::PAR::ScatraMatPoroECM* params)
    : ScatraReactionMat(params), params_(params), reaccoeff_(0.0)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatPoroECM::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // reaccoeff_
  add_to_pack(data, reaccoeff_);

  // add base class material
  ScatraReactionMat::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatPoroECM::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::ScatraMatPoroECM*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // reaccoeff_
  extract_from_pack(position, data, reaccoeff_);

  // extract base class material
  std::vector<char> basedata(0);
  extract_from_pack(position, data, basedata);
  ScatraReactionMat::Unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMatPoroECM::ComputeReacCoeff(double chempot)
{
  reaccoeff_ = params_->reaccoeff_ * exp(params_->reacscale_ * chempot);
  return;
}

FOUR_C_NAMESPACE_CLOSE
