/*----------------------------------------------------------------------*/
/*! \file
\brief Contains conductivity, permittivity and permeability of the medium for isotropic
       electromagetic field evolution.
       MAT 1 MAT_Electromagnetic CONDUCTIVITY 0.0 PERMITTIVITY 1.732 PERMEABILITY 1.732


\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_mat_electromagnetic.hpp"

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ElectromagneticMat::ElectromagneticMat(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata)
{
  Epetra_Map dummy_map(1, 1, 0, *(Global::Problem::Instance()->GetCommunicators()->LocalComm()));
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(sigma_)->PutScalar(matdata->Get<double>("CONDUCTIVITY"));
  matparams_.at(epsilon_)->PutScalar(matdata->Get<double>("PERMITTIVITY"));
  matparams_.at(mu_)->PutScalar(matdata->Get<double>("PERMEABILITY"));

  return;
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ElectromagneticMat::create_material()
{
  return Teuchos::rcp(new Mat::ElectromagneticMat(this));
}

Mat::ElectromagneticMatType Mat::ElectromagneticMatType::instance_;


Core::Communication::ParObject *Mat::ElectromagneticMatType::Create(const std::vector<char> &data)
{
  Mat::ElectromagneticMat *soundprop = new Mat::ElectromagneticMat();
  soundprop->Unpack(data);
  return soundprop;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
Mat::ElectromagneticMat::ElectromagneticMat() : params_(nullptr) {}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
Mat::ElectromagneticMat::ElectromagneticMat(Mat::PAR::ElectromagneticMat *params) : params_(params)
{
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void Mat::ElectromagneticMat::Pack(Core::Communication::PackBuffer &data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void Mat::ElectromagneticMat::Unpack(const std::vector<char> &data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter *mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::ElectromagneticMat *>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
