/*----------------------------------------------------------------------*/
/*! \file
\brief Contains conductivity, permittivity and permeability of the medium for isotropic
       electromagetic field evolution.
       MAT 1 MAT_Electromagnetic CONDUCTIVITY 0.0 PERMITTIVITY 1.732 PERMEABILITY 1.732


\level 2


*/
/*----------------------------------------------------------------------*/

#include "baci_mat_electromagnetic.hpp"

#include "baci_comm_utils.hpp"
#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElectromagneticMat::ElectromagneticMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  Epetra_Map dummy_map(1, 1, 0, *(GLOBAL::Problem::Instance()->GetCommunicators()->LocalComm()));
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(sigma_)->PutScalar(*matdata->Get<double>("CONDUCTIVITY"));
  matparams_.at(epsilon_)->PutScalar(*matdata->Get<double>("PERMITTIVITY"));
  matparams_.at(mu_)->PutScalar(*matdata->Get<double>("PERMEABILITY"));

  return;
}

Teuchos::RCP<MAT::Material> MAT::PAR::ElectromagneticMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ElectromagneticMat(this));
}

MAT::ElectromagneticMatType MAT::ElectromagneticMatType::instance_;


CORE::COMM::ParObject *MAT::ElectromagneticMatType::Create(const std::vector<char> &data)
{
  MAT::ElectromagneticMat *soundprop = new MAT::ElectromagneticMat();
  soundprop->Unpack(data);
  return soundprop;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ElectromagneticMat::ElectromagneticMat() : params_(nullptr) {}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ElectromagneticMat::ElectromagneticMat(MAT::PAR::ElectromagneticMat *params) : params_(params)
{
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::ElectromagneticMat::Pack(CORE::COMM::PackBuffer &data) const
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

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::ElectromagneticMat::Unpack(const std::vector<char> &data)
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
      MAT::PAR::Parameter *mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ElectromagneticMat *>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
