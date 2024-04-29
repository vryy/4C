/*----------------------------------------------------------------------*/
/*! \file
\brief scalar transport material

\level 1

*/
/*----------------------------------------------------------------------*/


#include "4C_mat_scatra.hpp"

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraMat::ScatraMat(Teuchos::RCP<MAT::PAR::Material> matdata) : Parameter(matdata)
{
  // extract relevant communicator
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem() == 0
                                ? *GLOBAL::Problem::Instance()->GetCommunicators()->LocalComm()
                                : *GLOBAL::Problem::Instance()->GetCommunicators()->SubComm();

  Epetra_Map dummy_map(1, 1, 0, comm);
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(diff)->PutScalar(matdata->Get<double>("DIFFUSIVITY"));
  matparams_.at(reac)->PutScalar(matdata->Get<double>("REACOEFF"));
  matparams_.at(densific)->PutScalar(matdata->Get<double>("DENSIFICATION"));
  matparams_.at(reacts_to_external_force)
      ->PutScalar(matdata->Get<bool>("REACTS_TO_EXTERNAL_FORCE"));
}


Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMat(this));
}


MAT::ScatraMatType MAT::ScatraMatType::instance_;


CORE::COMM::ParObject* MAT::ScatraMatType::Create(const std::vector<char>& data)
{
  MAT::ScatraMat* scatra_mat = new MAT::ScatraMat();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMat::ScatraMat() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMat::ScatraMat(MAT::PAR::ScatraMat* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMat::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::ScatraMat::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ScatraMat*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
