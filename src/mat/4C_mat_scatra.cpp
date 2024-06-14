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
Mat::PAR::ScatraMat::ScatraMat(const Core::Mat::PAR::Parameter::Data& matdata) : Parameter(matdata)
{
  // extract relevant communicator
  const Epetra_Comm& comm = Global::Problem::Instance()->Materials()->GetReadFromProblem() == 0
                                ? *Global::Problem::Instance()->GetCommunicators()->LocalComm()
                                : *Global::Problem::Instance()->GetCommunicators()->SubComm();

  Epetra_Map dummy_map(1, 1, 0, comm);
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(diff)->PutScalar(matdata.parameters.Get<double>("DIFFUSIVITY"));
  matparams_.at(reac)->PutScalar(matdata.parameters.Get<double>("REACOEFF"));
  matparams_.at(densific)->PutScalar(matdata.parameters.Get<double>("DENSIFICATION"));
  matparams_.at(reacts_to_external_force)
      ->PutScalar(matdata.parameters.Get<bool>("REACTS_TO_EXTERNAL_FORCE"));
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMat::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMat(this));
}


Mat::ScatraMatType Mat::ScatraMatType::instance_;


Core::Communication::ParObject* Mat::ScatraMatType::Create(const std::vector<char>& data)
{
  Mat::ScatraMat* scatra_mat = new Mat::ScatraMat();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMat::ScatraMat() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ScatraMat::ScatraMat(Mat::PAR::ScatraMat* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ScatraMat::Pack(Core::Communication::PackBuffer& data) const
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
void Mat::ScatraMat::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::ScatraMat*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
