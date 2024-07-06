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
  const Epetra_Comm& comm = Global::Problem::instance()->materials()->get_read_from_problem() == 0
                                ? *Global::Problem::instance()->get_communicators()->local_comm()
                                : *Global::Problem::instance()->get_communicators()->sub_comm();

  Epetra_Map dummy_map(1, 1, 0, comm);
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(diff)->PutScalar(matdata.parameters.get<double>("DIFFUSIVITY"));
  matparams_.at(reac)->PutScalar(matdata.parameters.get<double>("REACOEFF"));
  matparams_.at(densific)->PutScalar(matdata.parameters.get<double>("DENSIFICATION"));
  matparams_.at(reacts_to_external_force)
      ->PutScalar(matdata.parameters.get<bool>("REACTS_TO_EXTERNAL_FORCE"));
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::ScatraMat::create_material()
{
  return Teuchos::rcp(new Mat::ScatraMat(this));
}


Mat::ScatraMatType Mat::ScatraMatType::instance_;


Core::Communication::ParObject* Mat::ScatraMatType::create(const std::vector<char>& data)
{
  Mat::ScatraMat* scatra_mat = new Mat::ScatraMat();
  scatra_mat->unpack(data);
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
void Mat::ScatraMat::pack(Core::Communication::PackBuffer& data) const
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
void Mat::ScatraMat::unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::ScatraMat*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
