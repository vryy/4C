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
Mat::PAR::ElectromagneticMat::ElectromagneticMat(const Core::Mat::PAR::Parameter::Data &matdata)
    : Parameter(matdata)
{
  Epetra_Map dummy_map(1, 1, 0, *(Global::Problem::instance()->get_communicators()->local_comm()));
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(sigma_)->PutScalar(matdata.parameters.get<double>("CONDUCTIVITY"));
  matparams_.at(epsilon_)->PutScalar(matdata.parameters.get<double>("PERMITTIVITY"));
  matparams_.at(mu_)->PutScalar(matdata.parameters.get<double>("PERMEABILITY"));

  return;
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ElectromagneticMat::create_material()
{
  return Teuchos::rcp(new Mat::ElectromagneticMat(this));
}

Mat::ElectromagneticMatType Mat::ElectromagneticMatType::instance_;


Core::Communication::ParObject *Mat::ElectromagneticMatType::create(const std::vector<char> &data)
{
  Mat::ElectromagneticMat *soundprop = new Mat::ElectromagneticMat();
  soundprop->unpack(data);
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
void Mat::ElectromagneticMat::pack(Core::Communication::PackBuffer &data) const
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

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void Mat::ElectromagneticMat::unpack(const std::vector<char> &data)
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
      Core::Mat::PAR::Parameter *mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::ElectromagneticMat *>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
