// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_electromagnetic.hpp"

#include "4C_comm_pack_helpers.hpp"
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
    matparams_.push_back(Teuchos::make_rcp<Core::LinAlg::Vector<double>>(dummy_map, true));
  }
  matparams_.at(sigma_)->PutScalar(matdata.parameters.get<double>("CONDUCTIVITY"));
  matparams_.at(epsilon_)->PutScalar(matdata.parameters.get<double>("PERMITTIVITY"));
  matparams_.at(mu_)->PutScalar(matdata.parameters.get<double>("PERMEABILITY"));

  return;
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ElectromagneticMat::create_material()
{
  return Teuchos::make_rcp<Mat::ElectromagneticMat>(this);
}

Mat::ElectromagneticMatType Mat::ElectromagneticMatType::instance_;


Core::Communication::ParObject *Mat::ElectromagneticMatType::create(
    Core::Communication::UnpackBuffer &buffer)
{
  Mat::ElectromagneticMat *soundprop = new Mat::ElectromagneticMat();
  soundprop->unpack(buffer);
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
void Mat::ElectromagneticMat::unpack(Core::Communication::UnpackBuffer &buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // matid and recover params_
  int matid;
  extract_from_pack(buffer, matid);
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

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

FOUR_C_NAMESPACE_CLOSE
