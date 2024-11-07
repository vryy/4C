// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_particle_dem.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | define static class member                                 sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialDEMType Mat::ParticleMaterialDEMType::instance_;

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
Mat::PAR::ParticleMaterialDEM::ParticleMaterialDEM(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata), ParticleMaterialBase(matdata)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters  sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
std::shared_ptr<Core::Mat::Material> Mat::PAR::ParticleMaterialDEM::create_material()
{
  return std::make_shared<Mat::ParticleMaterialDEM>(this);
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::ParticleMaterialDEMType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::ParticleMaterialDEM* particlematdem = new Mat::ParticleMaterialDEM();
  particlematdem->unpack(buffer);
  return particlematdem;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                        sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialDEM::ParticleMaterialDEM() : params_(nullptr)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | constructor (with given material parameters)               sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialDEM::ParticleMaterialDEM(Mat::PAR::ParticleMaterialDEM* params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | pack                                                       sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void Mat::ParticleMaterialDEM::pack(Core::Communication::PackBuffer& data) const
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

/*---------------------------------------------------------------------------*
 | unpack                                                     sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void Mat::ParticleMaterialDEM::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // matid and recover params_
  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != nullptr)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      // note: dynamic_cast needed due diamond inheritance structure
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = dynamic_cast<Mat::PAR::ParticleMaterialDEM*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

FOUR_C_NAMESPACE_CLOSE
