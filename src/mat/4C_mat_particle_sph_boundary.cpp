/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material for SPH boundary

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "4C_mat_particle_sph_boundary.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | define static class member                                 sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialSPHBoundaryType Mat::ParticleMaterialSPHBoundaryType::instance_;

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::PAR::ParticleMaterialSPHBoundary::ParticleMaterialSPHBoundary(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata), ParticleMaterialBase(matdata), ParticleMaterialThermo(matdata)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters  sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::ParticleMaterialSPHBoundary::create_material()
{
  return Teuchos::rcp(new Mat::ParticleMaterialSPHBoundary(this));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::ParticleMaterialSPHBoundaryType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::ParticleMaterialSPHBoundary* particlematsphboundary = new Mat::ParticleMaterialSPHBoundary();
  particlematsphboundary->unpack(buffer);
  return particlematsphboundary;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                        sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialSPHBoundary::ParticleMaterialSPHBoundary() : params_(nullptr)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | constructor (with given material parameters)               sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Mat::ParticleMaterialSPHBoundary::ParticleMaterialSPHBoundary(
    Mat::PAR::ParticleMaterialSPHBoundary* params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | pack                                                       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void Mat::ParticleMaterialSPHBoundary::pack(Core::Communication::PackBuffer& data) const
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
 | unpack                                                     sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void Mat::ParticleMaterialSPHBoundary::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // matid and recover params_
  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      // note: dynamic_cast needed due diamond inheritance structure
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = dynamic_cast<Mat::PAR::ParticleMaterialSPHBoundary*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  FOUR_C_THROW_UNLESS(buffer.at_end(), "Buffer not fully consumed.");
}

FOUR_C_NAMESPACE_CLOSE
