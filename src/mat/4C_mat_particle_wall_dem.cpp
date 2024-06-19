/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall material for DEM

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
#include "4C_mat_particle_wall_dem.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | define static class member                                 sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
Mat::ParticleWallMaterialDEMType Mat::ParticleWallMaterialDEMType::instance_;

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
Mat::PAR::ParticleWallMaterialDEM::ParticleWallMaterialDEM(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      frictionTang_(matdata.parameters.get<double>("FRICT_COEFF_TANG")),
      frictionRoll_(matdata.parameters.get<double>("FRICT_COEFF_ROLL")),
      adhesionSurfaceEnergy_(matdata.parameters.get<double>("ADHESION_SURFACE_ENERGY"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters  sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::ParticleWallMaterialDEM::create_material()
{
  return Teuchos::rcp(new Mat::ParticleWallMaterialDEM(this));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::ParticleWallMaterialDEMType::Create(
    const std::vector<char>& data)
{
  Mat::ParticleWallMaterialDEM* particlewallmatdem = new Mat::ParticleWallMaterialDEM();
  particlewallmatdem->unpack(data);
  return particlewallmatdem;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                        sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
Mat::ParticleWallMaterialDEM::ParticleWallMaterialDEM() : params_(nullptr)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | constructor (with given material parameters)               sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
Mat::ParticleWallMaterialDEM::ParticleWallMaterialDEM(Mat::PAR::ParticleWallMaterialDEM* params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | pack                                                       sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void Mat::ParticleWallMaterialDEM::pack(Core::Communication::PackBuffer& data) const
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

/*---------------------------------------------------------------------------*
 | unpack                                                     sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void Mat::ParticleWallMaterialDEM::unpack(const std::vector<char>& data)
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
        params_ = static_cast<Mat::PAR::ParticleWallMaterialDEM*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
