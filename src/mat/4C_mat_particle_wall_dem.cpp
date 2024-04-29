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
MAT::ParticleWallMaterialDEMType MAT::ParticleWallMaterialDEMType::instance_;

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
MAT::PAR::ParticleWallMaterialDEM::ParticleWallMaterialDEM(
    Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata),
      frictionTang_(matdata->Get<double>("FRICT_COEFF_TANG")),
      frictionRoll_(matdata->Get<double>("FRICT_COEFF_ROLL")),
      adhesionSurfaceEnergy_(matdata->Get<double>("ADHESION_SURFACE_ENERGY"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters  sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<CORE::MAT::Material> MAT::PAR::ParticleWallMaterialDEM::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ParticleWallMaterialDEM(this));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
CORE::COMM::ParObject* MAT::ParticleWallMaterialDEMType::Create(const std::vector<char>& data)
{
  MAT::ParticleWallMaterialDEM* particlewallmatdem = new MAT::ParticleWallMaterialDEM();
  particlewallmatdem->Unpack(data);
  return particlewallmatdem;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                        sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
MAT::ParticleWallMaterialDEM::ParticleWallMaterialDEM() : params_(nullptr)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | constructor (with given material parameters)               sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
MAT::ParticleWallMaterialDEM::ParticleWallMaterialDEM(MAT::PAR::ParticleWallMaterialDEM* params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | pack                                                       sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void MAT::ParticleWallMaterialDEM::Pack(CORE::COMM::PackBuffer& data) const
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

/*---------------------------------------------------------------------------*
 | unpack                                                     sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void MAT::ParticleWallMaterialDEM::Unpack(const std::vector<char>& data)
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
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ParticleWallMaterialDEM*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
