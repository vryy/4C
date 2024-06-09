/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material for DEM

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "4C_mat_particle_dem.hpp"

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
Mat::PAR::ParticleMaterialDEM::ParticleMaterialDEM(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata), ParticleMaterialBase(matdata)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters  sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<Core::Mat::Material> Mat::PAR::ParticleMaterialDEM::create_material()
{
  return Teuchos::rcp(new Mat::ParticleMaterialDEM(this));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
Core::Communication::ParObject* Mat::ParticleMaterialDEMType::Create(const std::vector<char>& data)
{
  Mat::ParticleMaterialDEM* particlematdem = new Mat::ParticleMaterialDEM();
  particlematdem->Unpack(data);
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
void Mat::ParticleMaterialDEM::Pack(Core::Communication::PackBuffer& data) const
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
 | unpack                                                     sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void Mat::ParticleMaterialDEM::Unpack(const std::vector<char>& data)
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
      // note: dynamic_cast needed due diamond inheritance structure
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<Mat::PAR::ParticleMaterialDEM*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
