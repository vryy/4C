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

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | define static class member                                 sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::ParticleMaterialSPHBoundaryType MAT::ParticleMaterialSPHBoundaryType::instance_;

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::PAR::ParticleMaterialSPHBoundary::ParticleMaterialSPHBoundary(
    Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata), ParticleMaterialBase(matdata), ParticleMaterialThermo(matdata)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters  sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<CORE::MAT::Material> MAT::PAR::ParticleMaterialSPHBoundary::create_material()
{
  return Teuchos::rcp(new MAT::ParticleMaterialSPHBoundary(this));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
CORE::COMM::ParObject* MAT::ParticleMaterialSPHBoundaryType::Create(const std::vector<char>& data)
{
  MAT::ParticleMaterialSPHBoundary* particlematsphboundary = new MAT::ParticleMaterialSPHBoundary();
  particlematsphboundary->Unpack(data);
  return particlematsphboundary;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                        sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::ParticleMaterialSPHBoundary::ParticleMaterialSPHBoundary() : params_(nullptr)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | constructor (with given material parameters)               sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::ParticleMaterialSPHBoundary::ParticleMaterialSPHBoundary(
    MAT::PAR::ParticleMaterialSPHBoundary* params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | pack                                                       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void MAT::ParticleMaterialSPHBoundary::Pack(CORE::COMM::PackBuffer& data) const
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
 | unpack                                                     sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void MAT::ParticleMaterialSPHBoundary::Unpack(const std::vector<char>& data)
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
      // note: dynamic_cast needed due diamond inheritance structure
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<MAT::PAR::ParticleMaterialSPHBoundary*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
