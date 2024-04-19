/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material for SPH fluid

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "baci_mat_particle_sph_fluid.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | define static class member                                 sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::ParticleMaterialSPHFluidType MAT::ParticleMaterialSPHFluidType::instance_;

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::PAR::ParticleMaterialSPHFluid::ParticleMaterialSPHFluid(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      ParticleMaterialBase(matdata),
      ParticleMaterialThermo(matdata),
      refDensFac_(*matdata->Get<double>("REFDENSFAC")),
      exponent_(*matdata->Get<double>("EXPONENT")),
      backgroundPressure_(*matdata->Get<double>("BACKGROUNDPRESSURE")),
      bulkModulus_(*matdata->Get<double>("BULK_MODULUS")),
      dynamicViscosity_(*matdata->Get<double>("DYNAMIC_VISCOSITY")),
      bulkViscosity_(*matdata->Get<double>("BULK_VISCOSITY")),
      artificialViscosity_(*matdata->Get<double>("ARTIFICIAL_VISCOSITY"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters  sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ParticleMaterialSPHFluid::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ParticleMaterialSPHFluid(this));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
CORE::COMM::ParObject* MAT::ParticleMaterialSPHFluidType::Create(const std::vector<char>& data)
{
  MAT::ParticleMaterialSPHFluid* particlematsph = new MAT::ParticleMaterialSPHFluid();
  particlematsph->Unpack(data);
  return particlematsph;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                        sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::ParticleMaterialSPHFluid::ParticleMaterialSPHFluid() : params_(nullptr)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | constructor (with given material parameters)               sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::ParticleMaterialSPHFluid::ParticleMaterialSPHFluid(MAT::PAR::ParticleMaterialSPHFluid* params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | pack                                                       sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void MAT::ParticleMaterialSPHFluid::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::ParticleMaterialSPHFluid::Unpack(const std::vector<char>& data)
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
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<MAT::PAR::ParticleMaterialSPHFluid*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
