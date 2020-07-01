/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material for SPH fluid

\level 3


*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_material_sph_fluid.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

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
      refDensFac_(matdata->GetDouble("REFDENSFAC")),
      exponent_(matdata->GetDouble("EXPONENT")),
      backgroundPressure_(matdata->GetDouble("BACKGROUNDPRESSURE")),
      bulkModulus_(matdata->GetDouble("BULK_MODULUS")),
      dynamicViscosity_(matdata->GetDouble("DYNAMIC_VISCOSITY")),
      bulkViscosity_(matdata->GetDouble("BULK_VISCOSITY")),
      artificialViscosity_(matdata->GetDouble("ARTIFICIAL_VISCOSITY"))
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
DRT::ParObject* MAT::ParticleMaterialSPHFluidType::Create(const std::vector<char>& data)
{
  MAT::ParticleMaterialSPHFluid* particlematsph = new MAT::ParticleMaterialSPHFluid();
  particlematsph->Unpack(data);
  return particlematsph;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                        sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
MAT::ParticleMaterialSPHFluid::ParticleMaterialSPHFluid() : params_(NULL)
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
void MAT::ParticleMaterialSPHFluid::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*---------------------------------------------------------------------------*
 | unpack                                                     sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void MAT::ParticleMaterialSPHFluid::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      // note: dynamic_cast needed due diamond inheritance structure
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = dynamic_cast<MAT::PAR::ParticleMaterialSPHFluid*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
