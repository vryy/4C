/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall material for DEM

\level 3

\maintainer  Sebastian Fuchs

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
#include "particle_wall_material_dem.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*---------------------------------------------------------------------------*
 | define static class member                                 sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
MAT::ParticleWallMaterialDEMType MAT::ParticleWallMaterialDEMType::instance_;

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
MAT::PAR::ParticleWallMaterialDEM::ParticleWallMaterialDEM(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      frictionTang_(matdata->GetDouble("FRICT_COEFF_TANG")),
      frictionRoll_(matdata->GetDouble("FRICT_COEFF_ROLL")),
      adhesionSurfaceEnergy_(matdata->GetDouble("ADHESION_SURFACE_ENERGY"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters  sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ParticleWallMaterialDEM::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ParticleWallMaterialDEM(this));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
DRT::ParObject* MAT::ParticleWallMaterialDEMType::Create(const std::vector<char>& data)
{
  MAT::ParticleWallMaterialDEM* particlewallmatdem = new MAT::ParticleWallMaterialDEM();
  particlewallmatdem->Unpack(data);
  return particlewallmatdem;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                        sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
MAT::ParticleWallMaterialDEM::ParticleWallMaterialDEM() : params_(NULL)
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
void MAT::ParticleWallMaterialDEM::Pack(DRT::PackBuffer& data) const
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
 | unpack                                                     sfuchs 08/2019 |
 *---------------------------------------------------------------------------*/
void MAT::ParticleWallMaterialDEM::Unpack(const std::vector<char>& data)
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
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ParticleWallMaterialDEM*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
