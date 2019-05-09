/*---------------------------------------------------------------------------*/
/*!
\brief particle material for DEM

\level 3

\maintainer  Sebastian Fuchs

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_material_dem.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*---------------------------------------------------------------------------*
 | define static class member                                 sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
MAT::ParticleMaterialDEMType MAT::ParticleMaterialDEMType::instance_;

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
MAT::PAR::ParticleMaterialDEM::ParticleMaterialDEM(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), ParticleMaterialBase(matdata)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | create material instance of matching type with parameters  sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ParticleMaterialDEM::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ParticleMaterialDEM(this));
}

/*---------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
DRT::ParObject* MAT::ParticleMaterialDEMType::Create(const std::vector<char>& data)
{
  MAT::ParticleMaterialDEM* particlematsph = new MAT::ParticleMaterialDEM();
  particlematsph->Unpack(data);
  return particlematsph;
}

/*---------------------------------------------------------------------------*
 | constructor (empty material object)                        sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
MAT::ParticleMaterialDEM::ParticleMaterialDEM() : params_(NULL)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | constructor (with given material parameters)               sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
MAT::ParticleMaterialDEM::ParticleMaterialDEM(MAT::PAR::ParticleMaterialDEM* params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | pack                                                       sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void MAT::ParticleMaterialDEM::Pack(DRT::PackBuffer& data) const
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
 | unpack                                                     sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void MAT::ParticleMaterialDEM::Unpack(const std::vector<char>& data)
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
        params_ = dynamic_cast<MAT::PAR::ParticleMaterialDEM*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
