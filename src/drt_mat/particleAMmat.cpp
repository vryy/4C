/*!----------------------------------------------------------------------*/
/*!
\file particleAMmat.cpp

\brief Particle material with support for thermodynamics

\level 3

\maintainer Alessandro Cattabiani
*/
/*----------------------------------------------------------------------*/

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

#include "particleAMmat.H"

/*----------------------------------------------------------------------*
 | constructor                                              catta 06/16 |
 *----------------------------------------------------------------------*/
MAT::PAR::ParticleAMmat::ParticleAMmat(
    Teuchos::RCP<MAT::PAR::Material> matdata
    ) :
    ParticleMat(matdata),
    temperature_(matdata->GetDouble("INITTEMPERATURE")),
    CPS_(matdata->GetDouble("CPS")),
    CPL_(matdata->GetDouble("CPL")),
    SL_latent_heat_max_(matdata->GetDouble("SL_LATENT_HEAT")),
    SL_transitionTemperature_(matdata->GetDouble("SL_TRANSITION_TEMPERATURE"))
{
  return;
}


/*------------------------------------------------------------------*
 | create instance of ParticleAMmat material            catta 06/16 |
 *------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ParticleAMmat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ParticleAMmat(this));
}

MAT::ParticleAMmatType MAT::ParticleAMmatType::instance_;

DRT::ParObject* MAT::ParticleAMmatType::Create(const std::vector<char>& data)
{
  MAT::ParticleAMmat* particleAMmat = new MAT::ParticleAMmat();
  particleAMmat->Unpack(data);
  return particleAMmat;
}


/*------------------------------------------------------------------*
| construct empty ParticleAMmat material                catta 06/16 |
 *------------------------------------------------------------------*/
MAT::ParticleAMmat::ParticleAMmat() :
  params_(NULL)
{
  return;
}


/*--------------------------------------------------------------------------------*
 | construct ParticleAMmat material with specific material parameters catta 06/16 |
 *--------------------------------------------------------------------------------*/
MAT::ParticleAMmat::ParticleAMmat(MAT::PAR::ParticleAMmat* params) :
  ParticleMat(params),
  params_(params)
{
  return;
}


/*----------------------------------------------------------------------*
 | pack material for communication purposes                 catta 06/16 |
 *----------------------------------------------------------------------*/
void MAT::ParticleAMmat::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  int matid = -1;
  if(params_ != NULL)
    matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  // pack base class material
  ParticleMat::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 | unpack data from a char vector                           catta 06/16 |
 *----------------------------------------------------------------------*/
void MAT::ParticleAMmat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if(type != UniqueParObjectId())
    dserror("Wrong instance type data!");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if(DRT::Problem::Instance()->Materials() != Teuchos::null)
    if(DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if(mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ParticleAMmat*>(mat);
      else
        dserror("Type of parameter material %d does not match calling type %d!", mat->Type(), MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ParticleMat::Unpack(basedata);

  // final safety check
  if(position != data.size())
    dserror("Mismatch in size of data %d <-> %d!",data.size(),position);

  return;
}
