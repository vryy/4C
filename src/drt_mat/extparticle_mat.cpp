/*!----------------------------------------------------------------------*/
/*!
\file extparticle_mat.cpp

\brief Particle material with support for thermodynamics

\level 3

\maintainer Alessandro Cattabiani
*/
/*----------------------------------------------------------------------*/

#include "extparticle_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*
 | constructor                                              catta 06/16 |
 *----------------------------------------------------------------------*/
MAT::PAR::ExtParticleMat::ExtParticleMat(
    Teuchos::RCP<MAT::PAR::Material> matdata
    ) :
    ParticleMat(matdata),
    dismemberRadius_(matdata->GetDouble("DISMEMBER_RADIUS")),
    initTemperature_(matdata->GetDouble("INITTEMPERATURE")),
    CPS_(matdata->GetDouble("CPS")),
    CPL_(matdata->GetDouble("CPL")),
    latentHeatSL_(matdata->GetDouble("LATENT_HEAT_SL")),
    transitionTemperatureSL_(matdata->GetDouble("TRANSITION_TEMPERATURE_SL")),
    thermalExpansionS_(matdata->GetDouble("THERMAL_EXPANSION_S")),
    thermalExpansionL_(matdata->GetDouble("THERMAL_EXPANSION_L")),
    thermalExpansionSL_(matdata->GetDouble("THERMAL_EXPANSION_SL"))
{
  return;
}


/*------------------------------------------------------------------*
 | create instance of ExtParticleMat material            catta 06/16 |
 *------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ExtParticleMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ExtParticleMat(this));
}

MAT::ExtParticleMatType MAT::ExtParticleMatType::instance_;

DRT::ParObject* MAT::ExtParticleMatType::Create(const std::vector<char>& data)
{
  MAT::ExtParticleMat* ExtParticleMat = new MAT::ExtParticleMat();
  ExtParticleMat->Unpack(data);
  return ExtParticleMat;
}


/*------------------------------------------------------------------*
| construct empty ExtParticleMat material                catta 06/16 |
 *------------------------------------------------------------------*/
MAT::ExtParticleMat::ExtParticleMat() :
  params_(NULL)
{
  return;
}


/*--------------------------------------------------------------------------------*
 | construct ExtParticleMat material with specific material parameters catta 06/16 |
 *--------------------------------------------------------------------------------*/
MAT::ExtParticleMat::ExtParticleMat(MAT::PAR::ExtParticleMat* params) :
  ParticleMat(params),
  params_(params)
{
  return;
}


/*----------------------------------------------------------------------*
 | pack material for communication purposes                 catta 06/16 |
 *----------------------------------------------------------------------*/
void MAT::ExtParticleMat::Pack(DRT::PackBuffer& data) const
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
void MAT::ExtParticleMat::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ExtParticleMat*>(mat);
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
