/*----------------------------------------------------------------------*/
/*! \file
\brief material properties of particles

\level 2

\maintainer  Christoph Meier

*-----------------------------------------------------------------------*/


#include <vector>
#include "particle_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ParticleMat::ParticleMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      initDensity_(matdata->GetDouble("DENSITY")),
      initRadius_(matdata->GetDouble("INITRADIUS")),
      poissonRatio_(matdata->GetDouble("NUE")),
      youngModulus_(matdata->GetDouble("YOUNG")),
      yieldStrength_(matdata->GetDouble("YIELD"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::ParticleMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ParticleMat(this));
}

MAT::ParticleMatType MAT::ParticleMatType::instance_;

DRT::ParObject* MAT::ParticleMatType::Create(const std::vector<char>& data)
{
  MAT::ParticleMat* particlemat = new MAT::ParticleMat();
  particlemat->Unpack(data);
  return particlemat;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ParticleMat::ParticleMat() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ParticleMat::ParticleMat(MAT::PAR::ParticleMat* params) : params_(params) {}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ParticleMat::Pack(DRT::PackBuffer& data) const
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ParticleMat::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ParticleMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
