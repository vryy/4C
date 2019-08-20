/*----------------------------------------------------------------------*/
/*! \file
\brief material properties of ellipsoidal particles

\level 2

\maintainer Sebastian Fuchs
*/
/*----------------------------------------------------------------------*/
#include "particle_mat_ellipsoids.H"

#include "matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"

/*--------------------------------------------------------------------*
 | constructor                                             fang 09/17 |
 *--------------------------------------------------------------------*/
MAT::PAR::ParticleMatEllipsoids::ParticleMatEllipsoids(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ParticleMat(matdata), semiaxes_(true)
{
  // extract semi-axes
  const std::vector<double>& semiaxes = *matdata->Get<std::vector<double>>("SEMI-AXES");

  // safety check
  if (semiaxes.size() != 3) dserror("Must have exactly three semi-axes!");

  // store semi-axes
  semiaxes_.SetCopy(&semiaxes[0]);

  return;
}


/*--------------------------------------------------------------------*
 | create instance of particle material                    fang 09/17 |
 *--------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ParticleMatEllipsoids::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ParticleMatEllipsoids(this));
}


MAT::ParticleMatEllipsoidsType MAT::ParticleMatEllipsoidsType::instance_;


/*--------------------------------------------------------------------*
 | unpack instance of particle material                    fang 09/17 |
 *--------------------------------------------------------------------*/
DRT::ParObject* MAT::ParticleMatEllipsoidsType::Create(const std::vector<char>& data)
{
  MAT::ParticleMatEllipsoids* ParticleMatEllipsoids = new MAT::ParticleMatEllipsoids();
  ParticleMatEllipsoids->Unpack(data);
  return ParticleMatEllipsoids;
}


/*--------------------------------------------------------------------*
 | construct empty particle material                       fang 09/17 |
 *--------------------------------------------------------------------*/
MAT::ParticleMatEllipsoids::ParticleMatEllipsoids() : params_(NULL) { return; }


/*--------------------------------------------------------------------------------------*
 | construct particle material with specific material parameters             fang 09/17 |
 *--------------------------------------------------------------------------------------*/
MAT::ParticleMatEllipsoids::ParticleMatEllipsoids(MAT::PAR::ParticleMatEllipsoids* params)
    : ParticleMat(params), params_(params)
{
  return;
}


/*--------------------------------------------------------------------*
 | pack material for communication purposes                fang 09/17 |
 *--------------------------------------------------------------------*/
void MAT::ParticleMatEllipsoids::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  // pack base class material
  ParticleMat::Pack(data);

  return;
}


/*--------------------------------------------------------------------*
 | unpack data from a char vector                          fang 09/17 |
 *--------------------------------------------------------------------*/
void MAT::ParticleMatEllipsoids::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("Wrong instance type data!");

  // extract material ID
  int matid;
  ExtractfromPack(position, data, matid);

  // recover material parameters
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ParticleMatEllipsoids*>(mat);
      else
        dserror("Type of parameter material %d does not match calling type %d!", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ParticleMat::Unpack(basedata);

  // final safety check
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d!", data.size(), position);

  return;
}
