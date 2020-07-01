/*----------------------------------------------------------------------*/
/*! \file
\brief material according to mixture-fraction approach

\level 2

*----------------------------------------------------------------------*/


#include <vector>
#include "mixfrac.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::MixFrac::MixFrac(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      kinvisc_(matdata->GetDouble("KINVISC")),
      kindiff_(matdata->GetDouble("KINDIFF")),
      eosfaca_(matdata->GetDouble("EOSFACA")),
      eosfacb_(matdata->GetDouble("EOSFACB"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::MixFrac::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MixFrac(this));
}


MAT::MixFracType MAT::MixFracType::instance_;


DRT::ParObject* MAT::MixFracType::Create(const std::vector<char>& data)
{
  MAT::MixFrac* mixfrac = new MAT::MixFrac();
  mixfrac->Unpack(data);
  return mixfrac;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MixFrac::MixFrac() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MixFrac::MixFrac(MAT::PAR::MixFrac* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MixFrac::Pack(DRT::PackBuffer& data) const
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
void MAT::MixFrac::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::MixFrac*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::MixFrac::ComputeViscosity(const double mixfrac) const
{
  const double visc = KinVisc() / (EosFacA() * mixfrac + EosFacB());

  return visc;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::MixFrac::ComputeDiffusivity(const double mixfrac) const
{
  const double diffus = KinDiff() / (EosFacA() * mixfrac + EosFacB());

  return diffus;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::MixFrac::ComputeDensity(const double mixfrac) const
{
  const double density = 1.0 / (EosFacA() * mixfrac + EosFacB());

  return density;
}
