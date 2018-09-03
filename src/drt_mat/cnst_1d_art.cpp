/*----------------------------------------------------------------------*/
/*!
\file cnst_1d_art.cpp

\maintainer Lena Yoshihara

\level 3
*/
/*----------------------------------------------------------------------*/


#include <vector>
#include "cnst_1d_art.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Cnst_1d_art::Cnst_1d_art(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      viscosity_(matdata->GetDouble("VISCOSITY")),
      density_(matdata->GetDouble("DENS")),
      young_(matdata->GetDouble("YOUNG")),
      nue_(matdata->GetDouble("NUE")),
      diam_(matdata->GetDouble("DIAM")),
      th_(matdata->GetDouble("TH")),
      pext1_(matdata->GetDouble("PEXT1")),
      pext2_(matdata->GetDouble("PEXT2"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Cnst_1d_art::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Cnst_1d_art(this));
}


MAT::Cnst_1d_artType MAT::Cnst_1d_artType::instance_;


DRT::ParObject* MAT::Cnst_1d_artType::Create(const std::vector<char>& data)
{
  MAT::Cnst_1d_art* cnst_art = new MAT::Cnst_1d_art();
  cnst_art->Unpack(data);
  return cnst_art;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Cnst_1d_art::Cnst_1d_art() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Cnst_1d_art::Cnst_1d_art(MAT::PAR::Cnst_1d_art* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Cnst_1d_art::Pack(DRT::PackBuffer& data) const
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
void MAT::Cnst_1d_art::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
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
        params_ = static_cast<MAT::PAR::Cnst_1d_art*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
