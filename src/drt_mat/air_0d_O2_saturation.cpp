/*----------------------------------------------------------------------*/
/*! \file
\brief
Former file of Lena Yoshihara

\maintainer Martin Kronbichler

\level 3
*/
/*----------------------------------------------------------------------*/


#include <vector>
#include "air_0d_O2_saturation.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Air_0d_O2_saturation::Air_0d_O2_saturation(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      atmospheric_p_(matdata->GetDouble("AtmosphericPressure")),
      nO2_per_VO2_(matdata->GetDouble("NumberOfO2PerVO2"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Air_0d_O2_saturation::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Air_0d_O2_saturation(this));
}


MAT::Air_0d_O2_saturationType MAT::Air_0d_O2_saturationType::instance_;


DRT::ParObject* MAT::Air_0d_O2_saturationType::Create(const std::vector<char>& data)
{
  MAT::Air_0d_O2_saturation* air_0d_O2_sat = new MAT::Air_0d_O2_saturation();
  air_0d_O2_sat->Unpack(data);
  return air_0d_O2_sat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Air_0d_O2_saturation::Air_0d_O2_saturation() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Air_0d_O2_saturation::Air_0d_O2_saturation(MAT::PAR::Air_0d_O2_saturation* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Air_0d_O2_saturation::Pack(DRT::PackBuffer& data) const
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
void MAT::Air_0d_O2_saturation::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::Air_0d_O2_saturation*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
