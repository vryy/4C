/*----------------------------------------------------------------------*/
/*!
\brief material for macro-scale elements in multi-scale simulations of electrochemistry problems

\level 2

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/
#include "newman_multiscale.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/matpar_bundle.H"

/*--------------------------------------------------------------------*
 | constructor                                             fang 07/17 |
 *--------------------------------------------------------------------*/
MAT::PAR::NewmanMultiScale::NewmanMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Newman(matdata), ScatraMultiScale(matdata), sigma_(matdata->GetDouble("SIGMA"))
{
  return;
}


/*--------------------------------------------------------------------*
 | create instance of Newman multi-scale material          fang 07/17 |
 *--------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::NewmanMultiScale::CreateMaterial()
{
  return Teuchos::rcp(new MAT::NewmanMultiScale(this));
}


MAT::NewmanMultiScaleType MAT::NewmanMultiScaleType::instance_;


/*--------------------------------------------------------------------*
 | unpack instance of Newman multi-scale material          fang 07/17 |
 *--------------------------------------------------------------------*/
DRT::ParObject* MAT::NewmanMultiScaleType::Create(const std::vector<char>& data)
{
  MAT::NewmanMultiScale* NewmanMultiScale = new MAT::NewmanMultiScale();
  NewmanMultiScale->Unpack(data);
  return NewmanMultiScale;
}


/*--------------------------------------------------------------------*
 | construct empty Newman multi-scale material             fang 07/17 |
 *--------------------------------------------------------------------*/
MAT::NewmanMultiScale::NewmanMultiScale() : params_(NULL) { return; }


/*--------------------------------------------------------------------------------------*
 | construct Newman multi-scale material with specific material parameters   fang 07/17 |
 *--------------------------------------------------------------------------------------*/
MAT::NewmanMultiScale::NewmanMultiScale(MAT::PAR::NewmanMultiScale* params)
    : Newman(params), params_(params)
{
  return;
}


/*--------------------------------------------------------------------*
 | pack material for communication purposes                fang 07/17 |
 *--------------------------------------------------------------------*/
void MAT::NewmanMultiScale::Pack(DRT::PackBuffer& data) const
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
  Newman::Pack(data);

  return;
}


/*--------------------------------------------------------------------*
 | unpack data from a char vector                          fang 07/17 |
 *--------------------------------------------------------------------*/
void MAT::NewmanMultiScale::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("Wrong instance type data!");

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
        params_ = static_cast<MAT::PAR::NewmanMultiScale*>(mat);
      else
        dserror("Type of parameter material %d does not match calling type %d!", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  Newman::Unpack(basedata);

  // final safety check
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d!", data.size(), position);

  return;
}
