/*----------------------------------------------------------------------*/
/*! \file
\brief material for macro-scale elements in multi-scale simulations of scalar transport problems

\level 2

\maintainer Christoph Schmidt
*/
/*----------------------------------------------------------------------*/
#include "scatra_mat_multiscale.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/matpar_bundle.H"

/*--------------------------------------------------------------------*
 | constructor                                             fang 11/15 |
 *--------------------------------------------------------------------*/
MAT::PAR::ScatraMatMultiScale::ScatraMatMultiScale(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ScatraMat(matdata),
      ScatraMultiScale(matdata),
      porosity_(matdata->GetDouble("POROSITY")),
      tortuosity_(matdata->GetDouble("TORTUOSITY"))
{
  return;
}


/*--------------------------------------------------------------------*
 | create material                                         fang 11/15 |
 *--------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMatMultiScale::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMatMultiScale(this));
}


MAT::ScatraMatMultiScaleType MAT::ScatraMatMultiScaleType::instance_;


DRT::ParObject* MAT::ScatraMatMultiScaleType::Create(const std::vector<char>& data)
{
  MAT::ScatraMatMultiScale* ScatraMatMultiScale = new MAT::ScatraMatMultiScale();
  ScatraMatMultiScale->Unpack(data);
  return ScatraMatMultiScale;
}


/*--------------------------------------------------------------------*
 | construct empty material                                fang 11/15 |
 *--------------------------------------------------------------------*/
MAT::ScatraMatMultiScale::ScatraMatMultiScale() : params_(NULL) { return; }


/*--------------------------------------------------------------------*
 | construct material with specific material parameters    fang 11/15 |
 *--------------------------------------------------------------------*/
MAT::ScatraMatMultiScale::ScatraMatMultiScale(MAT::PAR::ScatraMatMultiScale* params)
    : ScatraMat(params), params_(params)
{
  return;
}


/*--------------------------------------------------------------------*
 | pack material for communication purposes                fang 11/15 |
 *--------------------------------------------------------------------*/
void MAT::ScatraMatMultiScale::Pack(DRT::PackBuffer& data) const
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
  ScatraMat::Pack(data);

  return;
}


/*--------------------------------------------------------------------*
 | unpack data from a char vector                          fang 11/15 |
 *--------------------------------------------------------------------*/
void MAT::ScatraMatMultiScale::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ScatraMatMultiScale*>(mat);
      else
        dserror("Type of parameter material %d does not match calling type %d!", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraMat::Unpack(basedata);

  // final safety check
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d!", data.size(), position);

  return;
}
