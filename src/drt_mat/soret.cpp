/*!----------------------------------------------------------------------*/
/*!
\brief material for heat transport due to Fourier-type thermal conduction and the Soret effect

\maintainer Christoph Schmidt

\level 2
*/
/*----------------------------------------------------------------------*/
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

#include "soret.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 06/15 |
 *----------------------------------------------------------------------*/
MAT::PAR::Soret::Soret(Teuchos::RCP<MAT::PAR::Material> matdata)
    : FourierIso(matdata), soretcoefficient_(matdata->GetDouble("SORET"))
{
  return;
}


/*------------------------------------------------------------------*
 | create instance of Soret material                     fang 06/15 |
 *------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Soret::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Soret(this));
}

MAT::SoretType MAT::SoretType::instance_;

DRT::ParObject* MAT::SoretType::Create(const std::vector<char>& data)
{
  MAT::Soret* soret = new MAT::Soret();
  soret->Unpack(data);
  return soret;
}


/*------------------------------------------------------------------*
 | construct empty Soret material                        fang 06/15 |
 *------------------------------------------------------------------*/
MAT::Soret::Soret() : params_(NULL) { return; }


/*-------------------------------------------------------------------------*
 | construct Soret material with specific material parameters   fang 06/15 |
 *-------------------------------------------------------------------------*/
MAT::Soret::Soret(MAT::PAR::Soret* params) : FourierIso(params), params_(params) { return; }


/*----------------------------------------------------------------------*
 | pack material for communication purposes                  fang 06/15 |
 *----------------------------------------------------------------------*/
void MAT::Soret::Pack(DRT::PackBuffer& data) const
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
  FourierIso::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 | unpack data from a char vector                            fang 06/15 |
 *----------------------------------------------------------------------*/
void MAT::Soret::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::Soret*>(mat);
      else
        dserror("Type of parameter material %d does not match calling type %d!", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  FourierIso::Unpack(basedata);

  // final safety check
  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d!", data.size(), position);

  return;
}
