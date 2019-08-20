/*----------------------------------------------------------------------*/
/*! \file
\brief class for handling of micro-macro transitions

\level 3

\maintainer Martin Kronbichler
*/
/*----------------------------------------------------------------------*/


#include "micromaterial.H"
#include "matpar_bundle.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"


// Be careful when adding new member functions of MicroMaterial that
// relate to MicroMaterialGP (which is NOT in the filter). See also
// comments in micromaterial_evaluate.cpp which is separated from
// this file for the very reason.


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::MicroMaterial::MicroMaterial(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      microfile_(*(matdata->Get<std::string>("MICROFILE"))),
      microdisnum_(matdata->GetInt("MICRODIS_NUM")),
      initvol_(matdata->GetDouble("INITVOL"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::MicroMaterial::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MicroMaterial(this));
}


MAT::MicroMaterialType MAT::MicroMaterialType::instance_;


DRT::ParObject* MAT::MicroMaterialType::Create(const std::vector<char>& data)
{
  MAT::MicroMaterial* micro = new MAT::MicroMaterial();
  micro->Unpack(data);
  return micro;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MicroMaterial::MicroMaterial() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MicroMaterial::MicroMaterial(MAT::PAR::MicroMaterial* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MicroMaterial::Pack(DRT::PackBuffer& data) const
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
void MAT::MicroMaterial::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::MicroMaterial*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
