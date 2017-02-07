/*----------------------------------------------------------------------*/
/*!
\file maxwell_0d_acinus.cpp

\maintainer Lena Yoshihara

\level 3
*/
/*----------------------------------------------------------------------*/


#include <vector>
#include "maxwell_0d_acinus.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Maxwell_0d_acinus::Maxwell_0d_acinus(
    Teuchos::RCP<MAT::PAR::Material> matdata) :
    Parameter(matdata), stiffness1_(matdata->GetDouble("Stiffness1")), stiffness2_(
        matdata->GetDouble("Stiffness2")), viscosity1_(
        matdata->GetDouble("Viscosity1")), viscosity2_(
        matdata->GetDouble("Viscosity2"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Maxwell_0d_acinus::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Maxwell_0d_acinus(this));
}


MAT::Maxwell_0d_acinusType MAT::Maxwell_0d_acinusType::instance_;


DRT::ParObject* MAT::Maxwell_0d_acinusType::Create( const std::vector<char> & data )
{
  MAT::Maxwell_0d_acinus* mxwll_0d_acin = new MAT::Maxwell_0d_acinus();
  mxwll_0d_acin->Unpack(data);
  return mxwll_0d_acin;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell_0d_acinus::Maxwell_0d_acinus()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Maxwell_0d_acinus::Maxwell_0d_acinus(MAT::PAR::Maxwell_0d_acinus* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Maxwell_0d_acinus*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Maxwell_0d_acinus::GetParams(std::string parametername)
{
  dserror("GetParams not implemented yet for this material!");
  return 0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Maxwell_0d_acinus::SetParams(std::string parametername, double new_value)
{
  dserror("SetParams not implemented yet for this material!");
}

