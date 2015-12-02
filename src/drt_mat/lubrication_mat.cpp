/*--------------------------------------------------------------------------*/
/*!
\file lubrication_mat.cpp

\brief Material model for the lubrication film

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/


#include <vector>
#include "lubrication_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::LubricationMat::LubricationMat(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata)
{
  Epetra_Map dummy_map(1,1,0,*(DRT::Problem::Instance()->GetNPGroup()->LocalComm()));
  for(int i=first ; i<=last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map,true)));
  }
  matparams_.at(viscosity)->PutScalar(matdata->GetDouble("VISCOSITY"));
  matparams_.at(density)->PutScalar(matdata->GetDouble("DENSITY"));

  return;
}


Teuchos::RCP<MAT::Material> MAT::PAR::LubricationMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::LubricationMat(this));
}


MAT::LubricationMatType MAT::LubricationMatType::instance_;


void MAT::PAR::LubricationMat::OptParams(std::map<std::string,int>* pnames)
{
  pnames->insert(std::pair<std::string,int>("VISC", viscosity));
  pnames->insert(std::pair<std::string,int>("DENS", density));
}


DRT::ParObject* MAT::LubricationMatType::Create( const std::vector<char> & data )
{
  MAT::LubricationMat* lubrication_mat = new MAT::LubricationMat();
  lubrication_mat->Unpack(data);
  return lubrication_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::LubricationMat::LubricationMat()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::LubricationMat::LubricationMat(MAT::PAR::LubricationMat* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::LubricationMat::Pack(DRT::PackBuffer& data) const
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::LubricationMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data. type = %d, UniqueParObjectId()=%d",type,UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::LubricationMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}
