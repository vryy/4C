/*-----------------------------------------------------------*/
/*!
\file crosslinkermat.cpp

\brief A class for a crosslinker material

\maintainer Jonas Eichinger

\date Oct, 2016

\level 3

*/
/*-----------------------------------------------------------*/

#include "crosslinkermat.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::CrosslinkerMat::CrosslinkerMat(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  poissonratio_(matdata->GetDouble("NUE")),
  density_(matdata->GetDouble("DENS"))
{
  if (youngs_<=0.)
    dserror("Young's modulus must be greater zero");
  if (poissonratio_>0.5 || poissonratio_<-1.)
    dserror("Poisson's ratio must be in [-1;0.5]");
}

Teuchos::RCP<MAT::Material> MAT::PAR::CrosslinkerMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::CrosslinkerMat(this));
}

MAT::CrosslinkerMatType MAT::CrosslinkerMatType::instance_;


DRT::ParObject* MAT::CrosslinkerMatType::Create( const std::vector<char> & data )
{
  MAT::CrosslinkerMat* linkermat = new MAT::CrosslinkerMat();
  linkermat->Unpack(data);
  return linkermat;
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::CrosslinkerMat::CrosslinkerMat()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::CrosslinkerMat::CrosslinkerMat(MAT::PAR::CrosslinkerMat* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::CrosslinkerMat::Pack(DRT::PackBuffer& data) const
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


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::CrosslinkerMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

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
        params_ = static_cast<MAT::PAR::CrosslinkerMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

