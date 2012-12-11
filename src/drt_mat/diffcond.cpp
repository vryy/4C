/*!----------------------------------------------------------------------
\file diffcond.cpp

<pre>
Maintainer: Andreas Bauer
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*----------------------------------------------------------------------*/


#include <vector>
#include "diffcond.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::DiffCond::DiffCond(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  valence_(matdata->GetDouble("VALENCE")),
  diffusivity_(matdata->GetDouble("DIFFUSIVITY")),
  transference_(matdata->GetDouble("TRANSFERENCE"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::DiffCond::CreateMaterial()
{
  return Teuchos::rcp(new MAT::DiffCond(this));
}

MAT::DiffCondType MAT::DiffCondType::instance_;


DRT::ParObject* MAT::DiffCondType::Create( const std::vector<char> & data )
{
  MAT::DiffCond* diffcond = new MAT::DiffCond();
  diffcond->Unpack(data);
  return diffcond;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::DiffCond::DiffCond()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::DiffCond::DiffCond(MAT::PAR::DiffCond* params)
  : params_(params)
{
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::DiffCond::Pack(DRT::PackBuffer& data) const
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
void MAT::DiffCond::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::DiffCond*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


