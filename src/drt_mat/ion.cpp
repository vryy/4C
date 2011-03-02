/*!----------------------------------------------------------------------
\file ion.cpp

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include "ion.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Ion::Ion(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  valence_(matdata->GetDouble("VALENCE")),
  diffusivity_(matdata->GetDouble("DIFFUSIVITY")),
  densification_(matdata->GetDouble("DENSIFICATION"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::Ion::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Ion(this));
}

MAT::IonType MAT::IonType::instance_;


DRT::ParObject* MAT::IonType::Create( const std::vector<char> & data )
{
  MAT::Ion* ion = new MAT::Ion();
  ion->Unpack(data);
  return ion;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Ion::Ion()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Ion::Ion(MAT::PAR::Ion* params)
  : params_(params)
{
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Ion::Pack(DRT::PackBuffer& data) const
{
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
void MAT::Ion::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::Ion*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

#endif
