/*----------------------------------------------------------------------*/
/*!
\file plasticneohooke.H
\brief Containing all plastic material parameters needed for
       the semi-smooth Newton method for plasticity.
       This is NOT a stand-alone material with proper evaluate
       routines. A hyperelastic material is needed to desribe
       the elastic part.


<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "plasticsemismooth.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PAR::PlasticSemiSmooth::PlasticSemiSmooth(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  yield_(matdata->GetDouble("YIELD")),
  isohard_(matdata->GetDouble("ISOHARD")),
  expisohard_(matdata->GetDouble("EXPISOHARD")),
  infyield_(matdata->GetDouble("INFYIELD")),
  kinhard_(matdata->GetDouble("KINHARD"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::PlasticSemiSmooth::CreateMaterial()
{
  return Teuchos::rcp(new MAT::PlasticSemiSmooth(this));
}

MAT::PlasticSemiSmoothType MAT::PlasticSemiSmoothType::instance_;


DRT::ParObject* MAT::PlasticSemiSmoothType::Create( const std::vector<char> & data )
{
  MAT::PlasticSemiSmooth* plastic = new MAT::PlasticSemiSmooth();
  plastic->Unpack(data);
  return plastic;
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PlasticSemiSmooth::PlasticSemiSmooth()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::PlasticSemiSmooth::PlasticSemiSmooth(MAT::PAR::PlasticSemiSmooth* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::PlasticSemiSmooth::Pack(DRT::PackBuffer& data) const
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
void MAT::PlasticSemiSmooth::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::PlasticSemiSmooth*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}
