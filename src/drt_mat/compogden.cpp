/*!----------------------------------------------------------------------
\file compogden.cpp

\level 3

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            089/289 15255
</pre>
\param  Epetra_SerialDenseVector* glstrain      (i) Green-Lagrange strains
\param  Epetra_SerialDenseVector* stress        (o) ele stress vector
\param  Epetra_SerialDenseMatrix* cmat          (o) constitutive matrix
*----------------------------------------------------------------------*/


#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "compogden.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::CompOgden::CompOgden(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  init_(-1),
  nue_(matdata->GetDouble("NUE")),
  beta_(matdata->GetDouble("BETA")),
  alfap_(),
  mup_(),
  density_(matdata->GetDouble("DENS")),
  lambda_(matdata->GetDouble("LAMBDA")),
  kappa_(matdata->GetDouble("KAPPA")),
  l_()
{
  alfap_[0] = matdata->GetDouble("ALFA1");
  alfap_[1] = matdata->GetDouble("ALFA2");
  alfap_[2] = matdata->GetDouble("ALFA3");

  mup_[0] = matdata->GetDouble("NU1");
  mup_[1] = matdata->GetDouble("NU2");
  mup_[2] = matdata->GetDouble("NU3");

  l_[0] = 0;
  l_[1] = 0;
  l_[2] = 0;
}

Teuchos::RCP<MAT::Material> MAT::PAR::CompOgden::CreateMaterial()
{
  return Teuchos::rcp(new MAT::CompOgden(this));
}


MAT::CompOgdenType MAT::CompOgdenType::instance_;


DRT::ParObject* MAT::CompOgdenType::Create( const std::vector<char> & data )
{
  return NULL;
}


/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
MAT::CompOgden::CompOgden()
  : params_(NULL)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::CompOgden::CompOgden(MAT::PAR::CompOgden* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
void MAT::CompOgden::Pack(DRT::PackBuffer& data) const
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
/*---------------------------------------------------------------------*/
void MAT::CompOgden::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::CompOgden*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}



