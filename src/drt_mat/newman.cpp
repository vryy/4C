/*!----------------------------------------------------------------------*/
/*!
\file newman.cpp

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/


#include <vector>
#include "newman.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

// TODO: math.H was included automatically


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Newman::Newman(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  valence_(matdata->GetDouble("VALENCE")),
  curvediff_(matdata->GetInt("CURVE_DIFF")),
  curvetrans_(matdata->GetInt("CURVE_TRANS")),
  cursolvar_((matdata->GetInt("CURSOLVAR")))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::Newman::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Newman(this));
}

MAT::NewmanType MAT::NewmanType::instance_;


DRT::ParObject* MAT::NewmanType::Create( const std::vector<char> & data )
{
  MAT::Newman* newman = new MAT::Newman();
  newman->Unpack(data);
  return newman;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Newman::Newman()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Newman::Newman(MAT::PAR::Newman* params)
  : params_(params)
{
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Newman::Pack(DRT::PackBuffer& data) const
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
void MAT::Newman::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::Newman*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeDiffusionCoefficient(const double cint) const
{
  const double diff = DRT::Problem::Instance()->Curve(CurveDiff()-1).f(cint);

  return diff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeFirstDerivDiffCoeff(const double cint) const
{
  double firstderiv = (DRT::Problem::Instance()->Curve(CurveDiff()-1).FctDer(cint,1))[1];

  return firstderiv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeTransferenceNumber(const double cint) const
{
  const double trans = DRT::Problem::Instance()->Curve(CurveTrans()-1).f(cint);

  return trans;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeFirstDerivTrans(const double cint) const
{
  double firstderiv = (DRT::Problem::Instance()->Curve(CurveTrans()-1).FctDer(cint,1))[1];

  return firstderiv;
}

