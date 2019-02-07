/*!----------------------------------------------------------------------*/
/*!
\file newman.cpp

\level 2

<pre>
\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
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
MAT::PAR::Newman::Newman(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ElchSingleMat(matdata),
      valence_(matdata->GetDouble("VALENCE")),
      transnrcurve_(matdata->GetInt("TRANSNR")),
      thermfaccurve_(matdata->GetInt("THERMFAC")),
      transnrparanum_(matdata->GetInt("TRANS_PARA_NUM")),
      transnrpara_(*matdata->Get<std::vector<double>>("TRANS_PARA")),
      thermfacparanum_(matdata->GetInt("THERM_PARA_NUM")),
      thermfacpara_(*matdata->Get<std::vector<double>>("THERM_PARA"))
{
  if (transnrparanum_ != (int)transnrpara_.size())
    dserror("number of materials %d does not fit to size of material vector %d", transnrparanum_,
        transnrpara_.size());
  if (thermfacparanum_ != (int)thermfacpara_.size())
    dserror("number of materials %d does not fit to size of material vector %d", thermfacparanum_,
        thermfacpara_.size());

  // check if number of provided parameter is valid for a the chosen predefined function
  CheckProvidedParams(transnrcurve_, transnrpara_);
  CheckProvidedParams(thermfaccurve_, thermfacpara_);
}


Teuchos::RCP<MAT::Material> MAT::PAR::Newman::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Newman(this));
}

MAT::NewmanType MAT::NewmanType::instance_;


DRT::ParObject* MAT::NewmanType::Create(const std::vector<char>& data)
{
  MAT::Newman* newman = new MAT::Newman();
  newman->Unpack(data);
  return newman;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Newman::Newman() : params_(NULL) { return; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Newman::Newman(MAT::PAR::Newman* params) : params_(params) { return; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Newman::Pack(DRT::PackBuffer& data) const
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

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Newman::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

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
        params_ = static_cast<MAT::PAR::Newman*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeTransferenceNumber(const double cint) const
{
  double trans = 0.0;

  if (TransNrCurve() < 0)
    trans = EvalFunctValue(TransNrCurve(), cint, TransNrParams());
  else if (TransNrCurve() == 0)
    trans = EvalFunctValue(-1, cint, TransNrParams());
  else
    trans = DRT::Problem::Instance()->Funct(TransNrCurve() - 1).EvaluateTime(cint);

  return trans;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeFirstDerivTrans(const double cint) const
{
  double firstderiv = 0.0;

  if (TransNrCurve() < 0)
    firstderiv = EvalFirstDerivFunctValue(TransNrCurve(), cint, TransNrParams());
  else if (TransNrCurve() == 0)
    firstderiv = EvalFirstDerivFunctValue(-1, cint, TransNrParams());
  else
    firstderiv =
        (DRT::Problem::Instance()->Funct(TransNrCurve() - 1).EvaluateTimeDerivative(cint, 1))[1];

  return firstderiv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeThermFac(const double cint) const
{
  double therm = 0.0;

  if (ThermFacCurve() < 0)
    therm = EvalFunctValue(ThermFacCurve(), cint, ThermFacParams());
  else if (ThermFacCurve() == 0)
    // thermodynamic factor has to be one if not defined
    therm = 1.0;
  else
    therm = DRT::Problem::Instance()->Funct(ThermFacCurve() - 1).EvaluateTime(cint);

  return therm;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeFirstDerivThermFac(const double cint) const
{
  double firstderiv = 0.0;

  if (ThermFacCurve() < 0)
    firstderiv = EvalFirstDerivFunctValue(ThermFacCurve(), cint, ThermFacParams());
  else if (ThermFacCurve() == 0)
    // thermodynamic factor has to be one if not defined
    // -> first derivative = 0.0
    firstderiv = 0.0;
  else
    firstderiv =
        (DRT::Problem::Instance()->Funct(ThermFacCurve() - 1).EvaluateTimeDerivative(cint, 1))[1];

  return firstderiv;
}
