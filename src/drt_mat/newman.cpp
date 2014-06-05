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
  diffcoefcurve_(matdata->GetInt("DIFFCOEF")),
  transnrcurve_(matdata->GetInt("TRANSNR")),
  thermfaccurve_(matdata->GetInt("THERMFAC")),
  condcurve_(matdata->GetInt("COND")),
  diffcoefparanum_(matdata->GetInt("DIFF_PARA_NUM")),
  diffcoefpara_(matdata->Get<std::vector<double> >("DIFF_PARA")),
  transnrparanum_(matdata->GetInt("TRANS_PARA_NUM")),
  transnrpara_(matdata->Get<std::vector<double> >("TRANS_PARA")),
  thermfacparanum_(matdata->GetInt("THERM_PARA_NUM")),
  thermfacpara_(matdata->Get<std::vector<double> >("THERM_PARA")),
  condparanum_(matdata->GetInt("COND_PARA_NUM")),
  condpara_(matdata->Get<std::vector<double> >("COND_PARA"))
{
  if (diffcoefparanum_ != (int)diffcoefpara_->size())
     dserror("number of materials %d does not fit to size of material vector %d", diffcoefparanum_, diffcoefpara_->size());
  if (transnrparanum_ != (int)transnrpara_->size())
     dserror("number of materials %d does not fit to size of material vector %d", transnrparanum_, transnrpara_->size());
  if (thermfacparanum_ != (int)thermfacpara_->size())
     dserror("number of materials %d does not fit to size of material vector %d", thermfacparanum_, thermfacpara_->size());
  if (condparanum_ != (int)condpara_->size())
     dserror("number of materials %d does not fit to size of material vector %d", condparanum_, condpara_->size());

  //check if number of provided parameter is valid for a the chosen predefined function
  CheckProvidedParams(diffcoefcurve_,diffcoefpara_->size());
  CheckProvidedParams(transnrcurve_,transnrpara_->size());
  CheckProvidedParams(thermfaccurve_,thermfacpara_->size());
  CheckProvidedParams(condcurve_,condpara_->size());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::PAR::Newman::CheckProvidedParams(
  const int           functnr,
  const unsigned int  numfunctparams
)
{
  bool error = false;
  std::string functionname;
  // needed parameter or a the predefined function
  unsigned int paraforfunction = 0;
  if(functnr<0)
  {
    switch (functnr)
    {
      case -1:
      {
        //constant value: functval=(*functparams)[0];
        functionname = "'constant value'";
        paraforfunction=1;
        if(numfunctparams != paraforfunction)
          error=true;
        break;
      }
      case -2:
      {
        // linear function: functval=(*functparams)[0]+(*functparams)[1]*cint;
        functionname = "'linear function'";
        paraforfunction=2;
        if(numfunctparams != paraforfunction)
          error=true;
        break;
      }
      case -3:
      {
        // quadratic function: functval=(*functparams)[0]+(*functparams)[1]*cint+(*functparams)[2]*cint*cint;
        functionname = "'quadratic function'";
        paraforfunction=3;
        if(numfunctparams != paraforfunction)
          error=true;
        break;
      }
      case -4:
      {
        // power function: functval=(*functparams)[0]*pow(cint,(*functparams)[1]);
        functionname = "'power function'";
        paraforfunction=2;
        if(numfunctparams != paraforfunction)
          error=true;
        break;
      }
      case -5:
      {
        // function 1 for conductivity;
        // const double nenner=(1.0+(*functparams)[2]*cint*cint-(*functparams)[3]*cint*cint*cint*cint);
        // functval=(*functparams)[0]*((*functparams)[1]*cint/nenner)+0.01;
        functionname = "'function 1 for conductivity'";
        paraforfunction=4;
        if(numfunctparams != paraforfunction)
          error=true;
        break;
      }
      case -6:
      {
        // a0*c + a1*c^1.5 + a2*c^3
        functionname = "'a0*c + a1*c^1.5 + a2*c^3'";
        paraforfunction=3;
        if(numfunctparams != paraforfunction)
          error=true;
        break;
      }
      case -7:
      {
        // a0 + a1*c + a2*c^2 + a3*c^3
        functionname = "'a0 + a1*c + a2*c^2 + a3*c^3'";
        paraforfunction=4;
        if(numfunctparams != paraforfunction)
          error=true;
        break;
      }
      case -8:
      {
        // thermodynamic factor Nyman 2008
        functionname = "'function thermodynamic factor (Nyman 2008)'";
        paraforfunction=7;
        if(numfunctparams != paraforfunction)
          error=true;
        break;
      }
      case -9:
      {
        // linear thermodynamic factor including Debye-Hückel theory
        functionname = "'function  linear thermodynamic factor (including Debye Hueckel theory)'";
        paraforfunction=2;
        if(numfunctparams != paraforfunction)
          error=true;
        break;
      }
      case -10:
      {
        // function 1 for conductivity
        functionname = "'function 1 for conductivity: own definition'";
        paraforfunction=6;
        if(numfunctparams != paraforfunction)
          error=true;
        break;
      }
      default: dserror("Curve number %i is not implemented",functnr); break;
    }

    if(error==true)
      dserror("number of %i provided parameter does not match the number of %i parameter "
              "which are needed for the predefined function with the number %i (%s)!!",numfunctparams,paraforfunction,functnr,functionname.c_str());
  }

  return;
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
  double diff=0.0;

  if(DiffCoefCurve()<0)
    diff = EvalFunctValue(DiffCoefCurve(),cint,DiffCoefParams());
  else if(DiffCoefCurve()==0)
    diff = EvalFunctValue(-1,cint,DiffCoefParams());
  else
    diff = DRT::Problem::Instance()->Curve(DiffCoefCurve()-1).f(cint);

  return diff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeFirstDerivDiffCoeff(const double cint) const
{
  double firstderiv=0.0;

  if(DiffCoefCurve()<0)
    firstderiv = EvalFirstDerivFunctValue(DiffCoefCurve(),cint,DiffCoefParams());
  else if(DiffCoefCurve()==0)
    firstderiv = EvalFirstDerivFunctValue(-1,cint,DiffCoefParams());
  else
    firstderiv = (DRT::Problem::Instance()->Curve(DiffCoefCurve()-1).FctDer(cint,1))[1];

  return firstderiv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeTransferenceNumber(const double cint) const
{
  double trans=0.0;

  if(TransNrCurve()<0)
    trans = EvalFunctValue(TransNrCurve(),cint,TransNrParams());
  else if(TransNrCurve()==0)
    trans = EvalFunctValue(-1,cint,TransNrParams());
  else trans = DRT::Problem::Instance()->Curve(TransNrCurve()-1).f(cint);

  return trans;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeFirstDerivTrans(const double cint) const
{
  double firstderiv=0.0;

  if(TransNrCurve()<0)
    firstderiv = EvalFirstDerivFunctValue(TransNrCurve(),cint,TransNrParams());
  else if(TransNrCurve()==0)
    firstderiv = EvalFirstDerivFunctValue(-1,cint,TransNrParams());
  else firstderiv=(DRT::Problem::Instance()->Curve(TransNrCurve()-1).FctDer(cint,1))[1];

  return firstderiv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeThermFac(const double cint) const
{
  double therm=0.0;

  if(ThermFacCurve()<0)
    therm = EvalFunctValue(ThermFacCurve(),cint,ThermFacParams());
  else if(ThermFacCurve()==0)
    // thermodynamic factor has to be one if not defined
    therm = 1.0;
  else therm = DRT::Problem::Instance()->Curve(ThermFacCurve()-1).f(cint);

  return therm;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeFirstDerivThermFac(const double cint) const
{
  double firstderiv=0.0;

  if(ThermFacCurve()<0)
    firstderiv = EvalFirstDerivFunctValue(ThermFacCurve(),cint,ThermFacParams());
  else if(ThermFacCurve()==0)
    // thermodynamic factor has to be one if not defined
    // -> first derivative = 0.0
    firstderiv = 0.0;
  else  firstderiv = (DRT::Problem::Instance()->Curve(ThermFacCurve()-1).FctDer(cint,1))[1];

  return firstderiv;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeConductivity(const double cint) const
{
  double cond = 0.0;

  if(CondCurve()<0)
    cond = EvalFunctValue(CondCurve(),cint,CondParams());
  else if(CondCurve()==0)
    cond = EvalFunctValue(-1,cint,CondParams());
  else cond = DRT::Problem::Instance()->Curve(CondCurve()-1).f(cint);


  return cond;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::ComputeFirstDerivCond(const double cint) const
{
  double firstderiv=0.0;

  if(CondCurve()<0)
    firstderiv = EvalFirstDerivFunctValue(CondCurve(),cint,CondParams());
  else if(CondCurve()==0)
    firstderiv = EvalFirstDerivFunctValue(-1,cint,CondParams());
  else firstderiv = (DRT::Problem::Instance()->Curve(CondCurve()-1).FctDer(cint,1))[1];

  return firstderiv;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::EvalFunctValue(
  const int             functnr,
  const double          cint,
  const std::vector<double>*   functparams) const
{
  double functval=0.0;

  switch (functnr)
  {
    // a0
    case -1: functval=(*functparams)[0]; break;
    // a0 + a1*c
    case -2: functval=(*functparams)[0]+(*functparams)[1]*cint; break;
    // a0 + a1*c + a2*c^2
    case -3: functval=(*functparams)[0]+(*functparams)[1]*cint+(*functparams)[2]*cint*cint; break;
    // a0 + c^a1
    case -4: functval=(*functparams)[0]*pow(cint,(*functparams)[1]); break;
    // conductivity
    case -5:
    {
      const double nenner=(1.0+(*functparams)[2]*cint*cint-(*functparams)[3]*cint*cint*cint*cint);
      // (*functparams)[0]*((*functparams)[1]*cint/nenner) + 0.01 -> constant level 0.01 deleted since it does not have a physical meaning (28.04.2014)
      functval=(*functparams)[0]*((*functparams)[1]*cint/nenner);
      break;
    }
    // a0*c + a1*c^1.5 + a2*c^3
    case -6: functval=(*functparams)[0]*cint+(*functparams)[1]*pow(cint,1.5)+(*functparams)[2]*cint*cint*cint; break;
    // a0 + a1*c + a2*c^2 + a3*c^3
    case -7: functval=(*functparams)[0]+(*functparams)[1]*cint+(*functparams)[2]*cint*cint+(*functparams)[3]*cint*cint*cint; break;
    // thermodynamic factor Nyman 2008
    case -8:
    {
      const double num=(*functparams)[0]+(*functparams)[1]*cint+(*functparams)[2]*cint*cint;
      const double denom=(*functparams)[3]+(*functparams)[4]*cint+(*functparams)[5]*cint*cint+(*functparams)[6]*cint*cint*cint;
      functval=num/denom;
      break;
    }
    // linear thermodynamic factor including Debye-Hückel theory
    // 1 + a1*0.5*c^0.5 + a2*c
    case -9: functval= 1.0 + (*functparams)[0]*0.5*pow(cint,0.5)+(*functparams)[1]*cint; break;
    // conductivity: own definition which also fulfills the Kohlrausches Square root law
    case -10:
    {
      const double num = (*functparams)[0]*cint+(*functparams)[1]*pow(cint,1.5)+(*functparams)[2]*cint*cint+(*functparams)[3]*cint*cint*cint;
      const double denom=(1.0+(*functparams)[4]*cint*cint+(*functparams)[5]*cint*cint*cint*cint);
      // (*functparams)[0]*((*functparams)[1]*cint/nenner) + 0.01 -> constant level 0.01 deleted since it does not have a physical meaning (28.04.2014)
      functval=num/denom;
      break;
    }
    default: dserror("Curve number %i is not implemented",functnr); break;
  }
  return functval;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Newman::EvalFirstDerivFunctValue(
    const int                   functnr,
    const double                cint,
    const std::vector<double>*  functparams) const
{
  double firstderivfunctval=0.0;

  switch (functnr)
  {
    // d/dc: a0
    case -1: firstderivfunctval=0.0; break;
    // d/dc: a0 + a1*c
    case -2: firstderivfunctval=(*functparams)[1]; break;
    // d/dc: a0 + a1*c + a2*c^2
    case -3: firstderivfunctval=(*functparams)[1]+2*(*functparams)[2]*cint; break;
    // d/dc: a0 + c^a1
    case -4: firstderivfunctval=(*functparams)[0]*(*functparams)[1]*pow(cint,(*functparams)[1]-1.0); break;
    // d/dc: conductivity
    case -5:
    {
      const double nenner=(1.0+(*functparams)[2]*cint*cint-(*functparams)[3]*cint*cint*cint*cint);
      const double nennernenner = nenner*nenner;
      firstderivfunctval=(*functparams)[0]*(((*functparams)[1]*nenner-(*functparams)[1]*cint*(2*(*functparams)[2]*cint-4*(*functparams)[3]*cint*cint*cint))/nennernenner);
      break;
    }
    // d/dc: a0*c + a1*c^1.5 + a2*c^3
    case -6: firstderivfunctval=(*functparams)[0]+1.5*(*functparams)[1]*pow(cint,0.5)+3*(*functparams)[2]*cint*cint; break;
    // d/dc: a0 + a1*c + a2*c^2 + a3*c^3
    case -7: firstderivfunctval=(*functparams)[1]+2*(*functparams)[2]*cint+3*(*functparams)[3]*cint*cint; break;
    // d/dc: thermodynamic factor Nyman 2008
    case -8:
    {
      const double num=(*functparams)[0]+(*functparams)[1]*cint+(*functparams)[2]*cint*cint;
      const double denom=(*functparams)[3]+(*functparams)[4]*cint+(*functparams)[5]*cint*cint+(*functparams)[6]*cint*cint*cint;
      const double denomdenom = denom*denom;
      const double derivnum=(*functparams)[1]+2*(*functparams)[2]*cint;
      const double derivdenom=(*functparams)[4]+2*(*functparams)[5]*cint+3*(*functparams)[6]*cint*cint;
      firstderivfunctval=(derivnum*denom-num*derivdenom)/denomdenom;
      break;
    }
    // linear thermodynamic factor including Debye-Hückel theory
    // d/dc: 1 + a1*0.5*c^0.5 + a2*c
    case -9: firstderivfunctval= (*functparams)[0]*0.5*0.5*pow(cint,-0.5)+(*functparams)[1]; break;
    // d/dc: conductivity: own definition which also fulfills the Kohlrausches Square root law
    case -10:
    {
      const double num = (*functparams)[0]*cint+(*functparams)[1]*pow(cint,1.5)+(*functparams)[2]*cint*cint+(*functparams)[3]*cint*cint*cint;
      const double denom=(1.0+(*functparams)[4]*cint*cint+(*functparams)[5]*cint*cint*cint*cint);
      const double denomdenom = denom*denom;
      const double derivnum=(*functparams)[0]+1.5*(*functparams)[1]*pow(cint,0.5)+2.0*(*functparams)[2]*cint+3.0*(*functparams)[3]*cint*cint;
      const double derivdenom=2.0*(*functparams)[4]*cint+4.0*(*functparams)[5]*cint*cint*cint;
      firstderivfunctval= ((derivnum*denom-num*derivdenom)/denomdenom);
      break;
    }
    default: dserror("Curve number %i is not implemented",functnr); break;
  }

  return firstderivfunctval;
}



