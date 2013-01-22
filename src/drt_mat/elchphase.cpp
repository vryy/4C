/*!----------------------------------------------------------------------*/
/*!
\file elchphase.cpp

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/


#include <vector>
#include "elchphase.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElchPhase::ElchPhase(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  epsilon_(matdata->GetDouble("EPSILON")),
  conductivity_(matdata->GetDouble("CONDUCTIVITY")),
  condcurvenr_(matdata->GetInt("NR"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::ElchPhase::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ElchPhase(this));
}

MAT::ElchPhaseType MAT::ElchPhaseType::instance_;


DRT::ParObject* MAT::ElchPhaseType::Create( const std::vector<char> & data )
{
  MAT::ElchPhase* elchphase = new MAT::ElchPhase();
  elchphase->Unpack(data);
  return elchphase;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElchPhase::ElchPhase()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElchPhase::ElchPhase(MAT::PAR::ElchPhase* params)
  : params_(params)
{
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchPhase::Pack(DRT::PackBuffer& data) const
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
void MAT::ElchPhase::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ElchPhase*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ElchPhase::ComputeConductivity(const double cint) const
{
  double cond = 0.0;

  if(CondCurveNr()==0)
  {
    cond = Conductivity();
  }
  else
  {
    cond = DRT::Problem::Instance()->Curve(CondCurveNr()-1).f(cint);
  }

  return cond;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::ElchPhase::ComputeFirstDerivCond(const double cint) const
{
  double firstderiv=0.0;

  if(CondCurveNr()==0)
  {
    firstderiv = 0.0;
  }
  else
  {
    firstderiv = (DRT::Problem::Instance()->Curve(CondCurveNr()-1).FctDer(cint,1))[1];
  }

  return firstderiv;
}

