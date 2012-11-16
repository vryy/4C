/*!----------------------------------------------------------------------
\file biofilm.cpp

<pre>
Maintainer: Mirella Coroneo
            coroneo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*----------------------------------------------------------------------*/


#include <vector>
#include "biofilm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Biofilm::Biofilm(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  diffusivity_(matdata->GetDouble("DIFFUSIVITY")),
  kinetics_(matdata->Get<string>("KINETICS")),
  rearate_(matdata->GetDouble("REARATE")),
  satcoeff_(matdata->GetDouble("SATCOEFF"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::Biofilm::CreateMaterial()
{
	return Teuchos::rcp(new MAT::Biofilm(this));
}


MAT::BiofilmType MAT::BiofilmType::instance_;


DRT::ParObject* MAT::BiofilmType::Create( const std::vector<char> & data )
{
  MAT::Biofilm* biofilm_mat = new MAT::Biofilm();
  biofilm_mat->Unpack(data);
  return biofilm_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Biofilm::Biofilm()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Biofilm::Biofilm(MAT::PAR::Biofilm* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Biofilm::Pack(DRT::PackBuffer& data) const
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
void MAT::Biofilm::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::Biofilm*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

double MAT::Biofilm::ComputeReactionCoeff(const double csnp) const
{
  const double reacoeff = ReaRate()/(SatCoeff()+csnp);

  return reacoeff;
}

double MAT::Biofilm::ComputeReactionCoeffDeriv(const double csnp) const
{
  const double reacoeffderiv = (ReaRate()*SatCoeff())/pow((SatCoeff()+csnp),2);

  return reacoeffderiv;
}

