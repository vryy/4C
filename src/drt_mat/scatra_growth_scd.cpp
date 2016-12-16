/*----------------------------------------------------------------------*/
/*!
 \file scatra_growth_scd.cpp

 \brief  scatra material used for scalar dependent growth modeling

\maintainer Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-10364

\level 3
 *----------------------------------------------------------------------*/

#include "scatra_growth_scd.H"

#include "../drt_mat/matpar_bundle.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraGrowthScd::ScatraGrowthScd(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  diffusivity_(matdata->GetDouble("DIFFUSIVITY")),
  strdensity_(matdata->GetDouble("STRDENSITY")),
  rearate_(matdata->GetDouble("REARATE")),
  satcoeff_(matdata->GetDouble("SATCOEFF")),
  sourcemass_(matdata->Get<std::string>("SOURCEMASS"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ScatraGrowthScd::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraGrowthScd(this));
}


MAT::ScatraGrowthScdType MAT::ScatraGrowthScdType::instance_;


DRT::ParObject* MAT::ScatraGrowthScdType::Create( const std::vector<char> & data )
{
  MAT::ScatraGrowthScd* ScatraGrowthScd_mat = new MAT::ScatraGrowthScd();
  ScatraGrowthScd_mat->Unpack(data);
  return ScatraGrowthScd_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraGrowthScd::ScatraGrowthScd()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraGrowthScd::ScatraGrowthScd(MAT::PAR::ScatraGrowthScd* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraGrowthScd::Pack(DRT::PackBuffer& data) const
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
void MAT::ScatraGrowthScd::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ScatraGrowthScd*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

double MAT::ScatraGrowthScd::ComputeReactionCoeff(const double csnp, const double theta, const double dtheta, const double detFe) const
{
  double strdensity = params_->strdensity_;
  const double rearate = params_->rearate_;
  const double satcoeff = params_->satcoeff_;
  const std::string* sourcemass =params_->sourcemass_;
  double reacoeff = 0.;

  if (*sourcemass == "Standard")
    reacoeff = strdensity*3.0*dtheta/theta/csnp/detFe;
  else if (*sourcemass == "MonodBased")
    reacoeff = rearate/(satcoeff+csnp);

  else
    dserror ("The chosen source mass term is not implemented");

  return reacoeff;
}

double MAT::ScatraGrowthScd::ComputeReactionCoeffDeriv(const double csnp, const double theta, const double thetaold, const double dt) const //Ableitung des gesamten reaction-term nach der conc
{
//  const double reacoeffderiv = 0.0; // no linearisation implemented
  const double rearate = params_->rearate_;
  const double satcoeff = params_->satcoeff_;
  const double reacoeffderiv = -1.0*rearate/pow((satcoeff+csnp),2);

  return reacoeffderiv;
}


