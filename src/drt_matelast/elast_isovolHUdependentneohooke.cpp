/*----------------------------------------------------------------------*/
/*!
\file elast_isovolHUdependentneohooke.cpp
\brief

This file contains the routines required to calculate the isochoric and volumetric contribution
of the cc neo hooke material

MAT 20 ELAST_IsoVolHUDependentNeoHooke ALPHA_MAX 8.929E6 CT_MIN 30.0 CT_MAX 600.0 NUE 0.49 BETA -2.0

<pre>
maintainer: Andreas Maier
</pre>

*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
#include "elast_isovolHUdependentneohooke.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke::IsoVolHUDependentNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  alphamax_(matdata->GetDouble("ALPHA_MAX")),
  ctmin_(matdata->GetDouble("CT_MIN")),
  ctmax_(matdata->GetDouble("CT_MAX")),
  nue_(matdata->GetDouble("NUE")),
  beta_(matdata->GetDouble("BETA"))
{
}

Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke::CreateMaterial()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public) AMaier 06/11 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVolHUDependentNeoHooke::IsoVolHUDependentNeoHooke()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                 (public)   AMaier 06/11 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVolHUDependentNeoHooke::IsoVolHUDependentNeoHooke(MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke* params)
  : params_(params)
{

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data,HU_);
  AddtoPack(data,HUlumen_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::UnpackSummand(const std::vector<char>& data,
																														std::vector<char>::size_type& position)
{
  ExtractfromPack(position,data,HU_);
  ExtractfromPack(position,data,HUlumen_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::Setup(DRT::INPUT::LineDefinition* linedef)
{
  if (linedef->HaveNamed("HU"))
  {
    linedef->ExtractDouble("HU",HU_);
  }
  else
  {
    HU_ = -999.0;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::SetupAAA(Teuchos::ParameterList& params)
{
  HUlumen_ = params.get("max hu lumen",-999);

  //if HU is smaller than the threshold for calcification or smaller than
  // the Lumen HU (+10 security factor), there is no contribution at all
  if (HU_ <= params_->ctmin_ || HU_ < (HUlumen_ + 10))
    alpha_ = 0.;
  else if (HU_ > params_->ctmax_)
    alpha_ = params_->alphamax_;
  else
    alpha_ = 0.5 * params_->alphamax_ * (sin(PI * (HU_ - params_->ctmin_)/(params_->ctmax_ - params_->ctmin_) -PI/2) + 1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::AddCoefficientsPrincipal(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{
  // principal coefficients
  gamma(0) += alpha_ * (2.*1.*pow(prinv(2),-1./3.));
  gamma(1) += 0.;
  gamma(2) += alpha_ * (-(2./3.)*1.*prinv(0)*pow(prinv(2),-1./3.) + (2./(1.-2.*params_->nue_))*1.*(1.-pow(prinv(2),-params_->beta_/2.))/params_->beta_);

  delta(0) += 0.;
  delta(1) += 0.;
  delta(2) += alpha_ * (-(4./3.)*1.*pow(prinv(2),-1./3.));
  delta(3) += 0.;
  delta(4) += 0.;
  delta(5) += alpha_ * ((4./9.)*1.*prinv(0)*pow(prinv(2),-1./3.) + (2./(1.-2.*params_->nue_))*1.*pow(prinv(2),-params_->beta_/2.));
  delta(6) += alpha_ * ((4./3.)*1.*prinv(0)*pow(prinv(2),-1./3.) + (2./(1.-2.*params_->nue_))*1.*2.*(pow(prinv(2),-params_->beta_/2.)-1.)/params_->beta_);
  delta(7) += 0.;

  return;
}

/*----------------------------------------------------------------------*/
