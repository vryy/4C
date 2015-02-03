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
#include "../drt_lib/drt_globalproblem.H"

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
  AddtoPack(data,alpha_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::UnpackSummand(
  const std::vector<char>& data,
  std::vector<char>::size_type& position
  )
{
  ExtractfromPack(position,data,alpha_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::Setup(DRT::INPUT::LineDefinition* linedef)
{
  double HU = 0.0;
  if (linedef->HaveNamed("HU"))
  {
    linedef->ExtractDouble("HU",HU);
  }
  else
  {
    HU = -999.0;
  }

  const Teuchos::ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  int HUlumen  = pslist.get<int>("MAXHULUMEN");

  //if HU is smaller than the threshold for calcification or smaller than
  // the Lumen HU (+10 security factor), there is no contribution at all
  if (HU <= params_->ctmin_ || HU < (HUlumen + 10))
    alpha_ = 0.;
  else if (HU > params_->ctmax_)
    alpha_ = params_->alphamax_;
  else
    alpha_ = 0.5 * params_->alphamax_ * (sin(PI * (HU - params_->ctmin_)/(params_->ctmax_ - params_->ctmin_) -PI/2) + 1.0);
}

/*----------------------------------------------------------------------
 *                                                      birzle 12/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::AddDerivativesModified(
    LINALG::Matrix<3,1>& dPmodI,
    LINALG::Matrix<6,1>& ddPmodII,
    const LINALG::Matrix<3,1>& modinv,
    const int eleGID
  )
{
  dPmodI(0) += alpha_;
  dPmodI(2) += (2.*alpha_*(1. - std::pow(modinv(2),-params_->beta_))) / ( (1.-2.*params_->nue_)*params_->beta_*modinv(2));


  ddPmodII(2) += (2.*alpha_*(-1. + std::pow(modinv(2),-params_->beta_)*(1.+params_->beta_))) / ( (1.-2.*params_->nue_)*params_->beta_*modinv(2)*modinv(2));

  return;
}

/*----------------------------------------------------------------------*/
