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
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
#include "elast_isovolHUdependentneohooke.H"

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
 |  Constructor                                   (public) AMaier 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVolHUDependentNeoHooke::IsoVolHUDependentNeoHooke()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                 (public)   AMaier 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVolHUDependentNeoHooke::IsoVolHUDependentNeoHooke(MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke* params)
  : params_(params)
{

}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::AddCoefficientsPrincCalcified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  LINALG::Matrix<3,1>& matparam,
  const LINALG::Matrix<3,1>& prinv
  )
{
  havecoefficients = havecoefficients or true;
 
  matparam(0) = params_->alphamax_;
  matparam(1) = params_->ctmin_;
  matparam(2) = params_->ctmax_;


  // principal coefficients
  gamma(0) += 2.*1.*pow(prinv(2),-1./3.);
  gamma(1) += 0.;
  gamma(2) += -(2./3.)*1.*prinv(0)*pow(prinv(2),-1./3.) + (2./(1.-2.*params_->nue_))*1.*(1.-pow(prinv(2),-params_->beta_/2.))/params_->beta_;

  delta(0) += 0.;
  delta(1) += 0.;
  delta(2) += -(4./3.)*1.*pow(prinv(2),-1./3.);
  delta(3) += 0.;
  delta(4) += 0.;
  delta(5) += (4./9.)*1.*prinv(0)*pow(prinv(2),-1./3.) + (2./(1.-2.*params_->nue_))*1.*pow(prinv(2),-params_->beta_/2.);
  delta(6) += (4./3.)*1.*prinv(0)*pow(prinv(2),-1./3.) + (2./(1.-2.*params_->nue_))*1.*2.*(pow(prinv(2),-params_->beta_/2.)-1.)/params_->beta_;
  delta(7) += 0.;

  return;
}

/*----------------------------------------------------------------------*/
#endif // CCADISCRET
