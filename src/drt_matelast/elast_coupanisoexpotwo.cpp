/*----------------------------------------------------------------------*/
/*!
\file elast_isoneohooke.cpp
\brief


the input line should read
  MAT 1 ELAST_CoupAnisoExpoTwo MUE 100 

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            rausch,kloeppel@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoexpotwo.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELAST::PAR::CoupAnisoExpoTwo::CoupAnisoExpoTwo(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  k1_(matdata->GetDouble("K1")),
  k2_(matdata->GetDouble("K2")),
  k3_(matdata->GetDouble("K3")),
  k4_(matdata->GetDouble("K4"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::CoupAnisoExpoTwo::CoupAnisoExpoTwo()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::CoupAnisoExpoTwo::CoupAnisoExpoTwo(MAT::ELAST::PAR::CoupAnisoExpoTwo* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELAST::CoupAnisoExpoTwo::AddCoefficientsPrincipalAniso(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<15,1>& delta,
  const LINALG::Matrix<6,1>& prinv
  )
{
  havecoefficients = havecoefficients or true;
  
  double k1=params_->k1_;
  double k2=params_->k2_;
  double k3=params_->k3_;
  double k4=params_->k4_;
  
  gamma(0) +=  2.*(k1*(prinv(3)-1.)* exp(k2*(prinv(3)-1.)*(prinv(3)-1.)));
  gamma(1) +=  2.*(k3*(prinv(4)-1.)* exp(k4*(prinv(4)-1.)*(prinv(4)-1.)));
  
  delta(0) += 2.*(1. + 2.*k2*(prinv(3)-1.)*(prinv(3)-1.))*2.*k1* exp(k2*(prinv(3)-1.)*(prinv(3)-1.));
  delta(1) += 2.*(1. + 2.*k4*(prinv(4)-1.)*(prinv(4)-1.))*2.*k3* exp(k4*(prinv(4)-1.)*(prinv(4)-1.));

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
