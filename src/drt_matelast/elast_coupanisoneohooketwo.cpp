/*----------------------------------------------------------------------*/
/*!
\file elast_coupanisoneohooketwo.cpp
\brief


the input line should read
  MAT 1 ELAST_CoupAnisoNeoHookeTwo C1 100 C2 100

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupanisoneohooketwo.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELAST::PAR::CoupAnisoNeoHookeTwo::CoupAnisoNeoHookeTwo(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c1_(matdata->GetDouble("C1")),
  c2_(matdata->GetDouble("C2"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::CoupAnisoNeoHookeTwo::CoupAnisoNeoHookeTwo()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::CoupAnisoNeoHookeTwo::CoupAnisoNeoHookeTwo(MAT::ELAST::PAR::CoupAnisoNeoHookeTwo* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELAST::CoupAnisoNeoHookeTwo::AddCoefficientsPrincipalAniso(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<15,1>& delta,
  const LINALG::Matrix<6,1>& prinv
  )
{
  havecoefficients = havecoefficients or true;

  double c1=params_->c1_;
  double c2=params_->c2_;

  gamma(0) += 2.*c1;
  gamma(1) += 2.*c2;

  delta(0) += 0.;
  delta(1) += 0.;

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
