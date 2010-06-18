/*----------------------------------------------------------------------*/
/*!
\file elast_varisocub.cpp
\brief


the input line should read
  MAT 1 ELAST_VarIsoCub FRAC 0.5 C 100

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
#include "elast_varisocub.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::VarIsoCub::VarIsoCub(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  frac_(matdata->GetDouble("FRAC")),
  c_(matdata->GetDouble("C"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VarIsoCub::VarIsoCub()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VarIsoCub::VarIsoCub(MAT::ELASTIC::PAR::VarIsoCub* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VarIsoCub::AddCoefficientsModified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  havecoefficients = havecoefficients or true;

  const double frac = params_ -> frac_;
  const double c = params_ -> c_;

  gamma(0) += 6*c*frac*(modinv(0)-3)*(modinv(0)-3);

  delta(0) +=  24*c*frac*(modinv(0)-3);

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
