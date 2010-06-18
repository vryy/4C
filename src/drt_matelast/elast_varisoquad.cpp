/*----------------------------------------------------------------------*/
/*!
\file elast_varisoquad.cpp
\brief


the input line should read
  MAT 1 ELAST_VarIsoQuad FRAC 0.5 C 100

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
#include "elast_varisoquad.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::VarIsoQuad::VarIsoQuad(
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
MAT::ELASTIC::VarIsoQuad::VarIsoQuad()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VarIsoQuad::VarIsoQuad(MAT::ELASTIC::PAR::VarIsoQuad* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VarIsoQuad::AddCoefficientsModified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{
  havecoefficients = havecoefficients or true;

  const double frac = params_ -> frac_;
  const double c = params_ -> c_;

  gamma(0) += 4*frac*c*(modinv(0)-3);

  delta(0) += 8*frac*c;

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
