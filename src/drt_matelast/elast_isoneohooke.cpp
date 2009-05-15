/*----------------------------------------------------------------------*/
/*!
\file elast_isoneohooke.cpp
\brief


the input line should read
  MAT 1 ELAST_IsoNeoHooke MUE 100 

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
#include "elast_isoneohooke.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELAST::PAR::IsoNeoHooke::IsoNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  mue_(matdata->GetDouble("MUE"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::IsoNeoHooke::IsoNeoHooke()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::IsoNeoHooke::IsoNeoHooke(MAT::ELAST::PAR::IsoNeoHooke* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELAST::IsoNeoHooke::AddCoefficientsModified(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{
  havecoefficients = havecoefficients or true;

  const double mue = params_ -> mue_;
  
  gamma(0) += mue;
  
  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
