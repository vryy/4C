/*----------------------------------------------------------------------*/
/*!
\file elast_couplogneohooke.cpp
\brief


the input line should read
  MAT 1 ELAST_CoupLogNeoHooke YOUNG 1.044E7 NUE 0.3

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_couplogneohooke.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELAST::PAR::CoupLogNeoHooke::CoupLogNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  nue_(matdata->GetDouble("NUE"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::CoupLogNeoHooke::CoupLogNeoHooke()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELAST::CoupLogNeoHooke::CoupLogNeoHooke(MAT::ELAST::PAR::CoupLogNeoHooke* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELAST::CoupLogNeoHooke::AddCoefficientsPrincipal(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{
  havecoefficients = havecoefficients or true;

  // material parameters for isochoric part
  const double youngs = params_->youngs_;  // Young's modulus
  const double nue = params_->nue_;  // Poisson's ratio
  const double lambda = (nue==0.5) ? 0.0 : youngs*nue/((1.0+nue)*(1.0-2.0*nue));  // Lame coeff.
  const double mue = youngs/(2.0*(1.0+nue));  // shear modulus

  // determinant of deformation gradient
  const double detf = std::sqrt(prinv(2));

  // gammas
  gamma(0) += mue;
  gamma(1) += 0.0;
  gamma(2) += -mue+lambda*std::log(detf);

  // deltas
  delta(5) += lambda;
  delta(6) += 2.0*(mue - lambda*std::log(detf));

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
