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
MAT::ELASTIC::PAR::CoupLogNeoHooke::CoupLogNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  mue_(matdata->GetDouble("MUE")),
  lambda_(matdata->GetDouble("LAMBDA")),
  parmode_(matdata->GetInt("MODE")),
  youngs_(matdata->GetDouble("YOUNG")),
  nue_(matdata->GetDouble("NUE"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupLogNeoHooke::CoupLogNeoHooke()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupLogNeoHooke::CoupLogNeoHooke(MAT::ELASTIC::PAR::CoupLogNeoHooke* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupLogNeoHooke::AddShearMod(
  bool& haveshearmod,  ///< non-zero shear modulus was added
  double& shearmod  ///< variable to add upon
  ) const
{
  haveshearmod = haveshearmod or true;

  // material parameters for isochoric part
  if (params_->parmode_ == 0) {
    shearmod += params_->mue_;
  }
  else if (params_->parmode_ == 1) {
    const double youngs = params_->youngs_;  // Young's modulus
    const double nue = params_->nue_;  // Poisson's ratio
    shearmod += youngs/(2.0*(1.0+nue));  // shear modulus
  }
  else {
    dserror("Cannot handle mode=%d", params_->parmode_);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupLogNeoHooke::AddCoefficientsPrincipal(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{
  havecoefficients = havecoefficients or true;

  // material parameters for isochoric part
  double lambda = 0.0;
  double mue = 0.0;
  if (params_->parmode_ == 0) {
    lambda = params_->lambda_;
    mue = params_->mue_;
  }
  else if (params_->parmode_ == 1) {
    const double youngs = params_->youngs_;  // Young's modulus
    const double nue = params_->nue_;  // Poisson's ratio
    lambda = (nue==0.5) ? 0.0 : youngs*nue/((1.0+nue)*(1.0-2.0*nue));  // Lame coeff.
    mue = youngs/(2.0*(1.0+nue));  // shear modulus
  }
  else {
    dserror("Cannot handle mode=%d", params_->parmode_);
  }

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
