/*----------------------------------------------------------------------*/
/*!
\file elast_couplogneohooke.cpp
\brief


the input line should read
  MAT 1 ELAST_CoupLogNeoHooke YOUNG 1.044E7 NUE 0.3 MODE YN
  MAT 1 ELAST_CoupLogNeoHooke MUE 1. LAMBDA 1. MODE Lame

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
  parmode_(*(matdata->Get<std::string>("MODE"))),
  youngs_(matdata->GetDouble("YOUNG")),
  nue_(matdata->GetDouble("NUE"))
{
  if (parmode_ == "YN")
  {
    if (nue_ < 0.5 and nue_ > -1.0)
    {
      lambda_ = youngs_*nue_/((1.0+nue_)*(1.0-2.0*nue_));  // Lame coeff.
      mue_ = youngs_/(2.0*(1.0+nue_));  // shear modulus
    }
    else
      dserror("Poisson's ratio must be between -1.0 and 0.5!");
  }
  else if (parmode_ != "Lame")
    dserror("unknown parameter set for NeoHooke material!\n Must be either YN (Young's modulus and Poisson's ratio) or Lame");
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::CoupLogNeoHooke::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::CoupLogNeoHooke( this ) );
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
  haveshearmod = true;

  shearmod += params_->mue_;

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
  havecoefficients = true;

   // determinant of deformation gradient
  const double logdetf = std::log(std::sqrt(prinv(2)));

  // gammas
  gamma(0) += params_->mue_;
  gamma(1) += 0.0;
  gamma(2) += - params_->mue_ + params_->lambda_ * logdetf;

  // deltas
  delta(5) += params_->lambda_;
  delta(6) += 2.0 * (params_->mue_ - params_->lambda_ * logdetf);

  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
