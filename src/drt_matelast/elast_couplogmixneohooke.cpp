/*----------------------------------------------------------------------*/
/*!
\file elast_couplogmixneohooke.cpp
\brief


the input line should read either
  MAT 1 ELAST_CoupLogMixNeoHooke YOUNG 1.044E7 NUE 0.3
or
  MAT 1 ELAST_CoupLogMixNeoHooke MUE 1. LAMBDA 1.

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/
/* headers */
#include "elast_couplogmixneohooke.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupLogMixNeoHooke::CoupLogMixNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata)
{
  std::string parmode = *(matdata->Get<std::string>("MODE"));
  double c1 = matdata->GetDouble("C1");
  double c2 = matdata->GetDouble("C2");

  if (parmode == "YN")
  {
    if (c2 <= 0.5 and c2 > -1.0)
    {
      lambda_ = (c2 == 0.5) ? 0.0 : c1*c2/((1.0+c2)*(1.0-2.0*c2));
      mue_ = c1/(2.0*(1.0+c2));  // shear modulus
    }
    else
      dserror("Poisson's ratio must be between -1.0 and 0.5!");
  }
  else if (parmode == "Lame")
  {
    mue_ = c1;
    lambda_ = c2;
  }
  else
    dserror("unknown parameter set for NeoHooke material!\n Must be either YN (Young's modulus and Poisson's ratio) or Lame");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupLogMixNeoHooke::CoupLogMixNeoHooke(MAT::ELASTIC::PAR::CoupLogMixNeoHooke* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupLogMixNeoHooke::AddShearMod(
  bool& haveshearmod,  ///< non-zero shear modulus was added
  double& shearmod     ///< variable to add upon
  ) const
{
  haveshearmod = true;

  shearmod += params_->mue_;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupLogMixNeoHooke::AddCoefficientsPrincipal(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{
  const double sq = std::sqrt(prinv(2));

  // gammas
  gamma(0) += params_->mue_;
  gamma(1) += 0.0;
  gamma(2) += - params_->mue_ + params_->lambda_ * (prinv(2)-sq);

  // deltas
  delta(5) += params_->lambda_*(2.0*prinv(2)-sq);
  delta(6) += 2.0 * (params_->mue_ - params_->lambda_ * (prinv(2)-sq));

  return;
}


/*----------------------------------------------------------------------*/
