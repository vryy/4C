/*----------------------------------------------------------------------*/
/*!
\file elast_couplogneohooke.cpp
\brief


the input line should read either
  MAT 1 ELAST_CoupLogNeoHooke YOUNG 1.044E7 NUE 0.3
or
  MAT 1 ELAST_CoupLogNeoHooke MUE 1. LAMBDA 1.

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_couplogneohooke.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupLogNeoHooke::CoupLogNeoHooke(
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
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{

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
