/*----------------------------------------------------------------------*/
/*!
\file elast_coupneohooke.cpp

\brief
This file contains the routines required to calculate the isochoric
contribution of a CoupNeoHookean material material.
The input line should read
  MAT 1 ELAST_CoupNeoHooke YOUNG 1 NUE 1

\level 1

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
<pre>

*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */

#include <limits>

#include "elast_coupneohooke.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupNeoHooke::CoupNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  nue_(matdata->GetDouble("NUE"))
{
  // Material Constants c and beta
  c_ = youngs_/(4.0*(1.0+nue_));
  beta_ = nue_/(1.0-2.0*nue_);

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupNeoHooke::CoupNeoHooke(MAT::ELASTIC::PAR::CoupNeoHooke* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupNeoHooke::AddShearMod(
  bool& haveshearmod,  ///< non-zero shear modulus was added
  double& shearmod     ///< variable to add upon
  ) const
{
  haveshearmod = true;

  shearmod += 2*params_->c_;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupNeoHooke::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1>& glstrain,
    const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const double beta = params_->beta_;

  // strain energy: psi = c / beta * (I3^{-beta} - 1) + c * (I1 - 3)
  double psiadd = c * (prinv(0) - 3.);
  if (beta != 0) // take care of possible division by zero in case or Poisson's ratio nu = 0.0
    psiadd += (c / beta) * (std::exp(std::log(prinv(2)) * (-beta)) - 1.);
  else
    psiadd -= c*std::log(prinv(2));

  // add to overall strain energy
  psi += psiadd;
}


/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupNeoHooke::AddDerivativesPrincipal(
    LINALG::Matrix<3,1>& dPI,
    LINALG::Matrix<6,1>& ddPII,
    const LINALG::Matrix<3,1>& prinv,
    const int eleGID
)
{
  const double beta  = params_->beta_;
  const double c     = params_->c_;

  dPI(0) += c;
  // computing exp(log(a)*b) is faster than pow(a,b)
  if (prinv(2) > 0)
    {
      const double prinv2_to_beta_m1 = std::exp(std::log(prinv(2))*(-beta-1.));
      dPI(2) -= c * prinv2_to_beta_m1;
      ddPII(2) += c*(beta+1.)*prinv2_to_beta_m1/prinv(2);
    }
  else
    dPI(2) = ddPII(2) = std::numeric_limits<double>::quiet_NaN();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupNeoHooke::AddThirdDerivativesPrincipalIso(LINALG::Matrix<10,1>& dddPIII_iso,
                                                                 const LINALG::Matrix<3,1>& prinv_iso,
                                                                 const int eleGID)
{
  const double beta  = params_->beta_;
  const double c     = params_->c_;

  dddPIII_iso(2) -= c*(beta + 1.0)*(beta + 2.0)*std::pow(prinv_iso(2),-beta-3.0);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupNeoHooke::AddCoupDerivVol(const double J, double& dPj1, double& dPj2, double& dPj3)
{
  const double beta  = params_->beta_;
  const double c     = params_->c_;

  dPj1+=2.*c*pow(J,-1./3.)-2.*c*pow(J*J,-beta)/J;
  dPj2+=-2./3.*c*pow(J,-4./3.)+4.*c*pow(J*J,-beta)*beta*pow(J,-2.)+2.*c*pow(J*J,-beta)*pow(J,-2.);
  dPj3+=0.8e1/0.9e1*c*pow(J,-0.7e1/3.)-0.8e1*c*pow(J*J,-beta)*beta*beta*pow(J,-3.)-0.12e2*c*pow(J*J,-beta)*beta*pow(J,-3.)-4.*c*pow(J*J,-beta)*pow(J,-3.);
}
