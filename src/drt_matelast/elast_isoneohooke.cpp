/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the isochoric routines required for logarithmic neo-Hooke material
according to Bonet and Wood, "Nonlinear continuum mechanics for finite
element analysis", Cambridge, 1997.
The input line should read
  MAT 1 ELAST_IsoNeoHooke MUE 100

\level 2

\maintainer Amadeus Gebauer
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isoneohooke.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoNeoHooke::IsoNeoHooke(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), mue_(matdata->GetDouble("MUE"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoNeoHooke::IsoNeoHooke(MAT::ELASTIC::PAR::IsoNeoHooke* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoNeoHooke::AddShearMod(bool& haveshearmod, double& shearmod) const
{
  haveshearmod = haveshearmod or true;

  shearmod += params_->mue_;

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoNeoHooke::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain, const int eleGID)
{
  const double mue = params_->mue_;

  // strain energy: Psi = frac{\mu}{2} (\overline{I}_{\boldsymbol{C}}-3) = \frac{\mu}{2}
  // (J^{-2/3}{I}_{\boldsymbol{C}}-3) add to overall strain energy
  psi += mue * 0.5 * (modinv(0) - 3.);
}


/*----------------------------------------------------------------------
 *                                                      birzle 11/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoNeoHooke::AddDerivativesModified(LINALG::Matrix<3, 1>& dPmodI,
    LINALG::Matrix<6, 1>& ddPmodII, const LINALG::Matrix<3, 1>& modinv, const int eleGID)
{
  const double mue = params_->mue_;

  dPmodI(0) += 0.5 * mue;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// void MAT::ELASTIC::IsoNeoHooke::AddCoefficientsStretchesPrincipal(
//   LINALG::Matrix<3,1>& gamma,  ///< see above, [gamma_1, gamma_2, gamma_3]
//   LINALG::Matrix<6,1>& delta,  ///< see above, [delta_11, delta_22, delta_33, delta_12, delta_23,
//   delta_31] const LINALG::Matrix<3,1>& prstr  ///< principal stretches, [lambda_1, lambda_2,
//   lambda_3]
//   )
// {
//   // NOTE
//   // This implementation was verified to deliver identical results as
//   // the description in terms of modified principal invariants
//   // of the right Cauchy-Green tensor (see above).
//
//   // energy density
//   //   Psi = mu/2 [ \bar{lam}_1^2 + \bar{lam}_2^2 + \bar{lam}_3^2 - 3 ]
//   //       = mu/2 [ (J^{-1/3} lam_1)^2 + (J^{-1/3} lam_2)^2 + (J^{-1/3} lam_3)^2 - 3 ]
//   //       = mu/2 [ J^{-2/3} (lam_1^2 + lam_2^2 + lam_3^2) - 3 ]
//   //       = mu/2 [ J^{-2/3} I_C - 3 ]
//   // with
//   //   mu ... shear modulus
//   //   lam_\alpha ... principal stretches, \alpha=1,2,3
//   //   \bar{lam}_\alpha = J^{-1/3} lam_\alpha ... modified principal stretches
//   //   I_C = lam_1^2 + lam_2^2 + lam_3^2 ... 1st invariant of right Cauchy-Green 2-tensor
//   //   J = lam_1*lam_2*lam_3 ... determinant of deformation gradient
//
//   // shear modulus
//   const double mue = params_->mue_;
//   // determinant of deformation gradient
//   const double jac = prstr(0)*prstr(1)*prstr(2);
//   // convenience
//   const double jac23 = std::pow(jac,-2.0/3.0);
//   const double jac53 = std::pow(jac,-5.0/3.0);
//   const double jac83 = std::pow(jac,-8.0/3.0);
//   // 1st invariant of right Cauchy-Green tensor
//   const double fst = prstr(0)*prstr(0) + prstr(1)*prstr(1) + prstr(2)*prstr(2);
//
//   // first derivatives \frac{\partial Psi}{\partial\lambda_alpha}
//   gamma(0)  // ,0
//     += mue*jac23*prstr(0)
//     + (-1.0/3.0)*mue*fst*jac53*prstr(1)*prstr(2);
//   gamma(1)  // ,1
//     += mue*jac23*prstr(1)
//     + (-1.0/3.0)*mue*fst*jac53*prstr(0)*prstr(2);
//   gamma(2)  // ,2
//     += mue*jac23*prstr(2)
//     + (-1.0/3.0)*mue*fst*jac53*prstr(0)*prstr(1);
//
//   // second derivatives \frac{\partial^2 Psi}{\partial\lambda_alpha \partial\lambda_\beta}
//   delta(0)  // ,00
//     += mue*jac23
//     + (-2.0/3.0)*mue*prstr(0)*jac53*prstr(1)*prstr(2)
//     + (-2.0/3.0)*mue*prstr(0)*jac53*prstr(1)*prstr(2)
//     + (5.0/9.0)*mue*fst*jac83*prstr(1)*prstr(2)*prstr(1)*prstr(2);
//   delta(1)  // ,11
//     += mue*jac23
//     + (-2.0/3.0)*mue*prstr(1)*jac53*prstr(0)*prstr(2)
//     + (-2.0/3.0)*mue*prstr(1)*jac53*prstr(0)*prstr(2)
//     + (5.0/9.0)*mue*fst*jac83*prstr(0)*prstr(2)*prstr(0)*prstr(2);
//   delta(2)  // ,22
//     += mue*jac23
//     + (-2.0/3.0)*mue*prstr(2)*jac53*prstr(0)*prstr(1)
//     + (-2.0/3.0)*mue*prstr(2)*jac53*prstr(0)*prstr(1)
//     + (5.0/9.0)*mue*fst*jac83*prstr(0)*prstr(1)*prstr(0)*prstr(1);
//   delta(3)  // ,01
//     += (-2.0/3.0)*mue*prstr(0)*jac53*prstr(0)*prstr(2)
//     + (-2.0/3.0)*mue*prstr(1)*jac53*prstr(1)*prstr(2)
//     + (5.0/9.0)*mue*fst*jac83*prstr(1)*prstr(2)*prstr(0)*prstr(2)
//     + (-1.0/3.0)*mue*fst*jac53*prstr(2);
//   delta(4)  // ,12
//     += (-2.0/3.0)*mue*prstr(1)*jac53*prstr(0)*prstr(1)
//     + (-2.0/3.0)*mue*prstr(2)*jac53*prstr(0)*prstr(2)
//     + (5.0/9.0)*mue*fst*jac83*prstr(0)*prstr(2)*prstr(0)*prstr(1)
//     + (-1.0/3.0)*mue*fst*jac53*prstr(0);
//   delta(5)  // ,20
//     += (-2.0/3.0)*mue*prstr(2)*jac53*prstr(1)*prstr(2)
//     + (-2.0/3.0)*mue*prstr(0)*jac53*prstr(0)*prstr(1)
//     + (5.0/9.0)*mue*fst*jac83*prstr(0)*prstr(1)*prstr(1)*prstr(2)
//     + (-1.0/3.0)*mue*fst*jac53*prstr(1);
//
//   // done
//   return;
// }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// void MAT::ELASTIC::IsoNeoHooke::AddCoefficientsStretchesModified(
//   LINALG::Matrix<3,1>& gamma,  ///< see above, [gamma_1, gamma_2, gamma_3]
//   LINALG::Matrix<6,1>& delta,  ///< see above, [delta_11, delta_22, delta_33, delta_12, delta_23,
//   delta_31] const LINALG::Matrix<3,1>& modstr  ///< modified principal stretches, [lambda_1,
//   lambda_2, lambda_3]
//   )
// {
//   // NOTE
//   // This implementation was verified to deliver identical results as
//   // the description in terms of modified principal invariants
//   // of the right Cauchy-Green tensor (see above).
//
//   // energy density
//   //   Psi = mu/2 [ \bar{lam}_1^2 + \bar{lam}_2^2 + \bar{lam}_3^2 - 3 ]
//   //       = mu/2 [ (J^{-1/3} lam_1)^2 + (J^{-1/3} lam_2)^2 + (J^{-1/3} lam_3)^2 - 3 ]
//   //       = mu/2 [ J^{-2/3} (lam_1^2 + lam_2^2 + lam_3^2) - 3 ]
//   //       = mu/2 [ J^{-2/3} I_C - 3 ]
//   // with
//   //   mu ... shear modulus
//   //   lam_\alpha ... principal stretches, \alpha=1,2,3
//   //   \bar{lam}_\alpha = J^{-1/3} lam_\alpha ... modified principal stretches
//   //   I_C ... 1st invariant of right Cauchy-Green 2-tensor
//   //   J = lam_1*lam_2*lam_3 ... determinant of deformation gradient
//
//   // shear modulus
//   const double mue = params_->mue_;
//
//   // first derivatives
//   gamma(0)  // ,0
//     += mue*modstr(0);
//   gamma(1)  // ,1
//     += mue*modstr(1);
//   gamma(2)  // ,2
//     += mue*modstr(2);
///
//   // second derivatives
//   delta(0)  // ,00
//     += mue;
//   delta(1)  // ,11
//     += mue;
//   delta(2)  // ,22
//     += mue;
//   delta(3)  // ,01
//     += 0.0;
//   delta(4)  // ,12
//     += 0.0;
//   delta(5)  // ,20
//     += 0.0;
//
//   // done
//   return;
// }

/*----------------------------------------------------------------------*/
