/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required to calculate the coupled
power-type material in invariant 1 multiplicative invariant 3^(-a).
The input line should read
  MAT 1 ELAST_Coup13aPow C 1 D 1 A 0.1

\level 1

*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coup13apow.H"

#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Coup13aPow::Coup13aPow(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      c_(matdata->GetDouble("C")),
      d_(matdata->GetInt("D")),
      a_(matdata->GetDouble("A"))
{
}

/*----------------------------------------------------------------------*
 *            Constructor Material Class                                *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Coup13aPow::Coup13aPow(MAT::ELASTIC::PAR::Coup13aPow* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup13aPow::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain, const int gp,
    const int eleGID)
{
  // material Constants
  const double c = params_->c_;
  const int d = params_->d_;
  const double a = params_->a_;

  // strain energy: Psi = c (I_{\boldsymbol{C}}*(III_{\boldsymbol{C}}^(-a))-3)^d
  // add to overall strain energy
  psi += c * pow((prinv(0) * (pow(prinv(2), -a)) - 3.), d);
}


/*----------------------------------------------------------------------
 *                                                       birzle 11/2016 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup13aPow::AddDerivativesPrincipal(LINALG::Matrix<3, 1>& dPI,
    LINALG::Matrix<6, 1>& ddPII, const LINALG::Matrix<3, 1>& prinv, const int gp, const int eleGID)
{
  const double c = params_->c_;
  const int d = params_->d_;
  const double a = params_->a_;

  // If d<2 the material model is not stress free in the reference configuration
  if (d < 2)
    dserror(
        "The Elast_Coup13aPow - material only works for positive integer exponents, which are "
        "larger than two.");
  // a==0 not necessary to implement extra as equal to coup1pow
  if (a == 0) dserror("Use Elast_Coup1Pow.");
  // Material not implemented for negative exponents a
  if (a < 0) dserror("Use positive values for the exponent a.");
  // Restriction for a (due to polyconvexity)
  if (a > (1. / 3. - 1. / (2. * d)))
    std::cout << "\nWARNING: A should be smaller then 1./3.- 1/(2*D). And D should be odd."
              << std::endl;

  const double I1I3a3 = prinv(0) * pow(prinv(2), -a) - 3.;

  dPI(0) += c * d * pow(prinv(2), -a) * pow(I1I3a3, d - 1.);
  dPI(2) += -a * c * d * prinv(0) * pow(prinv(2), -a - 1.) * pow(I1I3a3, d - 1.);

  if (d == 2)
  {
    ddPII(0) += c * d * (d - 1.) * pow(prinv(2), -2. * a);
    ddPII(2) += a * (a + 1.) * c * d * prinv(0) * pow(prinv(2), -a - 2.) * pow(I1I3a3, d - 1.) +
                a * a * c * d * (d - 1.) * prinv(0) * prinv(0) * pow(prinv(2), (-2. * a - 2.));
    ddPII(4) += -a * c * d * pow(prinv(2), -a - 1.) * pow(I1I3a3, d - 1.) -
                a * c * d * (d - 1.) * prinv(0) * pow(prinv(2), (-2. * a - 1.));
  }
  else
  {
    ddPII(0) += c * d * (d - 1.) * pow(prinv(2), -2. * a) * pow(I1I3a3, d - 2.);
    ddPII(2) += a * (a + 1.) * c * d * prinv(0) * pow(prinv(2), -a - 2.) * pow(I1I3a3, d - 1.) +
                a * a * c * d * (d - 1.) * prinv(0) * prinv(0) * pow(prinv(2), (-2. * a - 2.)) *
                    pow(I1I3a3, d - 2.);
    ddPII(4) +=
        -a * c * d * pow(prinv(2), -a - 1.) * pow(I1I3a3, d - 1.) -
        a * c * d * (d - 1.) * prinv(0) * pow(prinv(2), (-2. * a - 1.)) * pow(I1I3a3, d - 2.);
  }

  return;
}

/*----------------------------------------------------------------------*/
