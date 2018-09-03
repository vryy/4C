/*----------------------------------------------------------------------*/
/*!
\file elast_coupblatzko.cpp

This file contains the routines required for Blatz and Ko material model
according to Holzapfel, "Nonlinear solid mechanics", 2001.
The input line should read
  MAT 1 ELAST_CoupBlatzKo MUE 1.044E7 NUE 0.3 F 0.5

\level 1

<pre>
\maintainer Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupblatzko.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupBlatzKo::CoupBlatzKo(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      mue_(matdata->GetDouble("MUE")),
      nue_(matdata->GetDouble("NUE")),
      f_(matdata->GetDouble("F"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupBlatzKo::CoupBlatzKo(MAT::ELASTIC::PAR::CoupBlatzKo* params) : params_(params) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupBlatzKo::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain, const int eleGID)
{
  // material parameters for isochoric part
  const double mue = params_->mue_;  // Shear modulus
  const double nue = params_->nue_;  // Poisson's ratio
  const double f = params_->f_;      // interpolation parameter

  // introducing beta for consistency with Holzapfel and simplification
  const double beta = nue / (1. - 2. * nue);

  // strain energy: Psi= f \frac {\mu} 2 \left[ (I_{\boldsymbol C}-3)+\frac 1
  //         {\beta} ( III_{\boldsymbol C}^{-\beta} -1) \right]
  //         +(1-f) \frac {\mu} 2 \left[\left( \frac {II_{\boldsymbol
  //         C}}{III_{\boldsymbol C}}-3 \right) + \frac 1 {\beta}
  //         (III_{\boldsymbol C}^{\beta}-1)\right]

  double psiadd =
      f * mue * 0.5 * (prinv(0) - 3.) + (1. - f) * mue * 0.5 * (prinv(1) / prinv(2) - 3.);
  if (beta != 0)  // take care of possible division by zero in case of Poisson's ratio nu = 0.0

    // add to overall strain energy
    psi += psiadd;
}



/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupBlatzKo::AddDerivativesPrincipal(LINALG::Matrix<3, 1>& dPI,
    LINALG::Matrix<6, 1>& ddPII, const LINALG::Matrix<3, 1>& prinv, const int eleGID)
{
  // material parameters for isochoric part
  const double mue = params_->mue_;  // Shear modulus
  const double nue = params_->nue_;  // Poisson's ratio
  const double f = params_->f_;      // interpolation parameter

  // introducing beta for consistency with Holzapfel and simplification
  const double beta = nue / (1. - 2. * nue);


  dPI(0) += f * mue / 2.;
  dPI(1) += (1. - f) * mue / (2. * prinv(2));
  dPI(2) += -(f * mue) / (2. * pow(prinv(2), beta + 1.)) -
            (mue * (pow(prinv(2), beta - 1.) - prinv(1) / prinv(2) / prinv(2)) * (f - 1.)) * 0.5;

  ddPII(2) += (f * mue * (beta + 1.)) / (2. * pow(prinv(2), beta + 2.)) -
              (mue * (f - 1.) *
                  ((2. * prinv(1)) / prinv(2) / prinv(2) / prinv(2) +
                      pow(prinv(2), beta - 2.) * (beta - 1.))) *
                  0.5;
  ddPII(3) -= (1. - f) * 0.5 * mue / prinv(2) / prinv(2);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupBlatzKo::AddThirdDerivativesPrincipalIso(
    LINALG::Matrix<10, 1>& dddPIII_iso, const LINALG::Matrix<3, 1>& prinv, const int eleGID)
{
  // material parameters for isochoric part
  const double mu = params_->mue_;   // Shear modulus
  const double nue = params_->nue_;  // Poisson's ratio
  const double f = params_->f_;      // interpolation parameter

  // introducing beta for consistency with Holzapfel and simplification
  const double beta = nue / (1. - 2. * nue);

  dddPIII_iso(2) +=
      -f * mu * pow(prinv(2), -beta - 3.) * beta * beta * .5 -
      1.5 * f * mu * pow(prinv(2), -beta - 3.) * beta - f * mu * pow(prinv(2), -beta - 3.) +
      (1. - f) * mu *
          (-6. * prinv(1) * pow(prinv(2), -4.) + pow(prinv(2), beta - 3.) * beta * beta -
              3. * pow(prinv(2), beta - 3.) * beta + 2. * pow(prinv(2), beta - 3.)) *
          .5;

  dddPIII_iso(8) += -(-1. + f) * mu * pow(prinv(2), -3.);
}
/*----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupBlatzKo::AddCoupDerivVol(
    const double J, double* dPj1, double* dPj2, double* dPj3, double* dPj4)
{
  // material parameters for isochoric part
  const double mu = params_->mue_;   // Shear modulus
  const double nue = params_->nue_;  // Poisson's ratio
  const double f = params_->f_;      // interpolation parameter

  // introducing beta for consistency with Holzapfel and simplification
  const double beta = nue / (1. - 2. * nue);

  // generated with maple:
  if (dPj1)
    *dPj1 += f * mu * (2. * pow(J, -1. / 3.) - 2. * pow(J * J, -beta) / J) / 2. +
             (1. - f) * mu * (-2. * pow(J, -5. / 3.) + 2. * pow(J * J, beta) / J) / 2.;
  if (dPj2)
    *dPj2 += f * mu *
                 (-2. / 3. * pow(J, -4. / 3.) + 4. * pow(J * J, -beta) * beta * pow(J, -2.) +
                     2. * pow(J * J, -beta) * pow(J, -2.)) /
                 2. +
             (1. - f) * mu *
                 (10. / 3. * pow(J, -8. / 3.) + 4. * pow(J * J, beta) * beta * pow(J, -2.) -
                     2. * pow(J * J, beta) * pow(J, -2.)) /
                 2.;
  if (dPj3)
    *dPj3 +=
        f * mu *
            (8. / 9. * pow(J, -7. / 3.) - 8. * pow(J * J, -beta) * beta * beta * pow(J, -3.) -
                12. * pow(J * J, -beta) * beta * pow(J, -3.) -
                4. * pow(J * J, -beta) * pow(J, -3.)) /
            2. +
        (1. - f) * mu *
            (-80. / 9. * pow(J, -11. / 3.) + 8. * pow(J * J, beta) * beta * beta * pow(J, -3.) -
                12. * pow(J * J, beta) * beta * pow(J, -3.) + 4. * pow(J * J, beta) * pow(J, -3.)) /
            2.;
  if (dPj4)
    *dPj4 +=
        f * mu *
            (-56. / 27. * pow(J, -10. / 3.) +
                16. * pow(J * J, -beta) * pow(beta, 3.) * pow(J, -4.) +
                48. * pow(J * J, -beta) * beta * beta * pow(J, -4.) +
                44. * pow(J * J, -beta) * beta * pow(J, -4.) +
                12. * pow(J * J, -beta) * pow(J, -4.)) /
            2. +
        (1. - f) * mu *
            (880. / 27. * pow(J, -14. / 3.) + 16. * pow(J * J, beta) * pow(beta, 3.) * pow(J, -4.) -
                48. * pow(J * J, beta) * beta * beta * pow(J, -4.) +
                44. * pow(J * J, beta) * beta * pow(J, -4.) -
                12. * pow(J * J, beta) * pow(J, -4.)) /
            2.;
}
/*----------------------------------------------------------------------*/
