/*----------------------------------------------------------------------*/
/*!
\brief
This file contains the routines required to calculate the contribution
of a CoupMooneyRivlin material material.
The input line should read
  MAT 1 ELAST_CoupMooneyRivlin C1 1 C2 1 C3 1

\level 1

<pre>
\maintainer Fabian Braeu
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupmooneyrivlin.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupMooneyRivlin::CoupMooneyRivlin(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      c1_(matdata->GetDouble("C1")),
      c2_(matdata->GetDouble("C2")),
      c3_(matdata->GetDouble("C3"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupMooneyRivlin::CoupMooneyRivlin(MAT::ELASTIC::PAR::CoupMooneyRivlin* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupMooneyRivlin::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain, const int eleGID)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;
  const double c3 = params_->c3_;

  // strain energy: Psi = c_1 (I1 - 3)  +  c_2 (I2 - 3)  -  (2 c_1 + 4 c_2) ln(J) + c_3 * (J - 1)^2
  // add to overall strain energy

  psi += c1 * (prinv(0) - 3.) + c2 * (prinv(1) - 3.) - (2. * c1 + 4. * c2) * log(sqrt(prinv(2))) +
         c3 * pow((sqrt(prinv(2)) - 1.), 2.);
}

/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupMooneyRivlin::AddDerivativesPrincipal(LINALG::Matrix<3, 1>& dPI,
    LINALG::Matrix<6, 1>& ddPII, const LINALG::Matrix<3, 1>& prinv, const int eleGID)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;
  const double c3 = params_->c3_;

  dPI(0) += c1;
  dPI(1) += c2;
  dPI(2) += c3 * (1 - std::pow(prinv(2), -0.5)) - (c1 + 2. * c2) * std::pow(prinv(2), -1.);

  ddPII(2) += (c1 + 2 * c2) * std::pow(prinv(2), -2.) + 0.5 * c3 * std::pow(prinv(2), -1.5);


  return;
}

/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupMooneyRivlin::AddCoupDerivVol(
    const double J, double* dPj1, double* dPj2, double* dPj3, double* dPj4)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;
  const double c3 = params_->c3_;

  // generated with maple:
  if (dPj1)
    *dPj1 += 2. * c1 * pow(J, -1. / 3.) + 4. * c2 * pow(J, 1. / 3.) - (2. * c1 + 4. * c2) / J +
             2. * c3 * (sqrt(J * J) - 1.) * pow(J * J, -1. / 2.) * J;
  if (dPj2)
    *dPj2 += -2. / 3. * c1 * pow(J, -4. / 3.) + 4. / 3. * c2 * pow(J, -2. / 3.) +
             (2. * c1 + 4. * c2) * pow(J, -2.) + (2. * c3) -
             2. * c3 * (sqrt(J * J) - 1.) * pow(J * J, -3. / 2.) * J * J +
             2. * c3 * (sqrt(J * J) - 1.) * pow(J * J, -1. / 2.);
  if (dPj3)
    *dPj3 += 8. / 9. * c1 * pow(J, -0.7e1 / 3.) - 8. / 9. * c2 * pow(J, -5. / 3.) -
             2. * (2. * c1 + 4. * c2) * pow(J, -3.) +
             6. * c3 * (sqrt(J * J) - 1.) * pow(J * J, -5. / 2.) * pow(J, 3.) -
             6. * c3 * (sqrt(J * J) - 1.) * pow(J * J, -3. / 2.) * J;
  if (dPj4)
    *dPj4 += -56. / 27. * c1 * pow(J, -10. / 3.) + 40. / 27. * c2 * pow(J, -8. / 3.) +
             6. * (2. * c1 + 40 * c2) * pow(J, -4.) -
             30. * c3 * (sqrt(J * J) - 1.) * pow(J * J, -7. / 2.) * pow(J, 4.) +
             36. * c3 * (sqrt(J * J) - 1.) * pow(J * J, -5. / 2.) * J * J -
             6. * c3 * (sqrt(J * J) - 1.) * pow(J * J, -3. / 2.);
}
