/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required to calculate the contribution
of a general power-type material.
The input line should read
  MAT 1 ELAST_Coup2Pow C 1 D 1

\level 1

\maintainer Amadeus Gebauer
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coup2pow.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Coup2Pow::Coup2Pow(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), c_(matdata->GetDouble("C")), d_(matdata->GetInt("D"))
{
}


/*----------------------------------------------------------------------*
 *            Constructor Material Class                               *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Coup2Pow::Coup2Pow(MAT::ELASTIC::PAR::Coup2Pow* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup2Pow::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain, const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const int d = params_->d_;

  // strain energy: Psi = C (II_{\boldsymbol{C}}-3)^D
  // add to overall strain energy
  psi += c * pow((prinv(1) - 3.), d);
}



/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup2Pow::AddDerivativesPrincipal(LINALG::Matrix<3, 1>& dPI,
    LINALG::Matrix<6, 1>& ddPII, const LINALG::Matrix<3, 1>& prinv, const int eleGID)
{
  const double c = params_->c_;
  const int d = params_->d_;

  // If d<2 the material model is not stress free in the reference configuration
  if (d < 2)
    dserror(
        "The Elast_Coup2Pow - material only works for positive integer exponents, which are larger "
        "than two.");

  dPI(1) += c * d * pow((prinv(1) - 3.), d - 1.);

  if (d == 2)
    ddPII(1) += (c * d * d - c * d);
  else
    ddPII(1) += (c * d * d - c * d) * pow((prinv(1) - 3.), d - 2.);


  return;
}


/*----------------------------------------------------------------------*/
