/*----------------------------------------------------------------------*/
/*!
\file elast_volsussmanbathe.cpp

\brief
This file contains the routines required for logarithmic neo-Hooke material
according to Bonet and Wood, "Nonlinear continuum mechanics for finite
element analysis", Cambridge, 1997.
The input line should read
  MAT 1 ELAST_VolSussmanBathe KAPPA 100

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
#include "elast_volsussmanbathe.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::VolSussmanBathe::VolSussmanBathe(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  kappa_(matdata->GetDouble("KAPPA"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VolSussmanBathe::VolSussmanBathe(MAT::ELASTIC::PAR::VolSussmanBathe* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolSussmanBathe::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1>& glstrain,
    const int eleGID)
{
  const double kappa = params_ -> kappa_;

  // strain energy: Psi = \frac \kappa 2 (J-1)^2
  // add to overall strain energy
  psi += kappa*0.5 * (modinv(2)-1.)*(modinv(2)-1.);

}

/*----------------------------------------------------------------------
 *                                                      birzle 11/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolSussmanBathe::AddDerivativesModified(
    LINALG::Matrix<3,1>& dPmodI,
    LINALG::Matrix<6,1>& ddPmodII,
    const LINALG::Matrix<3,1>& modinv,
    const int eleGID
)
{
  const double kappa = params_ -> kappa_;

  dPmodI(2) += kappa*(modinv(2)-1.);

  ddPmodII(2) += kappa;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolSussmanBathe::Add3rdVolDeriv(const LINALG::Matrix<3,1>& modinv, double& d3PsiVolDJ3)
{
  d3PsiVolDJ3 += 0.;
  return;
}


/*----------------------------------------------------------------------*/
