/*----------------------------------------------------------------------*/
/*!
\file elast_isomooneyrivlin.cpp

\brief
This file contains the routines required to calculate the isochoric contribution
of a Mooney-Rivlin-type material.
The input line should read
MAT 1 ELAST_IsoMooneyRivlin C1 100 C2 50

\level 1

<pre>
\maintainer Fabian Braeu
</pre>

*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isomooneyrivlin.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoMooneyRivlin::IsoMooneyRivlin(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), c1_(matdata->GetDouble("C1")), c2_(matdata->GetDouble("C2"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoMooneyRivlin::IsoMooneyRivlin(MAT::ELASTIC::PAR::IsoMooneyRivlin* params)
    : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoMooneyRivlin::AddStrainEnergy(double& psi, const LINALG::Matrix<3, 1>& prinv,
    const LINALG::Matrix<3, 1>& modinv, const LINALG::Matrix<6, 1>& glstrain, const int eleGID)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;

  // strain energy: Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + C2
  // (\overline{II}_{\boldsymbol{C}}-3). add to overall strain energy
  psi += c1 * (modinv(0) - 3.) + c2 * (modinv(1) - 3.);
}

/*----------------------------------------------------------------------
 *                                                      birzle 11/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoMooneyRivlin::AddDerivativesModified(LINALG::Matrix<3, 1>& dPmodI,
    LINALG::Matrix<6, 1>& ddPmodII, const LINALG::Matrix<3, 1>& modinv, const int eleGID)
{
  const double c1 = params_->c1_;
  const double c2 = params_->c2_;

  dPmodI(0) += c1;
  dPmodI(1) += c2;

  return;
}


/*----------------------------------------------------------------------*/
