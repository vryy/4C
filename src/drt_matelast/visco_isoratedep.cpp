/*----------------------------------------------------------------------*/
/*! \file
\brief
This file contains the routines required to calculate the isochoric contribution
of a viscous rate dependent material law, modified from Pioletti,1997
The input line should read
  MAT 1 VISCO_IsoRateDep N 1

\level 2

*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "visco_isoratedep.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoRateDep::IsoRateDep(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata), n_(matdata->GetDouble("N"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoRateDep::IsoRateDep(MAT::ELASTIC::PAR::IsoRateDep* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoRateDep::AddCoefficientsViscoModified(const LINALG::Matrix<3, 1>& modinv,
    LINALG::Matrix<8, 1>& modmu, LINALG::Matrix<33, 1>& modxi, LINALG::Matrix<7, 1>& modrateinv,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  const double n = params_->n_;

  // get time algorithmic parameters.
  double dt = params.get<double>("delta time");

  modmu(1) += 2. * n * modrateinv(1);
  modmu(2) += (2. * n * (modinv(0) - 3.)) / dt;

  modxi(1) += (4. * n) / dt;
  modxi(2) += (4. * n * (modinv(0) - 3.)) / (dt * dt);

  return;
}


/*----------------------------------------------------------------------*/
