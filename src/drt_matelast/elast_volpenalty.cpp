/*----------------------------------------------------------------------*/
/*!
\file elast_volpenalty.cpp

\brief
This file contains the routines required for the volumetic function
of a penalty material
For more details see Roernbauer2008 (student thesis)
The input line should read
  MAT 1 ELAST_VolPenalty EPSILON 1. GAMMA 1.

\level 1

<pre>
\maintainer Anna Birzle
            birzle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_volpenalty.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::VolPenalty::VolPenalty(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  eps_(matdata->GetDouble("EPSILON")),
  gam_(matdata->GetDouble("GAMMA"))
{
  if (eps_<0. || gam_<=0.)
    dserror("VolPenalty parameters EPSILON and GAMMA have to be greater zero");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::VolPenalty::VolPenalty(MAT::ELASTIC::PAR::VolPenalty* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolPenalty::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1> glstrain,
    const int eleGID)
{
  const double eps = params_ -> eps_;
  const double gam = params_ -> gam_;

  // strain energy: Psi=\epsilon \left( J^{\gamma} + \frac 1 {J^{\gamma}} -2 \right)
  // add to overall strain energy
  psi += eps*(pow(modinv(2),gam) + pow(modinv(2),-gam) - 2.);

}

/*----------------------------------------------------------------------
 *                                                      birzle 11/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolPenalty::AddDerivativesModified(
    LINALG::Matrix<3,1>& dPmodI,
    LINALG::Matrix<6,1>& ddPmodII,
    const LINALG::Matrix<3,1>& modinv,
    const int eleGID

)
{
  const double eps = params_ -> eps_;
  const double gam = params_ -> gam_;

  dPmodI(2) += eps*gam*(pow(modinv(2),gam-1.)-pow(modinv(2),-gam-1.));

  ddPmodII(2) += eps*gam*((gam-1.)*pow(modinv(2),gam-2.)+(gam+1.)*pow(modinv(2),-gam-2.));

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::VolPenalty::Add3rdVolDeriv(const LINALG::Matrix<3,1>& modinv, double& d3PsiVolDJ3)
{
  const double eps = params_ -> eps_;
  const double gam = params_ -> gam_;
  const double J=modinv(2);
  d3PsiVolDJ3 += eps*(   gam *( gam-1.)*( gam-2.)*pow(J, gam-3.) +
                       (-gam)*(-gam-1.)*(-gam-2.)*pow(J,-gam-3.)
                     );
  return;
}


/*----------------------------------------------------------------------*/
