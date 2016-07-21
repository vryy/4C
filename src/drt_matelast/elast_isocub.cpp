/*----------------------------------------------------------------------*/
/*!
\file elast_isocub.cpp
\brief the input line should read
  MAT 1 ELAST_IsoCub C 1

\level 1

<pre>
\maintainer Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isocub.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoCub::IsoCub(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoCub::IsoCub(MAT::ELASTIC::PAR::IsoCub* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoCub::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1> glstrain,
    const int eleGID)
{
  // material Constant c
  const double c = params_->c_;

  // strain energy: Psi = C (\overline{I}_{\boldsymbol{C}}-3)^3.
  // add to overall strain energy
  psi += c * (modinv(0)-3.)*(modinv(0)-3.)*(modinv(0)-3.);
}


/*----------------------------------------------------------------------
 *                                                      birzle 11/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoCub::AddDerivativesModified(
    LINALG::Matrix<3,1>& dPmodI,
    LINALG::Matrix<6,1>& ddPmodII,
    const LINALG::Matrix<3,1>& modinv,
    const int eleGID
)
{

  const double c = params_ -> c_;

  dPmodI(0) += c*3.*(modinv(0)-3.)*(modinv(0)-3.);

  ddPmodII(0) += c*6.*(modinv(0)-3.);

  return;
}


/*----------------------------------------------------------------------*/
