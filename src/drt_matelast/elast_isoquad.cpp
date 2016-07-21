/*----------------------------------------------------------------------*/
/*!
\file elast_isoquad.cpp
\brief the input line should read
  MAT 1 ELAST_IsoQuad C 1

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
#include "elast_isoquad.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoQuad::IsoQuad(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoQuad::IsoQuad(MAT::ELASTIC::PAR::IsoQuad* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoQuad::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1> glstrain,
    const int eleGID)
{
  const double c = params_ -> c_;

  // strain energy: Psi = C (\overline{I}_{\boldsymbol{C}}-3)^2.
  // add to overall strain energy
  psi += c*(modinv(0)-3.)*(modinv(0)-3.);

}

/*----------------------------------------------------------------------
 *                                                      birzle 12/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoQuad::AddDerivativesModified(
    LINALG::Matrix<3,1>& dPmodI,
    LINALG::Matrix<6,1>& ddPmodII,
    const LINALG::Matrix<3,1>& modinv,
    const int eleGID
)
{
  const double c = params_ -> c_;

  dPmodI(0) += 2.*c*(modinv(0)-3.);

  ddPmodII(0) += 2.*c;

  return;
}


/*----------------------------------------------------------------------*/
