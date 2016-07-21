/*----------------------------------------------------------------------*/
/*!
\file elast_coup1pow.cpp
\brief the input line should read
  MAT 1 ELAST_Coup1Pow C 1 D 1

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
#include "elast_coup1pow.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Coup1Pow::Coup1Pow(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C")),
  d_(matdata->GetInt("D"))
{
}

/*----------------------------------------------------------------------*
 *            Constructor Material Class                                *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Coup1Pow::Coup1Pow(MAT::ELASTIC::PAR::Coup1Pow* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup1Pow::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1> glstrain,
    const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const int d = params_->d_;

  // strain energy: Psi = C (I_{\boldsymbol{C}}-3)^D
  // add to overall strain energy
  psi += c * pow((prinv(0) - 3.),d);

}


/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup1Pow::AddDerivativesPrincipal(
    LINALG::Matrix<3,1>& dPI,
    LINALG::Matrix<6,1>& ddPII,
    const LINALG::Matrix<3,1>& prinv,
    const int eleGID )
{
  const double c = params_ -> c_;
  const int    d = params_ -> d_;

  // If d<2 the material model is not stress free in the reference configuration
  if (d<2)
    dserror("The Elast_Coup1Pow - material only works for positive integer exponents, which are larger than two.");

  dPI(0) += c*d*pow((prinv(0)-3.),d-1.);

  if (d==2)
    ddPII(0) += (c*d*d-c*d);
  else
    ddPII(0) += (c*d*d-c*d)*pow((prinv(0)-3.),d-2.);

  return;
}

/*----------------------------------------------------------------------*/
