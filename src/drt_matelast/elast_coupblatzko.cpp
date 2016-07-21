/*----------------------------------------------------------------------*/
/*!
\file elast_coupblatzko.cpp
\brief the input line should read
  MAT 1 ELAST_CoupBlatzKo MUE 1.044E7 NUE 0.3 F 0.5

\level 1

<pre>
\maintainer Sophie Rausch
            rausch@lnm.mw.tum.de
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
MAT::ELASTIC::PAR::CoupBlatzKo::CoupBlatzKo(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  mue_(matdata->GetDouble("MUE")),
  nue_(matdata->GetDouble("NUE")),
  f_(matdata->GetDouble("F"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupBlatzKo::CoupBlatzKo(MAT::ELASTIC::PAR::CoupBlatzKo* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupBlatzKo::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1> glstrain,
    const int eleGID)
{
  // material parameters for isochoric part
  const double mue  = params_->mue_;  // Shear modulus
  const double nue  = params_->nue_;  // Poisson's ratio
  const double f    = params_->f_;    // interpolation parameter

  // introducing beta for consistency with Holzapfel and simplification
  const double beta = nue/(1. - 2.*nue);

  // strain energy: Psi= f \frac {\mu} 2 \left[ (I_{\boldsymbol C}-3)+\frac 1
  //         {\beta} ( III_{\boldsymbol C}^{-\beta} -1) \right]
  //         +(1-f) \frac {\mu} 2 \left[\left( \frac {II_{\boldsymbol
  //         C}}{III_{\boldsymbol C}}-3 \right) + \frac 1 {\beta}
  //         (III_{\boldsymbol C}^{\beta}-1)\right]

  double psiadd = f*mue*0.5*(prinv(0)-3.) + (1.-f)*mue*0.5*(prinv(1)/prinv(2)-3.);
  if (beta != 0) // take care of possible division by zero in case of Poisson's ratio nu = 0.0

  // add to overall strain energy
  psi += psiadd;
}




/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupBlatzKo::AddDerivativesPrincipal(
    LINALG::Matrix<3,1>& dPI,
    LINALG::Matrix<6,1>& ddPII,
    const LINALG::Matrix<3,1>& prinv,
    const int eleGID
)
{
  // material parameters for isochoric part
  const double mue  = params_->mue_;  // Shear modulus
  const double nue  = params_->nue_;  // Poisson's ratio
  const double f    = params_->f_;    // interpolation parameter

  // introducing beta for consistency with Holzapfel and simplification
  const double beta = nue/(1. - 2.*nue);


  dPI(0) += f*mue/2.;
  dPI(1) += (1.-f)*mue/(2.*prinv(2));
  dPI(2) += -(f*mue)/(2.*pow(prinv(2),beta+1.))-(mue*(pow(prinv(2), beta-1.)-prinv(1)/prinv(2)/prinv(2))*(f-1.))*0.5;

  ddPII(2) +=(f*mue*(beta+1.))/(2.*pow(prinv(2), beta+2.))-(mue*(f-1.)*((2.*prinv(1))/prinv(2)/prinv(2)/prinv(2)+pow(prinv(2), beta-2.)*(beta-1.)))*0.5;
  ddPII(3) -= (1.-f)*0.5*mue/prinv(2)/prinv(2);

  return;
}


/*----------------------------------------------------------------------*/
