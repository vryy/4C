/*----------------------------------------------------------------------*/
/*!
\file visco_coupmyocard.cpp
\brief


the input line should read
  MAT 1 VISCO_CoupMyocard N 1

<pre>
Maintainer: Martin Pfaller, pfaller@lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "visco_coupmyocard.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupMyocard::CoupMyocard(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  n_(matdata->GetDouble("N"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupMyocard::CoupMyocard(MAT::ELASTIC::PAR::CoupMyocard* params)
  : params_(params)
{
}

/*-----------------------------------------------------------------------/
 *                                                        pfaller Apr15 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupMyocard::AddCoefficientsViscoModified(
  const LINALG::Matrix<3,1>& modinv,
  LINALG::Matrix<8,1>& modmu,
  LINALG::Matrix<33,1>& modxi,
  LINALG::Matrix<7,1>& modrateinv,
  Teuchos::ParameterList& params,
  const int eleGID
  )
{
  // material parameter
  const double n = params_->n_;

  // get time algorithmic parameters.
  const double dt = params.get<double>("delta time");

  // HACK: pretend this is a vol/iso split material
  // (at the moment only split materials can be handled by viscoelasthyper)
  // \Psi = \eta/2 J_1^2 =  = \eta/2 \bar{J}_1^2 I_3^{2/3} = \eta/2 \bar{J}_1^2 \bar{I}_3^{4/3}

  // a_hack = \bar{I}_3^{4/3}
  const double a_hack = std::pow(modinv(2),4./3.);

  // my(1) += 2. * n/dt * rateinv(0);
  modmu(1) += 2. * n/dt * modrateinv(0) * a_hack;

  // xi(3) += 4. * n/(dt*dt);
  modxi(3) += 4. * n/(dt*dt) * a_hack;

  return;
}


/*----------------------------------------------------------------------*/
