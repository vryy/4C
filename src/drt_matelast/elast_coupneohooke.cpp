/*----------------------------------------------------------------------*/
/*!
\file elast_coupneohooke.cpp
\brief


the input line should read
  MAT 1 ELAST_CoupNeoHooke YOUNG 1 NUE 1

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupneohooke.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupNeoHooke::CoupNeoHooke(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  youngs_(matdata->GetDouble("YOUNG")),
  nue_(matdata->GetDouble("NUE"))
{
  // Material Constants c and beta
  c_ = youngs_/(4.0*(1.0+nue_));
  beta_ = nue_/(1.0-2.0*nue_);
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::CoupNeoHooke::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::CoupNeoHooke( this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupNeoHooke::CoupNeoHooke()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupNeoHooke::CoupNeoHooke(MAT::ELASTIC::PAR::CoupNeoHooke* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupNeoHooke::AddStrainEnergy(
  double& psi,
  const LINALG::Matrix<3,1>& prinv,
  const LINALG::Matrix<3,1>& modinv)
{
  // material Constants c and beta
  const double c = params_->c_;
  const double beta = params_->beta_;

  // strain energy psi = c / beta * (I3^{-beta} - 1) + c * (I1 - 3)
  double psiadd = (c/beta) * (pow(prinv(2),-beta)-1.) + c*(prinv(0)-3.);

  // add to overall strain energy
  psi += psiadd;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupNeoHooke::AddCoefficientsPrincipal(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{

  gamma(0) += 2.*params_->c_;
  gamma(2) -= 2.*params_->c_*pow(prinv(2),-params_->beta_);

  delta(5) += 4.*params_->beta_*params_->c_*pow(prinv(2),-params_->beta_);
  delta(6) += 4.*params_->c_*pow(prinv(2),-params_->beta_);


  return;
}


/*----------------------------------------------------------------------*/
