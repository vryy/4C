/*----------------------------------------------------------------------*/
/*!
\file elast_holzapfel_cardiac.cpp
\brief

the input line should read
    MAT 1 ELAST_Holzapfel_Cardiac A 57 B 8.023 A4 18472 B4 16.026 A6 2.481 B6 11.120 A8 216 B8 11.436

<pre>
Maintainer: Andreas Nagler
            nagler@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_holzapfel_cardiac.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Holzapfel_Cardiac::Holzapfel_Cardiac(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  A_(matdata->GetDouble("A")),
  B_(matdata->GetDouble("B")),
  A4_(matdata->GetDouble("A4")),
  B4_(matdata->GetDouble("B4")),
  A6_(matdata->GetDouble("A6")),
  B6_(matdata->GetDouble("B6")),
  A8_(matdata->GetDouble("A8")),
  B8_(matdata->GetDouble("B8"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::Holzapfel_Cardiac::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::CoupAnisoExpoTwo( this ) );
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)  bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Holzapfel_Cardiac::Holzapfel_Cardiac()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Holzapfel_Cardiac::Holzapfel_Cardiac(MAT::ELASTIC::PAR::Holzapfel_Cardiac* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Holzapfel_Cardiac::AddCoefficientsModified(
        LINALG::Matrix<3,1>& gamma,  ///< necessary coefficients for first derivative
        LINALG::Matrix<5,1>& delta,  ///< necessary coefficients for second derivative
        const LINALG::Matrix<3,1>& modinv  ///< modified invariants of right Cauchy-Green tensor
    )
{

  double A=params_->A_;
  double B=params_->B_;
  double Jmod = std::pow(modinv(2),2./3.);

  gamma(0) += Jmod * A * exp( B * (Jmod * modinv(0) - 3));

  delta(0) += 2.0 * A * B * std::pow(Jmod,2.0) * exp( B * (Jmod * modinv(0) - 3));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Holzapfel_Cardiac::AddCoefficientsPrincipalAniso(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<15,1>& delta,
  const LINALG::Matrix<6,1>& prinv
  )
{

  double A4=params_->A4_;
  double B4=params_->B4_;
  double A6=params_->A6_;
  double B6=params_->B6_;
  double A8=params_->A8_;
  double B8=params_->B8_;

  gamma(0) += 2.0 * A4 * (prinv(3) - 1) * exp( B4 * std::pow(prinv(3) - 1, 2.0));
  gamma(1) += 2.0 * A6 * (prinv(4) - 1) * exp( B6 * std::pow(prinv(4) - 1, 2.0));
  gamma(2) += A8 * prinv(5) * exp( B8 * std::pow(prinv(5) , 2.0));

  delta(0) +=  4.0 * A4 * exp( B4 * std::pow((prinv(3) - 1),2.0)) * (1 + 2 * B4 * std::pow((prinv(3) - 1),2.0));
  delta(1) +=  4.0 * A6 * exp( B6 * std::pow((prinv(4) - 1),2.0))*(1 + 2 * B6 * std::pow((prinv(4) - 1),2.0));
  delta(14)+=  A8 * exp( B8 * std::pow(prinv(5),2.0))*(1 + 2 * B8 * std::pow(prinv(5),2.0));
  return;
}

/*----------------------------------------------------------------------*/
#endif // CCADISCRET
