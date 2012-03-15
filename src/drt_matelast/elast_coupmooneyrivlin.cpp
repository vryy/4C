/*----------------------------------------------------------------------*/
/*!
\file elast_coupmooneyrivlin.cpp
\brief


the input line should read
  MAT 1 ELAST_CoupMooneyRivlin C1 1 C2 1 C3 1

<pre>
Maintainer: Sophie Rausch
            rausch@lnm.mw.tum.de
            089/289 15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupmooneyrivlin.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupMooneyRivlin::CoupMooneyRivlin(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c1_(matdata->GetDouble("C1")),
  c2_(matdata->GetDouble("C2")),
  c3_(matdata->GetDouble("C3"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::CoupMooneyRivlin::CreateMaterial()
{
  return Teuchos::null;
  //return Teuchos::rcp( new MAT::ELASTIC::CoupMooneyRivlin( this ) );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupMooneyRivlin::CoupMooneyRivlin()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupMooneyRivlin::CoupMooneyRivlin(MAT::ELASTIC::PAR::CoupMooneyRivlin* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupMooneyRivlin::AddCoefficientsPrincipal(
  bool& havecoefficients,
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{
  havecoefficients = havecoefficients or true;


   // determinant of deformation gradient
  const double logdetf = std::log(std::sqrt(prinv(2)));

  const double c1 = params_ -> c1_;
  const double c2 = params_ -> c2_;
  const double c3 = params_ -> c3_;




  gamma(0) += 2.*(c1+c2*prinv(0));
  gamma(1) -= 2.*c2;
  gamma(2) += 2.*c3*logdetf*(logdetf-1.)-2.*(c1+2.*c2);

  delta(0) += 4.*c2;
  delta(5) += 2.*c3*logdetf*(2.0*logdetf-1.);
  delta(6) -= 2.*(2.*c3*logdetf*(logdetf-1.)-2.*(c1+2.*c2));
  delta(7) -= 4.*c2;


  return;
}


/*----------------------------------------------------------------------*/
#endif // CCADISCRET
