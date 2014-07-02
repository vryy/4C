/*----------------------------------------------------------------------*/
/*!
\file elast_coupexppol.cpp
\brief


the input line should read
  MAT 1 ELAST_CoupExpPol A 600. B 2. C 5.

<pre>
Maintainer: Anna Birzle
            birzle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coupexppol.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupExpPol::CoupExpPol(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  a_(matdata->GetDouble("A")),
  b_(matdata->GetDouble("B")),
  c_(matdata->GetDouble("C"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                               (public)  |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupExpPol::CoupExpPol(MAT::ELASTIC::PAR::CoupExpPol* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupExpPol::AddCoefficientsPrincipal(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv
  )
{
  const double a = params_ -> a_;
  const double b = params_ -> b_;
  const double c = params_ -> c_;


  // ln of determinant of deformation gradient
  const double logdetf = log(sqrt(prinv(2)));

  // exponential function
  const double expfunc = exp(b*(prinv(0)-3) - (2*b+c)*logdetf + c*(sqrt(prinv(2))-1) );

  // gammas
  gamma(0) += 2*a*b*expfunc;
  gamma(2) += (-2*b -c + c*sqrt(prinv(2)))*a*expfunc;

  // deltas
  delta(0) += 4*a*b*b*expfunc;
  delta(2) += (-2*b -c +c*sqrt(prinv(2)))*2*a*b*expfunc;
  delta(5) += (4*b*b +c*c +4*b*c +(c-2*c*c-4*b*c)*sqrt(prinv(2)) + c*c*prinv(2))*a*expfunc;
  delta(6) += (4*b +2*c -2*c*sqrt(prinv(2)))*a*expfunc;


  return;
}


/*----------------------------------------------------------------------*/
