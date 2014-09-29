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
#include "../drt_lib/drt_globalproblem.H" //<-- just in this material, because of special use of inv ana

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
  // routine runs when calculation is done during the inverse analysis
  // inverse analysis runs with ln(b) and ln(c),
  // but material calculates with real parameters b and c=exp(ln(c))
  // if no inverse analysis: normal parameters are used
  std::string ia = DRT::Problem::Instance()->InverseAnalysisParams().get<std::string>("INV_ANALYSIS");
  if (ia != "none")
  {
    const double a = params_ -> a_;
    const double b = exp(params_->b_);
    const double c = exp(params_->c_);

    // ln of determinant of deformation gradient
    const double logdetf = std::log(std::sqrt(prinv(2)));

    // exponential function
    const double expfunc = std::exp(b*(prinv(0)-3.0) - (2.0*b + c)*logdetf + c*(std::sqrt(prinv(2)) - 1.0) );

    // gammas
    gamma(0) += 2.0*a*b*expfunc;
    gamma(2) += (-2.0*b - c + c*std::sqrt(prinv(2)))*a*expfunc;

    // deltas
    delta(0) += 4.0*a*b*b*expfunc;
    delta(2) += (-2.0*b - c + c*std::sqrt(prinv(2)))*2.0*a*b*expfunc;
    delta(5) += (4.0*b*b + c*c + 4.0*b*c + (c- 2.0*c*c - 4.0*b*c)*std::sqrt(prinv(2)) + c*c*prinv(2))*a*expfunc;
    delta(6) += (4.0*b + 2*c - 2.0*c*std::sqrt(prinv(2)))*a*expfunc;
  }
  else
    // normal routine
  {
    const double a = params_ -> a_;
    const double b = params_ -> b_;
    const double c = params_ -> c_;

    // ln of determinant of deformation gradient
    const double logdetf = std::log(std::sqrt(prinv(2)));

    // exponential function
    const double expfunc = std::exp(b*(prinv(0)-3.0) - (2.0*b + c)*logdetf + c*(std::sqrt(prinv(2)) - 1.0) );

    // gammas
    gamma(0) += 2.0*a*b*expfunc;
    gamma(2) += (-2.0*b - c + c*std::sqrt(prinv(2)))*a*expfunc;

    // deltas
    delta(0) += 4.0*a*b*b*expfunc;
    delta(2) += (-2.0*b - c + c*std::sqrt(prinv(2)))*2.0*a*b*expfunc;
    delta(5) += (4.0*b*b + c*c + 4.0*b*c + (c- 2.0*c*c - 4.0*b*c)*std::sqrt(prinv(2)) + c*c*prinv(2))*a*expfunc;
    delta(6) += (4.0*b + 2*c - 2.0*c*std::sqrt(prinv(2)))*a*expfunc;
  }


  return;
}


/*----------------------------------------------------------------------*/
