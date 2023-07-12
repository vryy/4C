/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an isochoric coupled viscous material with pseudo-potential representing
the collagen and elastin matrix surrounding the myocardial fiber (chappelle12)

\level 2
*/
/*----------------------------------------------------------------------*/

#include "matelast_visco_coupmyocard.H"
#include "mat_par_material.H"


MAT::ELASTIC::PAR::CoupMyocard::CoupMyocard(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata), n_(matdata->GetDouble("N"))
{
}

MAT::ELASTIC::CoupMyocard::CoupMyocard(MAT::ELASTIC::PAR::CoupMyocard* params) : params_(params) {}

void MAT::ELASTIC::CoupMyocard::AddCoefficientsViscoPrincipal(
    const CORE::LINALG::Matrix<3, 1>& prinv, CORE::LINALG::Matrix<8, 1>& mu,
    CORE::LINALG::Matrix<33, 1>& xi, CORE::LINALG::Matrix<7, 1>& rateinv,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  // material parameter
  const double eta = params_->n_;

  // get time algorithmic parameters.
  const double dt = params.get<double>("delta time");

  // contribution: \dot{C}
  mu(2) = .5 * eta;

  // contribution: id4sharp_{ijkl} = 1/2 (\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})
  xi(2) = eta / dt;
}