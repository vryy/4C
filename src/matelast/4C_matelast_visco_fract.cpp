/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a viscous material contribution, calculated according to an FSLS-model

\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_visco_fract.hpp"

#include "4C_mat_par_material.hpp"

FOUR_C_NAMESPACE_OPEN

MAT::ELASTIC::PAR::Fract::Fract(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      tau_(*matdata->Get<double>("TAU")),
      alpha_(*matdata->Get<double>("ALPHA")),
      beta_(*matdata->Get<double>("BETA"))
{
}

MAT::ELASTIC::Fract::Fract(MAT::ELASTIC::PAR::Fract* params) : params_(params) {}

void MAT::ELASTIC::Fract::ReadMaterialParametersVisco(
    double& tau, double& beta, double& alpha, std::string& solve)
{
  tau = params_->tau_;
  alpha = params_->alpha_;
  beta = params_->beta_;
}
FOUR_C_NAMESPACE_CLOSE
