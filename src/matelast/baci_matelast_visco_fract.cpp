/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a viscous material contribution, calculated according to an FSLS-model

\level 2
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_visco_fract.hpp"

#include "baci_mat_par_material.hpp"

BACI_NAMESPACE_OPEN

MAT::ELASTIC::PAR::Fract::Fract(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      tau_(matdata->GetDouble("TAU")),
      alpha_(matdata->GetDouble("ALPHA")),
      beta_(matdata->GetDouble("BETA"))
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
BACI_NAMESPACE_CLOSE
