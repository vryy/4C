/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for a viscous material contribution, calculated according to an
SLS-model


\level 1
*/
/*----------------------------------------------------------------------*/

#include "baci_matelast_visco_genmax.H"

#include "baci_mat_par_material.H"

BACI_NAMESPACE_OPEN

MAT::ELASTIC::PAR::GenMax::GenMax(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      tau_(matdata->GetDouble("TAU")),
      beta_(matdata->GetDouble("BETA")),
      solve_(*matdata->Get<std::string>("SOLVE"))
{
}

MAT::ELASTIC::GenMax::GenMax(MAT::ELASTIC::PAR::GenMax* params) : params_(params) {}

void MAT::ELASTIC::GenMax::ReadMaterialParametersVisco(
    double& tau, double& beta, double& alpha, std::string& solve)
{
  tau = params_->tau_;
  beta = params_->beta_;
  solve = params_->solve_;
}

BACI_NAMESPACE_CLOSE
