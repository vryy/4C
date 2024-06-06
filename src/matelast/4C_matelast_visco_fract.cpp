/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of a viscous material contribution, calculated according to an FSLS-model

\level 2
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_visco_fract.hpp"

#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

Mat::Elastic::PAR::Fract::Fract(const Teuchos::RCP<Core::Mat::PAR::Material>& matdata)
    : Parameter(matdata),
      tau_(matdata->Get<double>("TAU")),
      alpha_(matdata->Get<double>("ALPHA")),
      beta_(matdata->Get<double>("BETA"))
{
}

Mat::Elastic::Fract::Fract(Mat::Elastic::PAR::Fract* params) : params_(params) {}

void Mat::Elastic::Fract::read_material_parameters_visco(
    double& tau, double& beta, double& alpha, std::string& solve)
{
  tau = params_->tau_;
  alpha = params_->alpha_;
  beta = params_->beta_;
}
FOUR_C_NAMESPACE_CLOSE
