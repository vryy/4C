/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of classes for a viscous material contribution, calculated according to an
SLS-model


\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_visco_genmax.hpp"

#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

Mat::Elastic::PAR::GenMax::GenMax(const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      tau_(matdata.parameters.get<double>("TAU")),
      beta_(matdata.parameters.get<double>("BETA")),
      solve_(matdata.parameters.get<std::string>("SOLVE"))
{
}

Mat::Elastic::GenMax::GenMax(Mat::Elastic::PAR::GenMax* params) : params_(params) {}

void Mat::Elastic::GenMax::read_material_parameters_visco(
    double& tau, double& beta, double& alpha, std::string& solve)
{
  tau = params_->tau_;
  beta = params_->beta_;
  solve = params_->solve_;
}

FOUR_C_NAMESPACE_CLOSE
