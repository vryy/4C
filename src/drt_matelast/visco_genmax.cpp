/*----------------------------------------------------------------------*/
/*!
\brief
This file contains the routines required to read the material parameters
and add the contribution of the viscous part, calculated according to a
SLS-Model.
Depending on the other stresses this viscous part is added, always using
the same tau and beta.
The input line should read
  MAT 1 VISCO_GenMax TAU 0.1 BETA 1 SOLVE OST

\level 1

\maintainer Amadeus Gebauer

*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "visco_genmax.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::GenMax::GenMax(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      tau_(matdata->GetDouble("TAU")),
      beta_(matdata->GetDouble("BETA")),
      solve_(*matdata->Get<std::string>("SOLVE"))
{
}

/*----------------------------------------------------------------------*
 *            Constructor Material Class                                *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::GenMax::GenMax(MAT::ELASTIC::PAR::GenMax* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::GenMax::ReadMaterialParametersVisco(
    double& tau, double& beta, double& alpha, std::string& solve)
{
  tau = params_->tau_;
  beta = params_->beta_;
  solve = params_->solve_;

  return;
}
