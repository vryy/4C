/*----------------------------------------------------------------------*/
/*!
\brief
the input line should read
  MAT 1 VISCO_Fract TAU 0.1 ALPHA 0.5 BETA 1

\level 2

<pre>
\maintainer Fabian Braeu
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "visco_fract.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Fract::Fract(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      tau_(matdata->GetDouble("TAU")),
      alpha_(matdata->GetDouble("ALPHA")),
      beta_(matdata->GetDouble("BETA"))
{
}

/*----------------------------------------------------------------------*
 *            Constructor Material Class                                *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Fract::Fract(MAT::ELASTIC::PAR::Fract* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Fract::ReadMaterialParametersVisco(
    double& tau, double& beta, double& alpha, std::string& solve)
{
  tau = params_->tau_;
  alpha = params_->alpha_;
  beta = params_->beta_;

  return;
}


/*----------------------------------------------------------------------*/
