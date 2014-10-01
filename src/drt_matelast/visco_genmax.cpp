/*----------------------------------------------------------------------*/
/*!
\file visco_genmax.cpp
\brief


the input line should read
  MAT 1 VISCO_GenMax TAU 0.1 BETA 1

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
#include "visco_genmax.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::GenMax::GenMax(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  tau_(matdata->GetDouble("TAU")),
  beta_(matdata->GetDouble("BETA"))
{
}

/*----------------------------------------------------------------------*
 *            Constructor Material Class                                *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::GenMax::GenMax(MAT::ELASTIC::PAR::GenMax* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::GenMax::ReadMaterialParameters(
  double& tau,
  double& beta
  )
{
  tau = params_ -> tau_;
  beta = params_ -> beta_;

  return;
}


/*----------------------------------------------------------------------*/
