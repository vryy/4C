/*----------------------------------------------------------------------*/
/*!
\file visco_genmax.cpp

\brief
This file contains the routines required to read the material parameters
and add the contribution of the viscous part, calculated according to a
SLS-Model.
Depending on the other stresses this viscous part is added, always using
the same tau and beta.
The input line should read
  MAT 1 VISCO_GenMax TAU 0.1 BETA 1

\level 1

<pre>
\maintainer Anna Birzle
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
