/*----------------------------------------------------------------------*/
/*!
\file visco_coupmyocard.cpp

\brief
Isochoric coupled viscous material with pseudo-potential representing the collagen and
elastin matrix surrounding the myocardial fiber (chappelle12)
The input line should read
  MAT 1 VISCO_CoupMyocard N 1

\level 2

<pre>
\maintainer Martin Pfaller
            pfaller@lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "visco_coupmyocard.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupMyocard::CoupMyocard(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  n_(matdata->GetDouble("N"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupMyocard::CoupMyocard(MAT::ELASTIC::PAR::CoupMyocard* params)
  : params_(params)
{
}

/*-----------------------------------------------------------------------/
 *                                                        pfaller Apr15 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupMyocard::AddCoefficientsViscoPrincipal(
    const LINALG::Matrix<3,1>& prinv,
    LINALG::Matrix<8,1>& mu,
    LINALG::Matrix<33,1>& xi,
    LINALG::Matrix<7,1>& rateinv,
    Teuchos::ParameterList& params,
    const int eleGID
  )
{
  // material parameter
  const double eta = params_->n_;

  // get time algorithmic parameters.
  const double dt = params.get<double>("delta time");

  // contribution: \dot{C}
  mu(2) = .5 * eta / dt;

  // contribution: id4sharp_{ijkl} = 1/2 (\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})
  xi(2) = eta / (dt*dt);

  return;
}


/*----------------------------------------------------------------------*/
