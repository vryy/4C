/*----------------------------------------------------------------------*/
/*!
\file visco_isoratedep.cpp
\brief


the input line should read
  MAT 1 VISCO_IsoRateDep N 1

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
#include "visco_isoratedep.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoRateDep::IsoRateDep(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  n_(matdata->GetDouble("N"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoRateDep::IsoRateDep(MAT::ELASTIC::PAR::IsoRateDep* params)
  : params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoRateDep::AddCoefficientsViscoModified(
  const LINALG::Matrix<3,1>& modinv,
  LINALG::Matrix<8,1>& modmy,
  LINALG::Matrix<33,1>& modxi,
  LINALG::Matrix<7,1>& modrateinv,
  Teuchos::ParameterList& params,
  const int eleGID
  )
{

  const double n = params_ -> n_;

  // get time algorithmic parameters.
  double dt = params.get<double>("delta time");

  modmy(1) += 2.* n * modrateinv(1) ;
  modmy(2) += (2.* n *(modinv(0)-3.) ) / dt;

  modxi(1) += (4.* n) / dt;
  modxi(2) += (4.* n *(modinv(0)-3.) ) / (dt*dt);

  return;
}


/*----------------------------------------------------------------------*/
