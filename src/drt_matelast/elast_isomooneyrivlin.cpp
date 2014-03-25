/*----------------------------------------------------------------------*/
/*!
\file elast_isomooneyrivlin.cpp
\brief

This file contains the routines required to calculate the isochoric contribution
of a Mooney-Rivlin-type material

<pre>
Maintainer: Sophie Rausch & Thomas Kloeppel
            rausch,kloeppel@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isomooneyrivlin.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoMooneyRivlin::IsoMooneyRivlin(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c1_(matdata->GetDouble("C1")),
  c2_(matdata->GetDouble("C2"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                             (public)   bborn 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoMooneyRivlin::IsoMooneyRivlin(MAT::ELASTIC::PAR::IsoMooneyRivlin* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoMooneyRivlin::AddCoefficientsModified(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{

  const double c1 = params_ -> c1_;
  const double c2 = params_ -> c2_;

  gamma(0) += 2*c1+2*modinv(0)*c2;
  gamma(1) += - 2*c2;

  delta(0) += 4*c2;
  delta(3) += - 4*c2;

  return;
}


/*----------------------------------------------------------------------*/
