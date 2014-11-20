/*----------------------------------------------------------------------*/
/*!
\file elast_isotestmaterial.cpp
\brief

This file contains the routines required to calculate the isochoric contribution
of a Material, which is not realistic, but contains all possible derivations of invariants.
With this material in combination with volsussmannbathe, it is possible to test all
isochoric parts of the Elasthyper-Toolbox.


<pre>
Maintainer: Anna Birzle
            birzle@lnm.mw.tum.de
            089/289 15255
</pre>

*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isotestmaterial.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoTestMaterial::IsoTestMaterial(
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
MAT::ELASTIC::IsoTestMaterial::IsoTestMaterial(MAT::ELASTIC::PAR::IsoTestMaterial* params)
  : params_(params)
{
}

///*----------------------------------------------------------------------
// *                                                      birzle 11/2014  */
///*----------------------------------------------------------------------*/
//void MAT::ELASTIC::IsoTestMaterial::AddDerivationsModified(
//    LINALG::Matrix<3,1>& dPmodI,
//    LINALG::Matrix<6,1>& ddPmodII,
//    const LINALG::Matrix<3,1>& modinv
//  )
//{
//  const double c1 = params_ -> c1_;
//  const double c2 = params_ -> c2_;
//
//  const double d = c1+2.*c2;
//
//  dPmodI(0) += c1 + c1*(modinv(0)-3.) + d*(modinv(1)-3.);
//  dPmodI(1) += c2 + d*(modinv(0)-3.) + c2*(modinv(1)-3.) ;
//
//  ddPmodII(0) += c1;
//  ddPmodII(1) += c2;
//  ddPmodII(5) += d;
//
//  return;
//}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoTestMaterial::AddCoefficientsModified(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<5,1>& delta,
  const LINALG::Matrix<3,1>& modinv
  )
{

  const double c1 = params_ -> c1_;
  const double c2 = params_ -> c2_;

  const double d = c1+2.*c2;

  gamma(0) += 2.*(c1 + c1*(modinv(0)-3.) + d*(modinv(1)-3.) + modinv(0)*(c2 + c2*(modinv(1)-3.) + d*(modinv(0)-3.)));
  gamma(1) -= 2.*(c2 + c2*(modinv(1)-3.) + d*(modinv(0)-3.));

  delta(0) += 4.*(c1 + 2.*modinv(0)*d + c2 + c2*(modinv(1)-3.) + d*(modinv(0)-3.) + modinv(0)*modinv(0)*c2);
  delta(1) -= 4.*(d + modinv(0)*c2);
  delta(2) += 4.*c2;
  delta(3) -= 4.*(c2 + c2*(modinv(1)-3.) + d*(modinv(0)-3.));

  return;
}


/*----------------------------------------------------------------------*/
