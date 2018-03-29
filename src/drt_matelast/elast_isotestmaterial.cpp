/*----------------------------------------------------------------------*/
/*!
\file elast_isotestmaterial.cpp

\brief
This file contains the routines required to calculate the isochoric contribution
of a Material, which is not realistic, but contains all possible derivatives of invariants.
With this material in combination with volsussmannbathe, it is possible to test all
isochoric parts of the Elasthyper-Toolbox.

\level 1

<pre>
\maintainer Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoTestMaterial::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1>& glstrain,
    const int eleGID)
{
  const double c1 = params_ -> c1_;
  const double c2 = params_ -> c2_;
  const double d = c1+2.*c2;

  // strain energy: Psi = C (\overline{I}_{\boldsymbol{C}}-3)^2.
  ///   \Psi = C1 (\overline{I}_{\boldsymbol{C}}-3) + 0.5 C1 (\overline{I}_{\boldsymbol{C}}-3)^2
  ///        + C2 (\overline{II}_{\boldsymbol{C}}-3)  + 0.5 C2 (\overline{II}_{\boldsymbol{C}}-3)^2
  ///        + D (\overline{I}_{\boldsymbol{C}}-3) (\overline{II}_{\boldsymbol{C}}-3).

  //  // add to overall strain energy
  psi += c1*(modinv(0)-3) + 0.5*c1*(modinv(0)-3)*(modinv(0)-3) +c2*(modinv(1)-3) + 0.5*c2*(modinv(1)-3)*(modinv(1)-3)
          + d*(modinv(0)-3)*(modinv(1)-3);

}

/*----------------------------------------------------------------------
 *                                                      birzle 11/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoTestMaterial::AddDerivativesModified(
    LINALG::Matrix<3,1>& dPmodI,
    LINALG::Matrix<6,1>& ddPmodII,
    const LINALG::Matrix<3,1>& modinv,
    const int eleGID
  )
{
  const double c1 = params_ -> c1_;
  const double c2 = params_ -> c2_;

  const double d = c1+2.*c2;

  dPmodI(0) += c1 + c1*(modinv(0)-3.) + d*(modinv(1)-3.);
  dPmodI(1) += c2 + d*(modinv(0)-3.) + c2*(modinv(1)-3.);

  ddPmodII(0) += c1;
  ddPmodII(1) += c2;
  ddPmodII(5) += d;


  return;
}


/*----------------------------------------------------------------------*/



