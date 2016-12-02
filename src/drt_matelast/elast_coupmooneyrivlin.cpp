/*----------------------------------------------------------------------*/
/*!
\file elast_coupmooneyrivlin.cpp

\brief
This file contains the routines required to calculate the contribution
of a CoupMooneyRivlin material material.
The input line should read
  MAT 1 ELAST_CoupMooneyRivlin C1 1 C2 1 C3 1

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
#include "elast_coupmooneyrivlin.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::CoupMooneyRivlin::CoupMooneyRivlin(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c1_(matdata->GetDouble("C1")),
  c2_(matdata->GetDouble("C2")),
  c3_(matdata->GetDouble("C3"))
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
MAT::ELASTIC::CoupMooneyRivlin::CoupMooneyRivlin(MAT::ELASTIC::PAR::CoupMooneyRivlin* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupMooneyRivlin::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1> glstrain,
    const int eleGID)
{

  const double c1 = params_ -> c1_;
  const double c2 = params_ -> c2_;
  const double c3 = params_ -> c3_;

  // strain energy: Psi = c_1 (I1 - 3)  +  c_2 (I2 - 3)  -  (2 c_1 + 4 c_2) ln(J) + c_3 * (J - 1)^2
  // add to overall strain energy

  psi += c1*(prinv(0)-3.) + c2*(prinv(1)-3.) - (2.*c1+4.*c2)*log(sqrt(prinv(2))) + c3*pow((sqrt(prinv(2))-1.),2.);

}

/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::CoupMooneyRivlin::AddDerivativesPrincipal(
    LINALG::Matrix<3,1>& dPI,
    LINALG::Matrix<6,1>& ddPII,
    const LINALG::Matrix<3,1>& prinv,
    const int eleGID
)
{
  const double c1 = params_ -> c1_;
  const double c2 = params_ -> c2_;
  const double c3 = params_ -> c3_;

  dPI(0) += c1;
  dPI(1) += c2;
  dPI(2) += c3*(1-std::pow(prinv(2),-0.5)) - (c1+2.*c2)*std::pow(prinv(2),-1.);

  ddPII(2) += (c1+2*c2)*std::pow(prinv(2),-2.) + 0.5*c3*std::pow(prinv(2),-1.5);


  return;
}
