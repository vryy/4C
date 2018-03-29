/*----------------------------------------------------------------------*/
/*!
\file elast_coup1pow.cpp

\brief This file contains the routines required to calculate the contribution
of a general power-type material.
The input line should read
  MAT 1 ELAST_Coup1Pow C 1 D 1

\level 1

<pre>
\maintainer Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_coup1pow.H"
#include "../drt_mat/matpar_material.H"

/*----------------------------------------------------------------------*
 *         Constructor Material Parameter Class                         *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::Coup1Pow::Coup1Pow(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C")),
  d_(matdata->GetInt("D"))
{
}

/*----------------------------------------------------------------------*
 *            Constructor Material Class                                *
 *----------------------------------------------------------------------*/
MAT::ELASTIC::Coup1Pow::Coup1Pow(MAT::ELASTIC::PAR::Coup1Pow* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 * copy matparmas to summands
 *----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup1Pow::CopyStatInvAnaMatParams(
    std::vector< Teuchos::RCP<Epetra_Vector> > input)
{
  params_->ReturnMatparams()=input;
}

/*----------------------------------------------------------------------*
 * Add parameters of elasthyper-summand for stat inverse analysis to matparams_
 *----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup1Pow::SetStatInvAnaSummandMatParams()
{
  params_->ReturnMatparams().at(MAT::ELASTIC::PAR::coup1pow_c)->PutScalar(params_->c_);
  params_->ReturnMatparams().at(MAT::ELASTIC::PAR::coup1pow_d)->PutScalar(params_->d_);
}

/*----------------------------------------------------------------------*
 * Add parameters of elasthyper-summand for stat inverse analysis
 *----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup1Pow::AddElastOptParams(
    std::map<std::string,int>* pnames)
{
  pnames->insert(std::pair<std::string,int>("Coup1Pow_C", MAT::ELASTIC::PAR::coup1pow_c));
  pnames->insert(std::pair<std::string,int>("Coup1Pow_D", MAT::ELASTIC::PAR::coup1pow_d));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup1Pow::AddStrainEnergy(
    double& psi,
    const LINALG::Matrix<3,1>& prinv,
    const LINALG::Matrix<3,1>& modinv,
    const LINALG::Matrix<6,1>& glstrain,
    const int eleGID)
{
  // material Constants c and beta
  const double c = params_->c_;
  const int d = params_->d_;

  // strain energy: Psi = C (I_{\boldsymbol{C}}-3)^D
  // add to overall strain energy
  psi += c * pow((prinv(0) - 3.),d);

}


/*----------------------------------------------------------------------
 *                                                       birzle 12/2014 */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::Coup1Pow::AddDerivativesPrincipal(
    LINALG::Matrix<3,1>& dPI,
    LINALG::Matrix<6,1>& ddPII,
    const LINALG::Matrix<3,1>& prinv,
    const int eleGID )
{

  double c = 0.;
  int d = 0;

  // in case of stat inverse analysis use getparameter
  if (params_->ReturnMatparams().size()!=0)
  {
    c = params_->GetParameter(MAT::ELASTIC::PAR::coup1pow_c,eleGID);
    d = params_->GetParameter(MAT::ELASTIC::PAR::coup1pow_d,eleGID);
    dserror("Stat Inverse Analysis is not correct implemented with elasthyper-materials in the moment. "
        "See comments in elast_coup1pow.cpp -> AddDerivativesPrincipal.");
  }
  // in other cases (e.g. gen-inv-analysis) matparams_ does not exist
  else
  {
    c = params_ -> c_;
    d = params_ -> d_;
  }

  /*
   * This is the correct implementation for the use of stat inverse analysis
   * with elasthyper-materials.
   * However therefore the params-list params is necassary in this function.
   * This is invasive in the code; consequently I did not implement it
   * n the commited version.
   * To Do: Think about a version, that's not so invasive.
   *                                                                 abirzle 12/2017
  // in case of stat inverse analysis
  // calculate derivative of stress with respect to the parameters
  // and safe in stress
  int deriv = params.get<int>("matparderiv",-1);
  // c
  if (deriv == MAT::ELASTIC::PAR::coup1pow_c)
  {
    dPI(0) += d*pow((prinv(0)-3.),d-1.);
  }
  // normal case
  // e.g. forward problem
  else if(deriv == -1)
  {
    dPI(0) += c*d*pow((prinv(0)-3.),d-1.);

    if (d==2)
      ddPII(0) += (c*d*d-c*d);
    else
      ddPII(0) += (c*d*d-c*d)*pow((prinv(0)-3.),d-2.);
  }
  // other material --> do nothing
  */



  /* Correct implementation for stat inverse analysis replaces this part*/
  // If d<2 the material model is not stress free in the reference configuration
  if (d<2)
    dserror("The Elast_Coup1Pow - material only works for positive integer exponents, which are larger than two.");

  dPI(0) += c*d*pow((prinv(0)-3.),d-1.);

  if (d==2)
    ddPII(0) += (c*d*d-c*d);
  else
    ddPII(0) += (c*d*d-c*d)*pow((prinv(0)-3.),d-2.);

  return;
}

/*----------------------------------------------------------------------*/
