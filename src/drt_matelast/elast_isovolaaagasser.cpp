/*----------------------------------------------------------------------*/
/*!
\file elast_isovolaaagasser.cpp
\brief

This file contains the routines required to calculate the isochoric contribution
of the aaagasser material and the corresponding volumetric contribution.

MAT 20 ELAST_IsoVolAAAGasser CLUM 2.62E3 CMED 1.98E3 CABLUM 1.73E3 NUE 0.49 BETA -2.0

<pre>
maintainer: Andreas Maier
            a.maier@lnm.mw.tum.de
</pre>

*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include "elast_isovolaaagasser.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_lib/drt_linedefinition.H"

/*----------------------------------------------------------------------*
 |                                                         AMaier  06/11|
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoVolAAAGasser::IsoVolAAAGasser(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  clum_(matdata->GetDouble("CLUM")),
  cmed_(matdata->GetDouble("CMED")),
  cablum_(matdata->GetDouble("CABLUM")),
  nue_(matdata->GetDouble("NUE")),
  beta_(matdata->GetDouble("BETA"))
{
}


Teuchos::RCP<MAT::Material> MAT::ELASTIC::PAR::IsoVolAAAGasser::CreateMaterial()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 |  Constructor                                   (public) AMaier 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVolAAAGasser::IsoVolAAAGasser()
  : Summand(),
    params_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                 (public)   AMaier 04/09 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVolAAAGasser::IsoVolAAAGasser(MAT::ELASTIC::PAR::IsoVolAAAGasser* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolAAAGasser::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data,normdist_);
  AddtoPack(data,cele_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolAAAGasser::UnpackSummand(const std::vector<char>& data,
																									std::vector<char>::size_type& position)
{
  ExtractfromPack(position,data,normdist_);
  ExtractfromPack(position,data,cele_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolAAAGasser::SetupAAA(Teuchos::ParameterList& params)
{
  normdist_ = params.get("iltthick meanvalue",-999.0);

  if (normdist_==-999.0) dserror("Aneurysm mean ilt distance not found. Did you switch on 'PATSPEC'? (besides other possible errors of course)");

  if (0.0 <= normdist_ and normdist_ <= 0.5)
    cele_ = ((normdist_ - 0.5)/(-0.5))* params_->clum_ + (normdist_/0.5)*params_->cmed_;
  else if(0.5 < normdist_ and normdist_ <= 1.0)
    cele_ = ((normdist_ - 1.0)/(-0.5))*params_->cmed_ + ((normdist_ - 0.5)/0.5)*params_->cablum_;
  else
    dserror("Unable to calculate valid stiffness parameter in material AAAGasser");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolAAAGasser::AddStrainEnergy(
  double& psi,
  const LINALG::Matrix<3,1>& prinv,
  const LINALG::Matrix<3,1>& modinv)
{
  // An Ogden type material is chosen for the isochoric part
  // \f$\Psi=c\underset{i=1}{\overset{3}{\sum}}(\lambda_{i}^{4}-1)\f$
  // which is
  //Psi = c*(I_1^2*I_3^{-2/3} -2* I_2*I_3^{-2/3}-3)
  psi+=cele_*(pow(modinv(0),2.0)-2.0*modinv(1)-3.0 );
  // volumetric part
  // contribution is modeled by an Ogden-Simo_Miehe type SEF:
  // \f$\Psi=\frac {\kappa}{\beta^2}(\beta lnJ + J^{-\beta}-1)\f$
  // with kappa= 8*c/(1-2nu)
  // as Gasser paper states that referential stiffness E=24c and
  // K=24c/(3(1-2nu))

  double detF=sqrt(prinv(2));
  psi+= (8*cele_)/(1.0-2.0*params_->nue_)*1.0/(pow(params_->beta_,2.0))*(params_->beta_*log(detF)+pow(detF,-params_->beta_)-1.0);

}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolAAAGasser::AddCoefficientsPrincipal(
  LINALG::Matrix<3,1>& gamma,
  LINALG::Matrix<8,1>& delta,
  const LINALG::Matrix<3,1>& prinv)
{
  // principal coefficients
  gamma(0) += 0.;
  gamma(1) += cele_ * 4.*pow(prinv(2),-2./3.);
  gamma(2) += cele_ * (-(4./3.)*pow(prinv(2),-2./3.)*(prinv(0)*prinv(0) - 2.*prinv(1)) + (2./(1.-2.*params_->nue_))*4.*(1.-pow(prinv(2),-params_->beta_/2.))/params_->beta_);
  delta(0) += 0.0;
  delta(1) += 0.0;
  delta(2) += 0.0;
  delta(3) += 0.0;
  delta(4) += cele_ * (-16./3.)*pow(prinv(2),-2./3.);
  delta(5) += cele_ * ((16./9.)*pow(prinv(2),-2./3.)*(prinv(0)*prinv(0) - 2.*prinv(1)) + (2./(1.-2.*params_->nue_))*4.*pow(prinv(2),-params_->beta_/2.));
  delta(6) += cele_ * ((8./3.)*pow(prinv(2),-2./3.)*(prinv(0)*prinv(0) - 2.*prinv(1)) + (2./(1.-2.*params_->nue_))*4.*2.*(pow(prinv(2),-params_->beta_/2.)-1.)/params_->beta_);
  delta(7) += cele_ * 8.*pow(prinv(2),-2./3.);

  return;
}


/*----------------------------------------------------------------------*/
