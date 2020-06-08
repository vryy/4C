/*----------------------------------------------------------------------*/
/*! \file
\brief  This file contains the routines required to calculate the isochoric and volumetric
contribution of the cc neo hooke material MAT 20 ELAST_IsoVolHUDependentNeoHooke ALPHA_MAX 8.929E6
CT_MIN 30.0 CT_MAX 600.0 NUE 0.49 BETA -2.0

\level 2

\maintainer Amadeus Gebauer
*/
/*----------------------------------------------------------------------*/
#include "elast_isovolHUdependentneohooke.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke::IsoVolHUDependentNeoHooke(
    Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      alphamax_(matdata->GetDouble("ALPHA_MAX")),
      ctmin_(matdata->GetDouble("CT_MIN")),
      ctmax_(matdata->GetDouble("CT_MAX")),
      nue_(matdata->GetDouble("NUE")),
      beta_(matdata->GetDouble("BETA"))
{
}


/*----------------------------------------------------------------------*
 |  Constructor                                 (public)   AMaier 06/11 |
 *----------------------------------------------------------------------*/
MAT::ELASTIC::IsoVolHUDependentNeoHooke::IsoVolHUDependentNeoHooke(
    MAT::ELASTIC::PAR::IsoVolHUDependentNeoHooke* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::PackSummand(DRT::PackBuffer& data) const
{
  AddtoPack(data, alpha_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::UnpackSummand(
    const std::vector<char>& data, std::vector<char>::size_type& position)
{
  ExtractfromPack(position, data, alpha_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::Setup(int numgp, DRT::INPUT::LineDefinition* linedef)
{
  double HU = 0.0;
  if (linedef->HaveNamed("HU"))
  {
    linedef->ExtractDouble("HU", HU);
  }
  else
  {
    HU = -999.0;
  }

  const Teuchos::ParameterList& pslist = DRT::Problem::Instance()->PatSpecParams();
  int HUlumen = pslist.get<int>("MAXHULUMEN");

  // if HU is smaller than the threshold for calcification or smaller than
  // the Lumen HU (+10 security factor), there is no contribution at all
  if (HU <= params_->ctmin_ || HU < (HUlumen + 10))
    alpha_ = 0.;
  else if (HU > params_->ctmax_)
    alpha_ = params_->alphamax_;
  else
    alpha_ =
        0.5 * params_->alphamax_ *
        (sin(PI * (HU - params_->ctmin_) / (params_->ctmax_ - params_->ctmin_) - PI / 2) + 1.0);
}

/*----------------------------------------------------------------------
 *                                                    hemmler 09/2016  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::AddStrainEnergy(double& psi,
    const LINALG::Matrix<3, 1>& prinv, const LINALG::Matrix<3, 1>& modinv,
    const LINALG::Matrix<6, 1>& glstrain, const int gp, const int eleGID)
{
  // material Constants c and beta
  const double nue = params_->nue_;
  const double beta = params_->beta_;
  const double kappa = 2.0 * (alpha_ / (1 - 2.0 * nue));


  // isochoric part - NeoHookean
  // strain energy: psi = alpha * (I1_mod - 3.0)
  psi += alpha_ * (modinv(0) - 3.0);

  // volumetric part - Odgen-Simo-Miehe
  // strain energy: Psi = (kappa/beta^2)*(beta*ln(J)+J^(-beta)-1.0)
  if (beta != 0)
    psi += kappa / (beta * beta) * (beta * log(modinv(2)) + pow(modinv(2), -beta) - 1.);
  else
    psi += kappa / 2. * pow(std::log(modinv(2)), 2.);
}


/*----------------------------------------------------------------------
 *                                                      birzle 12/2014  */
/*----------------------------------------------------------------------*/
void MAT::ELASTIC::IsoVolHUDependentNeoHooke::AddDerivativesModified(LINALG::Matrix<3, 1>& dPmodI,
    LINALG::Matrix<6, 1>& ddPmodII, const LINALG::Matrix<3, 1>& modinv, const int gp,
    const int eleGID)
{
  dPmodI(0) += alpha_;
  dPmodI(2) += (2. * alpha_ * (1. - std::pow(modinv(2), -params_->beta_))) /
               ((1. - 2. * params_->nue_) * params_->beta_ * modinv(2));


  ddPmodII(2) +=
      (2. * alpha_ * (-1. + std::pow(modinv(2), -params_->beta_) * (1. + params_->beta_))) /
      ((1. - 2. * params_->nue_) * params_->beta_ * modinv(2) * modinv(2));

  return;
}

/*----------------------------------------------------------------------*/
