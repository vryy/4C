/*----------------------------------------------------------------------------*/
/*! \file
\brief material stores parameters for ion species in electrolyte solution. The newman material is
derived for a binary electrolyte using the electroneutrality condition to condense the non-reacting
species

\level 2

*/
/*----------------------------------------------------------------------------*/
#include "baci_mat_scl.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_utils_function_of_scalar.hpp"

#include <vector>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::Scl::Scl(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ElchSingleMat(matdata),
      valence_(matdata->GetDouble("VALENCE")),
      transnrcurve_(matdata->GetInt("TRANSNR")),
      transnrparanum_(matdata->GetInt("TRANS_PARA_NUM")),
      transnr_(*matdata->Get<std::vector<double>>("TRANS_PARA")),
      cmax_(matdata->GetDouble("MAX_CONC")),
      extrapolation_diffussion_coeff_strategy_(matdata->GetInt("EXTRAPOL_DIFF")),
      clim_(matdata->GetDouble("LIM_CONC")),
      cbulk_(matdata->GetDouble("BULK_CONC")),
      susceptibility_(matdata->GetDouble("SUSCEPT")),
      delta_nu_(matdata->GetDouble("DELTA_NU")),
      faraday_(GLOBAL::Problem::Instance()->ELCHControlParams().get<double>("FARADAY_CONSTANT")),
      epsilon_0_(GLOBAL::Problem::Instance()
                     ->ELCHControlParams()
                     .sublist("DIFFCOND")
                     .get<double>("PERMITTIVITY_VACUUM"))
{
  if (transnrparanum_ != static_cast<int>(transnr_.size()))
    dserror("number of materials %d does not fit to size of material vector %d", transnrparanum_,
        transnr_.size());

  // check if number of provided parameter is valid for a the chosen predefined function
  CheckProvidedParams(transnrcurve_, transnr_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<MAT::Material> MAT::PAR::Scl::CreateMaterial()
{
  return Teuchos::rcp(new MAT::Scl(this));
}

MAT::SclType MAT::SclType::instance_;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CORE::COMM::ParObject* MAT::SclType::Create(const std::vector<char>& data)
{
  auto* scl = new MAT::Scl();
  scl->Unpack(data);
  return scl;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Scl::Scl() : params_(nullptr) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::Scl::Scl(MAT::PAR::Scl* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Scl::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::Scl::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
  {
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::Scl*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }
  }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Scl::ComputeTransferenceNumber(const double cint) const
{
  if (TransNrCurve() < 0)
    return EvalPreDefinedFunct(TransNrCurve(), cint, TransNrParams());
  else if (TransNrCurve() == 0)
    return EvalPreDefinedFunct(-1, cint, TransNrParams());
  else
  {
    return GLOBAL::Problem::Instance()
        ->FunctionById<CORE::UTILS::FunctionOfScalar>(TransNrCurve() - 1)
        .Evaluate(cint);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Scl::ComputeFirstDerivTrans(const double cint) const
{
  if (TransNrCurve() < 0)
    return EvalFirstDerivPreDefinedFunct(TransNrCurve(), cint, TransNrParams());
  else if (TransNrCurve() == 0)
    return EvalFirstDerivPreDefinedFunct(-1, cint, TransNrParams());
  else
  {
    return GLOBAL::Problem::Instance()
        ->FunctionById<CORE::UTILS::FunctionOfScalar>(TransNrCurve() - 1)
        .EvaluateDerivative(cint, 1);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Scl::ComputeDiffusionCoefficient(
    const double concentration, const double temperature) const
{
  // L indicates the mobility factor corresponding to the linear onsager relation
  const double LRT = params_->R_ * InvValValenceFaradaySquared() *
                     ComputeConductivity(concentration, temperature) *
                     (1.0 / (1.0 - (concentration - params_->cbulk_) * params_->delta_nu_)) *
                     temperature;
  const double c_max = params_->cmax_;
  double diff_coeff = 0.0;

  // The diffusion coefficient within a SCL diverges when a fully depleted or fully occupied state
  // is reached. Extrapolation strategies are used (lin_diff_number != -1) to overcome numerical
  // instabilities and obtain a well converged solution for a wide scope of discretizations in time
  // and space. The model for the diffusion coefficient is based on derivations of Steinberger K. et
  // al. (2021).
  switch (params_->extrapolation_diffussion_coeff_strategy_)
  {
    case -1:
      // no extrapolation
      diff_coeff = LRT * c_max / (concentration * (c_max - concentration));
      break;
    case 0:
    {
      // extrapolation with 0th Taylor approximation
      const double c_lim_min = params_->clim_;
      const double c_lim_max = c_max - c_lim_min;

      // no extrapolation of D(c) is required between c_lim and c_lim_low
      //(L * R * T * c_max / (c_c * (c_max - c_c)));
      if (concentration < c_lim_max && concentration > c_lim_min)
      {
        diff_coeff = LRT * c_max / (concentration * (c_max - concentration));
      }
      // extrapolation in fully depleted region (concentration approaches 0)
      // constant function (dD/dc = 0)
      else if (concentration <= c_lim_min)
      {
        diff_coeff = LRT * c_max / (c_lim_min * (c_max - c_lim_min));
      }
      // extrapolation in fully occupied region (concentration approaches c_max)
      // constant function (dD/dc = 0)
      else
      {
        diff_coeff = LRT * c_max / (c_lim_max * (c_max - c_lim_max));
      }
      break;
    }
    default:
      dserror("Extrapolation strategy is not available");
  }
  return diff_coeff;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Scl::ComputeConcentrationDerivativeOfDiffusionCoefficient(
    const double concentration, const double temperature) const
{
  // L indicates the mobility factor corresponding to the linear onsager relation, which is
  // dependent on cation concentration
  const double LRT =
      params_->R_ * ComputeOnsagerCoefficient(concentration, temperature) * temperature;
  const double LRTderconc =
      params_->R_ * ComputeConcentrationDerivativeOfOnsagerCoefficient(concentration, temperature) *
      temperature;
  const double c_max = params_->cmax_;
  double diff_coeff_der = 0.0;

  switch (params_->extrapolation_diffussion_coeff_strategy_)
  {
    case -1:
      // no extrapolation
      diff_coeff_der = LRT * (-c_max * (c_max - 2.0 * concentration)) /
                           std::pow((concentration * (c_max - concentration)), 2) +
                       LRTderconc * c_max / (concentration * (c_max - concentration));
      break;
    case 0:
    {
      // extrapolation with 0th Taylor approximation ==> dD/dc = 0!
      const double c_lim_min = params_->clim_;
      const double c_lim_max = c_max - c_lim_min;
      // return (L * R * T * c_1 / (c_c * (c_1 - c_c)))
      // no extrapolation of D(c)/Dc is required within c_lim and c_lim_low
      if (concentration < c_lim_min && concentration > c_lim_max)
      {
        diff_coeff_der = LRT * (-c_max * (c_max - 2.0 * concentration)) /
                             std::pow((concentration * (c_max - concentration)), 2) +
                         LRTderconc * c_max / (concentration * (c_max - concentration));
      }
      break;
    }
    default:
      dserror("Linearization strategy is not available");
  }
  return diff_coeff_der;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Scl::InvValValenceFaradaySquared() const
{
  return std::pow(Valence() * params_->faraday_, -2);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Scl::ComputePermittivity() const
{
  const double susceptibility = ComputeSusceptibility();
  return ((1.0 + susceptibility) * params_->epsilon_0_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Scl::ComputeOnsagerCoefficient(
    const double concentration, const double temperature) const
{
  // onsager coefficient (mobility factor) is derived from the measurable ionic conductivity and is
  // also related to deltanu, the difference between the partial molar volumes of vacancies and
  // cations
  const double conductivity = ComputeConductivity(concentration, temperature);
  return InvValValenceFaradaySquared() * conductivity /
         (1.0 - (concentration - params_->cbulk_) * params_->delta_nu_);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double MAT::Scl::ComputeConcentrationDerivativeOfOnsagerCoefficient(
    const double concentration, const double temperature) const
{
  // derivative of mobility factor w.r.t concentration depends on the concentration dependence of
  // the conductivity and another factor deltanu, that takes volumetric effects into account
  // (usually, deltanu = 0.0)
  const double conductivity = ComputeConductivity(concentration, temperature);
  const double conductivityderconc =
      ComputeConcentrationDerivativeOfConductivity(concentration, temperature);
  const double cbulk = params_->cbulk_;
  const double delta_nu = params_->delta_nu_;

  const double onsagercoeffderconc =
      InvValValenceFaradaySquared() *
      (conductivity *
              (delta_nu / (std::pow(1.0 + cbulk * delta_nu - concentration * delta_nu, 2))) +
          conductivityderconc / (1.0 - (concentration - cbulk) * delta_nu));
  return onsagercoeffderconc;
}

BACI_NAMESPACE_CLOSE
