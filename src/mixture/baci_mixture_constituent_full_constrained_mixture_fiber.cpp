/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a quasi 1D full constrained mixture constituent
equations
\level 3
*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_constituent_full_constrained_mixture_fiber.H"

#include "baci_comm_parobject.H"
#include "baci_inpar_material.H"
#include "baci_lib_globalproblem.H"
#include "baci_mat_par_bundle.H"
#include "baci_matelast_aniso_structuraltensor_strategy.H"
#include "baci_mixture_constituent_remodelfiber_lib.H"
#include "baci_mixture_growth_evolution_linear_cauchy_poisson_turnover.H"
#include "baci_utils_function_of_time.H"

#include <algorithm>
#include <cstdlib>
#include <memory>
#include <numeric>

BACI_NAMESPACE_OPEN

// anonymous namespace for helper classes and functions
namespace
{
  [[nodiscard]] static inline CORE::LINALG::Matrix<3, 3> EvaluateC(
      const CORE::LINALG::Matrix<3, 3>& F)
  {
    CORE::LINALG::Matrix<3, 3> C(false);
    C.MultiplyTN(F, F);
    return C;
  }

  [[nodiscard]] static inline double GetTotalTime(const Teuchos::ParameterList& params)
  {
    double time = params.get<double>("total time");
    if (time < 0) return 0.0;  // Time has not been set by the time integrator during setup
    return time;
  }

  [[nodiscard]] static inline double GetDeltaTime(const Teuchos::ParameterList& params)
  {
    return params.get<double>("delta time");
  }

  MIXTURE::HistoryAdaptionStrategy GetHistoryAdaptionStrategyFromInput(const std::string& input)
  {
    if (input == "none")
    {
      return MIXTURE::HistoryAdaptionStrategy::none;
    }
    else if (input == "window")
    {
      return MIXTURE::HistoryAdaptionStrategy::window;
    }
    else if (input == "model_equation")
    {
      return MIXTURE::HistoryAdaptionStrategy::model_equation;
    }
    else if (input == "higher_order")
    {
      return MIXTURE::HistoryAdaptionStrategy::higher_order_integration;
    }
    else
    {
      dserror("Unknown history adaption strategy %s!", input.c_str());
    }
  }
}  // namespace

MIXTURE::PAR::MixtureConstituent_FullConstrainedMixtureFiber::
    MixtureConstituent_FullConstrainedMixtureFiber(const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureConstituent(matdata),
      fiber_id_(matdata->GetInt("FIBER_ID") - 1),
      init_(matdata->GetInt("INIT")),
      fiber_material_id_(matdata->GetInt("FIBER_MATERIAL_ID")),
      fiber_material_(FiberMaterialFactory(fiber_material_id_)),
      growth_enabled_(matdata->GetInt("GROWTH_ENABLED")),
      poisson_decay_time_(matdata->GetDouble("DECAY_TIME")),
      growth_constant_(matdata->GetDouble("GROWTH_CONSTANT")),
      deposition_stretch_(matdata->GetDouble("DEPOSITION_STRETCH")),
      initial_deposition_stretch_timefunc_num_(
          matdata->GetInt("INITIAL_DEPOSITION_STRETCH_TIMEFUNCT")),
      adaptive_history_strategy_(GetHistoryAdaptionStrategyFromInput(
          *matdata->Get<std::string>("ADAPTIVE_HISTORY_STRATEGY"))),
      adaptive_history_tolerance_(matdata->GetDouble("ADAPTIVE_HISTORY_TOLERANCE"))
{
}

std::unique_ptr<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituent_FullConstrainedMixtureFiber::CreateConstituent(int id)
{
  return std::make_unique<MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber>(this, id);
}

MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::
    MixtureConstituent_FullConstrainedMixtureFiber(
        MIXTURE::PAR::MixtureConstituent_FullConstrainedMixtureFiber* params, int id)
    : MixtureConstituent(params, id),
      params_(params),
      full_constrained_mixture_fiber_(),
      anisotropy_extension_(params_->init_, 0.0, false,
          Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {params_->fiber_id_})
{
  anisotropy_extension_.RegisterNeededTensors(
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS |
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

INPAR::MAT::MaterialType MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::MaterialType()
    const
{
  return INPAR::MAT::mix_full_constrained_mixture_fiber;
}

void MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::PackConstituent(
    CORE::COMM::PackBuffer& data) const
{
  MIXTURE::MixtureConstituent::PackConstituent(data);
  anisotropy_extension_.PackAnisotropy(data);

  for (const FullConstrainedMixtureFiber<double>& fiber : full_constrained_mixture_fiber_)
    fiber.Pack(data);

  CORE::COMM::ParObject::AddtoPack(data, last_lambda_f_);
}

void MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MIXTURE::MixtureConstituent::UnpackConstituent(position, data);
  Initialize();

  anisotropy_extension_.UnpackAnisotropy(data, position);

  for (FullConstrainedMixtureFiber<double>& fiber : full_constrained_mixture_fiber_)
    fiber.Unpack(position, data);

  CORE::COMM::ParObject::ExtractfromPack(position, data, last_lambda_f_);

  if (params_->growth_enabled_)
  {
    for (auto gp = 0; gp < NumGP(); ++gp)
    {
      full_constrained_mixture_fiber_[gp].ReinitializeHistory(
          last_lambda_f_[gp], full_constrained_mixture_fiber_[gp].GetLastTimeInHistory());
    }
  }
}

void MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::RegisterAnisotropyExtensions(
    MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropy_extension_);
}

void MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::Initialize()
{
  full_constrained_mixture_fiber_.clear();
  std::shared_ptr<const RemodelFiberMaterial<double>> material =
      params_->fiber_material_->CreateRemodelFiberMaterial();

  last_lambda_f_.resize(NumGP(), 1.0);

  for (int gp = 0; gp < NumGP(); ++gp)
  {
    LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<double> growth_evolution(
        params_->growth_constant_, params_->poisson_decay_time_);
    full_constrained_mixture_fiber_.emplace_back(material, growth_evolution,
        EvaluateInitialDepositionStretch(0.0), params_->adaptive_history_strategy_,
        params_->growth_enabled_);
    full_constrained_mixture_fiber_[gp].adaptive_tolerance_ = params_->adaptive_history_tolerance_;
  }
}

void MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::ReadElement(
    int numgp, INPUT::LineDefinition* linedef)
{
  MIXTURE::MixtureConstituent::ReadElement(numgp, linedef);
  Initialize();
}

void MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::Setup(
    Teuchos::ParameterList& params, int eleGID)
{
  MIXTURE::MixtureConstituent::Setup(params, eleGID);

  if (params_->growth_enabled_)
  {
    for (auto& fiber : full_constrained_mixture_fiber_)
    {
      // No deformation at t=0
      fiber.ReinitializeHistory(1.0, 0.0);
    }
  }
}

void MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::Update(
    const CORE::LINALG::Matrix<3, 3>& F, Teuchos::ParameterList& params, const int gp,
    const int eleGID)
{
  MixtureConstituent::Update(F, params, gp, eleGID);

  const double time = GetTotalTime(params);
  full_constrained_mixture_fiber_[gp].SetDepositionStretch(EvaluateInitialDepositionStretch(time));
  last_lambda_f_[gp] = EvaluateLambdaf(EvaluateC(F), gp, eleGID);

  // Update state
  full_constrained_mixture_fiber_[gp].Update();
}

void MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::RegisterOutputDataNames(
    std::unordered_map<std::string, int>& names_and_size) const
{
  MixtureConstituent::RegisterOutputDataNames(names_and_size);
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_sig_h"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_sig"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_growth_scalar"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_history_size"] = 1;
}

bool MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::EvaluateOutputData(
    const std::string& name, CORE::LINALG::SerialDenseMatrix& data) const
{
  if (name == "mixture_constituent_" + std::to_string(Id()) + "_sig_h")
  {
    for (int gp = 0; gp < NumGP(); ++gp)
    {
      data(gp, 0) = full_constrained_mixture_fiber_[gp].sig_h_;
    }
    return true;
  }
  else if (name == "mixture_constituent_" + std::to_string(Id()) + "_sig")
  {
    for (int gp = 0; gp < NumGP(); ++gp)
    {
      data(gp, 0) = full_constrained_mixture_fiber_[gp].computed_sigma_;
    }
    return true;
  }
  else if (name == "mixture_constituent_" + std::to_string(Id()) + "_growth_scalar")
  {
    for (int gp = 0; gp < NumGP(); ++gp)
    {
      data(gp, 0) = full_constrained_mixture_fiber_[gp].computed_growth_scalar_;
    }
    return true;
  }
  else if (name == "mixture_constituent_" + std::to_string(Id()) + "_history_size")
  {
    for (int gp = 0; gp < NumGP(); ++gp)
    {
      data(gp, 0) = std::accumulate(full_constrained_mixture_fiber_[gp].history_.begin(),
          full_constrained_mixture_fiber_[gp].history_.end(), 0,
          [](std::size_t number, const DepositionHistoryInterval<double>& item)
          { return number + item.timesteps.size(); });
    }
    return true;
  }
  return MixtureConstituent::EvaluateOutputData(name, data);
}

CORE::LINALG::Matrix<1, 6>
MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::EvaluateDLambdafsqDC(
    int gp, int eleGID) const
{
  CORE::LINALG::Matrix<1, 6> dLambdafDC(false);
  dLambdafDC.UpdateT(anisotropy_extension_.GetStructuralTensor_stress(gp, 0));
  return dLambdafDC;
}

CORE::LINALG::Matrix<6, 1>
MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::EvaluateCurrentPK2(
    int gp, int eleGID) const
{
  CORE::LINALG::Matrix<6, 1> S_stress(false);
  const double fiber_pk2 = full_constrained_mixture_fiber_[gp].EvaluateCurrentSecondPKStress();

  S_stress.Update(fiber_pk2, anisotropy_extension_.GetStructuralTensor_stress(gp, 0));

  return S_stress;
}

CORE::LINALG::Matrix<6, 6>
MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::EvaluateCurrentCmat(
    const int gp, const int eleGID) const
{
  const double dPK2dlambdafsq =
      full_constrained_mixture_fiber_[gp].EvaluateDCurrentFiberPK2StressDLambdafsq();

  CORE::LINALG::Matrix<6, 6> cmat(false);
  cmat.MultiplyNN(2.0 * dPK2dlambdafsq, anisotropy_extension_.GetStructuralTensor_stress(gp, 0),
      EvaluateDLambdafsqDC(gp, eleGID));

  return cmat;
}

void MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::Evaluate(
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<6, 1>& E_strain,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  const double time = GetTotalTime(params);
  const double delta_time = GetDeltaTime(params);

  CORE::LINALG::Matrix<3, 3> C = EvaluateC(F);

  const double lambda_f = EvaluateLambdaf(C, gp, eleGID);
  full_constrained_mixture_fiber_[gp].RecomputeState(lambda_f, time, delta_time);

  S_stress.Update(EvaluateCurrentPK2(gp, eleGID));
  cmat.Update(EvaluateCurrentCmat(gp, eleGID));
}

void MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::EvaluateElasticPart(
    const CORE::LINALG::Matrix<3, 3>& FM, const CORE::LINALG::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  dserror(
      "The full constrained mixture fiber cannot be evaluated with an additional inelastic "
      "deformation.");
}

double MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::GetGrowthScalar(int gp) const
{
  return full_constrained_mixture_fiber_[gp].computed_growth_scalar_;
}

CORE::LINALG::Matrix<1, 6>
MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::GetDGrowthScalarDC(
    int gp, int eleGID) const
{
  if (!params_->growth_enabled_) return CORE::LINALG::Matrix<1, 6>(true);
  CORE::LINALG::Matrix<1, 6> dGrowthScalarDE = EvaluateDLambdafsqDC(gp, eleGID);
  dGrowthScalarDE.Scale(
      2.0 * full_constrained_mixture_fiber_[gp].computed_dgrowth_scalar_dlambda_f_sq_);
  return dGrowthScalarDE;
}

double MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::EvaluateInitialDepositionStretch(
    const double time) const
{
  if (params_->initial_deposition_stretch_timefunc_num_ == 0)
  {
    return params_->deposition_stretch_;
  }

  return DRT::Problem::Instance()
      ->FunctionById<CORE::UTILS::FunctionOfTime>(
          params_->initial_deposition_stretch_timefunc_num_ - 1)
      .Evaluate(time);
}

double MIXTURE::MixtureConstituent_FullConstrainedMixtureFiber::EvaluateLambdaf(
    const CORE::LINALG::Matrix<3, 3>& C, const int gp, const int eleGID) const
{
  return std::sqrt(C.Dot(anisotropy_extension_.GetStructuralTensor(gp, 0)));
}
BACI_NAMESPACE_CLOSE
