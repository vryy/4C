/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a remodel constituent with explicit integration of the local evolution
equations
\level 3
*/
/*----------------------------------------------------------------------*/
#include "4C_mixture_constituent_remodelfiber_expl.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
#include "4C_mixture_constituent.hpp"
#include "4C_mixture_constituent_remodelfiber_lib.hpp"
#include "4C_mixture_constituent_remodelfiber_material_exponential.hpp"
#include "4C_mixture_constituent_remodelfiber_material_exponential_active.hpp"
#include "4C_mixture_growth_evolution_linear_cauchy_poisson_turnover.hpp"
#include "4C_utils_function_of_time.hpp"

#include <algorithm>
#include <cstdlib>
#include <memory>

FOUR_C_NAMESPACE_OPEN

// anonymous namespace for helper classes and functions
namespace
{
  [[nodiscard]] CORE::LINALG::Matrix<3, 3> EvaluateC(const CORE::LINALG::Matrix<3, 3>& F)
  {
    CORE::LINALG::Matrix<3, 3> C(false);
    C.MultiplyTN(F, F);
    return C;
  }

  [[nodiscard]] CORE::LINALG::Matrix<3, 3> EvaluateiCext(const CORE::LINALG::Matrix<3, 3>& iFext)
  {
    CORE::LINALG::Matrix<3, 3> iCext(false);
    iCext.MultiplyNT(iFext, iFext);
    return iCext;
  }
}  // namespace

MIXTURE::PAR::MixtureConstituentRemodelFiberExpl::MixtureConstituentRemodelFiberExpl(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureConstituent(matdata),
      fiber_id_(*matdata->Get<int>("FIBER_ID") - 1),
      init_(*matdata->Get<int>("INIT")),
      gamma_(*matdata->Get<double>("GAMMA")),
      fiber_material_id_(*matdata->Get<int>("FIBER_MATERIAL_ID")),
      fiber_material_(FiberMaterialFactory(fiber_material_id_)),
      growth_enabled_(*matdata->Get<bool>("GROWTH_ENABLED")),
      poisson_decay_time_(*matdata->Get<double>("DECAY_TIME")),
      growth_constant_(*matdata->Get<double>("GROWTH_CONSTANT")),
      deposition_stretch_(*matdata->Get<double>("DEPOSITION_STRETCH")),
      deposition_stretch_timefunc_num_(*matdata->Get<int>("DEPOSITION_STRETCH_TIMEFUNCT")),
      inelastic_external_deformation_(*matdata->Get<bool>("INELASTIC_GROWTH"))
{
}

std::unique_ptr<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituentRemodelFiberExpl::CreateConstituent(int id)
{
  return std::make_unique<MIXTURE::MixtureConstituentRemodelFiberExpl>(this, id);
}

MIXTURE::MixtureConstituentRemodelFiberExpl::MixtureConstituentRemodelFiberExpl(
    MIXTURE::PAR::MixtureConstituentRemodelFiberExpl* params, int id)
    : MixtureConstituent(params, id),
      params_(params),
      remodel_fiber_(),
      anisotropy_extension_(params_->init_, params_->gamma_, false,
          Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {params_->fiber_id_})
{
  anisotropy_extension_.RegisterNeededTensors(
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS |
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

INPAR::MAT::MaterialType MIXTURE::MixtureConstituentRemodelFiberExpl::MaterialType() const
{
  return INPAR::MAT::mix_remodelfiber_expl;
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::PackConstituent(
    CORE::COMM::PackBuffer& data) const
{
  MIXTURE::MixtureConstituent::PackConstituent(data);
  anisotropy_extension_.PackAnisotropy(data);

  for (const RemodelFiber<2>& fiber : remodel_fiber_) fiber.Pack(data);
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MIXTURE::MixtureConstituent::UnpackConstituent(position, data);
  Initialize();

  anisotropy_extension_.UnpackAnisotropy(data, position);
  for (RemodelFiber<2>& fiber : remodel_fiber_) fiber.Unpack(position, data);
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::RegisterAnisotropyExtensions(
    MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropy_extension_);
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::Initialize()
{
  std::shared_ptr<const RemodelFiberMaterial<double>> material =
      params_->fiber_material_->CreateRemodelFiberMaterial();

  remodel_fiber_.clear();
  for (int gp = 0; gp < NumGP(); ++gp)
  {
    LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<double> growth_evolution(
        params_->growth_constant_, params_->poisson_decay_time_);
    remodel_fiber_.emplace_back(material, growth_evolution, EvaluateDepositionStretch(0.0));
  }
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::ReadElement(
    int numgp, INPUT::LineDefinition* linedef)
{
  MIXTURE::MixtureConstituent::ReadElement(numgp, linedef);
  Initialize();
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::Setup(Teuchos::ParameterList& params, int eleGID)
{
  MIXTURE::MixtureConstituent::Setup(params, eleGID);
  UpdateHomeostaticValues(params, eleGID);
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::UpdateElasticPart(
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<3, 3>& iFext,
    Teuchos::ParameterList& params, const double dt, const int gp, const int eleGID)
{
  MixtureConstituent::UpdateElasticPart(F, iFext, params, dt, gp, eleGID);

  if (!params_->inelastic_external_deformation_)
  {
    FOUR_C_THROW(
        "You specified that there is no inelastic external deformation in the input file, but this "
        "method is only called if there is one. Probably, you are using a mixture rule with "
        "inelastic growth. You have to set INELASTIC_GROWTH to true or use a different growth "
        "rule.");
  }
  const double lambda_f = EvaluateLambdaf(EvaluateC(F), gp, eleGID);
  const double lambda_ext = EvaluateLambdaExt(iFext, gp, eleGID);
  remodel_fiber_[gp].SetState(lambda_f, lambda_ext);
  remodel_fiber_[gp].Update();

  if (params_->growth_enabled_)
  {
    const double dt = params.get<double>("delta time");

    remodel_fiber_[gp].IntegrateLocalEvolutionEquationsExplicit(dt);
  }
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::Update(const CORE::LINALG::Matrix<3, 3>& F,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  MixtureConstituent::Update(F, params, gp, eleGID);

  UpdateHomeostaticValues(params, eleGID);

  if (!params_->inelastic_external_deformation_)
  {
    // Update state
    const double lambda_f = EvaluateLambdaf(EvaluateC(F), gp, eleGID);
    remodel_fiber_[gp].SetState(lambda_f, 1.0);
    remodel_fiber_[gp].Update();

    UpdateHomeostaticValues(params, eleGID);
    if (params_->growth_enabled_)
    {
      const double dt = params.get<double>("delta time");

      remodel_fiber_[gp].IntegrateLocalEvolutionEquationsExplicit(dt);
    }
  }
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::RegisterOutputDataNames(
    std::unordered_map<std::string, int>& names_and_size) const
{
  MixtureConstituent::RegisterOutputDataNames(names_and_size);
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_sig_h"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_sig"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_growth_scalar"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_lambda_r"] = 1;
}

bool MIXTURE::MixtureConstituentRemodelFiberExpl::EvaluateOutputData(
    const std::string& name, CORE::LINALG::SerialDenseMatrix& data) const
{
  if (name == "mixture_constituent_" + std::to_string(Id()) + "_sig_h")
  {
    for (int gp = 0; gp < NumGP(); ++gp)
    {
      data(gp, 0) = remodel_fiber_[gp].EvaluateCurrentHomeostaticFiberCauchyStress();
    }
    return true;
  }
  else if (name == "mixture_constituent_" + std::to_string(Id()) + "_sig")
  {
    for (int gp = 0; gp < NumGP(); ++gp)
    {
      data(gp, 0) = remodel_fiber_[gp].EvaluateCurrentFiberCauchyStress();
    }
    return true;
  }
  else if (name == "mixture_constituent_" + std::to_string(Id()) + "_growth_scalar")
  {
    for (int gp = 0; gp < NumGP(); ++gp)
    {
      data(gp, 0) = remodel_fiber_[gp].EvaluateCurrentGrowthScalar();
    }
    return true;
  }
  else if (name == "mixture_constituent_" + std::to_string(Id()) + "_lambda_r")
  {
    for (int gp = 0; gp < NumGP(); ++gp)
    {
      data(gp, 0) = remodel_fiber_[gp].EvaluateCurrentLambdar();
    }
    return true;
  }
  return MixtureConstituent::EvaluateOutputData(name, data);
}

CORE::LINALG::Matrix<6, 1> MIXTURE::MixtureConstituentRemodelFiberExpl::EvaluateCurrentPK2(
    int gp, int eleGID) const
{
  CORE::LINALG::Matrix<6, 1> S_stress(false);
  const double fiber_pk2 = remodel_fiber_[gp].EvaluateCurrentFiberPK2Stress();

  S_stress.Update(fiber_pk2, anisotropy_extension_.GetStructuralTensor_stress(gp, 0));

  return S_stress;
}

CORE::LINALG::Matrix<6, 6> MIXTURE::MixtureConstituentRemodelFiberExpl::EvaluateCurrentCmat(
    const int gp, const int eleGID) const
{
  const double dPK2dlambdafsq = remodel_fiber_[gp].EvaluateDCurrentFiberPK2StressDLambdafsq();

  CORE::LINALG::Matrix<6, 6> cmat(false);
  cmat.MultiplyNT(2.0 * dPK2dlambdafsq, anisotropy_extension_.GetStructuralTensor_stress(gp, 0),
      anisotropy_extension_.GetStructuralTensor_stress(gp, 0));
  return cmat;
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::Evaluate(const CORE::LINALG::Matrix<3, 3>& F,
    const CORE::LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  if (params_->inelastic_external_deformation_)
  {
    FOUR_C_THROW(
        "You specified that there is inelastic external deformation in the input file, but this "
        "method is only called if there is none. Probably, you are using a mixture rule without "
        "inelastic growth. You have to set INELASTIC_GROWTH to false or use a different growth "
        "rule.");
  }

  CORE::LINALG::Matrix<3, 3> C = EvaluateC(F);

  const double lambda_f = EvaluateLambdaf(C, gp, eleGID);
  remodel_fiber_[gp].SetState(lambda_f, 1.0);

  S_stress.Update(EvaluateCurrentPK2(gp, eleGID));
  cmat.Update(EvaluateCurrentCmat(gp, eleGID));
}

void MIXTURE::MixtureConstituentRemodelFiberExpl::EvaluateElasticPart(
    const CORE::LINALG::Matrix<3, 3>& FM, const CORE::LINALG::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  if (!params_->inelastic_external_deformation_)
  {
    FOUR_C_THROW(
        "You specified that there is no inelastic external deformation in the input file, but this "
        "method is only called if there is one. Probably, you are using a mixture rule with "
        "inelastic growth. You have to set INELASTIC_GROWTH to true or use a different growth "
        "rule.");
  }

  CORE::LINALG::Matrix<3, 3> C = EvaluateC(FM);

  const double lambda_f = EvaluateLambdaf(C, gp, eleGID);
  const double lambda_ext = EvaluateLambdaExt(iFextin, gp, eleGID);
  remodel_fiber_[gp].SetState(lambda_f, lambda_ext);

  S_stress.Update(EvaluateCurrentPK2(gp, eleGID));
  cmat.Update(EvaluateCurrentCmat(gp, eleGID));
}

double MIXTURE::MixtureConstituentRemodelFiberExpl::GetGrowthScalar(int gp) const
{
  return remodel_fiber_[gp].EvaluateCurrentGrowthScalar();
}

double MIXTURE::MixtureConstituentRemodelFiberExpl::EvaluateDepositionStretch(
    const double time) const
{
  if (params_->deposition_stretch_timefunc_num_ == 0)
  {
    return params_->deposition_stretch_;
  }

  return GLOBAL::Problem::Instance()
      ->FunctionById<CORE::UTILS::FunctionOfTime>(params_->deposition_stretch_timefunc_num_ - 1)
      .Evaluate(time);
}
void MIXTURE::MixtureConstituentRemodelFiberExpl::UpdateHomeostaticValues(
    const Teuchos::ParameterList& params, const int eleGID)
{
  // Update deposition stretch / prestretch of fiber depending on time function
  const double time = std::invoke(
      [&]()
      {
        constexpr auto total_time_key = "total time";
        if (!params.isParameter(total_time_key)) return 0.0;

        const double total_time = params.get<double>(total_time_key);
        if (total_time < 0.0) return 0.0;

        return total_time;
      });

  double new_lambda_pre = EvaluateDepositionStretch(time);

  for (auto& fiber : remodel_fiber_)
  {
    fiber.UpdateDepositionStretch(new_lambda_pre);
  }
}

double MIXTURE::MixtureConstituentRemodelFiberExpl::EvaluateLambdaf(
    const CORE::LINALG::Matrix<3, 3>& C, const int gp, const int eleGID) const
{
  return std::sqrt(C.Dot(anisotropy_extension_.GetStructuralTensor(gp, 0)));
}

double MIXTURE::MixtureConstituentRemodelFiberExpl::EvaluateLambdaExt(
    const CORE::LINALG::Matrix<3, 3>& iFext, const int gp, const int eleGID) const
{
  return 1.0 /
         std::sqrt(EvaluateiCext(iFext).Dot(anisotropy_extension_.GetStructuralTensor(gp, 0)));
}
FOUR_C_NAMESPACE_CLOSE
