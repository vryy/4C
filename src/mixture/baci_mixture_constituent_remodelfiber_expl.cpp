/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a remodel constituent with explicit integration of the local evolution
equations
\level 3
*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_constituent_remodelfiber_expl.H"

#include "baci_global_data.H"
#include "baci_linalg_fixedsizematrix_voigt_notation.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_service.H"
#include "baci_matelast_aniso_structuraltensor_strategy.H"
#include "baci_mixture_constituent.H"
#include "baci_mixture_constituent_remodelfiber_lib.H"
#include "baci_mixture_constituent_remodelfiber_material_exponential.H"
#include "baci_mixture_constituent_remodelfiber_material_exponential_active.H"
#include "baci_mixture_growth_evolution_linear_cauchy_poisson_turnover.H"
#include "baci_utils_function_of_time.H"

#include <algorithm>
#include <cstdlib>
#include <memory>

BACI_NAMESPACE_OPEN

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

MIXTURE::PAR::MixtureConstituent_RemodelFiberExpl::MixtureConstituent_RemodelFiberExpl(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureConstituent(matdata),
      fiber_id_(matdata->GetInt("FIBER_ID") - 1),
      init_(matdata->GetInt("INIT")),
      gamma_(matdata->GetDouble("GAMMA")),
      fiber_material_id_(matdata->GetInt("FIBER_MATERIAL_ID")),
      fiber_material_(FiberMaterialFactory(fiber_material_id_)),
      growth_enabled_(matdata->GetInt("GROWTH_ENABLED")),
      poisson_decay_time_(matdata->GetDouble("DECAY_TIME")),
      growth_constant_(matdata->GetDouble("GROWTH_CONSTANT")),
      deposition_stretch_(matdata->GetDouble("DEPOSITION_STRETCH")),
      deposition_stretch_timefunc_num_(matdata->GetInt("DEPOSITION_STRETCH_TIMEFUNCT")),
      inelastic_external_deformation_(matdata->GetInt("INELASTIC_GROWTH"))
{
}

std::unique_ptr<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituent_RemodelFiberExpl::CreateConstituent(int id)
{
  return std::make_unique<MIXTURE::MixtureConstituent_RemodelFiberExpl>(this, id);
}

MIXTURE::MixtureConstituent_RemodelFiberExpl::MixtureConstituent_RemodelFiberExpl(
    MIXTURE::PAR::MixtureConstituent_RemodelFiberExpl* params, int id)
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

INPAR::MAT::MaterialType MIXTURE::MixtureConstituent_RemodelFiberExpl::MaterialType() const
{
  return INPAR::MAT::mix_remodelfiber_expl;
}

void MIXTURE::MixtureConstituent_RemodelFiberExpl::PackConstituent(
    CORE::COMM::PackBuffer& data) const
{
  MIXTURE::MixtureConstituent::PackConstituent(data);
  anisotropy_extension_.PackAnisotropy(data);

  for (const RemodelFiber<2>& fiber : remodel_fiber_) fiber.Pack(data);
}

void MIXTURE::MixtureConstituent_RemodelFiberExpl::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MIXTURE::MixtureConstituent::UnpackConstituent(position, data);
  Initialize();

  anisotropy_extension_.UnpackAnisotropy(data, position);
  for (RemodelFiber<2>& fiber : remodel_fiber_) fiber.Unpack(position, data);
}

void MIXTURE::MixtureConstituent_RemodelFiberExpl::RegisterAnisotropyExtensions(
    MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropy_extension_);
}

void MIXTURE::MixtureConstituent_RemodelFiberExpl::Initialize()
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

void MIXTURE::MixtureConstituent_RemodelFiberExpl::ReadElement(
    int numgp, INPUT::LineDefinition* linedef)
{
  MIXTURE::MixtureConstituent::ReadElement(numgp, linedef);
  Initialize();
}

void MIXTURE::MixtureConstituent_RemodelFiberExpl::Setup(Teuchos::ParameterList& params, int eleGID)
{
  MIXTURE::MixtureConstituent::Setup(params, eleGID);
  UpdateHomeostaticValues(params, eleGID);
}

void MIXTURE::MixtureConstituent_RemodelFiberExpl::UpdateElasticPart(
    const CORE::LINALG::Matrix<3, 3>& F, const CORE::LINALG::Matrix<3, 3>& iFext,
    Teuchos::ParameterList& params, const double dt, const int gp, const int eleGID)
{
  MixtureConstituent::UpdateElasticPart(F, iFext, params, dt, gp, eleGID);

  if (!params_->inelastic_external_deformation_)
  {
    dserror(
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

void MIXTURE::MixtureConstituent_RemodelFiberExpl::Update(const CORE::LINALG::Matrix<3, 3>& F,
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

void MIXTURE::MixtureConstituent_RemodelFiberExpl::RegisterOutputDataNames(
    std::unordered_map<std::string, int>& names_and_size) const
{
  MixtureConstituent::RegisterOutputDataNames(names_and_size);
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_sig_h"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_sig"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_growth_scalar"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_lambda_r"] = 1;
}

bool MIXTURE::MixtureConstituent_RemodelFiberExpl::EvaluateOutputData(
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

CORE::LINALG::Matrix<6, 1> MIXTURE::MixtureConstituent_RemodelFiberExpl::EvaluateCurrentPK2(
    int gp, int eleGID) const
{
  CORE::LINALG::Matrix<6, 1> S_stress(false);
  const double fiber_pk2 = remodel_fiber_[gp].EvaluateCurrentFiberPK2Stress();

  S_stress.Update(fiber_pk2, anisotropy_extension_.GetStructuralTensor_stress(gp, 0));

  return S_stress;
}

CORE::LINALG::Matrix<6, 6> MIXTURE::MixtureConstituent_RemodelFiberExpl::EvaluateCurrentCmat(
    const int gp, const int eleGID) const
{
  const double dPK2dlambdafsq = remodel_fiber_[gp].EvaluateDCurrentFiberPK2StressDLambdafsq();

  CORE::LINALG::Matrix<6, 6> cmat(false);
  cmat.MultiplyNT(2.0 * dPK2dlambdafsq, anisotropy_extension_.GetStructuralTensor_stress(gp, 0),
      anisotropy_extension_.GetStructuralTensor_stress(gp, 0));
  return cmat;
}

void MIXTURE::MixtureConstituent_RemodelFiberExpl::Evaluate(const CORE::LINALG::Matrix<3, 3>& F,
    const CORE::LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  if (params_->inelastic_external_deformation_)
  {
    dserror(
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

void MIXTURE::MixtureConstituent_RemodelFiberExpl::EvaluateElasticPart(
    const CORE::LINALG::Matrix<3, 3>& FM, const CORE::LINALG::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  if (!params_->inelastic_external_deformation_)
  {
    dserror(
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

double MIXTURE::MixtureConstituent_RemodelFiberExpl::GetGrowthScalar(int gp) const
{
  return remodel_fiber_[gp].EvaluateCurrentGrowthScalar();
}

double MIXTURE::MixtureConstituent_RemodelFiberExpl::EvaluateDepositionStretch(
    const double time) const
{
  if (params_->deposition_stretch_timefunc_num_ == 0)
  {
    return params_->deposition_stretch_;
  }

  return DRT::Problem::Instance()
      ->FunctionById<CORE::UTILS::FunctionOfTime>(params_->deposition_stretch_timefunc_num_ - 1)
      .Evaluate(time);
}
void MIXTURE::MixtureConstituent_RemodelFiberExpl::UpdateHomeostaticValues(
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

double MIXTURE::MixtureConstituent_RemodelFiberExpl::EvaluateLambdaf(
    const CORE::LINALG::Matrix<3, 3>& C, const int gp, const int eleGID) const
{
  return std::sqrt(C.Dot(anisotropy_extension_.GetStructuralTensor(gp, 0)));
}

double MIXTURE::MixtureConstituent_RemodelFiberExpl::EvaluateLambdaExt(
    const CORE::LINALG::Matrix<3, 3>& iFext, const int gp, const int eleGID) const
{
  return 1.0 /
         std::sqrt(EvaluateiCext(iFext).Dot(anisotropy_extension_.GetStructuralTensor(gp, 0)));
}
BACI_NAMESPACE_CLOSE
