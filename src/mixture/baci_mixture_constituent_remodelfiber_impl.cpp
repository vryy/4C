/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a remodel constituent with implicit integration of the local evolution
equations
\level 3
*/
/*----------------------------------------------------------------------*/
#include "baci_mixture_constituent_remodelfiber_impl.H"

#include "baci_lib_globalproblem.H"
#include "baci_linalg_fixedsizematrix_voigt_notation.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_service.H"
#include "baci_matelast_aniso_structuraltensor_strategy.H"
#include "baci_mixture_constituent_remodelfiber_lib.H"
#include "baci_mixture_constituent_remodelfiber_material_exponential.H"
#include "baci_mixture_constituent_remodelfiber_material_exponential_active.H"
#include "baci_mixture_growth_evolution_linear_cauchy_poisson_turnover.H"
#include "baci_utils_function_of_time.H"

#include <algorithm>
#include <cstdlib>
#include <memory>

// anonymous namespace for helper classes and functions
namespace
{
  [[nodiscard]] CORE::LINALG::Matrix<3, 3> EvaluateC(const CORE::LINALG::Matrix<3, 3>& F)
  {
    CORE::LINALG::Matrix<3, 3> C(false);
    C.MultiplyTN(F, F);
    return C;
  }
}  // namespace

MIXTURE::PAR::MixtureConstituent_RemodelFiberImpl::MixtureConstituent_RemodelFiberImpl(
    const Teuchos::RCP<MAT::PAR::Material>& matdata)
    : MixtureConstituent(matdata),
      fiber_id_(matdata->GetInt("FIBER_ID") - 1),
      init_(matdata->GetInt("INIT")),
      fiber_material_id_(matdata->GetInt("FIBER_MATERIAL_ID")),
      fiber_material_(FiberMaterialFactory(fiber_material_id_)),
      growth_enabled_(matdata->GetInt("GROWTH_ENABLED")),
      poisson_decay_time_(matdata->GetDouble("DECAY_TIME")),
      growth_constant_(matdata->GetDouble("GROWTH_CONSTANT")),
      deposition_stretch_(matdata->GetDouble("DEPOSITION_STRETCH")),
      deposition_stretch_timefunc_num_(matdata->GetInt("DEPOSITION_STRETCH_TIMEFUNCT"))
{
}

std::unique_ptr<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituent_RemodelFiberImpl::CreateConstituent(int id)
{
  return std::make_unique<MIXTURE::MixtureConstituent_RemodelFiberImpl>(this, id);
}

MIXTURE::MixtureConstituent_RemodelFiberImpl::MixtureConstituent_RemodelFiberImpl(
    MIXTURE::PAR::MixtureConstituent_RemodelFiberImpl* params, int id)
    : MixtureConstituent(params, id),
      params_(params),
      remodel_fiber_(),
      anisotropy_extension_(params_->init_, 0.0, false,
          Teuchos::rcp(new MAT::ELASTIC::StructuralTensorStrategyStandard(nullptr)),
          {params_->fiber_id_})
{
  anisotropy_extension_.RegisterNeededTensors(
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS |
      MAT::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

INPAR::MAT::MaterialType MIXTURE::MixtureConstituent_RemodelFiberImpl::MaterialType() const
{
  return INPAR::MAT::mix_remodelfiber_impl;
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::PackConstituent(
    CORE::COMM::PackBuffer& data) const
{
  MIXTURE::MixtureConstituent::PackConstituent(data);
  anisotropy_extension_.PackAnisotropy(data);

  for (const RemodelFiber<2>& fiber : remodel_fiber_) fiber.Pack(data);
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MIXTURE::MixtureConstituent::UnpackConstituent(position, data);
  Initialize();

  anisotropy_extension_.UnpackAnisotropy(data, position);
  for (RemodelFiber<2>& fiber : remodel_fiber_) fiber.Unpack(position, data);
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::RegisterAnisotropyExtensions(
    MAT::Anisotropy& anisotropy)
{
  anisotropy.RegisterAnisotropyExtension(anisotropy_extension_);
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::Initialize()
{
  dgrowthscalardC_.resize(NumGP());
  dlambdardC_.resize(NumGP());
  remodel_fiber_.clear();
  std::shared_ptr<const RemodelFiberMaterial<double>> material =
      params_->fiber_material_->CreateRemodelFiberMaterial();

  for (int gp = 0; gp < NumGP(); ++gp)
  {
    LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<double> growth_evolution(
        params_->growth_constant_, params_->poisson_decay_time_);
    remodel_fiber_.emplace_back(material, growth_evolution, EvaluateDepositionStretch(0.0));
  }
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::ReadElement(
    int numgp, DRT::INPUT::LineDefinition* linedef)
{
  MIXTURE::MixtureConstituent::ReadElement(numgp, linedef);
  Initialize();
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::Setup(Teuchos::ParameterList& params, int eleGID)
{
  MIXTURE::MixtureConstituent::Setup(params, eleGID);
  UpdateHomeostaticValues(params, eleGID);
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::Update(const CORE::LINALG::Matrix<3, 3>& F,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  MixtureConstituent::Update(F, params, gp, eleGID);

  // Update state
  remodel_fiber_[gp].Update();

  UpdateHomeostaticValues(params, eleGID);
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::RegisterOutputDataNames(
    std::unordered_map<std::string, int>& names_and_size) const
{
  MixtureConstituent::RegisterOutputDataNames(names_and_size);
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_sig_h"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_sig"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_growth_scalar"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_lambda_r"] = 1;
}

bool MIXTURE::MixtureConstituent_RemodelFiberImpl::EvaluateOutputData(
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

CORE::LINALG::Matrix<1, 6> MIXTURE::MixtureConstituent_RemodelFiberImpl::EvaluateDLambdafsqDC(
    int gp, int eleGID) const
{
  CORE::LINALG::Matrix<1, 6> dLambdafDC(false);
  dLambdafDC.UpdateT(anisotropy_extension_.GetStructuralTensor_stress(gp, 0));
  return dLambdafDC;
}

CORE::LINALG::Matrix<6, 1> MIXTURE::MixtureConstituent_RemodelFiberImpl::EvaluateCurrentPK2(
    int gp, int eleGID) const
{
  CORE::LINALG::Matrix<6, 1> S_stress(false);
  const double fiber_pk2 = remodel_fiber_[gp].EvaluateCurrentFiberPK2Stress();

  S_stress.Update(fiber_pk2, anisotropy_extension_.GetStructuralTensor_stress(gp, 0));

  return S_stress;
}

CORE::LINALG::Matrix<6, 6> MIXTURE::MixtureConstituent_RemodelFiberImpl::EvaluateCurrentCmat(
    const int gp, const int eleGID) const
{
  const double dPK2dlambdafsq = remodel_fiber_[gp].EvaluateDCurrentFiberPK2StressDLambdafsq();

  CORE::LINALG::Matrix<6, 6> cmat(false);
  cmat.MultiplyNT(2.0 * dPK2dlambdafsq, anisotropy_extension_.GetStructuralTensor_stress(gp, 0),
      anisotropy_extension_.GetStructuralTensor_stress(gp, 0));

  // additional linearization from implicit integration
  if (params_->growth_enabled_)
  {
    const double dpk2dlambdar = remodel_fiber_[gp].EvaluateDCurrentFiberPK2StressDLambdar();
    cmat.MultiplyNN(2.0 * dpk2dlambdar, anisotropy_extension_.GetStructuralTensor_stress(gp, 0),
        dlambdardC_[gp], 1.0);
  }

  return cmat;
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::IntegrateLocalEvolutionEquations(
    const double dt, int gp, int eleGID)
{
  dsassert(params_->growth_enabled_,
      "The integration of the local evolution equation should only be called if growth is "
      "enabled!");

  // Integrate local evolution equations
  CORE::LINALG::Matrix<2, 2> K = remodel_fiber_[gp].IntegrateLocalEvolutionEquationsImplicit(dt);

  // Compute increment w.r.t. C
  const CORE::LINALG::Matrix<2, 6> dRdC = std::invoke(
      [&]()
      {
        CORE::LINALG::Matrix<2, 6> dRdC;
        CORE::LINALG::Matrix<1, 6> dgrowthC;
        CORE::LINALG::Matrix<1, 6> dremodelC;

        const double dRgrowthDLambdafsq =
            remodel_fiber_[gp]
                .EvaluateDCurrentGrowthEvolutionImplicitTimeIntegrationResiduumDLambdafsq(dt);
        const double dRremodelDLambdafsq =
            remodel_fiber_[gp]
                .EvaluateDCurrentRemodelEvolutionImplicitTimeIntegrationResiduumDLambdafsq(dt);

        dgrowthC.Update(dRgrowthDLambdafsq, EvaluateDLambdafsqDC(gp, eleGID));
        dremodelC.Update(dRremodelDLambdafsq, EvaluateDLambdafsqDC(gp, eleGID));

        for (std::size_t i = 0; i < 6; ++i)
        {
          dRdC(0, i) = dgrowthC(i);
          dRdC(1, i) = dremodelC(i);
        }

        return dRdC;
      });

  K.Invert();
  CORE::LINALG::Matrix<1, 2> dgrowthscalardR(false);
  CORE::LINALG::Matrix<1, 2> dlambdardR(false);

  for (std::size_t i = 0; i < 2; ++i)
  {
    dgrowthscalardR(i) = K(0, i);
    dlambdardR(i) = K(1, i);
  }

  dgrowthscalardC_[gp].Multiply(-1.0, dgrowthscalardR, dRdC);
  dlambdardC_[gp].Multiply(-1.0, dlambdardR, dRdC);
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::Evaluate(const CORE::LINALG::Matrix<3, 3>& F,
    const CORE::LINALG::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    CORE::LINALG::Matrix<6, 1>& S_stress, CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  const double dt = params.get<double>("delta time");

  CORE::LINALG::Matrix<3, 3> C = EvaluateC(F);

  const double lambda_f = EvaluateLambdaf(C, gp, eleGID);
  remodel_fiber_[gp].SetState(lambda_f, 1.0);

  if (params_->growth_enabled_) IntegrateLocalEvolutionEquations(dt, gp, eleGID);

  S_stress.Update(EvaluateCurrentPK2(gp, eleGID));
  cmat.Update(EvaluateCurrentCmat(gp, eleGID));
}

void MIXTURE::MixtureConstituent_RemodelFiberImpl::EvaluateElasticPart(
    const CORE::LINALG::Matrix<3, 3>& FM, const CORE::LINALG::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, CORE::LINALG::Matrix<6, 1>& S_stress,
    CORE::LINALG::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  dserror(
      "The implicit remodel fiber cannot be evaluated with an additional inelastic deformation. "
      "You can either use the explicit remodel fiber or use a growth strategy without an inelastic "
      "external deformation.");
}

double MIXTURE::MixtureConstituent_RemodelFiberImpl::GetGrowthScalar(int gp) const
{
  return remodel_fiber_[gp].EvaluateCurrentGrowthScalar();
}

CORE::LINALG::Matrix<1, 6> MIXTURE::MixtureConstituent_RemodelFiberImpl::GetDGrowthScalarDC(
    int gp, int eleGID) const
{
  if (!params_->growth_enabled_) return CORE::LINALG::Matrix<1, 6>(true);
  return dgrowthscalardC_[gp];
}

double MIXTURE::MixtureConstituent_RemodelFiberImpl::EvaluateDepositionStretch(
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
void MIXTURE::MixtureConstituent_RemodelFiberImpl::UpdateHomeostaticValues(
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

  const double new_lambda_pre = EvaluateDepositionStretch(time);

  for (auto& fiber : remodel_fiber_)
  {
    fiber.UpdateDepositionStretch(new_lambda_pre);
  }
}

double MIXTURE::MixtureConstituent_RemodelFiberImpl::EvaluateLambdaf(
    const CORE::LINALG::Matrix<3, 3>& C, const int gp, const int eleGID) const
{
  return std::sqrt(C.Dot(anisotropy_extension_.GetStructuralTensor(gp, 0)));
}