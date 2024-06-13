/*----------------------------------------------------------------------*/
/*! \file
\brief Definition of a remodel constituent with implicit integration of the local evolution
equations
\level 3
*/
/*----------------------------------------------------------------------*/
#include "4C_mixture_constituent_remodelfiber_impl.hpp"

#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_matelast_aniso_structuraltensor_strategy.hpp"
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
  [[nodiscard]] Core::LinAlg::Matrix<3, 3> EvaluateC(const Core::LinAlg::Matrix<3, 3>& F)
  {
    Core::LinAlg::Matrix<3, 3> C(false);
    C.MultiplyTN(F, F);
    return C;
  }
}  // namespace

MIXTURE::PAR::MixtureConstituentRemodelFiberImpl::MixtureConstituentRemodelFiberImpl(
    const Teuchos::RCP<Core::Mat::PAR::Material>& matdata)
    : MixtureConstituent(matdata),
      fiber_id_(matdata->Get<int>("FIBER_ID") - 1),
      init_(matdata->Get<int>("INIT")),
      fiber_material_id_(matdata->Get<int>("FIBER_MATERIAL_ID")),
      fiber_material_(FiberMaterialFactory(fiber_material_id_)),
      growth_enabled_(matdata->Get<bool>("GROWTH_ENABLED")),
      poisson_decay_time_(matdata->Get<double>("DECAY_TIME")),
      growth_constant_(matdata->Get<double>("GROWTH_CONSTANT")),
      deposition_stretch_(matdata->Get<double>("DEPOSITION_STRETCH")),
      deposition_stretch_timefunc_num_(matdata->Get<int>("DEPOSITION_STRETCH_TIMEFUNCT"))
{
}

std::unique_ptr<MIXTURE::MixtureConstituent>
MIXTURE::PAR::MixtureConstituentRemodelFiberImpl::CreateConstituent(int id)
{
  return std::make_unique<MIXTURE::MixtureConstituentRemodelFiberImpl>(this, id);
}

MIXTURE::MixtureConstituentRemodelFiberImpl::MixtureConstituentRemodelFiberImpl(
    MIXTURE::PAR::MixtureConstituentRemodelFiberImpl* params, int id)
    : MixtureConstituent(params, id),
      params_(params),
      remodel_fiber_(),
      anisotropy_extension_(params_->init_, 0.0, false,
          Teuchos::rcp(new Mat::Elastic::StructuralTensorStrategyStandard(nullptr)),
          {params_->fiber_id_})
{
  anisotropy_extension_.register_needed_tensors(
      Mat::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR_STRESS |
      Mat::FiberAnisotropyExtension<1>::STRUCTURAL_TENSOR);
}

Core::Materials::MaterialType MIXTURE::MixtureConstituentRemodelFiberImpl::MaterialType() const
{
  return Core::Materials::mix_remodelfiber_impl;
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::PackConstituent(
    Core::Communication::PackBuffer& data) const
{
  MIXTURE::MixtureConstituent::PackConstituent(data);
  anisotropy_extension_.pack_anisotropy(data);

  for (const RemodelFiber<2>& fiber : remodel_fiber_) fiber.Pack(data);
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::UnpackConstituent(
    std::vector<char>::size_type& position, const std::vector<char>& data)
{
  MIXTURE::MixtureConstituent::UnpackConstituent(position, data);
  initialize();

  anisotropy_extension_.unpack_anisotropy(data, position);
  for (RemodelFiber<2>& fiber : remodel_fiber_) fiber.Unpack(position, data);
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::register_anisotropy_extensions(
    Mat::Anisotropy& anisotropy)
{
  anisotropy.register_anisotropy_extension(anisotropy_extension_);
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::initialize()
{
  dgrowthscalard_c_.resize(num_gp());
  dlambdard_c_.resize(num_gp());
  remodel_fiber_.clear();
  std::shared_ptr<const RemodelFiberMaterial<double>> material =
      params_->fiber_material_->create_remodel_fiber_material();

  for (int gp = 0; gp < num_gp(); ++gp)
  {
    LinearCauchyGrowthWithPoissonTurnoverGrowthEvolution<double> growth_evolution(
        params_->growth_constant_, params_->poisson_decay_time_);
    remodel_fiber_.emplace_back(material, growth_evolution, evaluate_deposition_stretch(0.0));
  }
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::ReadElement(
    int numgp, Input::LineDefinition* linedef)
{
  MIXTURE::MixtureConstituent::ReadElement(numgp, linedef);
  initialize();
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::Setup(Teuchos::ParameterList& params, int eleGID)
{
  MIXTURE::MixtureConstituent::Setup(params, eleGID);
  update_homeostatic_values(params, eleGID);
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::Update(const Core::LinAlg::Matrix<3, 3>& F,
    Teuchos::ParameterList& params, const int gp, const int eleGID)
{
  MixtureConstituent::Update(F, params, gp, eleGID);

  // Update state
  remodel_fiber_[gp].Update();

  update_homeostatic_values(params, eleGID);
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  MixtureConstituent::register_output_data_names(names_and_size);
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_sig_h"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_sig"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_growth_scalar"] = 1;
  names_and_size["mixture_constituent_" + std::to_string(Id()) + "_lambda_r"] = 1;
}

bool MIXTURE::MixtureConstituentRemodelFiberImpl::EvaluateOutputData(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  if (name == "mixture_constituent_" + std::to_string(Id()) + "_sig_h")
  {
    for (int gp = 0; gp < num_gp(); ++gp)
    {
      data(gp, 0) = remodel_fiber_[gp].evaluate_current_homeostatic_fiber_cauchy_stress();
    }
    return true;
  }
  else if (name == "mixture_constituent_" + std::to_string(Id()) + "_sig")
  {
    for (int gp = 0; gp < num_gp(); ++gp)
    {
      data(gp, 0) = remodel_fiber_[gp].evaluate_current_fiber_cauchy_stress();
    }
    return true;
  }
  else if (name == "mixture_constituent_" + std::to_string(Id()) + "_growth_scalar")
  {
    for (int gp = 0; gp < num_gp(); ++gp)
    {
      data(gp, 0) = remodel_fiber_[gp].evaluate_current_growth_scalar();
    }
    return true;
  }
  else if (name == "mixture_constituent_" + std::to_string(Id()) + "_lambda_r")
  {
    for (int gp = 0; gp < num_gp(); ++gp)
    {
      data(gp, 0) = remodel_fiber_[gp].evaluate_current_lambdar();
    }
    return true;
  }
  return MixtureConstituent::EvaluateOutputData(name, data);
}

Core::LinAlg::Matrix<1, 6> MIXTURE::MixtureConstituentRemodelFiberImpl::evaluate_d_lambdafsq_dc(
    int gp, int eleGID) const
{
  Core::LinAlg::Matrix<1, 6> dLambdafDC(false);
  dLambdafDC.UpdateT(anisotropy_extension_.get_structural_tensor_stress(gp, 0));
  return dLambdafDC;
}

Core::LinAlg::Matrix<6, 1> MIXTURE::MixtureConstituentRemodelFiberImpl::evaluate_current_p_k2(
    int gp, int eleGID) const
{
  Core::LinAlg::Matrix<6, 1> S_stress(false);
  const double fiber_pk2 = remodel_fiber_[gp].evaluate_current_fiber_p_k2_stress();

  S_stress.Update(fiber_pk2, anisotropy_extension_.get_structural_tensor_stress(gp, 0));

  return S_stress;
}

Core::LinAlg::Matrix<6, 6> MIXTURE::MixtureConstituentRemodelFiberImpl::evaluate_current_cmat(
    const int gp, const int eleGID) const
{
  const double dPK2dlambdafsq =
      remodel_fiber_[gp].evaluate_d_current_fiber_p_k2_stress_d_lambdafsq();

  Core::LinAlg::Matrix<6, 6> cmat(false);
  cmat.MultiplyNT(2.0 * dPK2dlambdafsq, anisotropy_extension_.get_structural_tensor_stress(gp, 0),
      anisotropy_extension_.get_structural_tensor_stress(gp, 0));

  // additional linearization from implicit integration
  if (params_->growth_enabled_)
  {
    const double dpk2dlambdar = remodel_fiber_[gp].evaluate_d_current_fiber_p_k2_stress_d_lambdar();
    cmat.MultiplyNN(2.0 * dpk2dlambdar, anisotropy_extension_.get_structural_tensor_stress(gp, 0),
        dlambdard_c_[gp], 1.0);
  }

  return cmat;
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::integrate_local_evolution_equations(
    const double dt, int gp, int eleGID)
{
  FOUR_C_ASSERT(params_->growth_enabled_,
      "The integration of the local evolution equation should only be called if growth is "
      "enabled!");

  // Integrate local evolution equations
  Core::LinAlg::Matrix<2, 2> K =
      remodel_fiber_[gp].integrate_local_evolution_equations_implicit(dt);

  // Compute increment w.r.t. C
  const Core::LinAlg::Matrix<2, 6> dRdC = std::invoke(
      [&]()
      {
        Core::LinAlg::Matrix<2, 6> dRdC;
        Core::LinAlg::Matrix<1, 6> dgrowthC;
        Core::LinAlg::Matrix<1, 6> dremodelC;

        const double dRgrowthDLambdafsq =
            remodel_fiber_[gp]
                .evaluate_d_current_growth_evolution_implicit_time_integration_residuum_d_lambdafsq(
                    dt);
        const double dRremodelDLambdafsq =
            remodel_fiber_[gp]
                .evaluate_d_current_remodel_evolution_implicit_time_integration_residuum_d_lambdafsq(
                    dt);

        dgrowthC.Update(dRgrowthDLambdafsq, evaluate_d_lambdafsq_dc(gp, eleGID));
        dremodelC.Update(dRremodelDLambdafsq, evaluate_d_lambdafsq_dc(gp, eleGID));

        for (std::size_t i = 0; i < 6; ++i)
        {
          dRdC(0, i) = dgrowthC(i);
          dRdC(1, i) = dremodelC(i);
        }

        return dRdC;
      });

  K.Invert();
  Core::LinAlg::Matrix<1, 2> dgrowthscalardR(false);
  Core::LinAlg::Matrix<1, 2> dlambdardR(false);

  for (std::size_t i = 0; i < 2; ++i)
  {
    dgrowthscalardR(i) = K(0, i);
    dlambdardR(i) = K(1, i);
  }

  dgrowthscalard_c_[gp].Multiply(-1.0, dgrowthscalardR, dRdC);
  dlambdard_c_[gp].Multiply(-1.0, dlambdardR, dRdC);
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::Evaluate(const Core::LinAlg::Matrix<3, 3>& F,
    const Core::LinAlg::Matrix<6, 1>& E_strain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>& S_stress, Core::LinAlg::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  const double dt = params.get<double>("delta time");

  Core::LinAlg::Matrix<3, 3> C = EvaluateC(F);

  const double lambda_f = evaluate_lambdaf(C, gp, eleGID);
  remodel_fiber_[gp].set_state(lambda_f, 1.0);

  if (params_->growth_enabled_) integrate_local_evolution_equations(dt, gp, eleGID);

  S_stress.Update(evaluate_current_p_k2(gp, eleGID));
  cmat.Update(evaluate_current_cmat(gp, eleGID));
}

void MIXTURE::MixtureConstituentRemodelFiberImpl::EvaluateElasticPart(
    const Core::LinAlg::Matrix<3, 3>& FM, const Core::LinAlg::Matrix<3, 3>& iFextin,
    Teuchos::ParameterList& params, Core::LinAlg::Matrix<6, 1>& S_stress,
    Core::LinAlg::Matrix<6, 6>& cmat, int gp, int eleGID)
{
  FOUR_C_THROW(
      "The implicit remodel fiber cannot be evaluated with an additional inelastic deformation. "
      "You can either use the explicit remodel fiber or use a growth strategy without an inelastic "
      "external deformation.");
}

double MIXTURE::MixtureConstituentRemodelFiberImpl::GetGrowthScalar(int gp) const
{
  return remodel_fiber_[gp].evaluate_current_growth_scalar();
}

Core::LinAlg::Matrix<1, 6> MIXTURE::MixtureConstituentRemodelFiberImpl::GetDGrowthScalarDC(
    int gp, int eleGID) const
{
  if (!params_->growth_enabled_) return Core::LinAlg::Matrix<1, 6>(true);
  return dgrowthscalard_c_[gp];
}

double MIXTURE::MixtureConstituentRemodelFiberImpl::evaluate_deposition_stretch(
    const double time) const
{
  if (params_->deposition_stretch_timefunc_num_ == 0)
  {
    return params_->deposition_stretch_;
  }

  return Global::Problem::Instance()
      ->FunctionById<Core::UTILS::FunctionOfTime>(params_->deposition_stretch_timefunc_num_ - 1)
      .Evaluate(time);
}
void MIXTURE::MixtureConstituentRemodelFiberImpl::update_homeostatic_values(
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

  const double new_lambda_pre = evaluate_deposition_stretch(time);

  for (auto& fiber : remodel_fiber_)
  {
    fiber.update_deposition_stretch(new_lambda_pre);
  }
}

double MIXTURE::MixtureConstituentRemodelFiberImpl::evaluate_lambdaf(
    const Core::LinAlg::Matrix<3, 3>& C, const int gp, const int eleGID) const
{
  return std::sqrt(C.Dot(anisotropy_extension_.get_structural_tensor(gp, 0)));
}
FOUR_C_NAMESPACE_CLOSE
