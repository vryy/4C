// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_viscoelasthyper.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_tensor_products.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_linalg_symmetric_tensor.hpp"
#include "4C_linalg_tensor_conversion.hpp"
#include "4C_linalg_utils_densematrix_inverse.hpp"
#include "4C_mat_elast_summand.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_mat_viscoelast_fsls.hpp"
#include "4C_mat_viscoelast_generalizedmaxwell.hpp"
#include "4C_utils_enum.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ViscoElastHyper::ViscoElastHyper(const Core::Mat::PAR::Parameter::Data& matdata)
    : Mat::PAR::ElastHyper(matdata)
{
  // polyconvexity check is just implemented for isotropic hyperlastic materials
  if (polyconvex_)
    FOUR_C_THROW(
        "This polyconvexity-check is just implemented for isotropic "
        "hyperelastic-materials (do not use for viscoelastic materials).");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::Mat::Material> Mat::PAR::ViscoElastHyper::create_material()
{
  return std::make_shared<Mat::ViscoElastHyper>(this);
}


Mat::ViscoElastHyperType Mat::ViscoElastHyperType::instance_;


Core::Communication::ParObject* Mat::ViscoElastHyperType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::ViscoElastHyper* elhy = new Mat::ViscoElastHyper();
  elhy->unpack(buffer);

  return elhy;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ViscoElastHyper::ViscoElastHyper() : Mat::ElastHyper()
{
  isovisco_ = false;
  visco_generalized_maxwell_ = false;
  visco_fsls_ = false;

  state_.clear();
  rebuild_active_model_sequence();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ViscoElastHyper::ViscoElastHyper(Mat::PAR::ViscoElastHyper* params)
    : Mat::ElastHyper(params),
      isovisco_(false),
      visco_generalized_maxwell_(false),
      visco_fsls_(false)
{
  rebuild_active_model_sequence();
}


Mat::ViscoElastHyper::EvaluateWorkspace::EvaluateWorkspace()
    : glstrain_mat(Core::LinAlg::Initialization::zero),
      c_stress(Core::LinAlg::Initialization::zero),
      i_c_stress(Core::LinAlg::Initialization::zero),
      c_strain(Core::LinAlg::Initialization::zero),
      mod_c_strain(Core::LinAlg::Initialization::zero),
      id2(Core::LinAlg::Initialization::zero),
      id4(Core::LinAlg::Initialization::zero),
      id4sharp(Core::LinAlg::Initialization::zero),
      prinv(Core::LinAlg::Initialization::zero),
      modinv(Core::LinAlg::Initialization::zero),
      rateinv(Core::LinAlg::Initialization::zero),
      modrateinv(Core::LinAlg::Initialization::zero),
      dPI(Core::LinAlg::Initialization::zero),
      ddPII(Core::LinAlg::Initialization::zero),
      scgrate(Core::LinAlg::Initialization::zero),
      modrcgrate(Core::LinAlg::Initialization::zero),
      mu(Core::LinAlg::Initialization::zero),
      modmu(Core::LinAlg::Initialization::zero),
      xi(Core::LinAlg::Initialization::zero),
      modxi(Core::LinAlg::Initialization::zero)
{
}


bool Mat::ViscoElastHyper::is_model_flag_enabled(const ViscoModelKind model_kind) const
{
  switch (model_kind)
  {
    case ViscoModelKind::iso_rate:
      return isovisco_;
    case ViscoModelKind::generalized_maxwell:
      return visco_generalized_maxwell_;
    case ViscoModelKind::fsls:
      return visco_fsls_;
  }

  FOUR_C_THROW("Unsupported visco model kind.");
}


bool Mat::ViscoElastHyper::is_model_active(const ViscoModelKind model_kind) const
{
  for (const ViscoModelKind active_model_kind : active_model_sequence_)
  {
    if (active_model_kind == model_kind) return true;
  }

  return false;
}


void Mat::ViscoElastHyper::rebuild_active_model_sequence()
{
  active_model_sequence_.clear();
  for (const ViscoModelKind model_kind : visco_model_registry())
  {
    if (is_model_flag_enabled(model_kind)) active_model_sequence_.push_back(model_kind);
  }
}


Mat::ViscoElastState::ActiveModels Mat::ViscoElastHyper::active_models_from_flags() const
{
  return Mat::ViscoElastState::ActiveModels{.iso_rate = isovisco_,
      .generalized_maxwell = visco_generalized_maxwell_,
      .fsls = visco_fsls_};
}


Mat::ViscoElastState::ActiveModels Mat::ViscoElastHyper::active_models_from_sequence() const
{
  Mat::ViscoElastState::ActiveModels models;
  for (const ViscoModelKind model_kind : active_model_sequence_)
  {
    switch (model_kind)
    {
      case ViscoModelKind::iso_rate:
        models.iso_rate = true;
        break;
      case ViscoModelKind::generalized_maxwell:
        models.generalized_maxwell = true;
        break;
      case ViscoModelKind::fsls:
        models.fsls = true;
        break;
    }
  }

  return models;
}


void Mat::ViscoElastHyper::ensure_model_activation_consistency(const char* context) const
{
  const Mat::ViscoElastState::ActiveModels models_from_flags = active_models_from_flags();
  const Mat::ViscoElastState::ActiveModels models_from_sequence = active_models_from_sequence();

  if (models_from_flags.iso_rate != models_from_sequence.iso_rate ||
      models_from_flags.generalized_maxwell != models_from_sequence.generalized_maxwell ||
      models_from_flags.fsls != models_from_sequence.fsls)
    FOUR_C_THROW(
        "Inconsistent visco model activation while {} in MAT_ViscoElastHyper (MAT {}): "
        "flags=(iso_rate={}, generalized_maxwell={}, fsls={}), "
        "sequence=(iso_rate={}, generalized_maxwell={}, fsls={}).",
        context, params_ != nullptr ? params_->id() : -1, models_from_flags.iso_rate,
        models_from_flags.generalized_maxwell, models_from_flags.fsls,
        models_from_sequence.iso_rate, models_from_sequence.generalized_maxwell,
        models_from_sequence.fsls);
}


Mat::ViscoElastState::ActiveModels Mat::ViscoElastHyper::active_models() const
{
  ensure_model_activation_consistency("deriving active model descriptor");
  return active_models_from_sequence();
}


Mat::ViscoElastHyper::GeneralizedMaxwellSolveKind
Mat::ViscoElastHyper::parse_generalized_maxwell_solve_kind(
    const std::string& solve, const int gp, const int eleGID) const
{
  if (solve == "OneStepTheta") return GeneralizedMaxwellSolveKind::one_step_theta;
  if (solve == "ExponentialTimeDiscretization")
    return GeneralizedMaxwellSolveKind::exponential_time_discretization;

  FOUR_C_THROW(
      "Invalid SOLVE='{}' in VISCO_GeneralizedMaxwell for MAT_ViscoElastHyper (MAT {}, GP {}, "
      "ELE {}). Use OneStepTheta or ExponentialTimeDiscretization.",
      solve, params_ != nullptr ? params_->id() : -1, gp, eleGID);

  return GeneralizedMaxwellSolveKind::exponential_time_discretization;
}


void Mat::ViscoElastHyper::clear_generalized_maxwell_metadata()
{
  generalized_maxwell_metadata_.reset();
}


void Mat::ViscoElastHyper::clear_runtime_context() { runtime_context_.reset(); }


void Mat::ViscoElastHyper::clear_fsls_metadata() { fsls_metadata_.reset(); }


void Mat::ViscoElastHyper::build_generalized_maxwell_metadata_for_setup(
    const int gp, const int eleGID)
{
  clear_generalized_maxwell_metadata();

  const int visco_mat_id = params_ != nullptr ? params_->id() : -1;
  int generalized_maxwell_model_count = 0;
  int generalized_maxwell_numbranch_value = -1;
  std::string generalized_maxwell_solve;
  const std::vector<int>* generalized_maxwell_matids = nullptr;

  std::shared_ptr<Mat::Elastic::GeneralizedMaxwell> generalized_maxwell = nullptr;
  int generalized_maxwell_summand_mat_id = -1;

  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    auto current_generalized_maxwell =
        std::dynamic_pointer_cast<Mat::Elastic::GeneralizedMaxwell>(potsum_[p]);
    if (current_generalized_maxwell != nullptr)
    {
      ++generalized_maxwell_model_count;
      generalized_maxwell = current_generalized_maxwell;
      generalized_maxwell_summand_mat_id = params_ != nullptr ? mat_id(p) : -1;
      generalized_maxwell->read_material_parameters(generalized_maxwell_numbranch_value,
          generalized_maxwell_matids, generalized_maxwell_solve);
    }
  }

  if (generalized_maxwell_model_count != 1)
    FOUR_C_THROW(
        "Invalid VISCO_GeneralizedMaxwell setup in MAT_ViscoElastHyper (MAT {}): expected "
        "exactly one VISCO_GeneralizedMaxwell summand but found {}.",
        visco_mat_id, generalized_maxwell_model_count);

  if (generalized_maxwell == nullptr)
    FOUR_C_THROW(
        "Failed to resolve VISCO_GeneralizedMaxwell summand in MAT_ViscoElastHyper (MAT {}).",
        visco_mat_id);

  if (generalized_maxwell_numbranch_value <= 0)
    FOUR_C_THROW(
        "Invalid VISCO_GeneralizedMaxwell setup in MAT_ViscoElastHyper (MAT {}): "
        "NUMBRANCH={} is not positive.",
        visco_mat_id, generalized_maxwell_numbranch_value);

  if (generalized_maxwell_matids == nullptr)
    FOUR_C_THROW(
        "Failed to read MATIDS for VISCO_GeneralizedMaxwell in MAT_ViscoElastHyper (MAT {}).",
        visco_mat_id);

  if (generalized_maxwell_matids->size() !=
      static_cast<unsigned int>(generalized_maxwell_numbranch_value))
    FOUR_C_THROW(
        "Invalid VISCO_GeneralizedMaxwell setup in MAT_ViscoElastHyper (MAT {}): "
        "NUMBRANCH={} but MATIDS has size {}.",
        visco_mat_id, generalized_maxwell_numbranch_value, generalized_maxwell_matids->size());

  const auto& branchespotsum = generalized_maxwell->get_branchespotsum();
  const auto& branchtau = generalized_maxwell->get_branchtaus();

  if (branchespotsum.size() != static_cast<unsigned int>(generalized_maxwell_numbranch_value) ||
      branchtau.size() != static_cast<unsigned int>(generalized_maxwell_numbranch_value))
    FOUR_C_THROW(
        "Failed to initialize VISCO_GeneralizedMaxwell branches in MAT_ViscoElastHyper (MAT {}, "
        "GP {}, ELE {}). Expected {} branches, got {} branch definitions and {} branch "
        "relaxation times.",
        visco_mat_id, gp, eleGID, generalized_maxwell_numbranch_value, branchespotsum.size(),
        branchtau.size());

  if (is_model_active(ViscoModelKind::iso_rate))
    FOUR_C_THROW(
        "Unsupported branch formulation in VISCO_GeneralizedMaxwell for MAT_ViscoElastHyper "
        "(MAT {}, GP {}, ELE {}): isovisco branch response is not implemented.",
        visco_mat_id, gp, eleGID);

  GeneralizedMaxwellMetadata generalized_maxwell_metadata;
  generalized_maxwell_metadata.summand_mat_id = generalized_maxwell_summand_mat_id;
  generalized_maxwell_metadata.solve_kind =
      parse_generalized_maxwell_solve_kind(generalized_maxwell_solve, gp, eleGID);
  generalized_maxwell_metadata.branches.reserve(generalized_maxwell_numbranch_value);

  for (int i = 0; i < generalized_maxwell_numbranch_value; ++i)
  {
    GeneralizedMaxwellBranchMetadata branch_metadata;
    branch_metadata.summands = branchespotsum.at(i);
    branch_metadata.tau = branchtau.at(i);

    if (branch_metadata.tau <= 0.0)
      FOUR_C_THROW(
          "Invalid branch relaxation time TAU={} in VISCO_GeneralizedMaxwell (MAT {}, branch {}, "
          "GP {}, ELE {}). Expected TAU > 0.",
          branch_metadata.tau, generalized_maxwell_summand_mat_id, i, gp, eleGID);

    branch_metadata.properties.clear();
    elast_hyper_properties(branch_metadata.summands, branch_metadata.properties);

    if (branch_metadata.properties.anisoprinc)
      FOUR_C_THROW(
          "Unsupported branch formulation in VISCO_GeneralizedMaxwell for MAT_ViscoElastHyper "
          "(MAT {}, GP {}, ELE {}): anisoprinc branch response is not implemented.",
          visco_mat_id, gp, eleGID);

    if (branch_metadata.properties.anisomod)
      FOUR_C_THROW(
          "Unsupported branch formulation in VISCO_GeneralizedMaxwell for MAT_ViscoElastHyper "
          "(MAT {}, GP {}, ELE {}): anisomod branch response is not implemented.",
          visco_mat_id, gp, eleGID);

    generalized_maxwell_metadata.branches.push_back(std::move(branch_metadata));
  }

  generalized_maxwell_metadata_ = std::move(generalized_maxwell_metadata);
}


const Mat::ViscoElastHyper::GeneralizedMaxwellMetadata&
Mat::ViscoElastHyper::require_generalized_maxwell_metadata(
    const char* context, const int gp, const int eleGID) const
{
  if (!generalized_maxwell_metadata_.has_value())
    FOUR_C_THROW(
        "Missing VISCO_GeneralizedMaxwell metadata cache while {} in MAT_ViscoElastHyper (MAT "
        "{}, GP {}, ELE {}). Run setup() before evaluation.",
        context, params_ != nullptr ? params_->id() : -1, gp, eleGID);

  if (generalized_maxwell_metadata_->branches.empty())
    FOUR_C_THROW(
        "Invalid VISCO_GeneralizedMaxwell metadata cache while {} in MAT_ViscoElastHyper (MAT "
        "{}, GP {}, ELE {}): no branch metadata available.",
        context, params_ != nullptr ? params_->id() : -1, gp, eleGID);

  return generalized_maxwell_metadata_.value();
}


std::size_t Mat::ViscoElastHyper::read_generalized_maxwell_branch_count_for_setup() const
{
  if (!is_model_active(ViscoModelKind::generalized_maxwell)) return 0;

  return require_generalized_maxwell_metadata("reading branch count", -1, -1).branches.size();
}


void Mat::ViscoElastHyper::build_runtime_context_for_setup(const int gp, const int eleGID)
{
  clear_runtime_context();

  ViscoRuntimeContext runtime_context;
  const bool requires_structural_dynamic_parameters =
      is_model_active(ViscoModelKind::generalized_maxwell) || is_model_active(ViscoModelKind::fsls);
  if (!requires_structural_dynamic_parameters)
  {
    runtime_context_ = std::move(runtime_context);
    return;
  }

  // Integration-boundary access: read structural-dynamic runtime settings once and cache them
  // for visco update/evaluation code paths.
  const Teuchos::ParameterList& structural_dynamic_parameters =
      Global::Problem::instance()->structural_dynamic_params();

  if (is_model_active(ViscoModelKind::generalized_maxwell))
  {
    if (!structural_dynamic_parameters.isParameter("DYNAMICTYPE"))
      FOUR_C_THROW(
          "Missing DYNAMICTYPE in STRUCTURAL DYNAMIC parameters while building generalized "
          "Maxwell runtime context for MAT_ViscoElastHyper (MAT {}, GP {}, ELE {}).",
          params_ != nullptr ? params_->id() : -1, gp, eleGID);

    GeneralizedMaxwellRuntimeContext generalized_maxwell_runtime_context;
    const auto dyntype = Teuchos::getIntegralValue<Inpar::Solid::DynamicType>(
        structural_dynamic_parameters, "DYNAMICTYPE");
    if (dyntype == Inpar::Solid::DynamicType::OneStepTheta)
      generalized_maxwell_runtime_context.one_step_theta =
          structural_dynamic_parameters.sublist("ONESTEPTHETA").get<double>("THETA");

    runtime_context.generalized_maxwell = generalized_maxwell_runtime_context;
  }

  if (is_model_active(ViscoModelKind::fsls))
  {
    if (!structural_dynamic_parameters.isParameter("NUMSTEP"))
      FOUR_C_THROW(
          "Missing NUMSTEP in STRUCTURAL DYNAMIC parameters while deriving FSLS history "
          "capacity for MAT_ViscoElastHyper (MAT {}, GP {}, ELE {}).",
          params_ != nullptr ? params_->id() : -1, gp, eleGID);

    const int numsteps = structural_dynamic_parameters.get<int>("NUMSTEP");
    if (numsteps < 0)
      FOUR_C_THROW(
          "Invalid NUMSTEP={} while deriving FSLS history capacity for MAT_ViscoElastHyper (MAT "
          "{}, GP {}, ELE {}). Expected NUMSTEP >= 0.",
          numsteps, params_ != nullptr ? params_->id() : -1, gp, eleGID);

    FslsRuntimeContext fsls_runtime_context;
    fsls_runtime_context.max_history_size = static_cast<unsigned int>(numsteps + 1);
    runtime_context.fsls = fsls_runtime_context;
  }

  runtime_context_ = std::move(runtime_context);
}


const Mat::ViscoElastHyper::ViscoRuntimeContext& Mat::ViscoElastHyper::require_runtime_context(
    const char* context, const int gp, const int eleGID) const
{
  if (!runtime_context_.has_value())
    FOUR_C_THROW(
        "Missing visco runtime context while {} in MAT_ViscoElastHyper (MAT {}, GP {}, ELE {}). "
        "Run setup() before evaluation or update.",
        context, params_ != nullptr ? params_->id() : -1, gp, eleGID);

  return runtime_context_.value();
}


void Mat::ViscoElastHyper::build_fsls_metadata_for_setup(const int gp, const int eleGID)
{
  clear_fsls_metadata();

  FslsMetadata fsls_metadata;
  const int visco_mat_id = params_ != nullptr ? params_->id() : -1;
  int fsls_model_count = 0;

  std::string solve;
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    std::shared_ptr<Mat::Elastic::Fsls> fsls =
        std::dynamic_pointer_cast<Mat::Elastic::Fsls>(potsum_[p]);
    if (fsls != nullptr)
    {
      ++fsls_model_count;
      fsls_metadata.summand_mat_id = params_ != nullptr ? mat_id(p) : -1;
      fsls->read_material_parameters_visco(
          fsls_metadata.tau, fsls_metadata.beta, fsls_metadata.alpha, solve);
    }
  }

  if (fsls_model_count != 1)
    FOUR_C_THROW(
        "Invalid VISCO_FSLS setup in MAT_ViscoElastHyper (MAT {}, GP {}, ELE {}): expected "
        "exactly one VISCO_FSLS summand but found {}.",
        visco_mat_id, gp, eleGID, fsls_model_count);

  if (fsls_metadata.tau <= 0.0)
    FOUR_C_THROW(
        "Invalid TAU={} in VISCO_FSLS (MAT {}, referenced by MAT_ViscoElastHyper MAT {}, GP {}, "
        "ELE {}). TAU has to be positive.",
        fsls_metadata.tau, fsls_metadata.summand_mat_id, visco_mat_id, gp, eleGID);

  if (fsls_metadata.alpha < 0.0 || fsls_metadata.alpha >= 1.0)
    FOUR_C_THROW(
        "Invalid ALPHA={} in VISCO_FSLS (MAT {}, referenced by MAT_ViscoElastHyper MAT {}, GP "
        "{}, ELE {}). Expected 0 <= ALPHA < 1.",
        fsls_metadata.alpha, fsls_metadata.summand_mat_id, visco_mat_id, gp, eleGID);

  fsls_metadata_ = std::move(fsls_metadata);
}


const Mat::ViscoElastHyper::FslsMetadata& Mat::ViscoElastHyper::require_fsls_metadata(
    const char* context, const int gp, const int eleGID) const
{
  if (!fsls_metadata_.has_value())
    FOUR_C_THROW(
        "Missing VISCO_FSLS metadata cache while {} in MAT_ViscoElastHyper (MAT {}, GP {}, ELE "
        "{}). Run setup() before evaluation.",
        context, params_ != nullptr ? params_->id() : -1, gp, eleGID);

  return fsls_metadata_.value();
}


Mat::ViscoElastHyper::FslsParameters Mat::ViscoElastHyper::read_fsls_parameters(
    const int gp, const int eleGID) const
{
  const FslsMetadata& fsls_metadata = require_fsls_metadata("reading FSLS parameters", gp, eleGID);

  FslsParameters fsls_parameters;
  fsls_parameters.tau = fsls_metadata.tau;
  fsls_parameters.alpha = fsls_metadata.alpha;
  fsls_parameters.beta = fsls_metadata.beta;
  fsls_parameters.summand_mat_id = fsls_metadata.summand_mat_id;

  return fsls_parameters;
}


unsigned int Mat::ViscoElastHyper::read_fsls_max_history_size_for_update() const
{
  if (!is_model_active(ViscoModelKind::fsls)) return 0;

  const ViscoRuntimeContext& runtime_context =
      require_runtime_context("reading FSLS history capacity", -1, -1);
  if (!runtime_context.fsls.has_value())
    FOUR_C_THROW(
        "Missing FSLS runtime context while reading FSLS history capacity in "
        "MAT_ViscoElastHyper (MAT {}).",
        params_ != nullptr ? params_->id() : -1);

  const unsigned int max_history_size = runtime_context.fsls->max_history_size;
  if (max_history_size == 0)
    FOUR_C_THROW(
        "Invalid FSLS runtime history capacity {} in MAT_ViscoElastHyper (MAT {}). Expected a "
        "positive history capacity.",
        max_history_size, params_ != nullptr ? params_->id() : -1);

  return max_history_size;
}


double Mat::ViscoElastHyper::read_visco_time_step_size(
    const EvaluationContext<3>& context, const int gp, const int eleGID) const
{
  if (context.time_step_size == nullptr)
    FOUR_C_THROW(
        "Missing EvaluationContext::time_step_size in MAT_ViscoElastHyper (MAT {}, GP {}, ELE "
        "{}).",
        params_ != nullptr ? params_->id() : -1, gp, eleGID);

  return *context.time_step_size;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);
  summandProperties_.pack(data);
  const Mat::ViscoElastState::ActiveModels models = active_models();
  add_to_pack(data, models.iso_rate);
  add_to_pack(data, models.generalized_maxwell);
  add_to_pack(data, models.fsls);

  anisotropy_.pack_anisotropy(data);

  Core::Communication::PotentiallyUnusedBufferScope summand_scope{data};
  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->pack_summand(data);
    }

    state_.serialize_state(data, active_models());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::unpack(Core::Communication::UnpackBuffer& buffer)
{
  // make sure we have a pristine material
  params_ = nullptr;
  potsum_.clear();

  summandProperties_.clear();
  isovisco_ = false;
  visco_generalized_maxwell_ = false;
  visco_fsls_ = false;
  state_.clear();
  active_model_sequence_.clear();
  clear_generalized_maxwell_metadata();
  clear_runtime_context();
  clear_fsls_metadata();

  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // Integration-boundary access: unpack resolves parameter material from the global
  // material bundle and then hands explicit data to visco internals.
  // matid and recover params_
  int matid;
  extract_from_pack(buffer, matid);
  if (Global::Problem::instance()->materials() != nullptr)
  {
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const unsigned int probinst =
          Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::ViscoElastHyper*>(mat);
      else
        FOUR_C_THROW("Type of parameter material {} does not fit to calling type {}", mat->type(),
            material_type());
    }
  }

  summandProperties_.unpack(buffer);
  extract_from_pack(buffer, isovisco_);
  extract_from_pack(buffer, visco_generalized_maxwell_);
  extract_from_pack(buffer, visco_fsls_);
  rebuild_active_model_sequence();
  ensure_model_activation_consistency("unpack model activation");

  anisotropy_.unpack_anisotropy(buffer);

  Core::Communication::PotentiallyUnusedBufferScope summand_scope{buffer};
  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = params_->matids_.begin(); m != params_->matids_.end(); ++m)
    {
      const int matid = *m;
      std::shared_ptr<Mat::Elastic::Summand> sum = Mat::Elastic::Summand::factory(matid);
      if (sum == nullptr) FOUR_C_THROW("Failed to allocate");
      potsum_.push_back(sum);
    }

    // loop map of associated potential summands
    for (auto& p : potsum_)
    {
      p->unpack_summand(buffer);
      p->register_anisotropy_extensions(anisotropy_);
    }

    if (is_model_active(ViscoModelKind::generalized_maxwell))
    {
      build_generalized_maxwell_metadata_for_setup(-1, -1);
    }

    build_runtime_context_for_setup(-1, -1);

    if (is_model_active(ViscoModelKind::fsls))
    {
      build_fsls_metadata_for_setup(-1, -1);
    }

    state_.deserialize_state(buffer, active_models());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::setup(int numgp, const Discret::Elements::Fibers& fibers,
    const std::optional<Discret::Elements::CoordinateSystem>& coord_system)
{
  // read anisotropy
  anisotropy_.set_number_of_gauss_points(numgp);
  anisotropy_.read_anisotropy_from_element(fibers, coord_system);

  // Setup summands
  for (auto& p : potsum_) p->setup(numgp, fibers, coord_system);

  // find out which formulations are used
  isovisco_ = false;
  visco_generalized_maxwell_ = false;
  visco_fsls_ = false;

  summandProperties_.clear();
  elast_hyper_properties(potsum_, summandProperties_);


  if (summandProperties_.viscoGeneral)
  {
    for (auto& p : potsum_)
    {
      p->specify_visco_formulation(isovisco_, visco_generalized_maxwell_, visco_fsls_);
    }
  }

  rebuild_active_model_sequence();
  ensure_model_activation_consistency("pre-loop setup orchestration");

  clear_generalized_maxwell_metadata();
  clear_runtime_context();
  clear_fsls_metadata();
  if (is_model_active(ViscoModelKind::generalized_maxwell))
  {
    build_generalized_maxwell_metadata_for_setup(-1, -1);
  }
  build_runtime_context_for_setup(-1, -1);
  if (is_model_active(ViscoModelKind::fsls))
  {
    build_fsls_metadata_for_setup(-1, -1);
  }

  const Mat::ViscoElastState::ActiveModels models = active_models();
  const std::size_t generalized_maxwell_numbranch =
      read_generalized_maxwell_branch_count_for_setup();

  state_.initialize_from_setup(numgp, models, generalized_maxwell_numbranch);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::update()
{
  Mat::ElastHyper::update();

  ensure_model_activation_consistency("pre-loop update orchestration");

  const unsigned int max_hist = read_fsls_max_history_size_for_update();

  const int visco_mat_id = params_ != nullptr ? params_->id() : -1;
  state_.advance_time_step(active_models(), max_hist, visco_mat_id);

  return;
}

void Mat::ViscoElastHyper::prepare_evaluate_kinematics(
    const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
    const EvaluationContext<3>& context, EvaluateWorkspace& workspace, const int gp,
    const int eleGID) const
{
  workspace.glstrain_mat = Core::LinAlg::make_strain_like_voigt_matrix(glstrain);

  evaluate_right_cauchy_green_strain_like_voigt(glstrain, workspace.c);
  workspace.i_c = Core::LinAlg::inv(workspace.c);

  workspace.c_stress = Core::LinAlg::make_stress_like_voigt_view(workspace.c);
  workspace.i_c_stress = Core::LinAlg::make_stress_like_voigt_view(workspace.i_c);

  Core::LinAlg::Voigt::Stresses::to_strain_like(workspace.c_stress, workspace.c_strain);
  Core::LinAlg::Voigt::Stresses::invariants_principal(workspace.prinv, workspace.c_stress);

  workspace.dt = 0.0;
  if (!active_model_sequence_.empty())
    workspace.dt = read_visco_time_step_size(context, gp, eleGID);

  Core::LinAlg::Voigt::identity_matrix(workspace.id2);

  using VoigtNotation = Core::LinAlg::Voigt::NotationType;
  Core::LinAlg::Voigt::fourth_order_identity_matrix<VoigtNotation::stress, VoigtNotation::stress>(
      workspace.id4sharp);
  Core::LinAlg::Voigt::fourth_order_identity_matrix<VoigtNotation::stress, VoigtNotation::strain>(
      workspace.id4);

  elast_hyper_evaluate_invariant_derivatives(
      workspace.prinv, workspace.dPI, workspace.ddPII, potsum_, summandProperties_, gp, eleGID);
}


void Mat::ViscoElastHyper::prepare_iso_rate_visco_inputs_if_active(
    const Teuchos::ParameterList& params, EvaluateWorkspace& workspace, const int gp,
    const int eleGID)
{
  for (const ViscoModelKind model_kind : active_model_sequence_)
  {
    if (model_kind != ViscoModelKind::iso_rate) continue;

    if (summandProperties_.isomod)
    {
      // calculate modified invariants
      invariants_modified(workspace.modinv, workspace.prinv);
    }

    // calculate viscous quantities
    evaluate_kin_quant_vis(workspace.c_strain, workspace.c_stress, workspace.i_c_stress,
        workspace.prinv, workspace.rateinv, workspace.mod_c_strain, workspace.dt, workspace.scgrate,
        workspace.modrcgrate, workspace.modrateinv, gp);
    evaluate_mu_xi(workspace.prinv, workspace.modinv, workspace.mu, workspace.modmu, workspace.xi,
        workspace.modxi, workspace.rateinv, workspace.modrateinv, params, workspace.dt, gp, eleGID);
    break;
  }
}


void Mat::ViscoElastHyper::initialize_elastic_response(const EvaluateWorkspace& workspace,
    Core::LinAlg::Matrix<6, 1>& stress_view, Core::LinAlg::Matrix<6, 6>& cmat_view,
    Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat) const
{
  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress_view.clear();
  cmat_view.clear();

  // add isotropic part
  elast_hyper_add_isotropic_stress_cmat(
      stress, cmat, workspace.c, workspace.i_c, workspace.prinv, workspace.dPI, workspace.ddPII);
}


void Mat::ViscoElastHyper::add_iso_rate_contribution(const EvaluateWorkspace& workspace,
    Core::LinAlg::Matrix<6, 1>& stress_view, Core::LinAlg::Matrix<6, 6>& cmat_view)
{
  if (summandProperties_.isomod)
  {
    // add viscous part decoupled
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> stressisomodisovisco(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisomodisovisco(
        Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> stressisomodvolvisco(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisomodvolvisco(
        Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 1> prinv(workspace.prinv);
    Core::LinAlg::Matrix<3, 1> modinv(workspace.modinv);
    Core::LinAlg::Matrix<8, 1> modmu(workspace.modmu);
    Core::LinAlg::Matrix<33, 1> modxi(workspace.modxi);
    Core::LinAlg::Matrix<6, 1> c_strain(workspace.c_strain);
    Core::LinAlg::Matrix<6, 1> id2(workspace.id2);
    Core::LinAlg::Matrix<6, 1> i_c_stress(workspace.i_c_stress);
    Core::LinAlg::Matrix<6, 6> id4(workspace.id4);
    Core::LinAlg::Matrix<6, 1> modrcgrate(workspace.modrcgrate);
    evaluate_iso_visco_modified(stressisomodisovisco, stressisomodvolvisco, cmatisomodisovisco,
        cmatisomodvolvisco, prinv, modinv, modmu, modxi, c_strain, id2, i_c_stress, id4,
        modrcgrate);
    stress_view.update(1.0, stressisomodisovisco, 1.0);
    stress_view.update(1.0, stressisomodvolvisco, 1.0);
    cmat_view.update(1.0, cmatisomodisovisco, 1.0);
    cmat_view.update(1.0, cmatisomodvolvisco, 1.0);
  }

  if (summandProperties_.isoprinc)
  {
    // add viscous part coupled
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> stressisovisco(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisovisco(
        Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<8, 1> mu(workspace.mu);
    Core::LinAlg::Matrix<33, 1> xi(workspace.xi);
    Core::LinAlg::Matrix<6, 6> id4sharp(workspace.id4sharp);
    Core::LinAlg::Matrix<6, 1> scgrate(workspace.scgrate);
    evaluate_iso_visco_principal(stressisovisco, cmatisovisco, mu, xi, id4sharp, scgrate);
    stress_view.update(1.0, stressisovisco, 1.0);
    cmat_view.update(1.0, cmatisovisco, 1.0);
  }
}


void Mat::ViscoElastHyper::add_generalized_maxwell_contribution(const EvaluateWorkspace& workspace,
    Core::LinAlg::Matrix<6, 1>& stress_view, Core::LinAlg::Matrix<6, 6>& cmat_view, const int gp,
    const int eleGID)
{
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> q(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatq(Core::LinAlg::Initialization::zero);
  evaluate_visco_generalized_maxwell(q, cmatq, workspace.dt, workspace.glstrain_mat, gp, eleGID);
  stress_view.update(1.0, q, 1.0);
  cmat_view.update(1.0, cmatq, 1.0);
}


void Mat::ViscoElastHyper::add_fsls_contribution(const EvaluateWorkspace& workspace,
    Core::LinAlg::Matrix<6, 1>& stress_view, Core::LinAlg::Matrix<6, 6>& cmat_view, const int gp,
    const int eleGID)
{
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> q(
      Core::LinAlg::Initialization::zero);  // artificial viscous stress
  Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatq(Core::LinAlg::Initialization::zero);
  evaluate_visco_fsls(stress_view, cmat_view, q, cmatq, workspace.dt, gp, eleGID);
  stress_view.update(1.0, q, 1.);
  cmat_view.update(1.0, cmatq, 1.);
}


void Mat::ViscoElastHyper::add_visco_contributions_in_sequence(const EvaluateWorkspace& workspace,
    Core::LinAlg::Matrix<6, 1>& stress_view, Core::LinAlg::Matrix<6, 6>& cmat_view, const int gp,
    const int eleGID)
{
  // add viscous part
  for (const ViscoModelKind model_kind : active_model_sequence_)
  {
    switch (model_kind)
    {
      case ViscoModelKind::iso_rate:
      {
        add_iso_rate_contribution(workspace, stress_view, cmat_view);
        break;
      }
      case ViscoModelKind::generalized_maxwell:
      {
        add_generalized_maxwell_contribution(workspace, stress_view, cmat_view, gp, eleGID);
        break;
      }
      case ViscoModelKind::fsls:
      {
        add_fsls_contribution(workspace, stress_view, cmat_view, gp, eleGID);
        break;
      }
    }
  }
}


void Mat::ViscoElastHyper::add_post_elastic_composition_hooks(const Teuchos::ParameterList& params,
    const EvaluationContext<3>& context, const EvaluateWorkspace& workspace,
    Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, const int gp, const int eleGID)
{
  /*----------------------------------------------------------------------*/
  // coefficients in principal stretches
  if (summandProperties_.coeffStretchesPrinc || summandProperties_.coeffStretchesMod)
  {
    elast_hyper_add_response_stretches(
        cmat, stress, workspace.c, potsum_, summandProperties_, gp, eleGID);
  }

  /*----------------------------------------------------------------------*/
  // Do all the anisotropic stuff!
  if (summandProperties_.anisoprinc)
  {
    elast_hyper_add_anisotropic_princ(stress, cmat, workspace.c, params, gp, eleGID, potsum_);
  }

  if (summandProperties_.anisomod)
  {
    elast_hyper_add_anisotropic_mod(
        stress, cmat, workspace.c, workspace.i_c, workspace.prinv, gp, eleGID, context, potsum_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::evaluate(const Core::LinAlg::Tensor<double, 3, 3>* defgrad,
    const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
    const Teuchos::ParameterList& params, const EvaluationContext<3>& context,
    Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp, int eleGID)
{
  ensure_model_activation_consistency("pre-loop evaluate orchestration");
  Core::LinAlg::Matrix<6, 1> stress_view = Core::LinAlg::make_stress_like_voigt_view(stress);
  Core::LinAlg::Matrix<6, 6> cmat_view = Core::LinAlg::make_stress_like_voigt_view(cmat);

  EvaluateWorkspace workspace;

  prepare_evaluate_kinematics(glstrain, context, workspace, gp, eleGID);
  prepare_iso_rate_visco_inputs_if_active(params, workspace, gp, eleGID);
  initialize_elastic_response(workspace, stress_view, cmat_view, stress, cmat);
  add_visco_contributions_in_sequence(workspace, stress_view, cmat_view, gp, eleGID);
  add_post_elastic_composition_hooks(params, context, workspace, stress, cmat, gp, eleGID);
}

/*----------------------------------------------------------------------*/
/* Evaluate Quantities for viscous Part                           09/13 */
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::evaluate_kin_quant_vis(Core::LinAlg::Matrix<6, 1>& rcg,
    Core::LinAlg::Matrix<6, 1>& scg, Core::LinAlg::Matrix<6, 1>& icg,
    Core::LinAlg::Matrix<3, 1>& prinv, Core::LinAlg::Matrix<7, 1>& rateinv,
    Core::LinAlg::Matrix<6, 1>& modrcg, const double dt, Core::LinAlg::Matrix<6, 1>& scgrate,
    Core::LinAlg::Matrix<6, 1>& modrcgrate, Core::LinAlg::Matrix<7, 1>& modrateinv, const int gp)
{
  if (dt <= 0.0)
    FOUR_C_THROW(
        "Invalid time step size dt={} in MAT_ViscoElastHyper (MAT {}, GP {}) for "
        "rate-dependent viscous update. Expected dt > 0.",
        dt, params_ != nullptr ? params_->id() : -1, gp);

  // modrcg : \overline{C} = J^{-\frac{2}{3}} C
  const double modscale = std::pow(prinv(2), -1. / 3.);
  modrcg.update(modscale, rcg);

  // read history
  const auto iso_rate_previous = state_.iso_rate_prev_point(gp);
  Core::LinAlg::Matrix<6, 1> scglast(iso_rate_previous.scg);
  Core::LinAlg::Matrix<6, 1> modrcglast(iso_rate_previous.modrcg);

  // Update history of Cauchy-Green Tensor
  state_.set_iso_rate_current_point(gp, scg, modrcg);  // store C^{n} and \overline{C}^{n}

  // rate of Cauchy-Green Tensor
  // REMARK: strain-like 6-Voigt vector
  scgrate.update(1.0, scg, 1.0);  // principal material: \dot{C} = \frac{C^n - C^{n-1}}{\Delta t}
  scgrate.update(-1.0, scglast, 1.0);
  scgrate.scale(1 / dt);

  modrcgrate.update(1.0, modrcg, 1.0);  // decoupled material: \overline{\dot{C}} =
                                        // \frac{\overline{C}^n - \overline{C}^{n-1}}{\Delta t}
  modrcgrate.update(-1.0, modrcglast, 1.0);
  modrcgrate.scale(1 / dt);

  // invariants
  // -------------------------------------------------------------------
  // Second Invariant of modrcgrate \bar{J}_2 = \frac{1}{2} \tr (\dot{\overline{C^2}}
  modrateinv(1) =
      0.5 * (modrcgrate(0) * modrcgrate(0) + modrcgrate(1) * modrcgrate(1) +
                modrcgrate(2) * modrcgrate(2) + .5 * modrcgrate(3) * modrcgrate(3) +
                .5 * modrcgrate(4) * modrcgrate(4) + .5 * modrcgrate(5) * modrcgrate(5));


  // For further extension of material law (not necessary at the moment)
  /*
  // necessary transfer variable: Core::LinAlg::Matrix<6,1>& modicgrate
  // \overline{J}_3 = determinant of modified rate of right Cauchy-Green-Tensor
  modrateinv(2) = modrcgrate(0)*modrcgrate(1)*modrcgrate(2)
      + 0.25 * modrcgrate(3)*modrcgrate(4)*modrcgrate(5)
      - 0.25 * modrcgrate(1)*modrcgrate(5)*modrcgrate(5)
      - 0.25 * modrcgrate(2)*modrcgrate(3)*modrcgrate(3)
      - 0.25 * modrcgrate(0)*modrcgrate(4)*modrcgrate(4);

  // invert modified rate of right Cauchy-Green tensor
  // REMARK: stress-like 6-Voigt vector
  {
    modicgrate(0) = ( modrcgrate(1)*modrcgrate(2) - 0.25*modrcgrate(4)*modrcgrate(4) ) /
  modrateinv(2); modicgrate(1) = ( modrcgrate(0)*modrcgrate(2) - 0.25*modrcgrate(5)*modrcgrate(5) )
  / modrateinv(2); modicgrate(2) = ( modrcgrate(0)*modrcgrate(1) - 0.25*modrcgrate(3)*modrcgrate(3)
  ) / modrateinv(2); modicgrate(3) = ( 0.25*modrcgrate(5)*modrcgrate(4) -
  0.5*modrcgrate(3)*modrcgrate(2) ) / modrateinv(2); modicgrate(4) = (
  0.25*modrcgrate(3)*modrcgrate(5) - 0.5*modrcgrate(0)*modrcgrate(4) ) / modrateinv(2);
    modicgrate(5) = ( 0.25*modrcgrate(3)*modrcgrate(4) - 0.5*modrcgrate(5)*modrcgrate(1) ) /
  modrateinv(2);
  }
   */
}

/*----------------------------------------------------------------------*/
/* Evaluate Factors for viscous Quantities                        09/13 */
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::evaluate_mu_xi(Core::LinAlg::Matrix<3, 1>& prinv,
    Core::LinAlg::Matrix<3, 1>& modinv, Core::LinAlg::Matrix<8, 1>& mu,
    Core::LinAlg::Matrix<8, 1>& modmu, Core::LinAlg::Matrix<33, 1>& xi,
    Core::LinAlg::Matrix<33, 1>& modxi, Core::LinAlg::Matrix<7, 1>& rateinv,
    Core::LinAlg::Matrix<7, 1>& modrateinv, const Teuchos::ParameterList& params, const double dt,
    const int gp, const int eleGID)
{
  // principal materials
  if (summandProperties_.isoprinc)
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->add_coefficients_visco_principal(prinv, mu, xi, rateinv, params, dt, gp, eleGID);
    }
  }

  // decoupled (volumetric or isochoric) materials
  if (summandProperties_.isomod)
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->add_coefficients_visco_modified(
          modinv, modmu, modxi, modrateinv, params, dt, gp, eleGID);
    }
  }
}

/*----------------------------------------------------------------------*/
/* stress and constitutive tensor for principal viscous materials       */
/*                                                        pfaller May15 */
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::evaluate_iso_visco_principal(Core::LinAlg::Matrix<6, 1>& stress,
    Core::LinAlg::Matrix<6, 6>& cmat, Core::LinAlg::Matrix<8, 1>& mu,
    Core::LinAlg::Matrix<33, 1>& xi, Core::LinAlg::Matrix<6, 6>& id4sharp,
    Core::LinAlg::Matrix<6, 1>& scgrate)
{
  // contribution: \dot{C}
  stress.update(mu(2), scgrate, 1.0);

  // contribution: id4sharp_{ijkl} = 1/2 (\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})
  cmat.update(xi(2), id4sharp, 1.0);

  return;
}

/*----------------------------------------------------------------------*/
/* stress and constitutive tensor for decoupled viscous materials 09/13 */
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::evaluate_iso_visco_modified(
    Core::LinAlg::Matrix<6, 1>& stressisomodisovisco,
    Core::LinAlg::Matrix<6, 1>& stressisomodvolvisco,
    Core::LinAlg::Matrix<6, 6>& cmatisomodisovisco, Core::LinAlg::Matrix<6, 6>& cmatisomodvolvisco,
    Core::LinAlg::Matrix<3, 1>& prinv, Core::LinAlg::Matrix<3, 1>& modinv,
    Core::LinAlg::Matrix<8, 1>& modmu, Core::LinAlg::Matrix<33, 1>& modxi,
    Core::LinAlg::Matrix<6, 1>& rcg, Core::LinAlg::Matrix<6, 1>& id2,
    Core::LinAlg::Matrix<6, 1>& icg, Core::LinAlg::Matrix<6, 6>& id4,
    Core::LinAlg::Matrix<6, 1>& modrcgrate)
{
  // define necessary variables
  const double modscale = std::pow(prinv(2), -1. / 3.);

  // 2nd Piola Kirchhoff stresses

  // isochoric contribution
  Core::LinAlg::Matrix<6, 1> modstress(Core::LinAlg::Initialization::zero);
  modstress.update(modmu(1), id2);
  modstress.update(modmu(2), modrcgrate, 1.0);
  // build 4-tensor for projection as 6x6 tensor
  Core::LinAlg::Matrix<6, 6> Projection;
  Projection.multiply_nt(1. / 3., icg, rcg);
  Projection.update(1.0, id4, -1.0);
  // isochoric stress
  stressisomodisovisco.multiply_nn(modscale, Projection, modstress, 1.0);

  // volumetric contribution:
  // with visco_isoratedep: no volumetric part added --> always 0


  // Constitutive Tensor

  // isochoric contribution
  // modified constitutive tensor
  Core::LinAlg::Matrix<6, 6> modcmat(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<6, 6> modcmat2(Core::LinAlg::Initialization::zero);
  // contribution:  Id \otimes \overline{\dot{C}} + \overline{\dot{C}} \otimes Id
  modcmat.multiply_nt(modxi(1), id2, modrcgrate);
  modcmat.multiply_nt(modxi(1), modrcgrate, id2, 1.0);
  // contribution: Id4
  modcmat.update(modxi(2), id4, 1.0);
  // scaling
  modcmat.scale(std::pow(modinv(2), -4. / 3.));
  // contribution: P:\overline{C}:P
  modcmat2.multiply_nn(Projection, modcmat);
  cmatisomodisovisco.multiply_nt(1.0, modcmat2, Projection, 1.0);
  // contribution: 2/3*Tr(J^(-2/3)modstress) (Cinv \odot Cinv - 1/3 Cinv \otimes Cinv)
  modcmat.clear();
  modcmat.multiply_nt(-1.0 / 3.0, icg, icg);
  Core::LinAlg::FourTensorOperations::add_holzapfel_product(modcmat, icg, 1.0);
  Core::LinAlg::Matrix<1, 1> tracemat;
  tracemat.multiply_tn(2. / 3. * std::pow(modinv(2), -2. / 3.), modstress, rcg);
  cmatisomodisovisco.update(tracemat(0, 0), modcmat, 1.0);
  // contribution: -2/3 (Cinv \otimes S_iso^v + S_iso^v \otimes Cinv)
  cmatisomodisovisco.multiply_nt(-2. / 3., icg, stressisomodisovisco, 1.0);
  cmatisomodisovisco.multiply_nt(-2. / 3., stressisomodisovisco, icg, 1.0);

  // volumetric contribution:
  // with visco_isoratedep: no volumetric part added --> always 0

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::evaluate_visco_generalized_maxwell(Core::LinAlg::Matrix<6, 1>& Q,
    Core::LinAlg::Matrix<6, 6>& cmatq, const double dt, const Core::LinAlg::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  const int visco_mat_id = params_ != nullptr ? params_->id() : -1;

  const auto& generalized_maxwell_metadata =
      require_generalized_maxwell_metadata("evaluating generalized Maxwell response", gp, eleGID);
  const int numbranch = static_cast<int>(generalized_maxwell_metadata.branches.size());
  const ViscoRuntimeContext& runtime_context =
      require_runtime_context("evaluating generalized Maxwell response", gp, eleGID);
  double one_step_theta = 0.5;
  if (generalized_maxwell_metadata.solve_kind == GeneralizedMaxwellSolveKind::one_step_theta)
  {
    if (!runtime_context.generalized_maxwell.has_value())
      FOUR_C_THROW(
          "Missing generalized Maxwell runtime context while evaluating generalized Maxwell "
          "response in MAT_ViscoElastHyper (MAT {}, GP {}, ELE {}).",
          visco_mat_id, gp, eleGID);

    one_step_theta = runtime_context.generalized_maxwell->one_step_theta;
  }

  Core::LinAlg::Matrix<6, 1> gl_stress(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Voigt::Strains::to_stress_like(glstrain, gl_stress);

  const Core::LinAlg::SymmetricTensor<double, 3, 3> gl =
      Core::LinAlg::make_symmetric_tensor_from_stress_like_voigt_matrix(gl_stress);

  Core::LinAlg::SymmetricTensor<double, 3, 3> C;
  evaluate_right_cauchy_green_strain_like_voigt(gl, C);
  const Core::LinAlg::SymmetricTensor<double, 3, 3> iC = Core::LinAlg::inv(C);

  Core::LinAlg::Matrix<6, 1> C_strain(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Voigt::Stresses::to_strain_like(
      Core::LinAlg::make_stress_like_voigt_view(C), C_strain);

  Core::LinAlg::Matrix<3, 1> prinv(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Voigt::Strains::invariants_principal(prinv, C_strain);

  Core::LinAlg::Matrix<6, 6> cmatqbranch(Core::LinAlg::Initialization::zero);
  std::vector<Core::LinAlg::Matrix<6, 1>> S(numbranch);
  std::vector<Core::LinAlg::Matrix<6, 1>> Qbranch(numbranch);

  // Read generalized Maxwell history
  const auto& S_n = state_.generalized_maxwell_prev_branch_elastic_stress(gp);
  const auto& Q_n = state_.generalized_maxwell_prev_branch_stress(gp);

  if (S_n.size() != static_cast<unsigned int>(numbranch))
    FOUR_C_THROW(
        "Invalid generalized Maxwell elastic branch history size in MAT_ViscoElastHyper (MAT {}, "
        "GP {}, ELE {}): expected {} entries but got {}.",
        visco_mat_id, gp, eleGID, numbranch, S_n.size());

  if (Q_n.size() != static_cast<unsigned int>(numbranch))
    FOUR_C_THROW(
        "Invalid generalized Maxwell viscous branch history size in MAT_ViscoElastHyper (MAT {}, "
        "GP {}, ELE {}): expected {} entries but got {}.",
        visco_mat_id, gp, eleGID, numbranch, Q_n.size());

  /////////////////////////////////////////////////
  // Loop over all viscoelastic Maxwell branches //
  /////////////////////////////////////////////////
  for (int i = 0; i < numbranch; ++i)
  {
    const auto& branch_metadata = generalized_maxwell_metadata.branches.at(i);
    const auto& branchpotsum = branch_metadata.summands;
    const double tau = branch_metadata.tau;

    Core::LinAlg::Matrix<3, 1> dPI(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<6, 1> ddPII(Core::LinAlg::Initialization::zero);

    elast_hyper_evaluate_invariant_derivatives(
        prinv, dPI, ddPII, branchpotsum, branch_metadata.properties, gp, eleGID);

    // blank resulting quantities
    // ... even if it is an implicit law that cmat is zero upon input
    S.at(i).clear();
    cmatqbranch.clear();

    // build stress response and elasticity tensor
    Core::LinAlg::SymmetricTensor<double, 3, 3> stressiso{};
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmatiso{};
    elast_hyper_add_isotropic_stress_cmat(stressiso, cmatiso, C, iC, prinv, dPI, ddPII);
    S.at(i).update(1.0, Core::LinAlg::make_stress_like_voigt_view(stressiso), 1.0);
    cmatqbranch.update(1.0, Core::LinAlg::make_stress_like_voigt_view(cmatiso), 1.0);

    // make sure Qbranch in this branch is empty
    Qbranch.at(i).clear();
    double deltascalar = 1.0;
    switch (generalized_maxwell_metadata.solve_kind)
    {
      case GeneralizedMaxwellSolveKind::one_step_theta:
      {
        const double theta = one_step_theta;

        // get time algorithmic parameters
        // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
        // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
        // evaluate scalars to compute
        // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + beta(S^(n+1) - S^n)]
        double lambdascalar1 = tau / (tau + theta * dt);
        double lambdascalar2 = (tau - dt + theta * dt) / tau;

        // same branch update factor as in the one-branch Maxwell case
        deltascalar = lambdascalar1;

        // calculate artificial viscous stresses Q
        // Q_(n+1) = lambdascalar1*[lamdascalar2* Q_n + (Sa_(n+1) - Sa_n)]
        Qbranch.at(i).update(lambdascalar2, Q_n.at(i), 1.0);
        Qbranch.at(i).update(1.0, S.at(i), 1.0);
        Qbranch.at(i).update(-1.0, S_n.at(i), 1.0);
        Qbranch.at(i).scale(lambdascalar1);
        break;
      }
      case GeneralizedMaxwellSolveKind::exponential_time_discretization:
      {
        double xiscalar1 = exp(-dt / tau);
        double xiscalar2 = exp(-dt / (2 * tau));

        deltascalar = xiscalar2;

        // calculate artificial stresses Q
        // Q_(n+1) = xiscalar1* Q_n + xiscalar2*(Sa_(n+1) - Sa_n)
        Qbranch.at(i).update(xiscalar1, Q_n.at(i), 1.0);
        Qbranch.at(i).update(xiscalar2, S.at(i), 1.0);
        Qbranch.at(i).update(-xiscalar2, S_n.at(i), 1.0);
        break;
      }
    }

    // sum up branches
    Q.update(1.0, Qbranch.at(i), 1.0);
    cmatq.update(deltascalar, cmatqbranch, 1.0);

  }  // end for loop over branches

  // update history
  state_.set_generalized_maxwell_current_point(gp, S, Qbranch);
}  // evaluate_visco_generalized_maxwell


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::evaluate_visco_fsls(Core::LinAlg::Matrix<6, 1> stress,
    Core::LinAlg::Matrix<6, 6> cmat, Core::LinAlg::Matrix<6, 1>& Q,
    Core::LinAlg::Matrix<6, 6>& cmatq, const double dt, const int gp, const int eleGID)
{
  const int visco_mat_id = params_ != nullptr ? params_->id() : -1;
  const FslsParameters fsls_parameters = read_fsls_parameters(gp, eleGID);
  const double tau = fsls_parameters.tau;
  const double alpha = fsls_parameters.alpha;
  const double beta = fsls_parameters.beta;

  if (dt < 0.0)
    FOUR_C_THROW(
        "Invalid time step size dt={} in VISCO_FSLS evaluation for MAT_ViscoElastHyper (MAT {}, "
        "GP {}, ELE {}). Expected dt >= 0.",
        dt, visco_mat_id, gp, eleGID);

  const auto& fsls_previous_history = state_.fsls_previous_history();
  if (fsls_previous_history.empty())
    FOUR_C_THROW("Missing FSLS history state in MAT_ViscoElastHyper (MAT {}, GP {}, ELE {}).",
        visco_mat_id, gp, eleGID);

  if (gp < 0 || gp >= static_cast<int>(fsls_previous_history.size()))
    FOUR_C_THROW(
        "Invalid Gauss point index GP={} for FSLS history in MAT_ViscoElastHyper (MAT {}, ELE "
        "{}). History container size is {}.",
        gp, visco_mat_id, eleGID, fsls_previous_history.size());

  const auto& fsls_history_at_gp = fsls_previous_history.at(gp);

  // read history of last time step at gp
  const int hs = fsls_history_at_gp.size();  // history size
  if (hs <= 0)
    FOUR_C_THROW(
        "Invalid FSLS history size {} at GP {} in MAT_ViscoElastHyper (MAT {}, ELE {}). "
        "Expected at least one entry.",
        hs, gp, visco_mat_id, eleGID);


  // calculate artificial history stress Qq with weights b_j
  // Qq = sum[j=1 up to j=n][b_j*Q_(n+1-j)] (short: b*Qj)

  // b_j = gamma(j-alpha)/[gamma(-alpha)*gamma(j+1)]
  // with recurstion formula for gamma functions b_j shortens to
  // b_j = (j-1-alpha)/j * b_(j-1)
  // with b_0 = 1 and b_1 = -alpha ...
  double bj = 1.;   // b_0=1
  double fac = 1.;  // pre-factor (j-1-alpha)/j  for calculation of b
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Qq(Core::LinAlg::Initialization::zero);

  // j=1...n, hs=n
  for (int j = 1; j <= hs; j++)
  {
    fac = (j - 1. - alpha) / j;
    bj = bj * fac;

    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Qj(fsls_history_at_gp.at(hs - j));
    Qq.update(bj, Qj, 1.0);
  }


  // calculate artificial stress Q

  // Version 1: As in Adolfson and Enelund (2003): Fractional Derivative Visocelasticity at Large
  // Deformations
  //  // initialize and evaluate scalars to compute
  //  // Q^(n+1) = [((dt/tau)^alpha)/(1+theta*(dt/tau)^alpha)]*[theta*S^(n+1)+(1-theta)(S^n-Q^n)]-
  //  //           [1/(1+theta*(dt/tau)^alpha)]*Qq^n


  // Version 2: Anna's Version of calculation
  // Difference:  1.) No one-step theta schema necessary
  //              2.) Introduce beta
  // Q^(n+1) = (dt^alpha / (dt^alpha + tau^alpha))*S^(n+1) - (tau^alpha / (dt^alpha +
  // tau^alpha))*Qq^n
  const double dtalpha = std::pow(dt, alpha);
  const double taualpha = std::pow(tau, alpha);
  const double denominator = dtalpha + taualpha;
  if (denominator <= 0.0)
    FOUR_C_THROW(
        "Invalid FSLS update denominator dt^alpha + tau^alpha = {} in MAT_ViscoElastHyper "
        "(MAT {}, GP {}, ELE {}): dt={}, tau={}, alpha={}. Expected a positive denominator.",
        denominator, visco_mat_id, gp, eleGID, dt, tau, alpha);

  const double lambdascalar1 = dtalpha / denominator;
  const double lambdascalar2 = -1. * taualpha / denominator;

  Q.update(lambdascalar1 * beta, stress, 0.);
  Q.update(lambdascalar2, Qq, 1.);


  // update history for next step
  state_.set_fsls_current_artificial_stress(gp, Q);  // Q_n+1


  // calculate final stress here and in Evaluate
  // S = elastic stress of Psi
  // S_2 = S ; S_1 = beta*S ; Q = Q(S1) = Q(beta*S)
  // S_final = S + beta*S - Q(beta*S)
  Q.update(beta, stress, -1.);

  // viscos constitutive tensor
  cmatq.update(lambdascalar1 * beta, cmat, 0.);  // contribution of Q
  cmatq.update(beta, cmat, -1.);


  return;
}

FOUR_C_NAMESPACE_CLOSE
