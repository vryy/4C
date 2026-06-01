// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_viscoelasthyper.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
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
#include "4C_mat_viscoelast_isoratedep.hpp"
#include "4C_mat_viscoelast_summand.hpp"
#include "4C_utils_enum.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

namespace
{
  bool is_visco_material_type_for_split_input(const Core::Materials::MaterialType material_type)
  {
    switch (material_type)
    {
      case Core::Materials::mes_isoratedep:
      case Core::Materials::mes_fsls:
      case Core::Materials::mes_generalizedmaxwell:
      case Core::Materials::mes_viscobranch:
      case Core::Materials::mes_coupmyocard:
        return true;
      default:
        return false;
    }
  }

}  // namespace

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ViscoElastHyper::ViscoElastHyper(const Core::Mat::PAR::Parameter::Data& matdata)
    : Mat::PAR::ElastHyper(matdata),
      numelast_(
          [&matdata, this]() -> int
          {
            const int* numelast = matdata.parameters.get_if<int>("NUMELAST");
            const std::vector<int>* elast_matids =
                matdata.parameters.get_if<std::vector<int>>("ELAST_MATIDS");
            const int* numvisco = matdata.parameters.get_if<int>("NUMVISCO");
            const std::vector<int>* visco_matids =
                matdata.parameters.get_if<std::vector<int>>("VISCO_MATIDS");

            const bool has_partial_split = numelast != nullptr || elast_matids != nullptr ||
                                           numvisco != nullptr || visco_matids != nullptr;
            const bool has_complete_split = numelast != nullptr && elast_matids != nullptr &&
                                            numvisco != nullptr && visco_matids != nullptr;
            const bool explicit_split =
                has_complete_split && !elast_matids->empty() && !visco_matids->empty();
            if (has_partial_split && !has_complete_split)
              FOUR_C_THROW(
                  "Split MAT_ViscoElastHyper input requires NUMELAST, ELAST_MATIDS, "
                  "NUMVISCO and VISCO_MATIDS together (MAT {}).",
                  id());
            if (explicit_split)
            {
              return *numelast;
            }

            int derived_numelast = 0;
            const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
            for (const int summand_mat_id : matids_)
            {
              auto* sumpar =
                  Global::Problem::instance(probinst)->materials()->parameter_by_id(summand_mat_id);
              if (!is_visco_material_type_for_split_input(sumpar->type())) ++derived_numelast;
            }
            return derived_numelast;
          }()),
      elast_matids_(
          [&matdata, this]() -> std::vector<int>
          {
            const int* numelast = matdata.parameters.get_if<int>("NUMELAST");
            const std::vector<int>* elast_matids =
                matdata.parameters.get_if<std::vector<int>>("ELAST_MATIDS");
            const int* numvisco = matdata.parameters.get_if<int>("NUMVISCO");
            const std::vector<int>* visco_matids =
                matdata.parameters.get_if<std::vector<int>>("VISCO_MATIDS");

            const bool has_complete_split = numelast != nullptr && elast_matids != nullptr &&
                                            numvisco != nullptr && visco_matids != nullptr;
            const bool explicit_split =
                has_complete_split && !elast_matids->empty() && !visco_matids->empty();
            if (explicit_split) return *elast_matids;

            std::vector<int> derived_elast;
            const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
            for (const int summand_mat_id : matids_)
            {
              auto* sumpar =
                  Global::Problem::instance(probinst)->materials()->parameter_by_id(summand_mat_id);
              if (!is_visco_material_type_for_split_input(sumpar->type()))
                derived_elast.push_back(summand_mat_id);
            }
            return derived_elast;
          }()),
      numvisco_(
          [&matdata, this]() -> int
          {
            const int* numelast = matdata.parameters.get_if<int>("NUMELAST");
            const std::vector<int>* elast_matids =
                matdata.parameters.get_if<std::vector<int>>("ELAST_MATIDS");
            const int* numvisco = matdata.parameters.get_if<int>("NUMVISCO");
            const std::vector<int>* visco_matids =
                matdata.parameters.get_if<std::vector<int>>("VISCO_MATIDS");

            const bool has_partial_split = numelast != nullptr || elast_matids != nullptr ||
                                           numvisco != nullptr || visco_matids != nullptr;
            const bool has_complete_split = numelast != nullptr && elast_matids != nullptr &&
                                            numvisco != nullptr && visco_matids != nullptr;
            const bool explicit_split =
                has_complete_split && !elast_matids->empty() && !visco_matids->empty();
            if (has_partial_split && !has_complete_split)
              FOUR_C_THROW(
                  "Split MAT_ViscoElastHyper input requires NUMELAST, ELAST_MATIDS, "
                  "NUMVISCO and VISCO_MATIDS together (MAT {}).",
                  id());
            if (explicit_split)
            {
              return *numvisco;
            }

            int derived_numvisco = 0;
            const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
            for (const int summand_mat_id : matids_)
            {
              auto* sumpar =
                  Global::Problem::instance(probinst)->materials()->parameter_by_id(summand_mat_id);
              if (is_visco_material_type_for_split_input(sumpar->type())) ++derived_numvisco;
            }
            return derived_numvisco;
          }()),
      visco_matids_(
          [&matdata, this]() -> std::vector<int>
          {
            const int* numelast = matdata.parameters.get_if<int>("NUMELAST");
            const std::vector<int>* elast_matids =
                matdata.parameters.get_if<std::vector<int>>("ELAST_MATIDS");
            const int* numvisco = matdata.parameters.get_if<int>("NUMVISCO");
            const std::vector<int>* visco_matids =
                matdata.parameters.get_if<std::vector<int>>("VISCO_MATIDS");

            const bool has_complete_split = numelast != nullptr && elast_matids != nullptr &&
                                            numvisco != nullptr && visco_matids != nullptr;
            const bool explicit_split =
                has_complete_split && !elast_matids->empty() && !visco_matids->empty();
            if (explicit_split) return *visco_matids;

            std::vector<int> derived_visco;
            const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
            for (const int summand_mat_id : matids_)
            {
              auto* sumpar =
                  Global::Problem::instance(probinst)->materials()->parameter_by_id(summand_mat_id);
              if (is_visco_material_type_for_split_input(sumpar->type()))
                derived_visco.push_back(summand_mat_id);
            }
            return derived_visco;
          }()),
      uses_legacy_matids_(
          [&matdata]() -> bool
          {
            const int* numelast = matdata.parameters.get_if<int>("NUMELAST");
            const std::vector<int>* elast_matids =
                matdata.parameters.get_if<std::vector<int>>("ELAST_MATIDS");
            const int* numvisco = matdata.parameters.get_if<int>("NUMVISCO");
            const std::vector<int>* visco_matids =
                matdata.parameters.get_if<std::vector<int>>("VISCO_MATIDS");
            const bool has_complete_split = numelast != nullptr && elast_matids != nullptr &&
                                            numvisco != nullptr && visco_matids != nullptr;
            const bool explicit_split =
                has_complete_split && !elast_matids->empty() && !visco_matids->empty();
            return !explicit_split;
          }())
{
  if (numelast_ != static_cast<int>(elast_matids_.size()))
    FOUR_C_THROW(
        "Invalid MAT_ViscoElastHyper setup (MAT {}): NUMELAST={} but ELAST_MATIDS has size {}.",
        id(), numelast_, elast_matids_.size());

  if (numvisco_ != static_cast<int>(visco_matids_.size()))
    FOUR_C_THROW(
        "Invalid MAT_ViscoElastHyper setup (MAT {}): NUMVISCO={} but VISCO_MATIDS has size {}.",
        id(), numvisco_, visco_matids_.size());

  if (uses_legacy_matids_ &&
      static_cast<int>(elast_matids_.size() + visco_matids_.size()) != nummat_)
    FOUR_C_THROW(
        "Invalid MAT_ViscoElastHyper legacy partitioning (MAT {}): NUMMAT={} but partitioned "
        "ELAST_MATIDS({}) + VISCO_MATIDS({}) do not match.",
        id(), nummat_, elast_matids_.size(), visco_matids_.size());

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
    : Mat::ElastHyper(), isovisco_(false), visco_generalized_maxwell_(false), visco_fsls_(false)
{
  params_ = params;
  rebuild_summand_sets();
  rebuild_active_model_sequence();
}


bool Mat::ViscoElastHyper::is_visco_material_type(const Core::Materials::MaterialType material_type)
{
  return is_visco_material_type_for_split_input(material_type);
}


const Mat::PAR::ViscoElastHyper* Mat::ViscoElastHyper::visco_params() const
{
  return static_cast<const Mat::PAR::ViscoElastHyper*>(params_);
}


int Mat::ViscoElastHyper::visco_mat_id(const unsigned index) const
{
  const Mat::PAR::ViscoElastHyper* visco_par = visco_params();
  if (visco_par == nullptr) return -1;
  if (index >= visco_par->visco_matids_.size())
    FOUR_C_THROW(
        "Invalid visco summand index {} for MAT_ViscoElastHyper (MAT {}).", index, visco_par->id());
  return visco_par->visco_matids_.at(index);
}


void Mat::ViscoElastHyper::rebuild_summand_sets()
{
  elast_potsum_.clear();
  visco_potsum_.clear();

  const Mat::PAR::ViscoElastHyper* visco_par = visco_params();
  if (visco_par == nullptr) return;

  for (const int summand_mat_id : visco_par->elast_matids_)
  {
    auto sum = Mat::Elastic::Summand::factory(summand_mat_id);
    if (sum == nullptr) FOUR_C_THROW("Failed to allocate elastic summand MAT {}.", summand_mat_id);
    elast_potsum_.push_back(sum);
    sum->register_anisotropy_extensions(anisotropy_);
  }

  for (const int summand_mat_id : visco_par->visco_matids_)
  {
    auto sum = Mat::ViscoElast::Summand::factory(summand_mat_id);
    if (sum == nullptr) FOUR_C_THROW("Failed to allocate visco summand MAT {}.", summand_mat_id);
    visco_potsum_.push_back(sum);
  }
}


void Mat::ViscoElastHyper::rebuild_effective_summand_properties()
{
  effective_summand_properties_.clear();
  effective_summand_properties_.merge(elast_summand_properties_);

  SummandProperties visco_formulation_properties;
  visco_formulation_properties.clear();

  for (const auto& p : visco_potsum_)
  {
    p->specify_formulation(visco_formulation_properties.isoprinc,
        visco_formulation_properties.isomod, visco_formulation_properties.anisoprinc,
        visco_formulation_properties.anisomod, visco_formulation_properties.viscoGeneral);
  }

  effective_summand_properties_.merge(visco_formulation_properties);

  // Keep stretch coefficient flags tied to elastic constitutive summands.
  effective_summand_properties_.coeffStretchesPrinc = elast_summand_properties_.coeffStretchesPrinc;
  effective_summand_properties_.coeffStretchesMod = elast_summand_properties_.coeffStretchesMod;
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

  std::shared_ptr<Mat::ViscoElast::GeneralizedMaxwell> generalized_maxwell = nullptr;
  int generalized_maxwell_summand_mat_id = -1;

  for (unsigned int p = 0; p < visco_potsum_.size(); ++p)
  {
    auto current_generalized_maxwell =
        std::dynamic_pointer_cast<Mat::ViscoElast::GeneralizedMaxwell>(visco_potsum_[p]);
    if (current_generalized_maxwell != nullptr)
    {
      ++generalized_maxwell_model_count;
      generalized_maxwell = current_generalized_maxwell;
      generalized_maxwell_summand_mat_id = this->visco_mat_id(p);
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
  for (unsigned int p = 0; p < visco_potsum_.size(); ++p)
  {
    std::shared_ptr<Mat::ViscoElast::Fsls> fsls =
        std::dynamic_pointer_cast<Mat::ViscoElast::Fsls>(visco_potsum_[p]);
    if (fsls != nullptr)
    {
      ++fsls_model_count;
      fsls_metadata.summand_mat_id = this->visco_mat_id(p);
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
  elast_summand_properties_.pack(data);
  const Mat::ViscoElastState::ActiveModels models = active_models();
  add_to_pack(data, models.iso_rate);
  add_to_pack(data, models.generalized_maxwell);
  add_to_pack(data, models.fsls);

  anisotropy_.pack_anisotropy(data);

  Core::Communication::PotentiallyUnusedBufferScope summand_scope{data};
  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    for (auto& p : elast_potsum_) p->pack_summand(data);
    for (auto& p : visco_potsum_) p->pack_summand(data);

    state_.serialize_state(data, active_models());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::unpack(Core::Communication::UnpackBuffer& buffer)
{
  // make sure we have a pristine material
  params_ = nullptr;
  elast_potsum_.clear();
  visco_potsum_.clear();
  elast_summand_properties_.clear();
  effective_summand_properties_.clear();
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

  elast_summand_properties_.unpack(buffer);
  extract_from_pack(buffer, isovisco_);
  extract_from_pack(buffer, visco_generalized_maxwell_);
  extract_from_pack(buffer, visco_fsls_);
  rebuild_active_model_sequence();
  ensure_model_activation_consistency("unpack model activation");

  anisotropy_.unpack_anisotropy(buffer);

  Core::Communication::PotentiallyUnusedBufferScope summand_scope{buffer};
  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    rebuild_summand_sets();

    for (auto& p : elast_potsum_) p->unpack_summand(buffer);
    for (auto& p : visco_potsum_) p->unpack_summand(buffer);

    if (is_model_active(ViscoModelKind::generalized_maxwell))
    {
      build_generalized_maxwell_metadata_for_setup(-1, -1);
    }

    build_runtime_context_for_setup(-1, -1);

    if (is_model_active(ViscoModelKind::fsls))
    {
      build_fsls_metadata_for_setup(-1, -1);
    }

    rebuild_effective_summand_properties();

    state_.deserialize_state(buffer, active_models());
  }
  else
  {
    effective_summand_properties_.update(elast_summand_properties_);
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

  for (auto& p : elast_potsum_)
  {
    p->setup(numgp, fibers, coord_system);
  }
  for (auto& p : visco_potsum_) p->setup(numgp, fibers, coord_system);

  // find out which formulations are used
  isovisco_ = false;
  visco_generalized_maxwell_ = false;
  visco_fsls_ = false;

  elast_summand_properties_.clear();
  elast_hyper_properties(elast_potsum_, elast_summand_properties_);
  rebuild_effective_summand_properties();

  for (auto& p : visco_potsum_)
    p->specify_visco_formulation(isovisco_, visco_generalized_maxwell_, visco_fsls_);

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
  for (auto& p : elast_potsum_) p->update();
  for (auto& p : visco_potsum_) p->update();

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

  elast_hyper_evaluate_invariant_derivatives(workspace.prinv, workspace.dPI, workspace.ddPII,
      elast_potsum_, elast_summand_properties_, gp, eleGID);
}


void Mat::ViscoElastHyper::prepare_iso_rate_visco_inputs_if_active(
    const Teuchos::ParameterList& params, EvaluateWorkspace& workspace, const int gp,
    const int eleGID)
{
  for (const ViscoModelKind model_kind : active_model_sequence_)
  {
    if (model_kind != ViscoModelKind::iso_rate) continue;

    if (effective_summand_properties_.isomod)
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
  if (effective_summand_properties_.isomod)
  {
    // add viscous part decoupled
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> stressisomodisovisco(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisomodisovisco(
        Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> stressisomodvolvisco(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisomodvolvisco(
        Core::LinAlg::Initialization::zero);
    Mat::ViscoElast::Kernels::evaluate_iso_visco_modified_kernel(stressisomodisovisco,
        stressisomodvolvisco, cmatisomodisovisco, cmatisomodvolvisco, workspace.prinv,
        workspace.modinv, workspace.modmu, workspace.modxi, workspace.c_strain, workspace.id2,
        workspace.i_c_stress, workspace.id4, workspace.modrcgrate);
    stress_view.update(1.0, stressisomodisovisco, 1.0);
    stress_view.update(1.0, stressisomodvolvisco, 1.0);
    cmat_view.update(1.0, cmatisomodisovisco, 1.0);
    cmat_view.update(1.0, cmatisomodvolvisco, 1.0);
  }

  if (effective_summand_properties_.isoprinc)
  {
    // add viscous part coupled
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> stressisovisco(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisovisco(
        Core::LinAlg::Initialization::zero);
    Mat::ViscoElast::Kernels::evaluate_iso_visco_principal_kernel(stressisovisco, cmatisovisco,
        workspace.mu, workspace.xi, workspace.id4sharp, workspace.scgrate);
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
  if (elast_summand_properties_.coeffStretchesPrinc || elast_summand_properties_.coeffStretchesMod)
  {
    elast_hyper_add_response_stretches(
        cmat, stress, workspace.c, elast_potsum_, elast_summand_properties_, gp, eleGID);
  }

  /*----------------------------------------------------------------------*/
  // Do all the anisotropic stuff!
  if (elast_summand_properties_.anisoprinc)
  {
    elast_hyper_add_anisotropic_princ(stress, cmat, workspace.c, params, gp, eleGID, elast_potsum_);
  }

  if (elast_summand_properties_.anisomod)
  {
    elast_hyper_add_anisotropic_mod(stress, cmat, workspace.c, workspace.i_c, workspace.prinv, gp,
        eleGID, context, elast_potsum_);
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
  const auto iso_rate_previous = state_.iso_rate_prev_point(gp);

  Mat::ViscoElast::Kernels::evaluate_kin_quant_vis_kernel(rcg, scg, prinv, iso_rate_previous.scg,
      iso_rate_previous.modrcg, dt, modrcg, scgrate, modrcgrate, modrateinv,
      params_ != nullptr ? params_->id() : -1, gp);

  state_.set_iso_rate_current_point(gp, scg, modrcg);
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
  Mat::ViscoElast::Kernels::evaluate_mu_xi_kernel(visco_potsum_,
      effective_summand_properties_.isoprinc, effective_summand_properties_.isomod, prinv, modinv,
      mu, modmu, xi, modxi, rateinv, modrateinv, params, dt, gp, eleGID);
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

  // Read generalized Maxwell history
  const auto& S_n = state_.generalized_maxwell_prev_branch_elastic_stress(gp);
  const auto& Q_n = state_.generalized_maxwell_prev_branch_stress(gp);
  const auto collect_branch_taus = [&]()
  {
    std::vector<double> branch_taus;
    branch_taus.reserve(numbranch);
    for (const auto& branch_metadata : generalized_maxwell_metadata.branches)
      branch_taus.push_back(branch_metadata.tau);
    return branch_taus;
  };

  const auto make_generalized_maxwell_kernel_input = [&]()
  {
    Mat::ViscoElast::Kernels::GeneralizedMaxwellKernelInput kernel_input;
    kernel_input.visco_mat_id = visco_mat_id;
    kernel_input.gp = gp;
    kernel_input.ele_gid = eleGID;
    kernel_input.dt = dt;
    kernel_input.one_step_theta = one_step_theta;
    kernel_input.previous_branch_elastic_stress = &S_n;
    kernel_input.previous_branch_stress = &Q_n;
    switch (generalized_maxwell_metadata.solve_kind)
    {
      case GeneralizedMaxwellSolveKind::one_step_theta:
        kernel_input.solve_kind =
            Mat::ViscoElast::Kernels::GeneralizedMaxwellSolveKind::one_step_theta;
        break;
      case GeneralizedMaxwellSolveKind::exponential_time_discretization:
        kernel_input.solve_kind =
            Mat::ViscoElast::Kernels::GeneralizedMaxwellSolveKind::exponential_time_discretization;
        break;
    }
    return kernel_input;
  };

  const std::vector<double> branch_taus = collect_branch_taus();
  const Mat::ViscoElast::Kernels::GeneralizedMaxwellKernelInput kernel_input =
      make_generalized_maxwell_kernel_input();

  Mat::ViscoElast::Kernels::PointHistory current_branch_elastic_stress;
  Mat::ViscoElast::Kernels::PointHistory current_branch_stress;

  const Mat::ViscoElast::Kernels::BranchResponseEvaluator evaluate_branch_response =
      [&](const int branch_index, Mat::ViscoElast::Kernels::StressVector& branch_elastic_stress,
          Mat::ViscoElast::Kernels::TangentMatrix& branch_cmat)
  {
    const auto& branch_metadata = generalized_maxwell_metadata.branches.at(branch_index);

    Core::LinAlg::Matrix<3, 1> dPI(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<6, 1> ddPII(Core::LinAlg::Initialization::zero);
    elast_hyper_evaluate_invariant_derivatives(
        prinv, dPI, ddPII, branch_metadata.summands, branch_metadata.properties, gp, eleGID);

    branch_elastic_stress.clear();
    branch_cmat.clear();

    Core::LinAlg::SymmetricTensor<double, 3, 3> stressiso{};
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3> cmatiso{};
    elast_hyper_add_isotropic_stress_cmat(stressiso, cmatiso, C, iC, prinv, dPI, ddPII);
    branch_elastic_stress.update(1.0, Core::LinAlg::make_stress_like_voigt_view(stressiso), 1.0);
    branch_cmat.update(1.0, Core::LinAlg::make_stress_like_voigt_view(cmatiso), 1.0);
  };

  Mat::ViscoElast::Kernels::evaluate_generalized_maxwell_kernel(Q, cmatq,
      current_branch_elastic_stress, current_branch_stress, branch_taus, kernel_input,
      evaluate_branch_response);

  state_.set_generalized_maxwell_current_point(
      gp, current_branch_elastic_stress, current_branch_stress);
}  // evaluate_visco_generalized_maxwell


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::evaluate_visco_fsls(Core::LinAlg::Matrix<6, 1> stress,
    Core::LinAlg::Matrix<6, 6> cmat, Core::LinAlg::Matrix<6, 1>& Q,
    Core::LinAlg::Matrix<6, 6>& cmatq, const double dt, const int gp, const int eleGID)
{
  const int visco_mat_id = params_ != nullptr ? params_->id() : -1;
  const FslsParameters fsls_parameters = read_fsls_parameters(gp, eleGID);
  const auto& fsls_previous_history = state_.fsls_previous_history();

  const auto make_fsls_kernel_input = [&]()
  {
    Mat::ViscoElast::Kernels::FslsKernelInput kernel_input;
    kernel_input.visco_mat_id = visco_mat_id;
    kernel_input.gp = gp;
    kernel_input.ele_gid = eleGID;
    kernel_input.dt = dt;
    kernel_input.tau = fsls_parameters.tau;
    kernel_input.alpha = fsls_parameters.alpha;
    kernel_input.beta = fsls_parameters.beta;
    kernel_input.previous_history = &fsls_previous_history;
    return kernel_input;
  };

  const Mat::ViscoElast::Kernels::FslsKernelInput kernel_input = make_fsls_kernel_input();

  Mat::ViscoElast::Kernels::FslsStressVector q_current_for_history(
      Core::LinAlg::Initialization::zero);
  Mat::ViscoElast::Kernels::evaluate_fsls_kernel(
      stress, cmat, q_current_for_history, Q, cmatq, kernel_input);

  state_.set_fsls_current_artificial_stress(gp, q_current_for_history);


  return;
}

FOUR_C_NAMESPACE_CLOSE
