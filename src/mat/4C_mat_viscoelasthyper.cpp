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
#include "4C_mat_elast_visco_fsls.hpp"
#include "4C_mat_elast_visco_generalizedmaxwell.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"
#include "4C_mat_so3_material.hpp"
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ViscoElastHyper::ViscoElastHyper(Mat::PAR::ViscoElastHyper* params)
    : Mat::ElastHyper(params),
      isovisco_(false),
      visco_generalized_maxwell_(false),
      visco_fsls_(false)
{
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
  add_to_pack(data, isovisco_);
  add_to_pack(data, visco_generalized_maxwell_);
  add_to_pack(data, visco_fsls_);

  anisotropy_.pack_anisotropy(data);

  Core::Communication::PotentiallyUnusedBufferScope summand_scope{data};
  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      potsum_[p]->pack_summand(data);
    }

    //  pack history data
    const int histsize = state_.packed_history_size();
    add_to_pack(data, histsize);  // Length of history vector(s)
    state_.pack_kinematic_history(data, histsize);

    if (visco_generalized_maxwell_) state_.pack_generalized_maxwell_history(data, histsize);

    // pack history of FSLS-model
    if (visco_fsls_) state_.pack_fsls_history(data);
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

  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

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

    // history data 09/13
    state_.set_state_initialized(true);
    int histsize;
    extract_from_pack(buffer, histsize);

    if (histsize == 0) state_.set_state_initialized(false);

    state_.unpack_kinematic_history(buffer, histsize);

    if (visco_generalized_maxwell_) state_.unpack_generalized_maxwell_history(buffer, histsize);

    // for FSLS-model
    if (visco_fsls_) state_.unpack_fsls_history(buffer, histsize);
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

  // Initialise/allocate history variables 09/13
  state_.initialize_kinematic_history(numgp);

  if (visco_generalized_maxwell_)
  {
    int generalized_maxwell_model_count = 0;
    int generalized_maxwell_numbranch = -1;
    std::string generalized_maxwell_solve;
    const std::vector<int>* generalized_maxwell_matids = nullptr;

    for (unsigned int p = 0; p < potsum_.size(); ++p)
    {
      std::shared_ptr<Mat::Elastic::GeneralizedMaxwell> generalized_maxwell =
          std::dynamic_pointer_cast<Mat::Elastic::GeneralizedMaxwell>(potsum_[p]);
      if (generalized_maxwell != nullptr)
      {
        ++generalized_maxwell_model_count;
        generalized_maxwell->read_material_parameters(
            generalized_maxwell_numbranch, generalized_maxwell_matids, generalized_maxwell_solve);
      }
    }

    if (generalized_maxwell_model_count != 1)
      FOUR_C_THROW(
          "Invalid VISCO_GeneralizedMaxwell setup in MAT_ViscoElastHyper (MAT {}): expected "
          "exactly one VISCO_GeneralizedMaxwell summand but found {}.",
          params_->id(), generalized_maxwell_model_count);

    if (generalized_maxwell_numbranch <= 0)
      FOUR_C_THROW(
          "Invalid VISCO_GeneralizedMaxwell setup in MAT_ViscoElastHyper (MAT {}): "
          "NUMBRANCH={} is not positive.",
          params_->id(), generalized_maxwell_numbranch);

    if (generalized_maxwell_matids == nullptr)
      FOUR_C_THROW(
          "Failed to read MATIDS for VISCO_GeneralizedMaxwell in MAT_ViscoElastHyper (MAT "
          "{}).",
          params_->id());

    if (generalized_maxwell_matids->size() !=
        static_cast<unsigned int>(generalized_maxwell_numbranch))
      FOUR_C_THROW(
          "Invalid VISCO_GeneralizedMaxwell setup in MAT_ViscoElastHyper (MAT {}): "
          "NUMBRANCH={} but MATIDS has size {}.",
          params_->id(), generalized_maxwell_numbranch, generalized_maxwell_matids->size());

    if (generalized_maxwell_solve != "OneStepTheta" &&
        generalized_maxwell_solve != "ExponentialTimeDiscretization")
      FOUR_C_THROW(
          "Invalid SOLVE='{}' in VISCO_GeneralizedMaxwell for MAT_ViscoElastHyper (MAT {}). "
          "Use OneStepTheta or ExponentialTimeDiscretization.",
          generalized_maxwell_solve, params_->id());

    state_.initialize_generalized_maxwell_history(numgp, generalized_maxwell_numbranch);
  }

  // in case of FSLS-model
  if (visco_fsls_)
  {
    state_.initialize_fsls_history(numgp);
  }

  state_.set_state_initialized(true);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::update()
{
  Mat::ElastHyper::update();

  // Update history values 09/13
  state_.commit_kinematic_history();

  // for FSLS-model
  if (visco_fsls_)
  {
    // To calculate the fractional derivative the history of all previous timesteps is saved in each
    // gauss-point

    // numsteps
    const Teuchos::ParameterList& sdyn = Global::Problem::instance()->structural_dynamic_params();
    const int numsteps = sdyn.get<int>("NUMSTEP");
    // maximal size of history (in time steps)
    const unsigned int max_hist = numsteps + 1;

    state_.append_fsls_history(max_hist);
  }

  const int numgp = state_.gauss_point_count();
  state_.reset_current_iteration(numgp);

  if (visco_generalized_maxwell_)
  {
    const int visco_mat_id = params_ != nullptr ? params_->id() : -1;
    const std::size_t numbranch = state_.rollover_generalized_maxwell_history(numgp, visco_mat_id);
    state_.reset_generalized_maxwell_current(numgp, numbranch);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::evaluate(const Core::LinAlg::Tensor<double, 3, 3>* defgrad,
    const Core::LinAlg::SymmetricTensor<double, 3, 3>& glstrain,
    const Teuchos::ParameterList& params, const EvaluationContext<3>& context,
    Core::LinAlg::SymmetricTensor<double, 3, 3>& stress,
    Core::LinAlg::SymmetricTensor<double, 3, 3, 3, 3>& cmat, int gp, int eleGID)
{
  const Core::LinAlg::Matrix<6, 1> glstrain_mat =
      Core::LinAlg::make_strain_like_voigt_matrix(glstrain);
  Core::LinAlg::Matrix<6, 1> stress_view = Core::LinAlg::make_stress_like_voigt_view(stress);
  Core::LinAlg::Matrix<6, 6> cmat_view = Core::LinAlg::make_stress_like_voigt_view(cmat);

  Core::LinAlg::Matrix<6, 1> modC_strain(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<6, 1> id2(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<6, 1> modrcg(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<6, 6> id4(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<6, 6> id4sharp(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<3, 1> prinv(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<3, 1> modinv(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<7, 1> rateinv(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<7, 1> modrateinv(Core::LinAlg::Initialization::zero);

  Core::LinAlg::Matrix<3, 1> dPI(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<6, 1> ddPII(Core::LinAlg::Initialization::zero);

  Core::LinAlg::Matrix<6, 1> scgrate(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<6, 1> modrcgrate(Core::LinAlg::Initialization::zero);

  // for extension: Core::LinAlg::Matrix<6,1> modicgrate(true);
  Core::LinAlg::Matrix<8, 1> mu(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<8, 1> modmu(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<33, 1> xi(Core::LinAlg::Initialization::zero);
  Core::LinAlg::Matrix<33, 1> modxi(Core::LinAlg::Initialization::zero);

  Core::LinAlg::SymmetricTensor<double, 3, 3> C{};
  evaluate_right_cauchy_green_strain_like_voigt(glstrain, C);
  Core::LinAlg::SymmetricTensor<double, 3, 3> iC = Core::LinAlg::inv(C);

  Core::LinAlg::Matrix<6, 1> C_stress = Core::LinAlg::make_stress_like_voigt_view(C);
  Core::LinAlg::Matrix<6, 1> iC_stress = Core::LinAlg::make_stress_like_voigt_view(iC);
  Core::LinAlg::Matrix<6, 1> C_strain;
  Core::LinAlg::Voigt::Stresses::to_strain_like(C_stress, C_strain);
  Core::LinAlg::Matrix<6, 1> iC_strain;
  Core::LinAlg::Voigt::Stresses::to_strain_like(iC_stress, iC_strain);
  Core::LinAlg::Voigt::Stresses::invariants_principal(prinv, C_stress);

  const bool has_visco_contribution = isovisco_ || visco_generalized_maxwell_ || visco_fsls_;
  double dt = 0.0;
  if (has_visco_contribution)
  {
    if (context.time_step_size == nullptr)
      FOUR_C_THROW(
          "Missing EvaluationContext::time_step_size in MAT_ViscoElastHyper (MAT {}, GP {}, "
          "ELE {}).",
          params_ != nullptr ? params_->id() : -1, gp, eleGID);
    dt = *context.time_step_size;
  }


  Core::LinAlg::Voigt::identity_matrix(id2);

  using VoigtNotation = Core::LinAlg::Voigt::NotationType;
  Core::LinAlg::Voigt::fourth_order_identity_matrix<VoigtNotation::stress, VoigtNotation::stress>(
      id4sharp);
  Core::LinAlg::Voigt::fourth_order_identity_matrix<VoigtNotation::stress, VoigtNotation::strain>(
      id4);

  elast_hyper_evaluate_invariant_derivatives(
      prinv, dPI, ddPII, potsum_, summandProperties_, gp, eleGID);

  if (isovisco_)
  {
    if (summandProperties_.isomod)
    {
      // calculate modified invariants
      invariants_modified(modinv, prinv);
    }
    // calculate viscous quantities
    evaluate_kin_quant_vis(C_strain, C_stress, iC_stress, prinv, rateinv, modC_strain, dt, scgrate,
        modrcgrate, modrateinv, gp);
    evaluate_mu_xi(
        prinv, modinv, mu, modmu, xi, modxi, rateinv, modrateinv, params, dt, gp, eleGID);
  }

  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress_view.clear();
  cmat_view.clear();

  // add isotropic part
  elast_hyper_add_isotropic_stress_cmat(stress, cmat, C, iC, prinv, dPI, ddPII);


  // add viscous part
  if (isovisco_)
  {
    if (summandProperties_.isomod)
    {
      // add viscous part decoupled
      Core::LinAlg::Matrix<NUM_STRESS_3D, 1> stressisomodisovisco(
          Core::LinAlg::Initialization::zero);
      Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisomodisovisco(
          Core::LinAlg::Initialization::zero);
      Core::LinAlg::Matrix<NUM_STRESS_3D, 1> stressisomodvolvisco(
          Core::LinAlg::Initialization::zero);
      Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatisomodvolvisco(
          Core::LinAlg::Initialization::zero);
      evaluate_iso_visco_modified(stressisomodisovisco, stressisomodvolvisco, cmatisomodisovisco,
          cmatisomodvolvisco, prinv, modinv, modmu, modxi, C_strain, id2, iC_stress, id4,
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
      evaluate_iso_visco_principal(stressisovisco, cmatisovisco, mu, xi, id4sharp, scgrate);
      stress_view.update(1.0, stressisovisco, 1.0);
      cmat_view.update(1.0, cmatisovisco, 1.0);
    }
  }

  // add contribution of generalized Maxwell model
  if (visco_generalized_maxwell_)
  {
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Q(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatq(Core::LinAlg::Initialization::zero);
    evaluate_visco_generalized_maxwell(Q, cmatq, dt, &glstrain_mat, gp, eleGID);
    stress_view.update(1.0, Q, 1.0);
    cmat_view.update(1.0, cmatq, 1.0);
  }

  // add contribution of FSLS material
  if (visco_fsls_)
  {
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Q(
        Core::LinAlg::Initialization::zero);  // artificial viscous stress
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> cmatq(Core::LinAlg::Initialization::zero);
    evaluate_visco_fsls(stress_view, cmat_view, Q, cmatq, dt, gp, eleGID);
    stress_view.update(1.0, Q, 1.);
    cmat_view.update(1.0, cmatq, 1.);
  }


  /*----------------------------------------------------------------------*/
  // coefficients in principal stretches
  if (summandProperties_.coeffStretchesPrinc || summandProperties_.coeffStretchesMod)
  {
    elast_hyper_add_response_stretches(cmat, stress, C, potsum_, summandProperties_, gp, eleGID);
  }

  /*----------------------------------------------------------------------*/
  // Do all the anisotropic stuff!
  if (summandProperties_.anisoprinc)
  {
    elast_hyper_add_anisotropic_princ(stress, cmat, C, params, gp, eleGID, potsum_);
  }

  if (summandProperties_.anisomod)
  {
    elast_hyper_add_anisotropic_mod(stress, cmat, C, iC, prinv, gp, eleGID, context, potsum_);
  }
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
  Core::LinAlg::Matrix<6, 1> scglast(state_.scg_last_at(gp));
  Core::LinAlg::Matrix<6, 1> modrcglast(state_.modrcg_last_at(gp));

  // Update history of Cauchy-Green Tensor
  state_.set_scg_current_at(gp, scg);        // principal material: store C^{n}
  state_.set_modrcg_current_at(gp, modrcg);  // decoupled material: store \overline{C}^{n}

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
    Core::LinAlg::Matrix<6, 6>& cmatq, const double dt, const Core::LinAlg::Matrix<6, 1>* glstrain,
    const int gp, const int eleGID)
{
  const int visco_mat_id = params_ != nullptr ? params_->id() : -1;

  int numbranch = -1;
  std::string solve = "";
  const std::vector<int>* matids = nullptr;
  std::vector<std::shared_ptr<Mat::Elastic::Summand>> branchpotsum(
      0);  // vector of summands in one branch
  std::vector<std::vector<std::shared_ptr<Mat::Elastic::Summand>>> branchespotsum(
      0);  // vector for each branch of vectors of summands in each branch
  std::vector<double> branchtau(0);
  int generalized_maxwell_model_count = 0;

  // get parameters of ViscoGeneralizedMaxwell
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    std::shared_ptr<Mat::Elastic::GeneralizedMaxwell> generalized_maxwell =
        std::dynamic_pointer_cast<Mat::Elastic::GeneralizedMaxwell>(potsum_[p]);

    if (generalized_maxwell != nullptr)
    {
      ++generalized_maxwell_model_count;
      generalized_maxwell->read_material_parameters(numbranch, matids, solve);
      branchespotsum = generalized_maxwell->get_branchespotsum();
      branchtau = generalized_maxwell->get_branchtaus();
    }
  }

  if (generalized_maxwell_model_count != 1)
    FOUR_C_THROW(
        "Invalid VISCO_GeneralizedMaxwell setup in MAT_ViscoElastHyper (MAT {}, GP {}, ELE {}): "
        "expected exactly one VISCO_GeneralizedMaxwell summand but found {}.",
        visco_mat_id, gp, eleGID, generalized_maxwell_model_count);

  if (matids == nullptr)
    FOUR_C_THROW(
        "Failed to read MATIDS for VISCO_GeneralizedMaxwell in MAT_ViscoElastHyper (MAT {}, "
        "GP {}, ELE {}).",
        visco_mat_id, gp, eleGID);

  if (numbranch != static_cast<int>(matids->size()))
    FOUR_C_THROW(
        "Invalid VISCO_GeneralizedMaxwell state in MAT_ViscoElastHyper (MAT {}, GP {}, ELE {}): "
        "NUMBRANCH={} but MATIDS has size {}.",
        visco_mat_id, gp, eleGID, numbranch, matids->size());

  if (numbranch < 0 || branchespotsum.size() != static_cast<unsigned int>(numbranch) ||
      branchtau.size() != static_cast<unsigned int>(numbranch))
    FOUR_C_THROW(
        "Failed to initialize VISCO_GeneralizedMaxwell branches in MAT_ViscoElastHyper (MAT {}, "
        "GP {}, ELE {}). Expected {} branches, got {} branch definitions and {} branch "
        "relaxation times.",
        visco_mat_id, gp, eleGID, numbranch, branchespotsum.size(), branchtau.size());

  Core::LinAlg::Matrix<6, 6> cmatqbranch(Core::LinAlg::Initialization::zero);
  std::vector<Core::LinAlg::Matrix<6, 1>> S(numbranch);
  std::vector<Core::LinAlg::Matrix<6, 1>> Qbranch(numbranch);
  std::vector<Core::LinAlg::Matrix<6, 1>> S_n(numbranch);
  std::vector<Core::LinAlg::Matrix<6, 1>> Q_n(numbranch);

  // read history
  S_n = state_.branch_elastic_stress_last_at(gp);
  Q_n = state_.branch_stress_last_at(gp);

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

  // save switches
  SummandProperties branchProperties;

  /////////////////////////////////////////////////
  // Loop over all viscoelastic Maxwell branches //
  /////////////////////////////////////////////////
  for (int i = 0; i < numbranch; ++i)
  {
    // get parameter of each visco branch
    branchpotsum = branchespotsum[i];
    const double tau = branchtau.at(i);

    branchProperties.clear();
    elast_hyper_properties(branchpotsum, branchProperties);

    if (isovisco_)
      FOUR_C_THROW(
          "Unsupported branch formulation in VISCO_GeneralizedMaxwell for MAT_ViscoElastHyper "
          "(MAT {}, GP {}, ELE {}): isovisco branch response is not implemented.",
          visco_mat_id, gp, eleGID);
    if (branchProperties.anisoprinc)
      FOUR_C_THROW(
          "Unsupported branch formulation in VISCO_GeneralizedMaxwell for MAT_ViscoElastHyper "
          "(MAT {}, GP {}, ELE {}): anisoprinc branch response is not implemented.",
          visco_mat_id, gp, eleGID);
    if (branchProperties.anisomod)
      FOUR_C_THROW(
          "Unsupported branch formulation in VISCO_GeneralizedMaxwell for MAT_ViscoElastHyper "
          "(MAT {}, GP {}, ELE {}): anisomod branch response is not implemented.",
          visco_mat_id, gp, eleGID);

    Core::LinAlg::Matrix<6, 1> modrcg(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 1> prinv(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<3, 1> dPI(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<6, 1> ddPII(Core::LinAlg::Initialization::zero);

    Core::LinAlg::Matrix<6, 1> gl_stress;
    Core::LinAlg::Voigt::Strains::to_stress_like(*glstrain, gl_stress);

    Core::LinAlg::SymmetricTensor<double, 3, 3> gl =
        Core::LinAlg::make_symmetric_tensor_from_stress_like_voigt_matrix(gl_stress);

    Core::LinAlg::SymmetricTensor<double, 3, 3> C;
    evaluate_right_cauchy_green_strain_like_voigt(gl, C);

    Core::LinAlg::SymmetricTensor<double, 3, 3> iC = Core::LinAlg::inv(C);

    Core::LinAlg::Matrix<6, 1> C_strain(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Matrix<6, 1> iC_strain(Core::LinAlg::Initialization::zero);
    Core::LinAlg::Voigt::Stresses::to_strain_like(
        Core::LinAlg::make_stress_like_voigt_view(C), C_strain);
    Core::LinAlg::Voigt::Stresses::to_strain_like(
        Core::LinAlg::make_stress_like_voigt_view(iC), iC_strain);

    Core::LinAlg::Voigt::Strains::invariants_principal(prinv, C_strain);
    elast_hyper_evaluate_invariant_derivatives(
        prinv, dPI, ddPII, branchpotsum, branchProperties, gp, eleGID);

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
    if (solve == "OneStepTheta")
    {
      // initialize scalars
      double lambdascalar1(true);
      double lambdascalar2(true);
      double theta = 0.5;

      // get theta of global time integration scheme to use it here
      // if global time integration scheme is not ONESTEPTHETA, theta is by default = 0.5 (abirzle
      // 09/14)
      const auto dyntype = Teuchos::getIntegralValue<Inpar::Solid::DynamicType>(
          Global::Problem::instance()->structural_dynamic_params(), "DYNAMICTYPE");
      if (dyntype == Inpar::Solid::DynamicType::OneStepTheta)
        theta = Global::Problem::instance()
                    ->structural_dynamic_params()
                    .sublist("ONESTEPTHETA")
                    .get<double>("THETA");

      // get time algorithmic parameters
      // NOTE: dt can be zero (in restart of STI) for Generalized Maxwell model
      // there is no special treatment required. Adaptation for Kelvin-Voigt were necessary.
      // evaluate scalars to compute
      // Q^(n+1) = tau/(tau+theta*dt) [(tau-dt+theta*dt)/tau Q + beta(S^(n+1) - S^n)]
      lambdascalar1 = tau / (tau + theta * dt);
      lambdascalar2 = (tau - dt + theta * dt) / tau;

      // same branch update factor as in the one-branch Maxwell case
      deltascalar = lambdascalar1;

      // calculate artificial viscous stresses Q
      // Q_(n+1) = lambdascalar1*[lamdascalar2* Q_n + (Sa_(n+1) - Sa_n)]
      Qbranch.at(i).update(lambdascalar2, Q_n.at(i), 1.0);
      Qbranch.at(i).update(1.0, S.at(i), 1.0);
      Qbranch.at(i).update(-1.0, S_n.at(i), 1.0);
      Qbranch.at(i).scale(lambdascalar1);
    }
    else if (solve == "ExponentialTimeDiscretization")
    {
      // initialize scalars
      double xiscalar1(true);
      double xiscalar2(true);
      xiscalar1 = exp(-dt / tau);
      xiscalar2 = exp(-dt / (2 * tau));

      deltascalar = xiscalar2;

      // calculate artificial stresses Q
      // Q_(n+1) = xiscalar1* Q_n + xiscalar2*(Sa_(n+1) - Sa_n)
      Qbranch.at(i).update(xiscalar1, Q_n.at(i), 1.0);
      Qbranch.at(i).update(xiscalar2, S.at(i), 1.0);
      Qbranch.at(i).update(-xiscalar2, S_n.at(i), 1.0);
    }
    else
      FOUR_C_THROW(
          "Invalid SOLVE='{}' in VISCO_GeneralizedMaxwell for MAT_ViscoElastHyper (MAT {}, GP "
          "{}, ELE {}). Use OneStepTheta or ExponentialTimeDiscretization.",
          solve, visco_mat_id, gp, eleGID);

    // sum up branches
    Q.update(1.0, Qbranch.at(i), 1.0);
    cmatq.update(deltascalar, cmatqbranch, 1.0);

  }  // end for loop over branches

  // update history
  state_.set_branch_elastic_stress_current_at(gp, S);
  state_.set_branch_stress_current_at(gp, Qbranch);
}  // evaluate_visco_generalized_maxwell


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ViscoElastHyper::evaluate_visco_fsls(Core::LinAlg::Matrix<6, 1> stress,
    Core::LinAlg::Matrix<6, 6> cmat, Core::LinAlg::Matrix<6, 1>& Q,
    Core::LinAlg::Matrix<6, 6>& cmatq, const double dt, const int gp, const int eleGID)
{
  const int visco_mat_id = params_ != nullptr ? params_->id() : -1;

  // initialize parameters
  double tau = 0.0;
  double alpha = 0.0;
  double beta = 0.0;
  int fsls_model_count = 0;
  int fsls_summand_mat_id = -1;
  // string not used in visco_fsls
  std::string solve = "";

  // read material parameters of FSLS material
  for (unsigned int p = 0; p < potsum_.size(); ++p)
  {
    std::shared_ptr<Mat::Elastic::Fsls> fsls =
        std::dynamic_pointer_cast<Mat::Elastic::Fsls>(potsum_[p]);

    if (fsls != nullptr)
    {
      ++fsls_model_count;
      fsls_summand_mat_id = params_ != nullptr ? mat_id(p) : -1;
      fsls->read_material_parameters_visco(tau, beta, alpha, solve);
    }
  }

  if (fsls_model_count != 1)
    FOUR_C_THROW(
        "Invalid VISCO_FSLS setup in MAT_ViscoElastHyper (MAT {}, GP {}, ELE {}): expected "
        "exactly one VISCO_FSLS summand but found {}.",
        visco_mat_id, gp, eleGID, fsls_model_count);

  if (tau <= 0.0)
    FOUR_C_THROW(
        "Invalid TAU={} in VISCO_FSLS (MAT {}, referenced by MAT_ViscoElastHyper MAT {}, GP {}, "
        "ELE {}). TAU has to be positive.",
        tau, fsls_summand_mat_id, visco_mat_id, gp, eleGID);

  if (alpha < 0.0 || alpha >= 1.0)
    FOUR_C_THROW(
        "Invalid ALPHA={} in VISCO_FSLS (MAT {}, referenced by MAT_ViscoElastHyper MAT {}, GP "
        "{}, ELE {}). Expected 0 <= ALPHA < 1.",
        alpha, fsls_summand_mat_id, visco_mat_id, gp, eleGID);

  if (dt < 0.0)
    FOUR_C_THROW(
        "Invalid time step size dt={} in VISCO_FSLS evaluation for MAT_ViscoElastHyper (MAT {}, "
        "GP {}, ELE {}). Expected dt >= 0.",
        dt, visco_mat_id, gp, eleGID);

  if (!state_.has_fsls_history())
    FOUR_C_THROW("Missing FSLS history state in MAT_ViscoElastHyper (MAT {}, GP {}, ELE {}).",
        visco_mat_id, gp, eleGID);

  if (gp < 0 || gp >= state_.fsls_num_gauss_points())
    FOUR_C_THROW(
        "Invalid Gauss point index GP={} for FSLS history in MAT_ViscoElastHyper (MAT {}, ELE "
        "{}). History container size is {}.",
        gp, visco_mat_id, eleGID, state_.fsls_num_gauss_points());

  // read history of last time step at gp
  // -> Q_n and history size
  const int hs = state_.fsls_history_size_at(gp);  // history size
  if (hs <= 0)
    FOUR_C_THROW(
        "Invalid FSLS history size {} at GP {} in MAT_ViscoElastHyper (MAT {}, ELE {}). "
        "Expected at least one entry.",
        hs, gp, visco_mat_id, eleGID);

  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Q_n(state_.fsls_history_at(gp, hs - 1));


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

    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> Qj(state_.fsls_history_at(gp, hs - j));
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
  double dtalpha = std::pow(dt, alpha);
  double taualpha = std::pow(tau, alpha);
  double lambdascalar1 = dtalpha / (dtalpha + taualpha);
  double lambdascalar2 = -1. * taualpha / (dtalpha + taualpha);

  Q.update(lambdascalar1 * beta, stress, 0.);
  Q.update(lambdascalar2, Qq, 1.);


  // update history for next step
  state_.set_fsls_current_at(gp, Q);  // Q_n+1


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
