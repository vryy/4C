#include "4C_scatra_timint_cardiac_monodomain.hpp"

#include "4C_fem_nurbs_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_krylov_projector.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mat_list.hpp"
#include "4C_material_base.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ScaTra::TimIntCardiacMonodomain::TimIntCardiacMonodomain(Teuchos::RCP<Core::FE::Discretization> dis,
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(dis, solver, sctratimintparams, extraparams, output),
      // Initialization of electrophysiology variables
      activation_time_np_(Teuchos::null),
      activation_threshold_(0.0),
      nb_max_mat_int_state_vars_(0),
      material_internal_state_np_(Teuchos::null),
      material_internal_state_np_component_(Teuchos::null),
      nb_max_mat_ionic_currents_(0),
      material_ionic_currents_np_(Teuchos::null),
      material_ionic_currents_np_component_(Teuchos::null),
      ep_params_(params)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomain::setup()
{
  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  // Activation time at time n+1
  activation_time_np_ = Core::LinAlg::create_vector(*dofrowmap, true);
  activation_threshold_ = ep_params_->get<double>("ACTTHRES");
  // Assumes that maximum nb_max_mat_int_state_vars_ internal state variables will be written
  nb_max_mat_int_state_vars_ = ep_params_->get<int>(
      "WRITEMAXINTSTATE");  // number of maximal internal state variables to be postprocessed
  if (nb_max_mat_int_state_vars_)
  {
    material_internal_state_np_ = Teuchos::make_rcp<Core::LinAlg::MultiVector<double>>(
        *(discret_->element_row_map()), nb_max_mat_int_state_vars_, true);
    material_internal_state_np_component_ =
        Core::LinAlg::create_vector(*(discret_->element_row_map()), true);
  }
  // Assumes that maximum nb_max_mat_ionic_currents_ ionic_currents variables will be written
  nb_max_mat_ionic_currents_ = ep_params_->get<int>(
      "WRITEMAXIONICCURRENTS");  // number of maximal internal state variables to be postprocessed
  if (nb_max_mat_ionic_currents_)
  {
    material_ionic_currents_np_ = Teuchos::make_rcp<Core::LinAlg::MultiVector<double>>(
        *(discret_->element_row_map()), nb_max_mat_ionic_currents_, true);
    material_ionic_currents_np_component_ =
        Core::LinAlg::create_vector(*(discret_->element_row_map()), true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomain::collect_runtime_output_data()
{
  // call base class first
  ScaTraTimIntImpl::collect_runtime_output_data();

  // electrophysiology
  // Compute and write activation time
  if (activation_time_np_ != Teuchos::null)
  {
    for (int k = 0; k < phinp_->MyLength(); k++)
    {
      if ((*phinp_)[k] >= activation_threshold_ && (*activation_time_np_)[k] <= dta_ * 0.9)
        (*activation_time_np_)[k] = time_;
    }
    std::vector<std::optional<std::string>> context(num_dof_per_node(), "activation_time");
    visualization_writer().append_result_data_vector_with_context(
        *activation_time_np_, Core::IO::OutputEntity::dof, context);
  }

  // Recover internal state of the material (for electrophysiology)
  if (material_internal_state_np_ != Teuchos::null and nb_max_mat_int_state_vars_)
  {
    Teuchos::ParameterList params;
    Core::Utils::add_enum_class_to_parameter_list<ScaTra::Action>(
        "action", ScaTra::Action::get_material_internal_state, params);
    params.set<Teuchos::RCP<Core::LinAlg::MultiVector<double>>>("material_internal_state",
        material_internal_state_np_);  // Probably do it once at the beginning
    discret_->evaluate(params);
    material_internal_state_np_ =
        params.get<Teuchos::RCP<Core::LinAlg::MultiVector<double>>>("material_internal_state");

    for (int k = 0; k < material_internal_state_np_->NumVectors(); ++k)
    {
      std::ostringstream temp;
      temp << k + 1;
      material_internal_state_np_component_ =
          Teuchos::make_rcp<Core::LinAlg::Vector<double>>((*material_internal_state_np_)(k));

      visualization_writer().append_result_data_vector_with_context(
          *material_internal_state_np_component_, Core::IO::OutputEntity::element,
          {"mat_int_state_" + temp.str()});
    }
  }

  // Recover internal ionic currents of the material (for electrophysiology)
  if (material_ionic_currents_np_ != Teuchos::null and nb_max_mat_ionic_currents_)
  {
    Teuchos::ParameterList params;
    Core::Utils::add_enum_class_to_parameter_list<ScaTra::Action>(
        "action", ScaTra::Action::get_material_ionic_currents, params);
    params.set<Teuchos::RCP<Core::LinAlg::MultiVector<double>>>("material_ionic_currents",
        material_ionic_currents_np_);  // Probably do it once at the beginning
    discret_->evaluate(params);
    material_internal_state_np_ =
        params.get<Teuchos::RCP<Core::LinAlg::MultiVector<double>>>("material_ionic_currents");

    for (int k = 0; k < material_ionic_currents_np_->NumVectors(); ++k)
    {
      std::ostringstream temp;
      temp << k + 1;
      material_ionic_currents_np_component_ =
          Teuchos::make_rcp<Core::LinAlg::Vector<double>>((*material_ionic_currents_np_)(k));

      visualization_writer().append_result_data_vector_with_context(
          *material_ionic_currents_np_component_, Core::IO::OutputEntity::element,
          {"mat_ionic_currents_" + temp.str()});
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomain::element_material_time_update()
{
  // create the parameters for the discretization
  Teuchos::ParameterList p;
  // action for elements
  Core::Utils::add_enum_class_to_parameter_list<ScaTra::Action>(
      "action", ScaTra::Action::time_update_material, p);
  // further required parameter
  p.set<double>("time-step length", dta_);

  // set vector values needed by elements
  discret_->clear_state();
  discret_->set_state("phinp", phinp_);

  // go to elements
  discret_->evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  discret_->clear_state();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomain::set_element_specific_scatra_parameters(
    Teuchos::ParameterList& eleparams) const
{
  // safety check
  if (params_->get<bool>("SEMIIMPLICIT"))
  {
    if (Inpar::ScaTra::timeint_gen_alpha ==
        Teuchos::getIntegralValue<Inpar::ScaTra::TimeIntegrationScheme>(*params_, "TIMEINTEGR"))
    {
      if (params_->get<double>("ALPHA_M") < 1.0 or params_->get<double>("ALPHA_F") < 1.0)
        FOUR_C_THROW(
            "EP calculation with semiimplicit timestepping scheme only tested for gen-alpha with "
            "alpha_f = alpha_m = 1!");
    }
  }

  eleparams.set<bool>("semiimplicit", params_->get<bool>("SEMIIMPLICIT"));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomain::write_restart() const
{
  // Compute and write activation time
  if (activation_time_np_ != Teuchos::null)
  {
    for (int k = 0; k < phinp_->MyLength(); k++)
    {
      if ((*phinp_)[k] >= activation_threshold_ && (*activation_time_np_)[k] <= dta_ * 0.9)
        (*activation_time_np_)[k] = time_;
    }
    output_->write_vector("activation_time_np", activation_time_np_);
  }
}

FOUR_C_NAMESPACE_CLOSE
