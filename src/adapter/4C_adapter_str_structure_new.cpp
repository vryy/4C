/*-----------------------------------------------------------*/
/*! \file

\brief Adapter for the new structural time integration framework.


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_adapter_str_structure_new.hpp"

#include "4C_adapter_str_constr_merged.hpp"
#include "4C_adapter_str_fbiwrapper.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_adapter_str_fsiwrapper_immersed.hpp"
#include "4C_adapter_str_pasiwrapper.hpp"
#include "4C_adapter_str_redairway.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_adapter_str_timeada.hpp"
#include "4C_adapter_str_timeloop.hpp"
#include "4C_adapter_str_wrapper.hpp"
#include "4C_beam3_kirchhoff.hpp"
#include "4C_beam3_reissner.hpp"
#include "4C_binstrategy.hpp"
#include "4C_comm_utils.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_beam_to_solid.hpp"
#include "4C_inpar_beamcontact.hpp"
#include "4C_inpar_beaminteraction.hpp"
#include "4C_inpar_contact.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_inpar_poroelast.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_rigidsphere.hpp"
#include "4C_shell7p_ele.hpp"
#include "4C_so3_hex8fbar.hpp"
#include "4C_so3_plast_ssn_eletypes.hpp"
#include "4C_so3_plast_ssn_sosh18.hpp"
#include "4C_so3_plast_ssn_sosh8.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_group_prepostoperator.hpp"
#include "4C_structure_new_model_evaluator_manager.hpp"
#include "4C_structure_new_solver_factory.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_factory.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Adapter::StructureBaseAlgorithmNew::StructureBaseAlgorithmNew()
    : str_wrapper_(Teuchos::null),
      prbdyn_(Teuchos::null),
      sdyn_(Teuchos::null),
      isinit_(false),
      issetup_(false)
{
  // empty
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::init(const Teuchos::ParameterList& prbdyn,
    Teuchos::ParameterList& sdyn, Teuchos::RCP<Core::FE::Discretization> actdis)
{
  issetup_ = false;

  // save a pointer of the input variables as class variables
  prbdyn_ = Teuchos::rcpFromRef<const Teuchos::ParameterList>(prbdyn);
  sdyn_ = Teuchos::rcpFromRef<Teuchos::ParameterList>(sdyn);
  actdis_ = actdis;

  isinit_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::setup()
{
  if (not is_init()) FOUR_C_THROW("You have to call init() first!");

  // major switch to different time integrators
  switch (Teuchos::getIntegralValue<Inpar::Solid::DynamicType>(*sdyn_, "DYNAMICTYP"))
  {
    case Inpar::Solid::dyna_statics:
    case Inpar::Solid::dyna_genalpha:
    case Inpar::Solid::dyna_genalpha_liegroup:
    case Inpar::Solid::dyna_onesteptheta:
    case Inpar::Solid::dyna_expleuler:
    case Inpar::Solid::dyna_centrdiff:
    case Inpar::Solid::dyna_ab2:
    case Inpar::Solid::dyna_ab4:
      setup_tim_int();  // <-- here is the show
      break;
    default:
      FOUR_C_THROW("Unknown time integration scheme '%s'",
          Teuchos::getStringValue<Inpar::Solid::DynamicType>(*sdyn_, "DYNAMICTYP").c_str());
      break;
  }

  issetup_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::register_model_evaluator(
    const std::string name, Teuchos::RCP<Solid::ModelEvaluator::Generic> me)
{
  // safety checks
  if (not is_init())
    FOUR_C_THROW("init(...) must be called before register_model_evaluator(...) !");
  if (is_setup()) FOUR_C_THROW("register_model_evaluator(...) must be called before setup() !");

  // set RCP ptr to model evaluator in problem dynamic parameter list
  const_cast<Teuchos::ParameterList&>(*prbdyn_).set<Teuchos::RCP<Solid::ModelEvaluator::Generic>>(
      name, me);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::setup_tim_int()
{
  if (not is_init()) FOUR_C_THROW("You have to call init() first!");

  // get the problem instance
  Global::Problem* problem = Global::Problem::instance();
  // get the restart step
  const int restart = problem->restart();

  // ---------------------------------------------------------------------------
  // Define, initialize and start the timer
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Teuchos::Time> t =
      Teuchos::TimeMonitor::getNewTimer("Adapter::StructureTimIntBaseAlgorithm::SetupStructure");
  Teuchos::TimeMonitor monitor(*t);

  // ---------------------------------------------------------------------------
  // Here we read the discretization at the current
  // time step from restart files
  // ---------------------------------------------------------------------------
  if (actdis_->get_condition("PointCoupling") != nullptr)
  {
    std::vector<Teuchos::RCP<Core::FE::Discretization>> actdis_vec(1, actdis_);
    Teuchos::ParameterList binning_params = Global::Problem::instance()->binning_strategy_params();
    Core::Utils::add_enum_class_to_parameter_list<Core::FE::ShapeFunctionType>(
        "spatial_approximation_type", Global::Problem::instance()->spatial_approximation_type(),
        binning_params);
    actdis_vec[0]->fill_complete(false, false, false);

    // Different types of structural elements may be present, so we need to help the binning
    // strategy understand their different shapes by providing the correct points to compute
    // bounding boxes.
    auto correct_node = [](const Core::Nodes::Node& node) -> decltype(auto)
    {
      const Core::Elements::Element* element = node.elements()[0];
      const auto* beamelement = dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(element);
      if (beamelement != nullptr && !beamelement->is_centerline_node(node))
        return *element->nodes()[0];
      else
        return node;
    };

    auto determine_relevant_points = [correct_node](const Core::FE::Discretization& discret,
                                         const Core::Elements::Element& ele,
                                         Teuchos::RCP<const Core::LinAlg::Vector<double>> disnp)
        -> std::vector<std::array<double, 3>>
    {
      if (dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(&ele))
      {
        return Core::Binstrategy::DefaultRelevantPoints{
            .correct_node = correct_node,
        }(discret, ele, disnp);
      }
      else if (ele.element_type() == Discret::ELEMENTS::RigidsphereType::instance())
      {
        double currpos[3] = {0.0, 0.0, 0.0};
        Core::Binstrategy::Utils::get_current_node_pos(discret, ele.nodes()[0], disnp, currpos);
        const double radius = dynamic_cast<const Discret::ELEMENTS::Rigidsphere&>(ele).radius();
        return {{currpos[0] - radius, currpos[1] - radius, currpos[2] - radius},
            {currpos[0] + radius, currpos[1] + radius, currpos[2] + radius}};
      }
      else
        return Core::Binstrategy::DefaultRelevantPoints{}(discret, ele, disnp);
    };

    Core::Rebalance::rebalance_discretizations_by_binning(binning_params,
        Global::Problem::instance()->output_control_file(), actdis_vec, correct_node,
        determine_relevant_points, true);
  }
  else if (not actdis_->filled() || not actdis_->have_dofs())
  {
    actdis_->fill_complete();
  }

  // ---------------------------------------------------------------------------
  // Setup a model type set by checking
  // the different conditions
  // ---------------------------------------------------------------------------
  // define and initial with default value
  Teuchos::RCP<std::set<enum Inpar::Solid::ModelType>> modeltypes =
      Teuchos::make_rcp<std::set<enum Inpar::Solid::ModelType>>();
  modeltypes->insert(Inpar::Solid::model_structure);
  set_model_types(*modeltypes);

  // ---------------------------------------------------------------------------
  // Setup a element technology set by checking
  // the elements of the discretization
  // ---------------------------------------------------------------------------
  Teuchos::RCP<std::set<enum Inpar::Solid::EleTech>> eletechs =
      Teuchos::make_rcp<std::set<enum Inpar::Solid::EleTech>>();
  detect_element_technologies(*eletechs);

  // ---------------------------------------------------------------------------
  // Setup the parameter lists for structural
  // time integration
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> ioflags =
      Teuchos::make_rcp<Teuchos::ParameterList>(problem->io_params());
  Teuchos::RCP<Teuchos::ParameterList> time_adaptivity_params =
      Teuchos::make_rcp<Teuchos::ParameterList>(sdyn_->sublist("TIMEADAPTIVITY"));
  Teuchos::RCP<Teuchos::ParameterList> xparams = Teuchos::make_rcp<Teuchos::ParameterList>();
  set_params(*ioflags, *xparams, *time_adaptivity_params);

  // ---------------------------------------------------------------------------
  // Setup and create model specific linear solvers
  // ---------------------------------------------------------------------------
  Teuchos::RCP<std::map<enum Inpar::Solid::ModelType, Teuchos::RCP<Core::LinAlg::Solver>>>
      linsolvers = Solid::SOLVER::build_lin_solvers(*modeltypes, *sdyn_, *actdis_);

  // ---------------------------------------------------------------------------
  // Checks in case of multi-scale simulations
  // ---------------------------------------------------------------------------
  {
    // make sure we IMR-like generalised-alpha requested for multi-scale
    // simulations
    Teuchos::RCP<Mat::PAR::Bundle> materials = problem->materials();
    for (const auto& [_, mat] : materials->map())
    {
      if (mat->type() == Core::Materials::m_struct_multiscale)
      {
        if (Teuchos::getIntegralValue<Inpar::Solid::DynamicType>(*sdyn_, "DYNAMICTYP") !=
            Inpar::Solid::dyna_genalpha)
          FOUR_C_THROW("In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha");
        else if (Teuchos::getIntegralValue<Inpar::Solid::MidAverageEnum>(
                     sdyn_->sublist("GENALPHA"), "GENAVG") != Inpar::Solid::midavg_trlike)
          FOUR_C_THROW(
              "In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha with GENAVG=TrLike");
        break;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Create context for output and restart
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = actdis_->writer();
  if (ioflags->get<bool>("OUTPUT_BIN"))
  {
    output->write_mesh(0, 0.0);
  }

  // ---------------------------------------------------------------------------
  // initialize/setup the input/output data container
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Solid::TimeInt::BaseDataIO> dataio = Teuchos::make_rcp<Solid::TimeInt::BaseDataIO>();
  dataio->init(*ioflags, *sdyn_, *xparams, output);
  dataio->setup();

  // ---------------------------------------------------------------------------
  // initialize/setup the structural dynamics data
  // container
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Solid::TimeInt::BaseDataSDyn> datasdyn = Solid::TimeInt::build_data_sdyn(*sdyn_);
  datasdyn->init(actdis_, *sdyn_, *xparams, modeltypes, eletechs, linsolvers);
  datasdyn->setup();

  // ---------------------------------------------------------------------------
  // initialize/setup the global state data container
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState> dataglobalstate = Teuchos::null;
  set_global_state(dataglobalstate, datasdyn);

  // ---------------------------------------------------------------------------
  // in case of non-additive rotation (pseudo-)vector DOFs:
  // ---------------------------------------------------------------------------
  if (eletechs->find(Inpar::Solid::EleTech::rotvec) != eletechs->end())
  {
    // -------------------------------------------------------------------------
    // set the RotVecUpdater as new precomputeX operator for the nox nln group
    // -------------------------------------------------------------------------
    Teuchos::ParameterList& p_grp_opt = datasdyn->get_nox_params().sublist("Group Options");
    // Get the current map. If there is no map, return a new empty one. (reference)
    NOX::Nln::GROUP::PrePostOperator::Map& prepostgroup_map =
        NOX::Nln::GROUP::PrePostOp::get_map(p_grp_opt);
    // create the new rotation vector update pre/post operator
    Teuchos::RCP<NOX::Nln::Abstract::PrePostOperator> prepostrotvec_ptr =
        Teuchos::make_rcp<NOX::Nln::GROUP::PrePostOp::TimeInt::RotVecUpdater>(dataglobalstate);
    // insert/replace the old pointer in the map
    prepostgroup_map[NOX::Nln::GROUP::prepost_rotvecupdate] = prepostrotvec_ptr;
  }

  // ---------------------------------------------------------------------------
  // Build time integrator
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Solid::TimeInt::Base> ti_strategy = Teuchos::null;
  set_time_integration_strategy(ti_strategy, dataio, datasdyn, dataglobalstate, restart);


  // ---------------------------------------------------------------------------
  // Create wrapper for the time integration strategy
  // ---------------------------------------------------------------------------
  set_structure_wrapper(*ioflags, *sdyn_, *xparams, *time_adaptivity_params, ti_strategy);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::set_model_types(
    std::set<enum Inpar::Solid::ModelType>& modeltypes) const
{
  if (not is_init()) FOUR_C_THROW("You have to call init() first!");
  // ---------------------------------------------------------------------------
  // check for meshtying and contact conditions
  // ---------------------------------------------------------------------------
  // --- contact conditions
  std::vector<Core::Conditions::Condition*> ccond(0);
  actdis_->get_condition("Contact", ccond);
  if (ccond.size())
  {
    // what's the current problem type?
    Core::ProblemType probtype = Global::Problem::instance()->get_problem_type();
    // ToDo: once the new structural time integration can handle
    //       condensed contact formulations, the model_evaluator
    //       can have its contact model. For now, the TSI Lagrange
    //       strategy resides in the TSI algorithm.
    if (probtype == Core::ProblemType::tsi)
    {
      const Teuchos::ParameterList& contact = Global::Problem::instance()->contact_dynamic_params();
      if (Teuchos::getIntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          Inpar::CONTACT::solution_nitsche)
        modeltypes.insert(Inpar::Solid::model_contact);
    }
    else
      modeltypes.insert(Inpar::Solid::model_contact);
  }
  // --- meshtying conditions
  std::vector<Core::Conditions::Condition*> mtcond(0);
  actdis_->get_condition("Mortar", mtcond);
  if (mtcond.size()) modeltypes.insert(Inpar::Solid::model_meshtying);

  // check for 0D cardiovascular conditions
  // ---------------------------------------------------------------------------
  std::vector<Core::Conditions::Condition*> cardiovasc0dcond_4elementwindkessel(0);
  std::vector<Core::Conditions::Condition*> cardiovasc0dcond_arterialproxdist(0);
  std::vector<Core::Conditions::Condition*> cardiovasc0dcond_syspulcirculation(0);
  std::vector<Core::Conditions::Condition*> cardiovascrespir0dcond_syspulperiphcirculation(0);
  actdis_->get_condition(
      "Cardiovascular0D4ElementWindkesselStructureCond", cardiovasc0dcond_4elementwindkessel);
  actdis_->get_condition(
      "Cardiovascular0DArterialProxDistStructureCond", cardiovasc0dcond_arterialproxdist);
  actdis_->get_condition("Cardiovascular0DSysPulCirculationStructureCond",
      cardiovascrespir0dcond_syspulperiphcirculation);
  actdis_->get_condition("CardiovascularRespiratory0DSysPulPeriphCirculationStructureCond",
      cardiovasc0dcond_syspulcirculation);
  if (cardiovasc0dcond_4elementwindkessel.size() or cardiovasc0dcond_arterialproxdist.size() or
      cardiovasc0dcond_syspulcirculation.size() or
      cardiovascrespir0dcond_syspulperiphcirculation.size())
    modeltypes.insert(Inpar::Solid::model_cardiovascular0d);

  // ---------------------------------------------------------------------------
  // check for constraint conditions
  // ---------------------------------------------------------------------------
  bool have_lag_constraint = false;
  bool have_pen_constraint = false;
  // --- enforcement by Lagrange multiplier
  std::vector<Core::Conditions::Condition*> lagcond_volconstr3d(0);
  std::vector<Core::Conditions::Condition*> lagcond_areaconstr3d(0);
  std::vector<Core::Conditions::Condition*> lagcond_areaconstr2d(0);
  std::vector<Core::Conditions::Condition*> lagcond_mpconline2d(0);
  std::vector<Core::Conditions::Condition*> lagcond_mpconplane3d(0);
  std::vector<Core::Conditions::Condition*> lagcond_mpcnormcomp3d(0);
  actdis_->get_condition("VolumeConstraint_3D", lagcond_volconstr3d);
  actdis_->get_condition("AreaConstraint_3D", lagcond_areaconstr3d);
  actdis_->get_condition("AreaConstraint_2D", lagcond_areaconstr2d);
  actdis_->get_condition("MPC_NodeOnLine_2D", lagcond_mpconline2d);
  actdis_->get_condition("MPC_NodeOnPlane_3D", lagcond_mpconplane3d);
  actdis_->get_condition("MPC_NormalComponent_3D", lagcond_mpcnormcomp3d);
  if (lagcond_volconstr3d.size() or lagcond_areaconstr3d.size() or lagcond_areaconstr2d.size() or
      lagcond_mpconline2d.size() or lagcond_mpconplane3d.size() or lagcond_mpcnormcomp3d.size())
    have_lag_constraint = true;
  // --- enforcement by penalty law
  std::vector<Core::Conditions::Condition*> pencond_volconstr3d(0);
  std::vector<Core::Conditions::Condition*> pencond_areaconstr3d(0);
  std::vector<Core::Conditions::Condition*> pencond_mpcnormcomp3d(0);
  actdis_->get_condition("VolumeConstraint_3D_Pen", pencond_volconstr3d);
  actdis_->get_condition("AreaConstraint_3D_Pen", pencond_areaconstr3d);
  actdis_->get_condition("MPC_NormalComponent_3D_Pen", pencond_mpcnormcomp3d);
  if (pencond_volconstr3d.size() or pencond_areaconstr3d.size() or pencond_mpcnormcomp3d.size())
    have_pen_constraint = true;
  if (have_lag_constraint or have_pen_constraint)
    modeltypes.insert(Inpar::Solid::model_lag_pen_constraint);

  // ---------------------------------------------------------------------------
  // check for spring dashpot conditions
  // ---------------------------------------------------------------------------
  std::vector<Core::Conditions::Condition*> sdp_cond(0);
  actdis_->get_condition("RobinSpringDashpot", sdp_cond);
  if (sdp_cond.size()) modeltypes.insert(Inpar::Solid::model_springdashpot);
  // ---------------------------------------------------------------------------
  // check for coupled problems
  // ---------------------------------------------------------------------------
  // get the problem instance
  Global::Problem* problem = Global::Problem::instance();
  // what's the current problem type?
  Core::ProblemType probtype = problem->get_problem_type();
  switch (probtype)
  {
    case Core::ProblemType::fsi:
    case Core::ProblemType::immersed_fsi:
    case Core::ProblemType::fbi:
    case Core::ProblemType::fsi_redmodels:
    case Core::ProblemType::gas_fsi:
    case Core::ProblemType::biofilm_fsi:
    case Core::ProblemType::thermo_fsi:
    case Core::ProblemType::fsi_xfem:
    case Core::ProblemType::pasi:
    case Core::ProblemType::ssi:
    case Core::ProblemType::ssti:
    {
      if (prbdyn_->INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<Solid::ModelEvaluator::Generic>>(
              "Partitioned Coupling Model"))
      {
        if (prbdyn_->INVALID_TEMPLATE_QUALIFIER
                isType<Teuchos::RCP<Solid::ModelEvaluator::Generic>>("Monolithic Coupling Model"))
          FOUR_C_THROW("Cannot have both partitioned and monolithic coupling at the same time!");
        const auto coupling_model_ptr =
            prbdyn_->INVALID_TEMPLATE_QUALIFIER get<Teuchos::RCP<Solid::ModelEvaluator::Generic>>(
                "Partitioned Coupling Model");
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW(
              "The partitioned coupling model pointer is not allowed to be Teuchos::null!");
        // set the model type
        modeltypes.insert(Inpar::Solid::model_partitioned_coupling);
        // copy the coupling model object pointer into the (temporal) sdyn parameter list
        sdyn_->set<Teuchos::RCP<Solid::ModelEvaluator::Generic>>(
            "Partitioned Coupling Model", coupling_model_ptr);
      }

      else if (prbdyn_->INVALID_TEMPLATE_QUALIFIER
                   isType<Teuchos::RCP<Solid::ModelEvaluator::Generic>>(
                       "Monolithic Coupling Model"))
      {
        const auto coupling_model_ptr =
            prbdyn_->INVALID_TEMPLATE_QUALIFIER get<Teuchos::RCP<Solid::ModelEvaluator::Generic>>(
                "Monolithic Coupling Model");
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW("The monolithic coupling model pointer is not allowed to be Teuchos::null!");
        // set the model type
        modeltypes.insert(Inpar::Solid::model_monolithic_coupling);
        // copy the coupling model object pointer into the (temporal) sdyn parameter list
        sdyn_->set<Teuchos::RCP<Solid::ModelEvaluator::Generic>>(
            "Monolithic Coupling Model", coupling_model_ptr);
      }

      else if (prbdyn_->INVALID_TEMPLATE_QUALIFIER
                   isType<Teuchos::RCP<Solid::ModelEvaluator::Generic>>("Basic Coupling Model"))
      {
        const auto coupling_model_ptr =
            prbdyn_->INVALID_TEMPLATE_QUALIFIER get<Teuchos::RCP<Solid::ModelEvaluator::Generic>>(
                "Basic Coupling Model");
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW("The basic coupling model pointer is not allowed to be Teuchos::null!");
        // set the model type
        modeltypes.insert(Inpar::Solid::model_basic_coupling);
        // copy the coupling model object pointer into the (temporal) sdyn parameter list
        sdyn_->set<Teuchos::RCP<Solid::ModelEvaluator::Generic>>(
            "Basic Coupling Model", coupling_model_ptr);
      }
      break;
    }
    default:
      // do nothing
      break;
  }  // switch (probtype)

  // ---------------------------------------------------------------------------
  // check for beam interactions (either contact or potential-based)
  // ---------------------------------------------------------------------------
  // get beam contact strategy
  const Teuchos::ParameterList& beamcontact = Global::Problem::instance()->beam_contact_params();
  auto strategy =
      Teuchos::getIntegralValue<Inpar::BEAMCONTACT::Strategy>(beamcontact, "BEAMS_STRATEGY");

  auto modelevaluator =
      Teuchos::getIntegralValue<Inpar::BEAMCONTACT::Modelevaluator>(beamcontact, "MODELEVALUATOR");

  // conditions for potential-based beam interaction
  std::vector<Core::Conditions::Condition*> beampotconditions(0);
  actdis_->get_condition("BeamPotentialLineCharge", beampotconditions);

  // conditions for beam penalty point coupling
  std::vector<Core::Conditions::Condition*> beampenaltycouplingconditions(0);
  actdis_->get_condition("PenaltyPointCouplingCondition", beampenaltycouplingconditions);


  if (strategy != Inpar::BEAMCONTACT::bstr_none and modelevaluator == Inpar::BEAMCONTACT::bstr_old)
    modeltypes.insert(Inpar::Solid::model_beam_interaction_old);

  // ---------------------------------------------------------------------------
  // check for brownian dynamics
  // ---------------------------------------------------------------------------
  if (Global::Problem::instance()->brownian_dynamics_params().get<bool>("BROWNDYNPROB"))
    modeltypes.insert(Inpar::Solid::model_browniandyn);

  // ---------------------------------------------------------------------------
  // check for beam interaction
  // ---------------------------------------------------------------------------
  if (Global::Problem::instance()
          ->beam_interaction_params()
          .sublist("CROSSLINKING")
          .get<bool>("CROSSLINKER") or
      Global::Problem::instance()
          ->beam_interaction_params()
          .sublist("SPHERE BEAM LINK")
          .get<bool>("SPHEREBEAMLINKING") or
      Teuchos::getIntegralValue<Inpar::BEAMINTERACTION::Strategy>(
          Global::Problem::instance()->beam_interaction_params().sublist("BEAM TO BEAM CONTACT"),
          "STRATEGY") != Inpar::BEAMINTERACTION::bstr_none or
      Teuchos::getIntegralValue<Inpar::BEAMINTERACTION::Strategy>(
          Global::Problem::instance()->beam_interaction_params().sublist("BEAM TO SPHERE CONTACT"),
          "STRATEGY") != Inpar::BEAMINTERACTION::bstr_none or
      Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          Global::Problem::instance()->beam_interaction_params().sublist(
              "BEAM TO SOLID VOLUME MESHTYING"),
          "CONTACT_DISCRETIZATION") != Inpar::BeamToSolid::BeamToSolidContactDiscretization::none or
      Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          Global::Problem::instance()->beam_interaction_params().sublist(
              "BEAM TO SOLID SURFACE MESHTYING"),
          "CONTACT_DISCRETIZATION") != Inpar::BeamToSolid::BeamToSolidContactDiscretization::none or
      Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          Global::Problem::instance()->beam_interaction_params().sublist(
              "BEAM TO SOLID SURFACE CONTACT"),
          "CONTACT_DISCRETIZATION") != Inpar::BeamToSolid::BeamToSolidContactDiscretization::none or
      beampotconditions.size() > 0 or beampenaltycouplingconditions.size() > 0)
  {
    modeltypes.insert(Inpar::Solid::model_beaminteraction);
  }

  // ---------------------------------------------------------------------------
  // check for constraints
  // ---------------------------------------------------------------------------
  std::vector<Teuchos::RCP<Core::Conditions::Condition>> linePeriodicRve, surfPeriodicRve,
      pointLinearCoupledEquation, embeddedMeshConditions;
  actdis_->get_condition("LinePeriodicRve", linePeriodicRve);
  actdis_->get_condition("SurfacePeriodicRve", surfPeriodicRve);
  actdis_->get_condition("PointLinearCoupledEquation", pointLinearCoupledEquation);
  actdis_->get_condition("EmbeddedMeshSolidSurfCoupling", embeddedMeshConditions);

  if (linePeriodicRve.size() > 0 || surfPeriodicRve.size() > 0 ||
      pointLinearCoupledEquation.size() > 0 || embeddedMeshConditions.size() > 0)
    modeltypes.insert(Inpar::Solid::model_constraints);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::detect_element_technologies(
    std::set<enum Inpar::Solid::EleTech>& eletechs) const
{
  int isplasticity_local = 0;
  int isplasticity_global = 0;

  int iseas_local = 0;
  int iseas_global = 0;

  int isfbar_local = 0;
  int isfbar_global = 0;

  int isrotvec_local = 0;
  int isrotvec_global = 0;

  for (int i = 0; i < actdis_->num_my_row_elements(); ++i)
  {
    Core::Elements::Element* actele = actdis_->l_row_element(i);
    // Detect plasticity -------------------------------------------------------
    if (actele->element_type() == Discret::ELEMENTS::SoHex8PlastType::instance() or
        actele->element_type() == Discret::ELEMENTS::SoHex27PlastType::instance() or
        actele->element_type() == Discret::ELEMENTS::SoSh8PlastType::instance() or
        actele->element_type() == Discret::ELEMENTS::SoHex18PlastType::instance() or
        actele->element_type() == Discret::ELEMENTS::SoSh18PlastType::instance())
    {
      if (actele->material()->material_type() == Core::Materials::m_plelasthyper)
        isplasticity_local = true;
    }

    // Detect EAS --------------------------------------------------------------
    Discret::ELEMENTS::SoBase* so_base_ele = dynamic_cast<Discret::ELEMENTS::SoBase*>(actele);
    if (so_base_ele != nullptr)
    {
      if (so_base_ele->have_eas()) iseas_local = 1;
    }

    Discret::ELEMENTS::Shell7p* shell7p = dynamic_cast<Discret::ELEMENTS::Shell7p*>(actele);
    if (shell7p)
      if (shell7p->get_ele_tech().find(Inpar::Solid::EleTech::eas) != shell7p->get_ele_tech().end())
        iseas_local = 1;

    Discret::ELEMENTS::Solid* solid = dynamic_cast<Discret::ELEMENTS::Solid*>(actele);
    if (solid != nullptr)
      if (solid->have_eas()) iseas_local = 1;

    // Detect fbar
    Discret::ELEMENTS::SoHex8fbar* so_hex8fbar_ele =
        dynamic_cast<Discret::ELEMENTS::SoHex8fbar*>(actele);
    if (so_hex8fbar_ele != nullptr) isfbar_local = 1;

    // Detect non-additive rotation-vector DOFs --------------------------------
    if (actele->element_type() == Discret::ELEMENTS::Beam3rType::instance() or
        actele->element_type() == Discret::ELEMENTS::Beam3kType::instance())
    {
      isrotvec_local = true;
      break;
    }
  }

  // plasticity - sum over all processors
  actdis_->get_comm().SumAll(&isplasticity_local, &isplasticity_global, 1);
  if (isplasticity_global > 0) eletechs.insert(Inpar::Solid::EleTech::plasticity);

  // eas - sum over all processors
  actdis_->get_comm().SumAll(&iseas_local, &iseas_global, 1);
  if (iseas_global > 0) eletechs.insert(Inpar::Solid::EleTech::eas);

  // fbar - sum over all processors
  actdis_->get_comm().SumAll(&isfbar_local, &isfbar_global, 1);
  if (isfbar_global > 0) eletechs.insert(Inpar::Solid::EleTech::fbar);

  // rotation vector DOFs - sum over all processors
  actdis_->get_comm().SumAll(&isrotvec_local, &isrotvec_global, 1);
  if (isrotvec_global > 0) eletechs.insert(Inpar::Solid::EleTech::rotvec);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::set_params(Teuchos::ParameterList& ioflags,
    Teuchos::ParameterList& xparams, Teuchos::ParameterList& time_adaptivity_params)
{
  // get the problem instance and the problem type
  Global::Problem* problem = Global::Problem::instance();
  Core::ProblemType probtype = problem->get_problem_type();

  // ---------------------------------------------------------------------------
  // get input parameter lists and copy them,
  // because a few parameters are overwritten
  // ---------------------------------------------------------------------------
  // nox parameter list
  Teuchos::ParameterList snox(problem->structural_nox_params());
  Teuchos::ParameterList& nox = xparams.sublist("NOX");
  nox = snox;

  /* overrule certain parameters
   *
   * These parameters are overwritten by the parameters of the current
   * problem type. This is done only once temporally, since we need the
   * structural dynamics parameter-list only for the setup routines of the
   * different data containers.                                         */
  sdyn_->set<double>("TIMESTEP", prbdyn_->get<double>("TIMESTEP"));
  sdyn_->set<double>("MAXTIME", prbdyn_->get<double>("MAXTIME"));

  // overrule certain parameters
  sdyn_->set<int>("NUMSTEP", prbdyn_->get<int>("NUMSTEP"));
  sdyn_->set<int>("RESTARTEVRY", prbdyn_->get<int>("RESTARTEVRY"));
  sdyn_->set<int>("RESULTSEVRY", prbdyn_->get<int>("RESULTSEVRY"));

  // Check if for chosen Rayleigh damping the regarding parameters are given explicitly in the .dat
  // file
  if (Teuchos::getIntegralValue<Inpar::Solid::DampKind>(*sdyn_, "DAMPING") ==
      Inpar::Solid::damp_rayleigh)
  {
    if (sdyn_->get<double>("K_DAMP") < 0.0)
    {
      FOUR_C_THROW("Rayleigh damping parameter K_DAMP not explicitly given.");
    }
    if (sdyn_->get<double>("M_DAMP") < 0.0)
    {
      FOUR_C_THROW("Rayleigh damping parameter M_DAMP not explicitly given.");
    }
  }

  // ---------------------------------------------------------------------------
  // ---------------------------------------------------------------------------
  /* Overwrite certain parameters in STRUCTURAL DYNAMIC/TIMADAPTIVITY by those from
   * FSI DYNAMIC/TIMEADAPTIVITY
   *
   * In case, that the structure field is part of an FSI simulation with time step
   * size adaptivity based on structure field error estimation, we have to provide
   * the following algorithmic control parameters:
   *
   * - ADAPTSTEPMAX
   * - STEPSIZEMAX
   * - STEPSIZEMIN
   * - SIZERATIOMAX
   * - SIZERATIOMIN
   * - SIZERATIOSCALE
   *
   * They are specified by the FSI algorithm since they have to be the same for
   * the structure and fluid field. Hence, we overwrite the corresponding
   * parameters in the structural parameter list in order to avoid redundant
   * parameter specification in the input file.
   *
   * Note: This is really ugly, but currently the only way to avoid that the user
   * has to specify these parameters twice in the input file.
   *
   * ToDO: Find something nicer here!
   *
   * \author mayr.mt \date 12/2013
   */
  // ---------------------------------------------------------------------------
  switch (probtype)
  {
    case Core::ProblemType::fsi:
    case Core::ProblemType::fsi_redmodels:
    {
      const Teuchos::ParameterList& fsidyn = problem->fsi_dynamic_params();
      const Teuchos::ParameterList& fsiada = fsidyn.sublist("TIMEADAPTIVITY");
      if (fsiada.get<bool>("TIMEADAPTON"))
      {
        // overrule time step size adaptivity control parameters
        if (time_adaptivity_params.get<Inpar::Solid::TimAdaKind>("KIND") !=
            Inpar::Solid::timada_kind_none)
        {
          time_adaptivity_params.set<int>("ADAPTSTEPMAX", fsiada.get<int>("ADAPTSTEPMAX"));
          time_adaptivity_params.set<double>("STEPSIZEMAX", fsiada.get<double>("DTMAX"));
          time_adaptivity_params.set<double>("STEPSIZEMIN", fsiada.get<double>("DTMIN"));
          time_adaptivity_params.set<double>("SIZERATIOMAX", fsiada.get<double>("SIZERATIOMAX"));
          time_adaptivity_params.set<double>("SIZERATIOMIN", fsiada.get<double>("SIZERATIOMIN"));
          time_adaptivity_params.set<double>("SIZERATIOSCALE", fsiada.get<double>("SAFETYFACTOR"));

          if (actdis_->get_comm().MyPID() == 0)
          {
            Core::IO::cout
                << "*** Due to FSI time step size adaptivity with structure based error "
                   "estimation,\n"
                   "algorithmic control parameters in STRUCTURAL DYNAMIC/TIMEADAPTIVITY have been\n"
                   "overwritten by those from FSI DYNAMIC/TIMEADAPTIVITY."
                << Core::IO::endl
                << Core::IO::endl;
          }
        }
      }
      break;
    }
    default:
    {
      // do nothing
      break;
    }
  }
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::set_global_state(
    Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& dataglobalstate,
    const Teuchos::RCP<const Solid::TimeInt::BaseDataSDyn>& datasdyn)
{
  dataglobalstate = Solid::TimeInt::build_data_global_state();
  dataglobalstate->init(actdis_, *sdyn_, datasdyn);
  dataglobalstate->setup();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::set_time_integration_strategy(
    Teuchos::RCP<Solid::TimeInt::Base>& ti_strategy,
    const Teuchos::RCP<Solid::TimeInt::BaseDataIO>& dataio,
    const Teuchos::RCP<Solid::TimeInt::BaseDataSDyn>& datasdyn,
    const Teuchos::RCP<Solid::TimeInt::BaseDataGlobalState>& dataglobalstate, const int& restart)
{
  ti_strategy = Solid::TimeInt::build_strategy(*sdyn_);
  ti_strategy->init(dataio, datasdyn, dataglobalstate);

  /* In the restart case, we Setup the structural time integration after the
   * discretization has been redistributed. See Solid::TimeInt::Base::read_restart()
   * for more information.                                     hiermeier 05/16*/
  if (not restart) ti_strategy->setup();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::set_structure_wrapper(
    const Teuchos::ParameterList& ioflags, const Teuchos::ParameterList& sdyn,
    const Teuchos::ParameterList& xparams, const Teuchos::ParameterList& time_adaptivity_params,
    Teuchos::RCP<Solid::TimeInt::Base> ti_strategy)
{
  // try to firstly create the adaptive wrapper
  if (str_wrapper_.is_null())
    str_wrapper_ = Adapter::StructureTimeAda::create(time_adaptivity_params, ti_strategy);

  // if no adaptive wrapper was found, we try to create a standard one
  if (str_wrapper_.is_null()) create_wrapper(ti_strategy);

  if (str_wrapper_.is_null()) FOUR_C_THROW("No proper time integration found!");
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::create_wrapper(
    Teuchos::RCP<Solid::TimeInt::Base> ti_strategy)
{
  // get the problem instance and the problem type
  Global::Problem* problem = Global::Problem::instance();
  Core::ProblemType probtype = problem->get_problem_type();

  switch (probtype)
  {
    case Core::ProblemType::fsi:
    case Core::ProblemType::fsi_redmodels:
    case Core::ProblemType::gas_fsi:
    case Core::ProblemType::biofilm_fsi:
    case Core::ProblemType::thermo_fsi:
    case Core::ProblemType::fsi_xfem:
    {
      // Are there any constraint conditions active?
      const std::set<Inpar::Solid::ModelType>& modeltypes =
          ti_strategy->get_data_sdyn().get_model_types();
      if (modeltypes.find(Inpar::Solid::model_lag_pen_constraint) != modeltypes.end())
      {
        if ((actdis_->get_comm()).MyPID() == 0)
          Core::IO::cout << "Using StructureNOXCorrectionWrapper()..." << Core::IO::endl;

        str_wrapper_ = Teuchos::make_rcp<StructureConstrMerged>(
            Teuchos::make_rcp<StructureNOXCorrectionWrapper>(ti_strategy));
      }
      else
      {
        // case of partitioned fsi
        str_wrapper_ = Teuchos::make_rcp<FSIStructureWrapper>(ti_strategy);
      }
      break;
    }
    case Core::ProblemType::fbi:
    {
      const Teuchos::ParameterList& fsidyn = problem->fsi_dynamic_params();
      if (Teuchos::getIntegralValue<Inpar::FSI::PartitionedCouplingMethod>(
              fsidyn.sublist("PARTITIONED SOLVER"), "PARTITIONED") == Inpar::FSI::DirichletNeumann)
        str_wrapper_ = Teuchos::make_rcp<FBIStructureWrapper>(ti_strategy);
      else
        FOUR_C_THROW("Only DirichletNeumann is implemented for FBI so far");
      break;
    }
    case Core::ProblemType::immersed_fsi:
    {
      str_wrapper_ = Teuchos::make_rcp<FSIStructureWrapperImmersed>(ti_strategy);
      break;
    }
    case Core::ProblemType::ssi:
    case Core::ProblemType::ssti:
    {
      str_wrapper_ = Teuchos::make_rcp<SSIStructureWrapper>(ti_strategy);
      break;
    }
    case Core::ProblemType::pasi:
    {
      str_wrapper_ = Teuchos::make_rcp<PASIStructureWrapper>(ti_strategy);
      break;
    }
    case Core::ProblemType::redairways_tissue:
      str_wrapper_ = Teuchos::make_rcp<StructureRedAirway>(ti_strategy);
      break;
    case Core::ProblemType::poroelast:
    case Core::ProblemType::poroscatra:
    case Core::ProblemType::fpsi:
    case Core::ProblemType::fps3i:
    case Core::ProblemType::fpsi_xfem:
    {
      const Teuchos::ParameterList& porodyn = problem->poroelast_dynamic_params();
      const auto coupling = Teuchos::getIntegralValue<Inpar::PoroElast::SolutionSchemeOverFields>(
          porodyn, "COUPALGO");
      // Are there any constraint conditions active?
      const std::set<Inpar::Solid::ModelType>& modeltypes =
          ti_strategy->get_data_sdyn().get_model_types();
      if (modeltypes.find(Inpar::Solid::model_lag_pen_constraint) != modeltypes.end())
      {
        if (coupling == Inpar::PoroElast::Monolithic_structuresplit or
            coupling == Inpar::PoroElast::Monolithic_fluidsplit or
            coupling == Inpar::PoroElast::Monolithic_nopenetrationsplit)
          str_wrapper_ = Teuchos::make_rcp<FPSIStructureWrapper>(ti_strategy);
        else
          str_wrapper_ = Teuchos::make_rcp<StructureConstrMerged>(ti_strategy);
      }
      else
      {
        str_wrapper_ = Teuchos::make_rcp<FPSIStructureWrapper>(ti_strategy);
      }
      break;
    }
    default:
      /// wrap time loop for pure structure problems
      str_wrapper_ = (Teuchos::make_rcp<StructureTimeLoop>(ti_strategy));
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
