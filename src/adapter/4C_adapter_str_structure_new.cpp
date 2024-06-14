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
#include "4C_adapter_str_fsi_timint_adaptive.hpp"
#include "4C_adapter_str_fsiwrapper_immersed.hpp"
#include "4C_adapter_str_lung.hpp"
#include "4C_adapter_str_pasiwrapper.hpp"
#include "4C_adapter_str_redairway.hpp"
#include "4C_adapter_str_ssiwrapper.hpp"
#include "4C_adapter_str_timeada.hpp"
#include "4C_adapter_str_timeloop.hpp"
#include "4C_adapter_str_timint_adaptive.hpp"
#include "4C_adapter_str_wrapper.hpp"
#include "4C_beam3_kirchhoff.hpp"
#include "4C_beam3_reissner.hpp"
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
#include "4C_inpar_validparameters.hpp"
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
#include "4C_so3_sh8p8.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_solver_nonlin_nox_group.hpp"
#include "4C_solver_nonlin_nox_group_prepostoperator.hpp"
#include "4C_structure_new_model_evaluator.hpp"
#include "4C_structure_new_solver_factory.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_structure_new_timint_factory.hpp"
#include "4C_structure_timada_create.hpp"
#include "4C_utils_exceptions.hpp"

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
  switch (Core::UTILS::IntegralValue<Inpar::STR::DynamicType>(*sdyn_, "DYNAMICTYP"))
  {
    case Inpar::STR::dyna_statics:
    case Inpar::STR::dyna_genalpha:
    case Inpar::STR::dyna_genalpha_liegroup:
    case Inpar::STR::dyna_onesteptheta:
    case Inpar::STR::dyna_expleuler:
    case Inpar::STR::dyna_centrdiff:
    case Inpar::STR::dyna_ab2:
    case Inpar::STR::dyna_ab4:
    case Inpar::STR::dyna_euma:
    case Inpar::STR::dyna_euimsto:
      setup_tim_int();  // <-- here is the show
      break;
    default:
      FOUR_C_THROW(
          "Unknown time integration scheme '%s'", sdyn_->get<std::string>("DYNAMICTYP").c_str());
      break;
  }

  issetup_ = true;
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::register_model_evaluator(
    const std::string name, Teuchos::RCP<STR::MODELEVALUATOR::Generic> me)
{
  // safety checks
  if (not is_init())
    FOUR_C_THROW("init(...) must be called before register_model_evaluator(...) !");
  if (is_setup()) FOUR_C_THROW("register_model_evaluator(...) must be called before setup() !");

  // set RCP ptr to model evaluator in problem dynamic parameter list
  const_cast<Teuchos::ParameterList&>(*prbdyn_).set<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>(
      name, me);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::setup_tim_int()
{
  if (not is_init()) FOUR_C_THROW("You have to call init() first!");

  // get the problem instance
  Global::Problem* problem = Global::Problem::Instance();
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
  if (actdis_->GetCondition("PointCoupling") != nullptr)
  {
    std::vector<Teuchos::RCP<Core::FE::Discretization>> actdis_vec(1, actdis_);
    Teuchos::ParameterList binning_params = Global::Problem::Instance()->binning_strategy_params();
    Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
        "spatial_approximation_type", Global::Problem::Instance()->spatial_approximation_type(),
        binning_params);
    actdis_vec[0]->fill_complete(false, false, false);
    auto element_filter = [](const Core::Elements::Element* element)
    {
      if (dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(element))
        return BINSTRATEGY::UTILS::SpecialElement::beam;
      else if (element->ElementType() == Discret::ELEMENTS::RigidsphereType::Instance())
        return BINSTRATEGY::UTILS::SpecialElement::rigid_sphere;
      else
        return BINSTRATEGY::UTILS::SpecialElement::none;
    };

    auto rigid_sphere_radius = [](const Core::Elements::Element* element)
    {
      if (element->ElementType() == Discret::ELEMENTS::RigidsphereType::Instance())
        return dynamic_cast<const Discret::ELEMENTS::Rigidsphere*>(element)->Radius();
      else
        return 0.0;
    };
    auto correct_beam_center_node = [](const Core::Nodes::Node* node)
    {
      const Core::Elements::Element* element = node->Elements()[0];
      const auto* beamelement = dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(element);
      if (beamelement != nullptr and not beamelement->IsCenterlineNode(*node))
        return element->Nodes()[0];
      else
        return node;
    };
    Core::Rebalance::RebalanceDiscretizationsByBinning(binning_params,
        Global::Problem::Instance()->OutputControlFile(), actdis_vec, element_filter,
        rigid_sphere_radius, correct_beam_center_node, true);
  }
  else if (not actdis_->Filled() || not actdis_->HaveDofs())
  {
    actdis_->fill_complete();
  }

  // ---------------------------------------------------------------------------
  // Setup a model type set by checking
  // the different conditions
  // ---------------------------------------------------------------------------
  // define and initial with default value
  Teuchos::RCP<std::set<enum Inpar::STR::ModelType>> modeltypes =
      Teuchos::rcp(new std::set<enum Inpar::STR::ModelType>());
  modeltypes->insert(Inpar::STR::model_structure);
  set_model_types(*modeltypes);

  // ---------------------------------------------------------------------------
  // Setup a element technology set by checking
  // the elements of the discretization
  // ---------------------------------------------------------------------------
  Teuchos::RCP<std::set<enum Inpar::STR::EleTech>> eletechs =
      Teuchos::rcp(new std::set<enum Inpar::STR::EleTech>());
  detect_element_technologies(*eletechs);

  // ---------------------------------------------------------------------------
  // Setup the parameter lists for structural
  // time integration
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Teuchos::ParameterList> ioflags =
      Teuchos::rcp(new Teuchos::ParameterList(problem->IOParams()));
  Teuchos::RCP<Teuchos::ParameterList> time_adaptivity_params =
      Teuchos::rcp(new Teuchos::ParameterList(sdyn_->sublist("TIMEADAPTIVITY")));
  Teuchos::RCP<Teuchos::ParameterList> xparams = Teuchos::rcp(new Teuchos::ParameterList());
  set_params(*ioflags, *xparams, *time_adaptivity_params);

  // ---------------------------------------------------------------------------
  // Setup and create model specific linear solvers
  // ---------------------------------------------------------------------------
  Teuchos::RCP<std::map<enum Inpar::STR::ModelType, Teuchos::RCP<Core::LinAlg::Solver>>>
      linsolvers = STR::SOLVER::build_lin_solvers(*modeltypes, *sdyn_, *actdis_);

  // ---------------------------------------------------------------------------
  // Checks in case of multi-scale simulations
  // ---------------------------------------------------------------------------
  {
    // make sure we IMR-like generalised-alpha requested for multi-scale
    // simulations
    Teuchos::RCP<Mat::PAR::Bundle> materials = problem->Materials();
    for (const auto& [_, mat] : materials->Map())
    {
      if (mat->Type() == Core::Materials::m_struct_multiscale)
      {
        if (Core::UTILS::IntegralValue<Inpar::STR::DynamicType>(*sdyn_, "DYNAMICTYP") !=
            Inpar::STR::dyna_genalpha)
          FOUR_C_THROW("In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha");
        else if (Core::UTILS::IntegralValue<Inpar::STR::MidAverageEnum>(
                     sdyn_->sublist("GENALPHA"), "GENAVG") != Inpar::STR::midavg_trlike)
          FOUR_C_THROW(
              "In multi-scale simulations, you have to use DYNAMICTYP=GenAlpha with GENAVG=TrLike");
        break;
      }
    }
  }

  // ---------------------------------------------------------------------------
  // Create context for output and restart
  // ---------------------------------------------------------------------------
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = actdis_->Writer();
  if (Core::UTILS::IntegralValue<int>(*ioflags, "OUTPUT_BIN"))
  {
    output->write_mesh(0, 0.0);
  }

  // ---------------------------------------------------------------------------
  // initialize/setup the input/output data container
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TimeInt::BaseDataIO> dataio = Teuchos::rcp(new STR::TimeInt::BaseDataIO());
  dataio->init(*ioflags, *sdyn_, *xparams, output);
  dataio->setup();

  // ---------------------------------------------------------------------------
  // initialize/setup the structural dynamics data
  // container
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TimeInt::BaseDataSDyn> datasdyn = STR::TimeInt::build_data_sdyn(*sdyn_);
  datasdyn->init(actdis_, *sdyn_, *xparams, modeltypes, eletechs, linsolvers);
  datasdyn->setup();

  // ---------------------------------------------------------------------------
  // initialize/setup the global state data container
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TimeInt::BaseDataGlobalState> dataglobalstate = Teuchos::null;
  set_global_state(dataglobalstate, datasdyn);

  // ---------------------------------------------------------------------------
  // in case of non-additive rotation (pseudo-)vector DOFs:
  // ---------------------------------------------------------------------------
  if (eletechs->find(Inpar::STR::EleTech::rotvec) != eletechs->end())
  {
    // -------------------------------------------------------------------------
    // set the RotVecUpdater as new precomputeX operator for the nox nln group
    // -------------------------------------------------------------------------
    Teuchos::ParameterList& p_grp_opt = datasdyn->get_nox_params().sublist("Group Options");
    // Get the current map. If there is no map, return a new empty one. (reference)
    NOX::Nln::GROUP::PrePostOperator::Map& prepostgroup_map =
        NOX::Nln::GROUP::PrePostOp::GetMap(p_grp_opt);
    // create the new rotation vector update pre/post operator
    Teuchos::RCP<NOX::Nln::Abstract::PrePostOperator> prepostrotvec_ptr =
        Teuchos::rcp(new NOX::Nln::GROUP::PrePostOp::TimeInt::RotVecUpdater(dataglobalstate));
    // insert/replace the old pointer in the map
    prepostgroup_map[NOX::Nln::GROUP::prepost_rotvecupdate] = prepostrotvec_ptr;
  }

  // ---------------------------------------------------------------------------
  // Build time integrator
  // ---------------------------------------------------------------------------
  Teuchos::RCP<STR::TimeInt::Base> ti_strategy = Teuchos::null;
  set_time_integration_strategy(ti_strategy, dataio, datasdyn, dataglobalstate, restart);


  // ---------------------------------------------------------------------------
  // Create wrapper for the time integration strategy
  // ---------------------------------------------------------------------------
  set_structure_wrapper(*ioflags, *sdyn_, *xparams, *time_adaptivity_params, ti_strategy);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::set_model_types(
    std::set<enum Inpar::STR::ModelType>& modeltypes) const
{
  if (not is_init()) FOUR_C_THROW("You have to call init() first!");
  // ---------------------------------------------------------------------------
  // check for meshtying and contact conditions
  // ---------------------------------------------------------------------------
  // --- contact conditions
  std::vector<Core::Conditions::Condition*> ccond(0);
  actdis_->GetCondition("Contact", ccond);
  if (ccond.size())
  {
    // what's the current problem type?
    Core::ProblemType probtype = Global::Problem::Instance()->GetProblemType();
    // ToDo: once the new structural time integration can handle
    //       condensed contact formulations, the model_evaluator
    //       can have its contact model. For now, the TSI Lagrange
    //       strategy resides in the TSI algorithm.
    if (probtype == Core::ProblemType::tsi)
    {
      const Teuchos::ParameterList& contact = Global::Problem::Instance()->contact_dynamic_params();
      if (Core::UTILS::IntegralValue<Inpar::CONTACT::SolvingStrategy>(contact, "STRATEGY") ==
          Inpar::CONTACT::solution_nitsche)
        modeltypes.insert(Inpar::STR::model_contact);
    }
    else
      modeltypes.insert(Inpar::STR::model_contact);
  }
  // --- meshtying conditions
  std::vector<Core::Conditions::Condition*> mtcond(0);
  actdis_->GetCondition("Mortar", mtcond);
  if (mtcond.size()) modeltypes.insert(Inpar::STR::model_meshtying);

  // check for 0D cardiovascular conditions
  // ---------------------------------------------------------------------------
  std::vector<Core::Conditions::Condition*> cardiovasc0dcond_4elementwindkessel(0);
  std::vector<Core::Conditions::Condition*> cardiovasc0dcond_arterialproxdist(0);
  std::vector<Core::Conditions::Condition*> cardiovasc0dcond_syspulcirculation(0);
  std::vector<Core::Conditions::Condition*> cardiovascrespir0dcond_syspulperiphcirculation(0);
  actdis_->GetCondition(
      "Cardiovascular0D4ElementWindkesselStructureCond", cardiovasc0dcond_4elementwindkessel);
  actdis_->GetCondition(
      "Cardiovascular0DArterialProxDistStructureCond", cardiovasc0dcond_arterialproxdist);
  actdis_->GetCondition("Cardiovascular0DSysPulCirculationStructureCond",
      cardiovascrespir0dcond_syspulperiphcirculation);
  actdis_->GetCondition("CardiovascularRespiratory0DSysPulPeriphCirculationStructureCond",
      cardiovasc0dcond_syspulcirculation);
  if (cardiovasc0dcond_4elementwindkessel.size() or cardiovasc0dcond_arterialproxdist.size() or
      cardiovasc0dcond_syspulcirculation.size() or
      cardiovascrespir0dcond_syspulperiphcirculation.size())
    modeltypes.insert(Inpar::STR::model_cardiovascular0d);

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
  actdis_->GetCondition("VolumeConstraint_3D", lagcond_volconstr3d);
  actdis_->GetCondition("AreaConstraint_3D", lagcond_areaconstr3d);
  actdis_->GetCondition("AreaConstraint_2D", lagcond_areaconstr2d);
  actdis_->GetCondition("MPC_NodeOnLine_2D", lagcond_mpconline2d);
  actdis_->GetCondition("MPC_NodeOnPlane_3D", lagcond_mpconplane3d);
  actdis_->GetCondition("MPC_NormalComponent_3D", lagcond_mpcnormcomp3d);
  if (lagcond_volconstr3d.size() or lagcond_areaconstr3d.size() or lagcond_areaconstr2d.size() or
      lagcond_mpconline2d.size() or lagcond_mpconplane3d.size() or lagcond_mpcnormcomp3d.size())
    have_lag_constraint = true;
  // --- enforcement by penalty law
  std::vector<Core::Conditions::Condition*> pencond_volconstr3d(0);
  std::vector<Core::Conditions::Condition*> pencond_areaconstr3d(0);
  std::vector<Core::Conditions::Condition*> pencond_mpcnormcomp3d(0);
  actdis_->GetCondition("VolumeConstraint_3D_Pen", pencond_volconstr3d);
  actdis_->GetCondition("AreaConstraint_3D_Pen", pencond_areaconstr3d);
  actdis_->GetCondition("MPC_NormalComponent_3D_Pen", pencond_mpcnormcomp3d);
  if (pencond_volconstr3d.size() or pencond_areaconstr3d.size() or pencond_mpcnormcomp3d.size())
    have_pen_constraint = true;
  if (have_lag_constraint or have_pen_constraint)
    modeltypes.insert(Inpar::STR::model_lag_pen_constraint);

  // ---------------------------------------------------------------------------
  // check for spring dashpot conditions
  // ---------------------------------------------------------------------------
  std::vector<Core::Conditions::Condition*> sdp_cond(0);
  actdis_->GetCondition("RobinSpringDashpot", sdp_cond);
  if (sdp_cond.size()) modeltypes.insert(Inpar::STR::model_springdashpot);
  // ---------------------------------------------------------------------------
  // check for coupled problems
  // ---------------------------------------------------------------------------
  // get the problem instance
  Global::Problem* problem = Global::Problem::Instance();
  // what's the current problem type?
  Core::ProblemType probtype = problem->GetProblemType();
  switch (probtype)
  {
    case Core::ProblemType::fsi:
    case Core::ProblemType::immersed_fsi:
    case Core::ProblemType::fbi:
    case Core::ProblemType::fsi_redmodels:
    case Core::ProblemType::fsi_lung:
    case Core::ProblemType::gas_fsi:
    case Core::ProblemType::ac_fsi:
    case Core::ProblemType::biofilm_fsi:
    case Core::ProblemType::thermo_fsi:
    case Core::ProblemType::fsi_xfem:
    case Core::ProblemType::pasi:
    case Core::ProblemType::ssi:
    case Core::ProblemType::ssti:
    {
      if (prbdyn_->INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>(
              "Partitioned Coupling Model"))
      {
        if (prbdyn_->INVALID_TEMPLATE_QUALIFIER isType<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>(
                "Monolithic Coupling Model"))
          FOUR_C_THROW("Cannot have both partitioned and monolithic coupling at the same time!");
        const auto coupling_model_ptr =
            prbdyn_->INVALID_TEMPLATE_QUALIFIER get<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>(
                "Partitioned Coupling Model");
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW(
              "The partitioned coupling model pointer is not allowed to be Teuchos::null!");
        // set the model type
        modeltypes.insert(Inpar::STR::model_partitioned_coupling);
        // copy the coupling model object pointer into the (temporal) sdyn parameter list
        sdyn_->set<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>(
            "Partitioned Coupling Model", coupling_model_ptr);
      }

      else if (prbdyn_->INVALID_TEMPLATE_QUALIFIER
                   isType<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>("Monolithic Coupling Model"))
      {
        const auto coupling_model_ptr =
            prbdyn_->INVALID_TEMPLATE_QUALIFIER get<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>(
                "Monolithic Coupling Model");
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW("The monolithic coupling model pointer is not allowed to be Teuchos::null!");
        // set the model type
        modeltypes.insert(Inpar::STR::model_monolithic_coupling);
        // copy the coupling model object pointer into the (temporal) sdyn parameter list
        sdyn_->set<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>(
            "Monolithic Coupling Model", coupling_model_ptr);
      }

      else if (prbdyn_->INVALID_TEMPLATE_QUALIFIER
                   isType<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>("Basic Coupling Model"))
      {
        const auto coupling_model_ptr =
            prbdyn_->INVALID_TEMPLATE_QUALIFIER get<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>(
                "Basic Coupling Model");
        if (coupling_model_ptr.is_null())
          FOUR_C_THROW("The basic coupling model pointer is not allowed to be Teuchos::null!");
        // set the model type
        modeltypes.insert(Inpar::STR::model_basic_coupling);
        // copy the coupling model object pointer into the (temporal) sdyn parameter list
        sdyn_->set<Teuchos::RCP<STR::MODELEVALUATOR::Generic>>(
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
  const Teuchos::ParameterList& beamcontact = Global::Problem::Instance()->beam_contact_params();
  Inpar::BEAMCONTACT::Strategy strategy =
      Core::UTILS::IntegralValue<Inpar::BEAMCONTACT::Strategy>(beamcontact, "BEAMS_STRATEGY");

  Inpar::BEAMCONTACT::Modelevaluator modelevaluator =
      Core::UTILS::IntegralValue<Inpar::BEAMCONTACT::Modelevaluator>(beamcontact, "MODELEVALUATOR");

  // conditions for potential-based beam interaction
  std::vector<Core::Conditions::Condition*> beampotconditions(0);
  actdis_->GetCondition("BeamPotentialLineCharge", beampotconditions);

  // conditions for beam penalty point coupling
  std::vector<Core::Conditions::Condition*> beampenaltycouplingconditions(0);
  actdis_->GetCondition("PenaltyPointCouplingCondition", beampenaltycouplingconditions);


  if (strategy != Inpar::BEAMCONTACT::bstr_none and modelevaluator == Inpar::BEAMCONTACT::bstr_old)
    modeltypes.insert(Inpar::STR::model_beam_interaction_old);

  // ---------------------------------------------------------------------------
  // check for brownian dynamics
  // ---------------------------------------------------------------------------
  if (Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->brownian_dynamics_params(), "BROWNDYNPROB"))
    modeltypes.insert(Inpar::STR::model_browniandyn);

  // ---------------------------------------------------------------------------
  // check for beam interaction
  // ---------------------------------------------------------------------------
  if (Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->beam_interaction_params().sublist("CROSSLINKING"),
          "CROSSLINKER") or
      Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->beam_interaction_params().sublist("SPHERE BEAM LINK"),
          "SPHEREBEAMLINKING") or
      Core::UTILS::IntegralValue<Inpar::BEAMINTERACTION::Strategy>(
          Global::Problem::Instance()->beam_interaction_params().sublist("BEAM TO BEAM CONTACT"),
          "STRATEGY") != Inpar::BEAMINTERACTION::bstr_none or
      Core::UTILS::IntegralValue<Inpar::BEAMINTERACTION::Strategy>(
          Global::Problem::Instance()->beam_interaction_params().sublist("BEAM TO SPHERE CONTACT"),
          "STRATEGY") != Inpar::BEAMINTERACTION::bstr_none or
      Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          Global::Problem::Instance()->beam_interaction_params().sublist(
              "BEAM TO SOLID VOLUME MESHTYING"),
          "CONTACT_DISCRETIZATION") != Inpar::BeamToSolid::BeamToSolidContactDiscretization::none or
      Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          Global::Problem::Instance()->beam_interaction_params().sublist(
              "BEAM TO SOLID SURFACE MESHTYING"),
          "CONTACT_DISCRETIZATION") != Inpar::BeamToSolid::BeamToSolidContactDiscretization::none or
      Teuchos::getIntegralValue<Inpar::BeamToSolid::BeamToSolidContactDiscretization>(
          Global::Problem::Instance()->beam_interaction_params().sublist(
              "BEAM TO SOLID SURFACE CONTACT"),
          "CONTACT_DISCRETIZATION") != Inpar::BeamToSolid::BeamToSolidContactDiscretization::none or
      beampotconditions.size() > 0 or beampenaltycouplingconditions.size() > 0)
    modeltypes.insert(Inpar::STR::model_beaminteraction);

  // ---------------------------------------------------------------------------
  // check for constraints
  // ---------------------------------------------------------------------------
  std::vector<Teuchos::RCP<Core::Conditions::Condition>> linePeriodicRve, surfPeriodicRve,
      pointLinearCoupledEquation;
  actdis_->GetCondition("LinePeriodicRve", linePeriodicRve);
  actdis_->GetCondition("SurfacePeriodicRve", surfPeriodicRve);
  actdis_->GetCondition("PointLinearCoupledEquation", pointLinearCoupledEquation);

  if (linePeriodicRve.size() > 0 || surfPeriodicRve.size() > 0 ||
      pointLinearCoupledEquation.size() > 0)
    modeltypes.insert(Inpar::STR::model_constraints);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::detect_element_technologies(
    std::set<enum Inpar::STR::EleTech>& eletechs) const
{
  int isplasticity_local = 0;
  int isplasticity_global = 0;

  int iseas_local = 0;
  int iseas_global = 0;

  int isfbar_local = 0;
  int isfbar_global = 0;

  int ispressure_local = 0;
  int ispressure_global = 0;

  int isrotvec_local = 0;
  int isrotvec_global = 0;

  for (int i = 0; i < actdis_->NumMyRowElements(); ++i)
  {
    Core::Elements::Element* actele = actdis_->lRowElement(i);
    // Detect plasticity -------------------------------------------------------
    if (actele->ElementType() == Discret::ELEMENTS::SoHex8PlastType::Instance() or
        actele->ElementType() == Discret::ELEMENTS::SoHex27PlastType::Instance() or
        actele->ElementType() == Discret::ELEMENTS::SoSh8PlastType::Instance() or
        actele->ElementType() == Discret::ELEMENTS::SoHex18PlastType::Instance() or
        actele->ElementType() == Discret::ELEMENTS::SoSh18PlastType::Instance())
    {
      if (actele->Material()->MaterialType() == Core::Materials::m_plelasthyper)
        isplasticity_local = true;
    }

    // Detect EAS --------------------------------------------------------------
    Discret::ELEMENTS::SoBase* so_base_ele = dynamic_cast<Discret::ELEMENTS::SoBase*>(actele);
    if (so_base_ele != nullptr)
    {
      if (so_base_ele->HaveEAS()) iseas_local = 1;
    }

    Discret::ELEMENTS::Shell7p* shell7p = dynamic_cast<Discret::ELEMENTS::Shell7p*>(actele);
    if (shell7p)
      if (shell7p->GetEleTech().find(Inpar::STR::EleTech::eas) != shell7p->GetEleTech().end())
        iseas_local = 1;

    Discret::ELEMENTS::Solid* solid = dynamic_cast<Discret::ELEMENTS::Solid*>(actele);
    if (solid != nullptr)
      if (solid->HaveEAS()) iseas_local = 1;

    // Detect additional pressure dofs -----------------------------------------
    if (actele->ElementType() == Discret::ELEMENTS::SoSh8p8Type::Instance()) ispressure_local = 1;

    // Detect fbar
    Discret::ELEMENTS::SoHex8fbar* so_hex8fbar_ele =
        dynamic_cast<Discret::ELEMENTS::SoHex8fbar*>(actele);
    if (so_hex8fbar_ele != nullptr) isfbar_local = 1;

    // Detect non-additive rotation-vector DOFs --------------------------------
    if (actele->ElementType() == Discret::ELEMENTS::Beam3rType::Instance() or
        actele->ElementType() == Discret::ELEMENTS::Beam3kType::Instance())
    {
      isrotvec_local = true;
      break;
    }
  }

  // plasticity - sum over all processors
  actdis_->Comm().SumAll(&isplasticity_local, &isplasticity_global, 1);
  if (isplasticity_global > 0) eletechs.insert(Inpar::STR::EleTech::plasticity);

  // eas - sum over all processors
  actdis_->Comm().SumAll(&iseas_local, &iseas_global, 1);
  if (iseas_global > 0) eletechs.insert(Inpar::STR::EleTech::eas);

  // pressure - sum over all processors
  actdis_->Comm().SumAll(&ispressure_local, &ispressure_global, 1);
  if (ispressure_global > 0) eletechs.insert(Inpar::STR::EleTech::pressure);

  // fbar - sum over all processors
  actdis_->Comm().SumAll(&isfbar_local, &isfbar_global, 1);
  if (isfbar_global > 0) eletechs.insert(Inpar::STR::EleTech::fbar);

  // rotation vector DOFs - sum over all processors
  actdis_->Comm().SumAll(&isrotvec_local, &isrotvec_global, 1);
  if (isrotvec_global > 0) eletechs.insert(Inpar::STR::EleTech::rotvec);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::set_params(Teuchos::ParameterList& ioflags,
    Teuchos::ParameterList& xparams, Teuchos::ParameterList& time_adaptivity_params)
{
  // get the problem instance and the problem type
  Global::Problem* problem = Global::Problem::Instance();
  Core::ProblemType probtype = problem->GetProblemType();

  // ---------------------------------------------------------------------------
  // show default parameters
  // ---------------------------------------------------------------------------
  if ((actdis_->Comm()).MyPID() == 0) Input::PrintDefaultParameters(Core::IO::cout, *sdyn_);

  // ---------------------------------------------------------------------------
  // get input parameter lists and copy them,
  // because a few parameters are overwritten
  // ---------------------------------------------------------------------------
  // nox parameter list
  Teuchos::RCP<Teuchos::ParameterList> snox =
      Teuchos::rcp(new Teuchos::ParameterList(problem->StructuralNoxParams()));
  Teuchos::ParameterList& nox = xparams.sublist("NOX");
  nox = *snox;

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
  if (Core::UTILS::IntegralValue<Inpar::STR::DampKind>(*sdyn_, "DAMPING") ==
      Inpar::STR::damp_rayleigh)
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
      const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
      const Teuchos::ParameterList& fsiada = fsidyn.sublist("TIMEADAPTIVITY");
      if (Core::UTILS::IntegralValue<bool>(fsiada, "TIMEADAPTON"))
      {
        // overrule time step size adaptivity control parameters
        if (time_adaptivity_params.get<std::string>("KIND") != "NONE")
        {
          time_adaptivity_params.set<int>("ADAPTSTEPMAX", fsiada.get<int>("ADAPTSTEPMAX"));
          time_adaptivity_params.set<double>("STEPSIZEMAX", fsiada.get<double>("DTMAX"));
          time_adaptivity_params.set<double>("STEPSIZEMIN", fsiada.get<double>("DTMIN"));
          time_adaptivity_params.set<double>("SIZERATIOMAX", fsiada.get<double>("SIZERATIOMAX"));
          time_adaptivity_params.set<double>("SIZERATIOMIN", fsiada.get<double>("SIZERATIOMIN"));
          time_adaptivity_params.set<double>("SIZERATIOSCALE", fsiada.get<double>("SAFETYFACTOR"));

          if (actdis_->Comm().MyPID() == 0)
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
    Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& dataglobalstate,
    const Teuchos::RCP<const STR::TimeInt::BaseDataSDyn>& datasdyn)
{
  dataglobalstate = STR::TimeInt::build_data_global_state();
  dataglobalstate->init(actdis_, *sdyn_, datasdyn);
  dataglobalstate->setup();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::set_time_integration_strategy(
    Teuchos::RCP<STR::TimeInt::Base>& ti_strategy,
    const Teuchos::RCP<STR::TimeInt::BaseDataIO>& dataio,
    const Teuchos::RCP<STR::TimeInt::BaseDataSDyn>& datasdyn,
    const Teuchos::RCP<STR::TimeInt::BaseDataGlobalState>& dataglobalstate, const int& restart)
{
  ti_strategy = STR::TimeInt::build_strategy(*sdyn_);
  ti_strategy->init(dataio, datasdyn, dataglobalstate);

  /* In the restart case, we Setup the structural time integration after the
   * discretization has been redistributed. See STR::TimeInt::Base::read_restart()
   * for more information.                                     hiermeier 05/16*/
  if (not restart) ti_strategy->setup();
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::set_structure_wrapper(
    const Teuchos::ParameterList& ioflags, const Teuchos::ParameterList& sdyn,
    const Teuchos::ParameterList& xparams, const Teuchos::ParameterList& time_adaptivity_params,
    Teuchos::RCP<STR::TimeInt::Base> ti_strategy)
{
  // try to firstly create the adaptive wrapper
  if (str_wrapper_.is_null())
    str_wrapper_ = Adapter::StructureTimeAda::Create(time_adaptivity_params, ti_strategy);

  // if no adaptive wrapper was found, we try to create a standard one
  if (str_wrapper_.is_null()) create_wrapper(ti_strategy);

  if (str_wrapper_.is_null()) FOUR_C_THROW("No proper time integration found!");
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void Adapter::StructureBaseAlgorithmNew::create_wrapper(
    Teuchos::RCP<STR::TimeInt::Base> ti_strategy)
{
  // get the problem instance and the problem type
  Global::Problem* problem = Global::Problem::Instance();
  Core::ProblemType probtype = problem->GetProblemType();

  switch (probtype)
  {
    case Core::ProblemType::fsi:
    case Core::ProblemType::fsi_redmodels:
    case Core::ProblemType::fsi_lung:
    case Core::ProblemType::gas_fsi:
    case Core::ProblemType::ac_fsi:
    case Core::ProblemType::biofilm_fsi:
    case Core::ProblemType::thermo_fsi:
    case Core::ProblemType::fsi_xfem:
    {
      const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
      const int coupling = Core::UTILS::IntegralValue<int>(fsidyn, "COUPALGO");

      // Are there any constraint conditions active?
      const std::set<Inpar::STR::ModelType>& modeltypes =
          ti_strategy->get_data_sdyn().GetModelTypes();
      if (modeltypes.find(Inpar::STR::model_lag_pen_constraint) != modeltypes.end())
      {
        if ((actdis_->Comm()).MyPID() == 0)
          Core::IO::cout << "Using StructureNOXCorrectionWrapper()..." << Core::IO::endl;

        if (coupling == fsi_iter_constr_monolithicstructuresplit or
            coupling == fsi_iter_constr_monolithicfluidsplit)
          str_wrapper_ = Teuchos::rcp(new FSIStructureWrapper(
              Teuchos::rcp(new StructureNOXCorrectionWrapper(ti_strategy))));
        else
          str_wrapper_ = Teuchos::rcp(new StructureConstrMerged(
              Teuchos::rcp(new StructureNOXCorrectionWrapper(ti_strategy))));
      }
      else
      {
        if (coupling == fsi_iter_lung_monolithicstructuresplit or
            coupling == fsi_iter_lung_monolithicfluidsplit)
        {
          if ((actdis_->Comm()).MyPID() == 0)
            Core::IO::cout << "Using StructureNOXCorrectionWrapper()..." << Core::IO::endl;
          str_wrapper_ = Teuchos::rcp(
              new StructureLung(Teuchos::rcp(new StructureNOXCorrectionWrapper(ti_strategy))));
        }
        else
        {
          str_wrapper_ =
              Teuchos::rcp(new FSIStructureWrapper(ti_strategy));  // case of partitioned fsi
        }
      }
      break;
    }
    case Core::ProblemType::fbi:
    {
      const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
      if (Core::UTILS::IntegralValue<Inpar::FSI::PartitionedCouplingMethod>(
              fsidyn.sublist("PARTITIONED SOLVER"), "PARTITIONED") == Inpar::FSI::DirichletNeumann)
        str_wrapper_ = Teuchos::rcp(new FBIStructureWrapper(ti_strategy));
      else
        FOUR_C_THROW("Only DirichletNeumann is implemented for FBI so far");
      break;
    }
    case Core::ProblemType::immersed_fsi:
    {
      str_wrapper_ = Teuchos::rcp(new FSIStructureWrapperImmersed(ti_strategy));
      break;
    }
    case Core::ProblemType::ssi:
    case Core::ProblemType::ssti:
    {
      str_wrapper_ = Teuchos::rcp(new SSIStructureWrapper(ti_strategy));
      break;
    }
    case Core::ProblemType::pasi:
    {
      str_wrapper_ = Teuchos::rcp(new PASIStructureWrapper(ti_strategy));
      break;
    }
    case Core::ProblemType::redairways_tissue:
      str_wrapper_ = Teuchos::rcp(new StructureRedAirway(ti_strategy));
      break;
    case Core::ProblemType::poroelast:
    case Core::ProblemType::poroscatra:
    case Core::ProblemType::fpsi:
    case Core::ProblemType::fps3i:
    case Core::ProblemType::fpsi_xfem:
    {
      const Teuchos::ParameterList& porodyn = problem->poroelast_dynamic_params();
      const Inpar::PoroElast::SolutionSchemeOverFields coupling =
          Core::UTILS::IntegralValue<Inpar::PoroElast::SolutionSchemeOverFields>(
              porodyn, "COUPALGO");
      // Are there any constraint conditions active?
      const std::set<Inpar::STR::ModelType>& modeltypes =
          ti_strategy->get_data_sdyn().GetModelTypes();
      if (modeltypes.find(Inpar::STR::model_lag_pen_constraint) != modeltypes.end())
      {
        if (coupling == Inpar::PoroElast::Monolithic_structuresplit or
            coupling == Inpar::PoroElast::Monolithic_fluidsplit or
            coupling == Inpar::PoroElast::Monolithic_nopenetrationsplit)
          str_wrapper_ = Teuchos::rcp(new FPSIStructureWrapper(ti_strategy));
        else
          str_wrapper_ = Teuchos::rcp(new StructureConstrMerged(ti_strategy));
      }
      else
      {
        str_wrapper_ = Teuchos::rcp(new FPSIStructureWrapper(ti_strategy));
      }
      break;
    }
    default:
      /// wrap time loop for pure structure problems
      str_wrapper_ = (Teuchos::rcp(new StructureTimeLoop(ti_strategy)));
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE
