/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for biomedical simulations

\level 3


*/
/*----------------------------------------------------------------------*/



#include "4C_inpar_bio.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_linecomponent.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


void Inpar::ArtDyn::set_valid_parameters(Teuchos::ParameterList& list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  Teuchos::ParameterList& andyn = list.sublist("ARTERIAL DYNAMIC", false, "");

  setStringToIntegralParameter<TimeIntegrationScheme>("DYNAMICTYP", "ExpTaylorGalerkin",
      "Explicit Taylor Galerkin Scheme", tuple<std::string>("ExpTaylorGalerkin", "Stationary"),
      tuple<TimeIntegrationScheme>(tay_gal, stationary), &andyn);

  Core::UTILS::double_parameter("TIMESTEP", 0.01, "Time increment dt", &andyn);
  Core::UTILS::int_parameter("NUMSTEP", 0, "Number of Time Steps", &andyn);
  Core::UTILS::double_parameter("MAXTIME", 1000.0, "total simulation time", &andyn);
  Core::UTILS::int_parameter("RESTARTEVRY", 1, "Increment for writing restart", &andyn);
  Core::UTILS::int_parameter("RESULTSEVRY", 1, "Increment for writing solution", &andyn);

  Core::UTILS::bool_parameter(
      "SOLVESCATRA", "no", "Flag to (de)activate solving scalar transport in blood", &andyn);

  // number of linear solver used for arterial dynamics
  Core::UTILS::int_parameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for arterial dynamics", &andyn);

  // initial function number
  Core::UTILS::int_parameter("INITFUNCNO", -1, "function number for artery initial field", &andyn);

  // type of initial field
  setStringToIntegralParameter<InitialField>("INITIALFIELD", "zero_field",
      "Initial Field for artery problem",
      tuple<std::string>("zero_field", "field_by_function", "field_by_condition"),
      tuple<InitialField>(
          initfield_zero_field, initfield_field_by_function, initfield_field_by_condition),
      &andyn);
}



void Inpar::ArteryNetwork::set_valid_parameters(Teuchos::ParameterList& list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& redtisdyn =
      list.sublist("COUPLED REDUCED-D AIRWAYS AND TISSUE DYNAMIC", false, "");
  Core::UTILS::double_parameter("CONVTOL_P", 1E-6,
      "Coupled red_airway and tissue iteration convergence for pressure", &redtisdyn);
  Core::UTILS::double_parameter("CONVTOL_Q", 1E-6,
      "Coupled red_airway and tissue iteration convergence for flux", &redtisdyn);
  Core::UTILS::int_parameter("MAXITER", 5, "Maximum coupling iterations", &redtisdyn);
  setStringToIntegralParameter<Relaxtype3D0D>("RELAXTYPE", "norelaxation",
      "Dynamic Relaxation Type",
      tuple<std::string>("norelaxation", "fixedrelaxation", "Aitken", "SD"),
      tuple<Relaxtype3D0D>(norelaxation, fixedrelaxation, Aitken, SD), &redtisdyn);
  Core::UTILS::double_parameter("TIMESTEP", 0.01, "Time increment dt", &redtisdyn);
  Core::UTILS::int_parameter("NUMSTEP", 1, "Number of Time Steps", &redtisdyn);
  Core::UTILS::double_parameter("MAXTIME", 4.0, "", &redtisdyn);
  Core::UTILS::double_parameter("NORMAL", 1.0, "", &redtisdyn);
}



void Inpar::ArteryNetwork::set_valid_conditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // 1D-Artery connector condition

  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_connection_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE 1D ARTERY JUNCTION CONDITIONS", "ArtJunctionCond",
          "Artery junction boundary condition", Core::Conditions::ArtJunctionCond, true,
          Core::Conditions::geometry_type_point);

  art_connection_bc->add_component(Teuchos::make_rcp<Input::IntComponent>("ConditionID"));
  art_connection_bc->add_component(Teuchos::make_rcp<Input::RealComponent>("Kr"));

  condlist.push_back(art_connection_bc);

  /*--------------------------------------------------------------------*/
  // Export 1D-Arterial nefrk in gnuplot format

  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_write_gnuplot_c =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN LINE EXPORT 1D-ARTERIAL NETWORK GNUPLOT FORMAT", "ArtWriteGnuplotCond",
          "Artery write gnuplot format condition", Core::Conditions::ArtWriteGnuplotCond, false,
          Core::Conditions::geometry_type_line);

  art_write_gnuplot_c->add_component(Teuchos::make_rcp<Input::IntComponent>("ArteryNumber"));

  condlist.push_back(art_write_gnuplot_c);

  /*--------------------------------------------------------------------*/
  // 1D artery prescribed BC

  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_in_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE 1D ARTERY PRESCRIBED CONDITIONS", "ArtPrescribedCond",
          "Artery prescribed boundary condition", Core::Conditions::ArtPrescribedCond, true,
          Core::Conditions::geometry_type_point);

  art_in_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("boundarycond", "flow",
      Teuchos::tuple<std::string>("flow", "pressure", "velocity", "area", "characteristicWave"),
      Teuchos::tuple<std::string>("flow", "pressure", "velocity", "area", "characteristicWave"),
      true));
  art_in_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("type", "forced",
      Teuchos::tuple<std::string>("forced", "absorbing"),
      Teuchos::tuple<std::string>("forced", "absorbing"), true));

  std::vector<Teuchos::RCP<Input::LineComponent>> artinletcomponents;
  artinletcomponents.push_back(Teuchos::make_rcp<Input::RealVectorComponent>("VAL", 2));
  artinletcomponents.push_back(
      Teuchos::make_rcp<Input::IntVectorComponent>("curve", 2, IntComponentData{0, true, true}));
  for (unsigned i = 0; i < artinletcomponents.size(); ++i)
    art_in_bc->add_component(artinletcomponents[i]);

  condlist.push_back(art_in_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery reflective BC
  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_rf_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE 1D ARTERY REFLECTIVE CONDITIONS", "ArtRfCond", "Artery reflection condition",
          Core::Conditions::ArtRfCond, true, Core::Conditions::geometry_type_point);

  std::vector<Teuchos::RCP<Input::LineComponent>> artrfcomponents;
  artrfcomponents.push_back(Teuchos::make_rcp<Input::RealVectorComponent>("VAL", 1));
  artrfcomponents.push_back(
      Teuchos::make_rcp<Input::IntVectorComponent>("curve", 1, IntComponentData{0, true, true}));
  for (unsigned i = 0; i < artrfcomponents.size(); ++i)
    art_rf_bc->add_component(artrfcomponents[i]);

  condlist.push_back(art_rf_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery windkessel BC
  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_wk_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE 1D ARTERY WINDKESSEL CONDITIONS", "ArtWkCond", "Artery windkessel condition",
          Core::Conditions::ArtWkCond, true, Core::Conditions::geometry_type_point);

  std::vector<Teuchos::RCP<Input::LineComponent>> artwkcomponents;

  art_wk_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("intigrationType",
      "ExplicitWindkessel", Teuchos::tuple<std::string>("ExplicitWindkessel", "ImpedaceWindkessel"),
      Teuchos::tuple<std::string>("ExplicitWindkessel", "ImpedaceWindkessel"), true));

  art_wk_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("windkesselType", "RCR",
      Teuchos::tuple<std::string>("R", "RC", "RCR", "RCRL", "none"),
      Teuchos::tuple<std::string>("R", "RC", "RCR", "RCRL", "none"), true));

  artwkcomponents.push_back(Teuchos::make_rcp<Input::RealVectorComponent>("VAL", 5));
  artwkcomponents.push_back(
      Teuchos::make_rcp<Input::IntVectorComponent>("curve", 5, IntComponentData{0, true, true}));
  for (unsigned i = 0; i < artwkcomponents.size(); ++i)
    art_wk_bc->add_component(artwkcomponents[i]);

  condlist.push_back(art_wk_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery in/out condition

  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_in_outlet_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE 1D ARTERY IN_OUTLET CONDITIONS", "ArtInOutCond",
          "Artery terminal in_outlet condition", Core::Conditions::ArtInOutletCond, true,
          Core::Conditions::geometry_type_point);

  art_in_outlet_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("terminaltype",
      "inlet", Teuchos::tuple<std::string>("inlet", "outlet"),
      Teuchos::tuple<std::string>("inlet", "outlet"), true));

  condlist.push_back(art_in_outlet_bc);
  /*--------------------------------------------------------------------*/
  // 1D artery scalar transport condition
  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_scatra_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE 1D ARTERY SCATRA PRESCRIBED CONDITIONS", "ArtPrescribedScatraCond",
          "Artery prescribed scatra boundary condition", Core::Conditions::ArtPrescribedScatraCond,
          true, Core::Conditions::geometry_type_point);

  std::vector<Teuchos::RCP<Input::LineComponent>> artscatracomponents;
  artscatracomponents.push_back(Teuchos::make_rcp<Input::RealComponent>("VAL"));
  artscatracomponents.push_back(
      Teuchos::make_rcp<Input::IntComponent>("curve", IntComponentData{0, true, true, false}));
  for (unsigned i = 0; i < artscatracomponents.size(); ++i)
    art_scatra_bc->add_component(artscatracomponents[i]);

  condlist.push_back(art_scatra_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery-to-porofluid coupling BC

  Teuchos::RCP<Core::Conditions::ConditionDefinition> artcoup =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE 1D ARTERY TO POROFLUID COUPLING CONDITIONS", "ArtPorofluidCouplConNodebased",
          "Artery coupling with porofluid", Core::Conditions::ArtPorofluidCouplingCondNodebased,
          true, Core::Conditions::geometry_type_point);

  add_named_int(artcoup, "COUPID");

  condlist.push_back(artcoup);

  /*--------------------------------------------------------------------*/
  // 1D artery-to-scatra coupling BC

  Teuchos::RCP<Core::Conditions::ConditionDefinition> artscatracoup =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE 1D ARTERY TO SCATRA COUPLING CONDITIONS", "ArtScatraCouplConNodebased",
          "Artery coupling with porofluid", Core::Conditions::ArtScatraCouplingCondNodebased, true,
          Core::Conditions::geometry_type_point);

  add_named_int(artscatracoup, "COUPID");

  condlist.push_back(artscatracoup);

  /*--------------------------------------------------------------------*/
  // 1D artery-to-porofluid coupling BC Node-To-Point
  Teuchos::RCP<Core::Conditions::ConditionDefinition> artcoup_ntp =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN 1D ARTERY/AIRWAY TO POROFLUID NONCONF COUPLING CONDITIONS",
          "ArtPorofluidCouplConNodeToPoint", "Artery coupling with porofluid nonconf",
          Core::Conditions::ArtPorofluidCouplingCondNodeToPoint, true,
          Core::Conditions::geometry_type_point);

  artcoup_ntp->add_component(Teuchos::make_rcp<Input::SelectionComponent>("coupling_type", "ARTERY",
      Teuchos::tuple<std::string>("ARTERY", "AIRWAY"),
      Teuchos::tuple<std::string>("ARTERY", "AIRWAY"), true));
  Input::add_named_int(artcoup_ntp, "NUMDOF");
  Input::add_named_int_vector(
      artcoup_ntp, "COUPLEDDOF_REDUCED", "coupling dofs of reduced airways or arteries", "NUMDOF");
  Input::add_named_int_vector(
      artcoup_ntp, "COUPLEDDOF_PORO", "coupling dofs in porous domain", "NUMDOF");
  Input::add_named_real(artcoup_ntp, "PENALTY");

  condlist.push_back(artcoup_ntp);

  /*--------------------------------------------------------------------*/
  // 1D artery-to-scatra coupling BC Node-To-Point
  Teuchos::RCP<Core::Conditions::ConditionDefinition> artscatracoup_ntp =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN 1D ARTERY/AIRWAY TO SCATRA NONCONF COUPLING CONDITIONS",
          "ArtScatraCouplConNodeToPoint", "Artery coupling with scatra nonconf",
          Core::Conditions::ArtScatraCouplingCondNodeToPoint, true,
          Core::Conditions::geometry_type_point);

  artscatracoup_ntp->add_component(Teuchos::make_rcp<Input::SelectionComponent>("coupling_type",
      "ARTERY", Teuchos::tuple<std::string>("ARTERY", "AIRWAY"),
      Teuchos::tuple<std::string>("ARTERY", "AIRWAY"), true));
  Input::add_named_int(artscatracoup_ntp, "NUMDOF");
  Input::add_named_int_vector(artscatracoup_ntp, "COUPLEDDOF_REDUCED",
      "coupling dofs of reduced airways or arteries", "NUMDOF");
  Input::add_named_int_vector(
      artscatracoup_ntp, "COUPLEDDOF_PORO", "coupling dofs in porous domain", "NUMDOF");
  Input::add_named_real(artscatracoup_ntp, "PENALTY");

  condlist.push_back(artscatracoup_ntp);
}



void Inpar::BioFilm::set_valid_parameters(Teuchos::ParameterList& list)
{
  Teuchos::ParameterList& biofilmcontrol =
      list.sublist("BIOFILM CONTROL", false, "control parameters for biofilm problems\n");

  Core::UTILS::bool_parameter(
      "BIOFILMGROWTH", "No", "Scatra algorithm for biofilm growth", &biofilmcontrol);
  Core::UTILS::bool_parameter("AVGROWTH", "No",
      "The calculation of growth parameters is based on averaged values", &biofilmcontrol);
  Core::UTILS::double_parameter(
      "FLUXCOEF", 0.0, "Coefficient for growth due to scalar flux", &biofilmcontrol);
  Core::UTILS::double_parameter("NORMFORCEPOSCOEF", 0.0,
      "Coefficient for erosion due to traction normal surface forces", &biofilmcontrol);
  Core::UTILS::double_parameter("NORMFORCENEGCOEF", 0.0,
      "Coefficient for erosion due to compression normal surface forces", &biofilmcontrol);
  Core::UTILS::double_parameter("TANGONEFORCECOEF", 0.0,
      "Coefficient for erosion due to the first tangential surface force", &biofilmcontrol);
  Core::UTILS::double_parameter("TANGTWOFORCECOEF", 0.0,
      "Coefficient for erosion due to the second tangential surface force", &biofilmcontrol);
  Core::UTILS::double_parameter(
      "BIOTIMESTEP", 0.05, "Time step size for biofilm growth", &biofilmcontrol);
  Core::UTILS::int_parameter(
      "BIONUMSTEP", 0, "Maximum number of steps for biofilm growth", &biofilmcontrol);
  Core::UTILS::bool_parameter(
      "OUTPUT_GMSH", "No", "Do you want to write Gmsh postprocessing files?", &biofilmcontrol);
}



void Inpar::BioFilm::set_valid_conditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // Additional coupling for biofilm growth

  std::vector<Teuchos::RCP<Input::LineComponent>> biogrcomponents;

  biogrcomponents.push_back(Teuchos::make_rcp<Input::IntComponent>("coupling id"));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linebiogr =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN BIOFILM GROWTH COUPLING LINE CONDITIONS", "BioGrCoupling", "BioGrCoupling",
          Core::Conditions::BioGrCoupling, true, Core::Conditions::geometry_type_line);
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfbiogr =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN BIOFILM GROWTH COUPLING SURF CONDITIONS", "BioGrCoupling", "BioGrCoupling",
          Core::Conditions::BioGrCoupling, true, Core::Conditions::geometry_type_surface);

  for (unsigned i = 0; i < biogrcomponents.size(); ++i)
  {
    linebiogr->add_component(biogrcomponents[i]);
    surfbiogr->add_component(biogrcomponents[i]);
  }

  condlist.push_back(linebiogr);
  condlist.push_back(surfbiogr);
}


void Inpar::ReducedLung::set_valid_parameters(Teuchos::ParameterList& list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& redawdyn = list.sublist("REDUCED DIMENSIONAL AIRWAYS DYNAMIC", false, "");

  setStringToIntegralParameter<RedAirwaysDyntype>("DYNAMICTYP", "OneStepTheta",
      "OneStepTheta Scheme", tuple<std::string>("OneStepTheta"),
      tuple<RedAirwaysDyntype>(one_step_theta), &redawdyn);

  setStringToIntegralParameter<RedAirwaysDyntype>("SOLVERTYPE", "Linear", "Solver type",
      tuple<std::string>("Linear", "Nonlinear"), tuple<RedAirwaysDyntype>(linear, nonlinear),
      &redawdyn);

  Core::UTILS::double_parameter("TIMESTEP", 0.01, "Time increment dt", &redawdyn);
  Core::UTILS::int_parameter("NUMSTEP", 0, "Number of Time Steps", &redawdyn);
  Core::UTILS::int_parameter("RESTARTEVRY", 1, "Increment for writing restart", &redawdyn);
  Core::UTILS::int_parameter("RESULTSEVRY", 1, "Increment for writing solution", &redawdyn);
  Core::UTILS::double_parameter("THETA", 1.0, "One-step-theta time integration factor", &redawdyn);

  Core::UTILS::int_parameter("MAXITERATIONS", 1, "maximum iteration steps", &redawdyn);
  Core::UTILS::double_parameter("TOLERANCE", 1.0E-6, "tolerance", &redawdyn);

  // number of linear solver used for reduced dimensional airways dynamic
  Core::UTILS::int_parameter("LINEAR_SOLVER", -1,
      "number of linear solver used for reduced dim arterial dynamics", &redawdyn);

  Core::UTILS::bool_parameter(
      "SOLVESCATRA", "no", "Flag to (de)activate solving scalar transport in blood", &redawdyn);

  Core::UTILS::bool_parameter("COMPAWACINTER", "no",
      "Flag to (de)activate computation of airway-acinus interdependency", &redawdyn);

  Core::UTILS::bool_parameter("CALCV0PRESTRESS", "no",
      "Flag to (de)activate initial acini volume adjustment with pre-stress condition", &redawdyn);

  Core::UTILS::double_parameter("TRANSPULMPRESS", 800.0,
      "Transpulmonary pressure needed for recalculation of acini volumes", &redawdyn);
}



void Inpar::ReducedLung::set_valid_conditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // 3-D/reduced-D coupling boundary condition
  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_red_to_3d_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE REDUCED D To 3D FLOW COUPLING CONDITIONS", "Art_redD_3D_CouplingCond",
          "Artery reduced D 3D coupling condition", Core::Conditions::ArtRedTo3DCouplingCond, true,
          Core::Conditions::geometry_type_point);

  art_red_to_3d_bc->add_component(Teuchos::make_rcp<Input::IntComponent>("ConditionID"));

  art_red_to_3d_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("CouplingType",
      "forced", Teuchos::tuple<std::string>("forced", "absorbing"),
      Teuchos::tuple<std::string>("forced", "absorbing"), true));

  art_red_to_3d_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("ReturnedVariable",
      "pressure", Teuchos::tuple<std::string>("pressure", "flow"),
      Teuchos::tuple<std::string>("pressure", "flow"), true));
  Input::add_named_real(art_red_to_3d_bc, "Tolerance");
  Input::add_named_int(art_red_to_3d_bc, "MaximumIterations");

  condlist.push_back(art_red_to_3d_bc);

  /*--------------------------------------------------------------------*/
  // 3-D/reduced-D coupling boundary condition
  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_3d_to_red_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN SURF 3D To REDUCED D FLOW COUPLING CONDITIONS", "Art_3D_redD_CouplingCond",
          "Artery 3D reduced D coupling condition", Core::Conditions::Art3DToRedCouplingCond, true,
          Core::Conditions::geometry_type_surface);

  art_3d_to_red_bc->add_component(Teuchos::make_rcp<Input::IntComponent>("ConditionID"));

  art_3d_to_red_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("ReturnedVariable",
      "flow", Teuchos::tuple<std::string>("pressure", "flow"),
      Teuchos::tuple<std::string>("pressure", "flow"), true));
  Input::add_named_real(art_3d_to_red_bc, "Tolerance");
  Input::add_named_int(art_3d_to_red_bc, "MaximumIterations");

  condlist.push_back(art_3d_to_red_bc);

  /*--------------------------------------------------------------------*/
  // Coupling of 3D tissue models and reduced-D airway tree

  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfredairtis =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN SURF TISSUE REDAIRWAY CONDITIONS", "SurfaceNeumann",
          "tissue RedAirway coupling surface condition", Core::Conditions::RedAirwayTissue, true,
          Core::Conditions::geometry_type_surface);

  surfredairtis->add_component(Teuchos::make_rcp<Input::IntComponent>("coupling id"));

  condlist.push_back(surfredairtis);


  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> noderedairtis =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE TISSUE REDAIRWAY CONDITIONS", "RedAirwayPrescribedCond",
          "tissue RedAirway coupling node condition", Core::Conditions::RedAirwayNodeTissue, true,
          Core::Conditions::geometry_type_point);


  noderedairtis->add_component(Teuchos::make_rcp<Input::IntComponent>("coupling id"));

  condlist.push_back(noderedairtis);



  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_in_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE Reduced D AIRWAYS PRESCRIBED CONDITIONS", "RedAirwayPrescribedCond",
          "Reduced d airway prescribed boundary condition",
          Core::Conditions::RedAirwayPrescribedCond, true, Core::Conditions::geometry_type_point);

  raw_in_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("boundarycond", "flow",
      Teuchos::tuple<std::string>(
          "flow", "pressure", "switchFlowPressure", "VolumeDependentPleuralPressure"),
      Teuchos::tuple<std::string>(
          "flow", "pressure", "switchFlowPressure", "VolumeDependentPleuralPressure"),
      true));

  // reduced airway inlet components
  raw_in_bc->add_component(Teuchos::make_rcp<Input::RealVectorComponent>("VAL", 1));
  raw_in_bc->add_component(
      Teuchos::make_rcp<Input::IntVectorComponent>("curve", 2, IntComponentData{0, true, true}));
  raw_in_bc->add_component(Teuchos::make_rcp<Input::IntVectorComponent>(
      "funct", 1, IntComponentData{0, false, true, true}));

  condlist.push_back(raw_in_bc);

  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways switching between different types of boundary
  // conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_in_switch_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE Reduced D AIRWAYS SWITCH FLOW PRESSURE CONDITIONS",
          "RedAirwaySwitchFlowPressureCond",
          "Reduced d airway switch flow pressure boundary condition",
          Core::Conditions::RedAirwayPrescribedSwitchCond, true,
          Core::Conditions::geometry_type_point);

  Input::add_named_int(raw_in_switch_bc, "FUNCT_ID_FLOW");
  Input::add_named_int(raw_in_switch_bc, "FUNCT_ID_PRESSURE");
  Input::add_named_int(raw_in_switch_bc, "FUNCT_ID_PRESSURE_ACTIVE");

  condlist.push_back(raw_in_switch_bc);

  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways external pressure

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_pext_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN LINE Reduced D AIRWAYS PRESCRIBED EXTERNAL PRESSURE CONDITIONS",
          "RedAirwayPrescribedExternalPressure",
          "Reduced d airway prescribed external pressure boundary condition",
          Core::Conditions::RedAirwayPrescribedExternalPressure, true,
          Core::Conditions::geometry_type_line);

  raw_pext_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("boundarycond",
      "ExternalPressure", Teuchos::tuple<std::string>("ExternalPressure"),
      Teuchos::tuple<std::string>("ExternalPressure"), true));

  // reduced airway pext components
  raw_pext_bc->add_component(Teuchos::make_rcp<Input::RealVectorComponent>("VAL", 1));
  raw_pext_bc->add_component(
      Teuchos::make_rcp<Input::IntVectorComponent>("curve", 2, IntComponentData{0, true, true}));
  raw_pext_bc->add_component(Teuchos::make_rcp<Input::IntVectorComponent>(
      "funct", 1, IntComponentData{0, false, false, true}));

  condlist.push_back(raw_pext_bc);


  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional scalar transport in airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_in_scatra_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE Reduced D AIRWAYS PRESCRIBED SCATRA CONDITIONS",
          "RedAirwayPrescribedScatraCond", "Reduced d airway prescribed scatra boundary condition",
          Core::Conditions::RedAirwayPrescribedScatraCond, true,
          Core::Conditions::geometry_type_point);

  // reduced airway inlet scatra components
  raw_in_scatra_bc->add_component(Teuchos::make_rcp<Input::RealVectorComponent>("VAL", 1));
  raw_in_scatra_bc->add_component(
      Teuchos::make_rcp<Input::IntVectorComponent>("curve", 1, IntComponentData{0, true, true}));
  raw_in_scatra_bc->add_component(Teuchos::make_rcp<Input::IntVectorComponent>(
      "funct", 1, IntComponentData{0, false, false, true}));

  condlist.push_back(raw_in_scatra_bc);

  /*--------------------------------------------------------------------*/
  // Prescribed BC for initial values of the scalar transport in reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_int_scatra_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN LINE Reduced D AIRWAYS INITIAL SCATRA CONDITIONS", "RedAirwayInitialScatraCond",
          "Reduced d airway initial scatra boundary condition",
          Core::Conditions::RedAirwayInitialScatraCond, true, Core::Conditions::geometry_type_line);

  raw_int_scatra_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("scalar", "O2",
      Teuchos::tuple<std::string>("O2", "CO2"), Teuchos::tuple<std::string>("O2", "CO2"), true));
  Input::add_named_real(raw_int_scatra_bc, "CONCENTRATION");

  condlist.push_back(raw_int_scatra_bc);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions of scatra exchange
  Teuchos::RCP<Core::Conditions::ConditionDefinition> scatra_exchange_cond =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN LINE Reduced D AIRWAYS SCATRA EXCHANGE CONDITIONS", "RedAirwayScatraExchangeCond",
          "scatra exchange condition", Core::Conditions::RedAirwayScatraExchangeCond, true,
          Core::Conditions::geometry_type_line);

  scatra_exchange_cond->add_component(Teuchos::make_rcp<Input::IntComponent>("ConditionID"));

  condlist.push_back(scatra_exchange_cond);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions with hemoglobin
  Teuchos::RCP<Core::Conditions::ConditionDefinition> scatra_hemoglobin_cond =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN LINE Reduced D AIRWAYS HEMOGLOBIN CONDITIONS", "RedAirwayScatraHemoglobinCond",
          "scatra hemoglobin condition", Core::Conditions::RedAirwayScatraHemoglobinCond, false,
          Core::Conditions::geometry_type_line);

  Input::add_named_real(scatra_hemoglobin_cond, "INITIAL_CONCENTRATION");

  condlist.push_back(scatra_hemoglobin_cond);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions with hemoglobin
  Teuchos::RCP<Core::Conditions::ConditionDefinition> scatra_air_cond =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN LINE Reduced D AIRWAYS AIR CONDITIONS", "RedAirwayScatraAirCond",
          "scatra air condition", Core::Conditions::RedAirwayScatraAirCond, false,
          Core::Conditions::geometry_type_line);

  Input::add_named_real(scatra_air_cond, "INITIAL_CONCENTRATION");

  condlist.push_back(scatra_air_cond);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions with hemoglobin
  Teuchos::RCP<Core::Conditions::ConditionDefinition> scatra_capillary_cond =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN LINE Reduced D AIRWAYS CAPILLARY CONDITIONS", "RedAirwayScatraCapillaryCond",
          "scatra capillary condition", Core::Conditions::RedAirwayScatraCapillaryCond, false,
          Core::Conditions::geometry_type_line);

  condlist.push_back(scatra_capillary_cond);

  /*--------------------------------------------------------------------*/
  // Prescribed Ventilator BC for reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_vent_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN NODE Reduced D AIRWAYS VENTILATOR CONDITIONS", "RedAirwayVentilatorCond",
          "Reduced d airway prescribed ventilator condition",
          Core::Conditions::RedAirwayVentilatorCond, true, Core::Conditions::geometry_type_point);

  raw_vent_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("phase1", "flow",
      Teuchos::tuple<std::string>("flow", "volume", "pressure"),
      Teuchos::tuple<std::string>("flow", "volume", "pressure"), true));

  raw_vent_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("Phase1Smoothness",
      "smooth", Teuchos::tuple<std::string>("smooth", "discontinous"),
      Teuchos::tuple<std::string>("smooth", "discontinous"), true));

  raw_vent_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("phase2", "pressure",
      Teuchos::tuple<std::string>("pressure", "flow", "volume"),
      Teuchos::tuple<std::string>("pressure", "flow", "volume"), true));

  raw_vent_bc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("Phase2Smoothness",
      "smooth", Teuchos::tuple<std::string>("smooth", "discontinous"),
      Teuchos::tuple<std::string>("smooth", "discontinous"), true));

  Input::add_named_real(raw_vent_bc, "period");
  Input::add_named_real(raw_vent_bc, "phase1_period");
  Input::add_named_real(raw_vent_bc, "smoothness_period1");
  Input::add_named_real(raw_vent_bc, "smoothness_period2");

  // reduced airway ventilation components
  raw_vent_bc->add_component(Teuchos::make_rcp<Input::RealVectorComponent>("VAL", 2));
  raw_vent_bc->add_component(
      Teuchos::make_rcp<Input::IntVectorComponent>("curve", 2, IntComponentData{0, true, true}));

  condlist.push_back(raw_vent_bc);



  /*--------------------------------------------------------------------*/
  // Prescribed volume dependent pleural pressure for reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_volPpl_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN LINE REDUCED D AIRWAYS VOL DEPENDENT PLEURAL PRESSURE CONDITIONS",
          "RedAirwayVolDependentPleuralPressureCond",
          "Reduced D airways volume-dependent peural pressure condition",
          Core::Conditions::RedAirwayVolDependentPleuralPressureCond, true,
          Core::Conditions::geometry_type_line);

  raw_volPpl_bc->add_component(
      Teuchos::make_rcp<Input::SelectionComponent>("TYPE", "Linear_Exponential",
          Teuchos::tuple<std::string>("Linear_Polynomial", "Linear_Exponential", "Linear_Ogden",
              "Nonlinear_Polynomial", "Nonlinear_Exponential", "Nonlinear_Ogden"),
          Teuchos::tuple<std::string>("Linear_Polynomial", "Linear_Exponential", "Linear_Ogden",
              "Nonlinear_Polynomial", "Nonlinear_Exponential", "Nonlinear_Ogden"),
          true));

  Input::add_named_real(raw_volPpl_bc, "TLC");
  Input::add_named_real(raw_volPpl_bc, "RV");

  Input::add_named_real(raw_volPpl_bc, "P_PLEURAL_0");
  Input::add_named_real(raw_volPpl_bc, "P_PLEURAL_LIN");
  Input::add_named_real(raw_volPpl_bc, "P_PLEURAL_NONLIN");
  Input::add_named_real(raw_volPpl_bc, "TAU");

  // raw_volPpl_bc_components
  raw_volPpl_bc->add_component(Teuchos::make_rcp<Input::RealVectorComponent>("VAL", 1));
  raw_volPpl_bc->add_component(
      Teuchos::make_rcp<Input::IntVectorComponent>("curve", 1, IntComponentData{0, true, true}));

  condlist.push_back(raw_volPpl_bc);

  /*--------------------------------------------------------------------*/
  // Evaluate lung volume condition for reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_eval_lungV_bc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>(
          "DESIGN LINE REDUCED D AIRWAYS EVALUATE LUNG VOLUME CONDITIONS",
          "RedAirwayEvalLungVolCond", "Reduced D airways evaluate lung volume condition",
          Core::Conditions::RedAirwayEvalLungVolCond, true, Core::Conditions::geometry_type_line);


  condlist.push_back(raw_eval_lungV_bc);


  /*--------------------------------------------------------------------*/
  // Impedance condition

  Teuchos::RCP<Core::Conditions::ConditionDefinition> impedancebc =
      Teuchos::make_rcp<Core::Conditions::ConditionDefinition>("DESIGN SURF IMPEDANCE CONDITIONS",
          "ImpedanceCond", "Impedance boundary condition", Core::Conditions::ImpedanceCond, true,
          Core::Conditions::geometry_type_surface);

  impedancebc->add_component(Teuchos::make_rcp<Input::IntComponent>("ConditionID"));

  impedancebc->add_component(Teuchos::make_rcp<Input::SelectionComponent>("TYPE", "windkessel",
      Teuchos::tuple<std::string>("windkessel", "resistive", "pressure_by_funct"),
      Teuchos::tuple<std::string>("windkessel", "resistive", "pressure_by_funct"), true));
  Input::add_named_real(impedancebc, "R1");
  Input::add_named_real(impedancebc, "R2");
  Input::add_named_real(impedancebc, "C");
  Input::add_named_real(impedancebc, "TIMEPERIOD");
  Input::add_named_int(impedancebc, "FUNCT");

  condlist.push_back(impedancebc);
}

FOUR_C_NAMESPACE_CLOSE
