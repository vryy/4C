/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for biomedical simulations

\level 3


*/
/*----------------------------------------------------------------------*/



#include "4C_inpar_bio.hpp"

#include "4C_discretization_condition_definition.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_linecomponent.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


void Inpar::ArtDyn::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;
  Teuchos::ParameterList& andyn = list->sublist("ARTERIAL DYNAMIC", false, "");

  setStringToIntegralParameter<int>("DYNAMICTYP", "ExpTaylorGalerkin",
      "Explicit Taylor Galerkin Scheme", tuple<std::string>("ExpTaylorGalerkin", "Stationary"),
      tuple<int>(tay_gal, stationary), &andyn);

  Core::UTILS::DoubleParameter("TIMESTEP", 0.01, "Time increment dt", &andyn);
  Core::UTILS::IntParameter("NUMSTEP", 0, "Number of Time Steps", &andyn);
  Core::UTILS::DoubleParameter("MAXTIME", 1000.0, "total simulation time", &andyn);
  Core::UTILS::IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &andyn);
  Core::UTILS::IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &andyn);
  setStringToIntegralParameter<int>("SOLVESCATRA", "no",
      "Flag to (de)activate solving scalar transport in blood", tuple<std::string>("no", "yes"),
      tuple<std::string>("do not solve scatra", "solve scatra"), tuple<int>(0, 1), &andyn);

  // number of linear solver used for arterial dynamics
  Core::UTILS::IntParameter(
      "LINEAR_SOLVER", -1, "number of linear solver used for arterial dynamics", &andyn);

  // initial function number
  Core::UTILS::IntParameter("INITFUNCNO", -1, "function number for artery initial field", &andyn);

  // type of initial field
  setStringToIntegralParameter<int>("INITIALFIELD", "zero_field",
      "Initial Field for artery problem",
      tuple<std::string>("zero_field", "field_by_function", "field_by_condition"),
      tuple<int>(initfield_zero_field, initfield_field_by_function, initfield_field_by_condition),
      &andyn);
}



void Inpar::ArteryNetwork::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& redtisdyn =
      list->sublist("COUPLED REDUCED-D AIRWAYS AND TISSUE DYNAMIC", false, "");
  Core::UTILS::DoubleParameter("CONVTOL_P", 1E-6,
      "Coupled red_airway and tissue iteration convergence for pressure", &redtisdyn);
  Core::UTILS::DoubleParameter("CONVTOL_Q", 1E-6,
      "Coupled red_airway and tissue iteration convergence for flux", &redtisdyn);
  Core::UTILS::IntParameter("MAXITER", 5, "Maximum coupling iterations", &redtisdyn);
  setStringToIntegralParameter<int>("RELAXTYPE", "norelaxation", "Dynamic Relaxation Type",
      tuple<std::string>("norelaxation", "fixedrelaxation", "Aitken", "SD"),
      tuple<int>(norelaxation, fixedrelaxation, Aitken, SD), &redtisdyn);
  Core::UTILS::DoubleParameter("TIMESTEP", 0.01, "Time increment dt", &redtisdyn);
  Core::UTILS::IntParameter("NUMSTEP", 1, "Number of Time Steps", &redtisdyn);
  Core::UTILS::DoubleParameter("MAXTIME", 4.0, "", &redtisdyn);
  Core::UTILS::DoubleParameter("NORMAL", 1.0, "", &redtisdyn);
}



void Inpar::ArteryNetwork::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // 1D-Artery connector condition

  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_connection_bc = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN NODE 1D ARTERY JUNCTION CONDITIONS",
          "ArtJunctionCond", "Artery junction boundary condition",
          Core::Conditions::ArtJunctionCond, true, Core::Conditions::geometry_type_point));

  art_connection_bc->AddComponent(Teuchos::rcp(new Input::IntComponent("ConditionID")));
  art_connection_bc->AddComponent(Teuchos::rcp(new Input::RealComponent("Kr")));

  condlist.push_back(art_connection_bc);

  /*--------------------------------------------------------------------*/
  // Export 1D-Arterial nefrk in gnuplot format

  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_write_gnuplot_c =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE EXPORT 1D-ARTERIAL NETWORK GNUPLOT FORMAT", "ArtWriteGnuplotCond",
          "Artery write gnuplot format condition", Core::Conditions::ArtWriteGnuplotCond, false,
          Core::Conditions::geometry_type_line));

  art_write_gnuplot_c->AddComponent(Teuchos::rcp(new Input::IntComponent("ArteryNumber")));

  condlist.push_back(art_write_gnuplot_c);

  /*--------------------------------------------------------------------*/
  // 1D artery prescribed BC

  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_in_bc = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN NODE 1D ARTERY PRESCRIBED CONDITIONS",
          "ArtPrescribedCond", "Artery prescribed boundary condition",
          Core::Conditions::ArtPrescribedCond, true, Core::Conditions::geometry_type_point));

  art_in_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("boundarycond", "flow",
      Teuchos::tuple<std::string>("flow", "pressure", "velocity", "area", "characteristicWave"),
      Teuchos::tuple<std::string>("flow", "pressure", "velocity", "area", "characteristicWave"),
      true)));
  art_in_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("type", "forced",
      Teuchos::tuple<std::string>("forced", "absorbing"),
      Teuchos::tuple<std::string>("forced", "absorbing"), true)));

  std::vector<Teuchos::RCP<Input::LineComponent>> artinletcomponents;
  artinletcomponents.push_back(Teuchos::rcp(new Input::RealVectorComponent("val", 2)));
  artinletcomponents.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("curve", 2, {0, true, true})));
  for (unsigned i = 0; i < artinletcomponents.size(); ++i)
    art_in_bc->AddComponent(artinletcomponents[i]);

  condlist.push_back(art_in_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery reflective BC
  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_rf_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN NODE 1D ARTERY REFLECTIVE CONDITIONS", "ArtRfCond", "Artery reflection condition",
          Core::Conditions::ArtRfCond, true, Core::Conditions::geometry_type_point));

  std::vector<Teuchos::RCP<Input::LineComponent>> artrfcomponents;
  artrfcomponents.push_back(Teuchos::rcp(new Input::RealVectorComponent("val", 1)));
  artrfcomponents.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("curve", 1, {0, true, true})));
  for (unsigned i = 0; i < artrfcomponents.size(); ++i) art_rf_bc->AddComponent(artrfcomponents[i]);

  condlist.push_back(art_rf_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery windkessel BC
  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_wk_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN NODE 1D ARTERY WINDKESSEL CONDITIONS", "ArtWkCond", "Artery windkessel condition",
          Core::Conditions::ArtWkCond, true, Core::Conditions::geometry_type_point));

  std::vector<Teuchos::RCP<Input::LineComponent>> artwkcomponents;

  art_wk_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("intigrationType",
      "ExplicitWindkessel", Teuchos::tuple<std::string>("ExplicitWindkessel", "ImpedaceWindkessel"),
      Teuchos::tuple<std::string>("ExplicitWindkessel", "ImpedaceWindkessel"), true)));

  art_wk_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("windkesselType", "RCR",
      Teuchos::tuple<std::string>("R", "RC", "RCR", "RCRL", "none"),
      Teuchos::tuple<std::string>("R", "RC", "RCR", "RCRL", "none"), true)));

  artwkcomponents.push_back(Teuchos::rcp(new Input::RealVectorComponent("val", 5)));
  artwkcomponents.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("curve", 5, {0, true, true})));
  for (unsigned i = 0; i < artwkcomponents.size(); ++i) art_wk_bc->AddComponent(artwkcomponents[i]);

  condlist.push_back(art_wk_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery in/out condition

  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_in_outlet_bc = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN NODE 1D ARTERY IN_OUTLET CONDITIONS",
          "ArtInOutCond", "Artery terminal in_outlet condition", Core::Conditions::ArtInOutletCond,
          true, Core::Conditions::geometry_type_point));

  art_in_outlet_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("terminaltype", "inlet",
      Teuchos::tuple<std::string>("inlet", "outlet"),
      Teuchos::tuple<std::string>("inlet", "outlet"), true)));

  condlist.push_back(art_in_outlet_bc);
  /*--------------------------------------------------------------------*/
  // 1D artery scalar transport condition
  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_scatra_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN NODE 1D ARTERY SCATRA PRESCRIBED CONDITIONS", "ArtPrescribedScatraCond",
          "Artery prescribed scatra boundary condition", Core::Conditions::ArtPrescribedScatraCond,
          true, Core::Conditions::geometry_type_point));

  std::vector<Teuchos::RCP<Input::LineComponent>> artscatracomponents;
  artscatracomponents.push_back(Teuchos::rcp(new Input::RealComponent("val")));
  artscatracomponents.push_back(Teuchos::rcp(new Input::IntComponent("curve", {0, true, true})));
  for (unsigned i = 0; i < artscatracomponents.size(); ++i)
    art_scatra_bc->AddComponent(artscatracomponents[i]);

  condlist.push_back(art_scatra_bc);

  /*--------------------------------------------------------------------*/
  // 1D artery-to-porofluid coupling BC
  std::vector<Teuchos::RCP<Input::LineComponent>> artcoupcomponents;

  Teuchos::RCP<Core::Conditions::ConditionDefinition> artcoup =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN NODE 1D ARTERY TO POROFLUID COUPLING CONDITIONS", "ArtPorofluidCouplConNodebased",
          "Artery coupling with porofluid", Core::Conditions::ArtPorofluidCouplingCondNodebased,
          true, Core::Conditions::geometry_type_point));

  artcoupcomponents.push_back(Teuchos::rcp(new Input::SeparatorComponent("COUPID")));
  artcoupcomponents.push_back(Teuchos::rcp(new Input::IntComponent("coupling id")));

  for (unsigned i = 0; i < artcoupcomponents.size(); ++i)
  {
    artcoup->AddComponent(artcoupcomponents[i]);
  }
  condlist.push_back(artcoup);

  /*--------------------------------------------------------------------*/
  // 1D artery-to-scatra coupling BC
  std::vector<Teuchos::RCP<Input::LineComponent>> artscatracoupcomponents;

  Teuchos::RCP<Core::Conditions::ConditionDefinition> artscatracoup =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN NODE 1D ARTERY TO SCATRA COUPLING CONDITIONS", "ArtScatraCouplConNodebased",
          "Artery coupling with porofluid", Core::Conditions::ArtScatraCouplingCondNodebased, true,
          Core::Conditions::geometry_type_point));

  artscatracoupcomponents.push_back(Teuchos::rcp(new Input::SeparatorComponent("COUPID")));
  artscatracoupcomponents.push_back(Teuchos::rcp(new Input::IntComponent("coupling id")));

  for (unsigned i = 0; i < artscatracoupcomponents.size(); ++i)
  {
    artscatracoup->AddComponent(artscatracoupcomponents[i]);
  }

  condlist.push_back(artscatracoup);

  /*--------------------------------------------------------------------*/
  // 1D artery-to-porofluid coupling BC Node-To-Point
  Teuchos::RCP<Core::Conditions::ConditionDefinition> artcoup_ntp =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN 1D ARTERY/AIRWAY TO POROFLUID NONCONF COUPLING CONDITIONS",
          "ArtPorofluidCouplConNodeToPoint", "Artery coupling with porofluid nonconf",
          Core::Conditions::ArtPorofluidCouplingCondNodeToPoint, true,
          Core::Conditions::geometry_type_point));

  artcoup_ntp->AddComponent(Teuchos::rcp(new Input::SelectionComponent("coupling_type", "ARTERY",
      Teuchos::tuple<std::string>("ARTERY", "AIRWAY"),
      Teuchos::tuple<std::string>("ARTERY", "AIRWAY"), true)));
  Input::AddNamedInt(artcoup_ntp, "NUMDOF");
  Input::AddNamedIntVector(
      artcoup_ntp, "COUPLEDDOF_REDUCED", "coupling dofs of reduced airways or arteries", "NUMDOF");
  Input::AddNamedIntVector(
      artcoup_ntp, "COUPLEDDOF_PORO", "coupling dofs in porous domain", "NUMDOF");
  Input::AddNamedReal(artcoup_ntp, "PENALTY");

  condlist.push_back(artcoup_ntp);

  /*--------------------------------------------------------------------*/
  // 1D artery-to-scatra coupling BC Node-To-Point
  Teuchos::RCP<Core::Conditions::ConditionDefinition> artscatracoup_ntp =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN 1D ARTERY/AIRWAY TO SCATRA NONCONF COUPLING CONDITIONS",
          "ArtScatraCouplConNodeToPoint", "Artery coupling with scatra nonconf",
          Core::Conditions::ArtScatraCouplingCondNodeToPoint, true,
          Core::Conditions::geometry_type_point));

  artscatracoup_ntp->AddComponent(Teuchos::rcp(new Input::SelectionComponent("coupling_type",
      "ARTERY", Teuchos::tuple<std::string>("ARTERY", "AIRWAY"),
      Teuchos::tuple<std::string>("ARTERY", "AIRWAY"), true)));
  Input::AddNamedInt(artscatracoup_ntp, "NUMDOF");
  Input::AddNamedIntVector(artscatracoup_ntp, "COUPLEDDOF_REDUCED",
      "coupling dofs of reduced airways or arteries", "NUMDOF");
  Input::AddNamedIntVector(
      artscatracoup_ntp, "COUPLEDDOF_PORO", "coupling dofs in porous domain", "NUMDOF");
  Input::AddNamedReal(artscatracoup_ntp, "PENALTY");

  condlist.push_back(artscatracoup_ntp);
}



void Inpar::BioFilm::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  Teuchos::ParameterList& biofilmcontrol =
      list->sublist("BIOFILM CONTROL", false, "control parameters for biofilm problems\n");

  Core::UTILS::BoolParameter(
      "BIOFILMGROWTH", "No", "Scatra algorithm for biofilm growth", &biofilmcontrol);
  Core::UTILS::BoolParameter("AVGROWTH", "No",
      "The calculation of growth parameters is based on averaged values", &biofilmcontrol);
  Core::UTILS::DoubleParameter(
      "FLUXCOEF", 0.0, "Coefficient for growth due to scalar flux", &biofilmcontrol);
  Core::UTILS::DoubleParameter("NORMFORCEPOSCOEF", 0.0,
      "Coefficient for erosion due to traction normal surface forces", &biofilmcontrol);
  Core::UTILS::DoubleParameter("NORMFORCENEGCOEF", 0.0,
      "Coefficient for erosion due to compression normal surface forces", &biofilmcontrol);
  Core::UTILS::DoubleParameter("TANGONEFORCECOEF", 0.0,
      "Coefficient for erosion due to the first tangential surface force", &biofilmcontrol);
  Core::UTILS::DoubleParameter("TANGTWOFORCECOEF", 0.0,
      "Coefficient for erosion due to the second tangential surface force", &biofilmcontrol);
  Core::UTILS::DoubleParameter(
      "BIOTIMESTEP", 0.05, "Time step size for biofilm growth", &biofilmcontrol);
  Core::UTILS::IntParameter(
      "BIONUMSTEP", 0, "Maximum number of steps for biofilm growth", &biofilmcontrol);
  Core::UTILS::BoolParameter(
      "OUTPUT_GMSH", "No", "Do you want to write Gmsh postprocessing files?", &biofilmcontrol);
}



void Inpar::BioFilm::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // Additional coupling for biofilm growth

  std::vector<Teuchos::RCP<Input::LineComponent>> biogrcomponents;

  biogrcomponents.push_back(Teuchos::rcp(new Input::IntComponent("coupling id")));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linebiogr =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN BIOFILM GROWTH COUPLING LINE CONDITIONS", "BioGrCoupling", "BioGrCoupling",
          Core::Conditions::BioGrCoupling, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfbiogr =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN BIOFILM GROWTH COUPLING SURF CONDITIONS", "BioGrCoupling", "BioGrCoupling",
          Core::Conditions::BioGrCoupling, true, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < biogrcomponents.size(); ++i)
  {
    linebiogr->AddComponent(biogrcomponents[i]);
    surfbiogr->AddComponent(biogrcomponents[i]);
  }

  condlist.push_back(linebiogr);
  condlist.push_back(surfbiogr);
}


void Inpar::ReducedLung::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace Input;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& redawdyn =
      list->sublist("REDUCED DIMENSIONAL AIRWAYS DYNAMIC", false, "");

  setStringToIntegralParameter<int>("DYNAMICTYP", "OneStepTheta", "OneStepTheta Scheme",
      tuple<std::string>("OneStepTheta"), tuple<int>(one_step_theta), &redawdyn);

  setStringToIntegralParameter<int>("SOLVERTYPE", "Linear", "Solver type",
      tuple<std::string>("Linear", "Nonlinear"), tuple<int>(linear, nonlinear), &redawdyn);

  Core::UTILS::DoubleParameter("TIMESTEP", 0.01, "Time increment dt", &redawdyn);
  Core::UTILS::IntParameter("NUMSTEP", 0, "Number of Time Steps", &redawdyn);
  Core::UTILS::IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &redawdyn);
  Core::UTILS::IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &redawdyn);
  Core::UTILS::DoubleParameter("THETA", 1.0, "One-step-theta time integration factor", &redawdyn);

  Core::UTILS::IntParameter("MAXITERATIONS", 1, "maximum iteration steps", &redawdyn);
  Core::UTILS::DoubleParameter("TOLERANCE", 1.0E-6, "tolerance", &redawdyn);

  // number of linear solver used for reduced dimensional airways dynamic
  Core::UTILS::IntParameter("LINEAR_SOLVER", -1,
      "number of linear solver used for reduced dim arterial dynamics", &redawdyn);

  // Solve scatra flag
  setStringToIntegralParameter<int>("SOLVESCATRA", "no",
      "Flag to (de)activate solving scalar transport in blood", tuple<std::string>("no", "yes"),
      tuple<std::string>("do not solve scatra", "solve scatra"), tuple<int>(0, 1), &redawdyn);
  // Compute airway-acinus interdependency flag
  setStringToIntegralParameter<int>("COMPAWACINTER", "no",
      "Flag to (de)activate computation of airway-acinus interdependency ",
      tuple<std::string>("no", "yes"),
      tuple<std::string>("do not compute interdependency", "compute interdependency"),
      tuple<int>(0, 1), &redawdyn);
  // Re-calculate initial acini volume flag
  setStringToIntegralParameter<int>("CALCV0PRESTRESS", "no",
      "Flag to (de)activate initial acini volume adjustment with pre-stress condition ",
      tuple<std::string>("no", "yes"), tuple<std::string>("do not adjust", "adjust volumes"),
      tuple<int>(0, 1), &redawdyn);
  Core::UTILS::DoubleParameter("TRANSPULMPRESS", 800.0,
      "Transpulmonary pressure needed for recalculation of acini volumes", &redawdyn);
}



void Inpar::ReducedLung::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  /*--------------------------------------------------------------------*/
  // 3-D/reduced-D coupling boundary condition
  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_red_to_3d_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN NODE REDUCED D To 3D FLOW COUPLING CONDITIONS", "Art_redD_3D_CouplingCond",
          "Artery reduced D 3D coupling condition", Core::Conditions::ArtRedTo3DCouplingCond, true,
          Core::Conditions::geometry_type_point));

  art_red_to_3d_bc->AddComponent(Teuchos::rcp(new Input::IntComponent("ConditionID")));

  art_red_to_3d_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("CouplingType",
      "forced", Teuchos::tuple<std::string>("forced", "absorbing"),
      Teuchos::tuple<std::string>("forced", "absorbing"), true)));

  art_red_to_3d_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("ReturnedVariable",
      "pressure", Teuchos::tuple<std::string>("pressure", "flow"),
      Teuchos::tuple<std::string>("pressure", "flow"), true)));
  Input::AddNamedReal(art_red_to_3d_bc, "Tolerance");
  Input::AddNamedInt(art_red_to_3d_bc, "MaximumIterations");

  condlist.push_back(art_red_to_3d_bc);

  /*--------------------------------------------------------------------*/
  // 3-D/reduced-D coupling boundary condition
  Teuchos::RCP<Core::Conditions::ConditionDefinition> art_3d_to_red_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN SURF 3D To REDUCED D FLOW COUPLING CONDITIONS", "Art_3D_redD_CouplingCond",
          "Artery 3D reduced D coupling condition", Core::Conditions::Art3DToRedCouplingCond, true,
          Core::Conditions::geometry_type_surface));

  art_3d_to_red_bc->AddComponent(Teuchos::rcp(new Input::IntComponent("ConditionID")));

  art_3d_to_red_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("ReturnedVariable",
      "flow", Teuchos::tuple<std::string>("pressure", "flow"),
      Teuchos::tuple<std::string>("pressure", "flow"), true)));
  Input::AddNamedReal(art_3d_to_red_bc, "Tolerance");
  Input::AddNamedInt(art_3d_to_red_bc, "MaximumIterations");

  condlist.push_back(art_3d_to_red_bc);

  /*--------------------------------------------------------------------*/
  // Coupling of 3D tissue models and reduced-D airway tree

  std::vector<Teuchos::RCP<Input::LineComponent>> redairtiscomponents;

  redairtiscomponents.push_back(Teuchos::rcp(new Input::IntComponent("coupling id")));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfredairtis = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN SURF TISSUE REDAIRWAY CONDITIONS",
          "SurfaceNeumann", "tissue RedAirway coupling surface condition",
          Core::Conditions::RedAirwayTissue, true, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < redairtiscomponents.size(); ++i)
  {
    surfredairtis->AddComponent(redairtiscomponents[i]);
  }

  condlist.push_back(surfredairtis);


  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways

  std::vector<Teuchos::RCP<Input::LineComponent>> noderedairtiscomponents;

  noderedairtiscomponents.push_back(Teuchos::rcp(new Input::IntComponent("coupling id")));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> noderedairtis = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN NODE TISSUE REDAIRWAY CONDITIONS",
          "RedAirwayPrescribedCond", "tissue RedAirway coupling node condition",
          Core::Conditions::RedAirwayNodeTissue, true, Core::Conditions::geometry_type_point));


  for (unsigned i = 0; i < noderedairtiscomponents.size(); ++i)
  {
    noderedairtis->AddComponent(noderedairtiscomponents[i]);
  }

  condlist.push_back(noderedairtis);



  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_in_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN NODE Reduced D AIRWAYS PRESCRIBED CONDITIONS", "RedAirwayPrescribedCond",
          "Reduced d airway prescribed boundary condition",
          Core::Conditions::RedAirwayPrescribedCond, true, Core::Conditions::geometry_type_point));

  raw_in_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("boundarycond", "flow",
      Teuchos::tuple<std::string>(
          "flow", "pressure", "switchFlowPressure", "VolumeDependentPleuralPressure"),
      Teuchos::tuple<std::string>(
          "flow", "pressure", "switchFlowPressure", "VolumeDependentPleuralPressure"),
      true)));

  std::vector<Teuchos::RCP<Input::LineComponent>> redairwayinletcomponents;
  redairwayinletcomponents.push_back(Teuchos::rcp(new Input::RealVectorComponent("val", 1)));
  redairwayinletcomponents.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("curve", 2, {0, true, true})));
  redairwayinletcomponents.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("funct", 1, {0, false, true, true})));
  for (unsigned i = 0; i < redairwayinletcomponents.size(); ++i)
    raw_in_bc->AddComponent(redairwayinletcomponents[i]);

  condlist.push_back(raw_in_bc);

  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways switching between different types of boundary
  // conditions

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_in_switch_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN NODE Reduced D AIRWAYS SWITCH FLOW PRESSURE CONDITIONS",
          "RedAirwaySwitchFlowPressureCond",
          "Reduced d airway switch flow pressure boundary condition",
          Core::Conditions::RedAirwayPrescribedSwitchCond, true,
          Core::Conditions::geometry_type_point));

  Input::AddNamedInt(raw_in_switch_bc, "FUNCT_ID_FLOW");
  Input::AddNamedInt(raw_in_switch_bc, "FUNCT_ID_PRESSURE");
  Input::AddNamedInt(raw_in_switch_bc, "FUNCT_ID_PRESSURE_ACTIVE");

  condlist.push_back(raw_in_switch_bc);

  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional airways external pressure

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_pext_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE Reduced D AIRWAYS PRESCRIBED EXTERNAL PRESSURE CONDITIONS",
          "RedAirwayPrescribedExternalPressure",
          "Reduced d airway prescribed external pressure boundary condition",
          Core::Conditions::RedAirwayPrescribedExternalPressure, true,
          Core::Conditions::geometry_type_line));

  raw_pext_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("boundarycond",
      "ExternalPressure", Teuchos::tuple<std::string>("ExternalPressure"),
      Teuchos::tuple<std::string>("ExternalPressure"), true)));

  std::vector<Teuchos::RCP<Input::LineComponent>> redairwaypextcomponents;
  redairwaypextcomponents.push_back(Teuchos::rcp(new Input::RealVectorComponent("val", 1)));
  redairwaypextcomponents.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("curve", 2, {0, true, true})));
  redairwaypextcomponents.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("funct", 1, {0, false, false, true})));
  for (unsigned i = 0; i < redairwaypextcomponents.size(); ++i)
    raw_pext_bc->AddComponent(redairwaypextcomponents[i]);

  condlist.push_back(raw_pext_bc);


  /*--------------------------------------------------------------------*/
  // Prescribed BC for reduced dimensional scalar transport in airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_in_scatra_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN NODE Reduced D AIRWAYS PRESCRIBED SCATRA CONDITIONS",
          "RedAirwayPrescribedScatraCond", "Reduced d airway prescribed scatra boundary condition",
          Core::Conditions::RedAirwayPrescribedScatraCond, true,
          Core::Conditions::geometry_type_point));

  std::vector<Teuchos::RCP<Input::LineComponent>> redairwayinletscatracomponents;
  redairwayinletscatracomponents.push_back(Teuchos::rcp(new Input::RealVectorComponent("val", 1)));
  redairwayinletscatracomponents.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("curve", 1, {0, true, true})));
  redairwayinletscatracomponents.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("funct", 1, {0, false, false, true})));
  for (unsigned i = 0; i < redairwayinletscatracomponents.size(); ++i)
    raw_in_scatra_bc->AddComponent(redairwayinletscatracomponents[i]);

  condlist.push_back(raw_in_scatra_bc);

  /*--------------------------------------------------------------------*/
  // Prescribed BC for initial values of the scalar transport in reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_int_scatra_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE Reduced D AIRWAYS INITIAL SCATRA CONDITIONS", "RedAirwayInitialScatraCond",
          "Reduced d airway initial scatra boundary condition",
          Core::Conditions::RedAirwayInitialScatraCond, true,
          Core::Conditions::geometry_type_line));

  raw_int_scatra_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("scalar", "O2",
      Teuchos::tuple<std::string>("O2", "CO2"), Teuchos::tuple<std::string>("O2", "CO2"), true)));

  Input::AddNamedReal(raw_int_scatra_bc, "CONCENTRATION");
  condlist.push_back(raw_int_scatra_bc);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions of scatra exchange
  Teuchos::RCP<Core::Conditions::ConditionDefinition> scatra_exchange_cond =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE Reduced D AIRWAYS SCATRA EXCHANGE CONDITIONS", "RedAirwayScatraExchangeCond",
          "scatra exchange condition", Core::Conditions::RedAirwayScatraExchangeCond, true,
          Core::Conditions::geometry_type_line));

  scatra_exchange_cond->AddComponent(Teuchos::rcp(new Input::IntComponent("ConditionID")));

  condlist.push_back(scatra_exchange_cond);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions with hemoglobin
  Teuchos::RCP<Core::Conditions::ConditionDefinition> scatra_hemoglobin_cond =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE Reduced D AIRWAYS HEMOGLOBIN CONDITIONS", "RedAirwayScatraHemoglobinCond",
          "scatra hemoglobin condition", Core::Conditions::RedAirwayScatraHemoglobinCond, false,
          Core::Conditions::geometry_type_line));

  Input::AddNamedReal(scatra_hemoglobin_cond, "INITIAL_CONCENTRATION");
  condlist.push_back(scatra_hemoglobin_cond);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions with hemoglobin
  Teuchos::RCP<Core::Conditions::ConditionDefinition> scatra_air_cond = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN LINE Reduced D AIRWAYS AIR CONDITIONS",
          "RedAirwayScatraAirCond", "scatra air condition",
          Core::Conditions::RedAirwayScatraAirCond, false, Core::Conditions::geometry_type_line));

  Input::AddNamedReal(scatra_air_cond, "INITIAL_CONCENTRATION");
  condlist.push_back(scatra_air_cond);

  /*--------------------------------------------------------------------*/
  // Reduced D airway Scatra condition for regions with hemoglobin
  Teuchos::RCP<Core::Conditions::ConditionDefinition> scatra_capillary_cond =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE Reduced D AIRWAYS CAPILLARY CONDITIONS", "RedAirwayScatraCapillaryCond",
          "scatra capillary condition", Core::Conditions::RedAirwayScatraCapillaryCond, false,
          Core::Conditions::geometry_type_line));

  condlist.push_back(scatra_capillary_cond);

  /*--------------------------------------------------------------------*/
  // Prescribed Ventilator BC for reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_vent_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN NODE Reduced D AIRWAYS VENTILATOR CONDITIONS", "RedAirwayVentilatorCond",
          "Reduced d airway prescribed ventilator condition",
          Core::Conditions::RedAirwayVentilatorCond, true, Core::Conditions::geometry_type_point));

  raw_vent_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("phase1", "flow",
      Teuchos::tuple<std::string>("flow", "volume", "pressure"),
      Teuchos::tuple<std::string>("flow", "volume", "pressure"), true)));

  raw_vent_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("Phase1Smoothness", "smooth",
      Teuchos::tuple<std::string>("smooth", "discontinous"),
      Teuchos::tuple<std::string>("smooth", "discontinous"), true)));

  raw_vent_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("phase2", "pressure",
      Teuchos::tuple<std::string>("pressure", "flow", "volume"),
      Teuchos::tuple<std::string>("pressure", "flow", "volume"), true)));

  raw_vent_bc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("Phase2Smoothness", "smooth",
      Teuchos::tuple<std::string>("smooth", "discontinous"),
      Teuchos::tuple<std::string>("smooth", "discontinous"), true)));

  Input::AddNamedReal(raw_vent_bc, "period");
  Input::AddNamedReal(raw_vent_bc, "phase1_period");
  Input::AddNamedReal(raw_vent_bc, "smoothness_period1");
  Input::AddNamedReal(raw_vent_bc, "smoothness_period2");

  std::vector<Teuchos::RCP<Input::LineComponent>> redairwayventcomponents;
  redairwayventcomponents.push_back(Teuchos::rcp(new Input::RealVectorComponent("val", 2)));
  redairwayventcomponents.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("curve", 2, {0, true, true})));
  for (unsigned i = 0; i < redairwayventcomponents.size(); ++i)
    raw_vent_bc->AddComponent(redairwayventcomponents[i]);

  condlist.push_back(raw_vent_bc);



  /*--------------------------------------------------------------------*/
  // Prescribed volume dependent pleural pressure for reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_volPpl_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE REDUCED D AIRWAYS VOL DEPENDENT PLEURAL PRESSURE CONDITIONS",
          "RedAirwayVolDependentPleuralPressureCond",
          "Reduced D airways volume-dependent peural pressure condition",
          Core::Conditions::RedAirwayVolDependentPleuralPressureCond, true,
          Core::Conditions::geometry_type_line));

  raw_volPpl_bc->AddComponent(
      Teuchos::rcp(new Input::SelectionComponent("TYPE", "Linear_Exponential",
          Teuchos::tuple<std::string>("Linear_Polynomial", "Linear_Exponential", "Linear_Ogden",
              "Nonlinear_Polynomial", "Nonlinear_Exponential", "Nonlinear_Ogden"),
          Teuchos::tuple<std::string>("Linear_Polynomial", "Linear_Exponential", "Linear_Ogden",
              "Nonlinear_Polynomial", "Nonlinear_Exponential", "Nonlinear_Ogden"),
          true)));

  Input::AddNamedReal(raw_volPpl_bc, "TLC");
  Input::AddNamedReal(raw_volPpl_bc, "RV");

  Input::AddNamedReal(raw_volPpl_bc, "P_PLEURAL_0");
  Input::AddNamedReal(raw_volPpl_bc, "P_PLEURAL_LIN");
  Input::AddNamedReal(raw_volPpl_bc, "P_PLEURAL_NONLIN");
  Input::AddNamedReal(raw_volPpl_bc, "TAU");


  std::vector<Teuchos::RCP<Input::LineComponent>> raw_volPpl_bc_components;
  raw_volPpl_bc_components.push_back(Teuchos::rcp(new Input::RealVectorComponent("val", 1)));
  raw_volPpl_bc_components.push_back(
      Teuchos::rcp(new Input::IntVectorComponent("curve", 1, {0, true, true})));
  for (unsigned i = 0; i < raw_volPpl_bc_components.size(); ++i)
    raw_volPpl_bc->AddComponent(raw_volPpl_bc_components[i]);

  condlist.push_back(raw_volPpl_bc);

  /*--------------------------------------------------------------------*/
  // Evaluate lung volume condition for reduced dimensional airways

  Teuchos::RCP<Core::Conditions::ConditionDefinition> raw_eval_lungV_bc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN LINE REDUCED D AIRWAYS EVALUATE LUNG VOLUME CONDITIONS",
          "RedAirwayEvalLungVolCond", "Reduced D airways evaluate lung volume condition",
          Core::Conditions::RedAirwayEvalLungVolCond, true, Core::Conditions::geometry_type_line));


  condlist.push_back(raw_eval_lungV_bc);


  /*--------------------------------------------------------------------*/
  // Impedance condition

  Teuchos::RCP<Core::Conditions::ConditionDefinition> impedancebc =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN SURF IMPEDANCE CONDITIONS",
          "ImpedanceCond", "Impedance boundary condition", Core::Conditions::ImpedanceCond, true,
          Core::Conditions::geometry_type_surface));

  impedancebc->AddComponent(Teuchos::rcp(new Input::IntComponent("ConditionID")));

  impedancebc->AddComponent(Teuchos::rcp(new Input::SelectionComponent("TYPE", "windkessel",
      Teuchos::tuple<std::string>("windkessel", "resistive", "pressure_by_funct"),
      Teuchos::tuple<std::string>("windkessel", "resistive", "pressure_by_funct"), true)));
  Input::AddNamedReal(impedancebc, "R1");
  Input::AddNamedReal(impedancebc, "R2");
  Input::AddNamedReal(impedancebc, "C");
  Input::AddNamedReal(impedancebc, "TIMEPERIOD");
  Input::AddNamedInt(impedancebc, "FUNCT");

  condlist.push_back(impedancebc);
}

FOUR_C_NAMESPACE_CLOSE
