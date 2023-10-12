/*----------------------------------------------------------------------*/
/*! \file
 \brief input parameters for porous multiphase problem with scalar transport

   \level 3

 *----------------------------------------------------------------------*/

#ifndef INPAR_INPAR_POROMULTIPHASE_SCATRA_CPP
#define INPAR_INPAR_POROMULTIPHASE_SCATRA_CPP


#include "baci_inpar_poromultiphase_scatra.H"

#include "baci_inpar_poroelast.H"
#include "baci_inpar_scatra.H"
#include "baci_inpar_validparameters.H"
#include "baci_lib_conditiondefinition.H"
#include "baci_linalg_equilibrate.H"


void INPAR::POROMULTIPHASESCATRA::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  // ----------------------------------------------------------------------
  // (1) general control parameters
  Teuchos::ParameterList& poromultiphasescatradyn = list->sublist("POROMULTIPHASESCATRA DYNAMIC",
      false, "Control paramters for scatra porous multiphase media coupling");

  // Output type
  IntParameter("RESTARTEVRY", 1, "write restart possibility every RESTARTEVRY steps",
      &poromultiphasescatradyn);
  // Time loop control
  IntParameter("NUMSTEP", 200, "maximum number of Timesteps", &poromultiphasescatradyn);
  DoubleParameter("MAXTIME", 1000.0, "total simulation time", &poromultiphasescatradyn);
  DoubleParameter("TIMESTEP", 0.05, "time step size dt", &poromultiphasescatradyn);
  IntParameter("RESULTSEVRY", 1, "increment for writing solution", &poromultiphasescatradyn);
  IntParameter("ITEMAX", 10, "maximum number of iterations over fields", &poromultiphasescatradyn);
  IntParameter("ITEMIN", 1, "minimal number of iterations over fields", &poromultiphasescatradyn);

  // Coupling strategy for poroscatra solvers
  setStringToIntegralParameter<int>("COUPALGO", "twoway_partitioned_nested",
      "Coupling strategies for poroscatra solvers",
      tuple<std::string>(
          "twoway_partitioned_nested", "twoway_partitioned_sequential", "twoway_monolithic"),
      tuple<int>(solscheme_twoway_partitioned_nested, solscheme_twoway_partitioned_sequential,
          solscheme_twoway_monolithic),
      &poromultiphasescatradyn);

  // coupling with 1D artery network active
  BoolParameter(
      "ARTERY_COUPLING", "No", "Coupling with 1D blood vessels.", &poromultiphasescatradyn);

  // no convergence of coupling scheme
  setStringToIntegralParameter<int>("DIVERCONT", "stop",
      "What to do with time integration when Poromultiphase-Scatra iteration failed",
      tuple<std::string>("stop", "continue"), tuple<int>(divcont_stop, divcont_continue),
      &poromultiphasescatradyn);

  // ----------------------------------------------------------------------
  // (2) monolithic parameters
  Teuchos::ParameterList& poromultiphasescatradynmono = poromultiphasescatradyn.sublist(
      "MONOLITHIC", false, "Parameters for monolithic Poro-Multiphase-Scatra Interaction");

  setStringToIntegralParameter<int>("VECTORNORM_RESF", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<int>(INPAR::POROMULTIPHASESCATRA::norm_l1, INPAR::POROMULTIPHASESCATRA::norm_l1_scaled,
          INPAR::POROMULTIPHASESCATRA::norm_l2, INPAR::POROMULTIPHASESCATRA::norm_rms,
          INPAR::POROMULTIPHASESCATRA::norm_inf),
      &poromultiphasescatradynmono);

  setStringToIntegralParameter<int>("VECTORNORM_INC", "L2",
      "type of norm to be applied to residuals",
      tuple<std::string>("L1", "L1_Scaled", "L2", "Rms", "Inf"),
      tuple<int>(INPAR::POROMULTIPHASESCATRA::norm_l1, INPAR::POROMULTIPHASESCATRA::norm_l1_scaled,
          INPAR::POROMULTIPHASESCATRA::norm_l2, INPAR::POROMULTIPHASESCATRA::norm_rms,
          INPAR::POROMULTIPHASESCATRA::norm_inf),
      &poromultiphasescatradynmono);

  // convergence criteria adaptivity --> note ADAPTCONV_BETTER set pretty small
  BoolParameter("ADAPTCONV", "yes",
      "Switch on adaptive control of linear solver tolerance for nonlinear solution",
      &poromultiphasescatradynmono);
  DoubleParameter("ADAPTCONV_BETTER", 0.001,
      "The linear solver shall be this much better "
      "than the current nonlinear residual in the nonlinear convergence limit",
      &poromultiphasescatradynmono);

  // Iterationparameters
  DoubleParameter("TOLRES_GLOBAL", 1e-8, "tolerance in the residual norm for the Newton iteration",
      &poromultiphasescatradynmono);
  DoubleParameter("TOLINC_GLOBAL", 1e-8, "tolerance in the increment norm for the Newton iteration",
      &poromultiphasescatradynmono);

  // number of linear solver used for poroelasticity
  IntParameter("LINEAR_SOLVER", -1,
      "number of linear solver used for monolithic poroscatra problems",
      &poromultiphasescatradynmono);

  // parameters for finite difference check
  setStringToIntegralParameter<int>("FDCHECK", "none",
      "flag for finite difference check: none or global",
      tuple<std::string>("none",
          "global"),  // perform finite difference check on time integrator level
      tuple<int>(fdcheck_none, fdcheck_global), &poromultiphasescatradynmono);

  // flag for equilibration of global system of equations
  setStringToIntegralParameter<CORE::LINALG::EquilibrationMethod>("EQUILIBRATION", "none",
      "flag for equilibration of global system of equations",
      tuple<std::string>("none", "rows_full", "rows_maindiag", "columns_full", "columns_maindiag",
          "rowsandcolumns_full", "rowsandcolumns_maindiag"),
      tuple<CORE::LINALG::EquilibrationMethod>(CORE::LINALG::EquilibrationMethod::none,
          CORE::LINALG::EquilibrationMethod::rows_full,
          CORE::LINALG::EquilibrationMethod::rows_maindiag,
          CORE::LINALG::EquilibrationMethod::columns_full,
          CORE::LINALG::EquilibrationMethod::columns_maindiag,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_full,
          CORE::LINALG::EquilibrationMethod::rowsandcolumns_maindiag),
      &poromultiphasescatradynmono);

  // ----------------------------------------------------------------------
  // (3) partitioned parameters
  Teuchos::ParameterList& poromultiphasescatradynpart = poromultiphasescatradyn.sublist(
      "PARTITIONED", false, "Parameters for partitioned Poro-Multiphase-Scatra Interaction");

  // convergence tolerance of outer iteration loop
  DoubleParameter("CONVTOL", 1e-6, "tolerance for convergence check of outer iteration",
      &poromultiphasescatradynpart);
}

void INPAR::POROMULTIPHASESCATRA::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // oxygen partial pressure calculation condition
  {
    // definition of oxygen partial pressure calculation condition
    Teuchos::RCP<ConditionDefinition> oxypartpressline = Teuchos::rcp(
        new ConditionDefinition("DESIGN OXYGEN PARTIAL PRESSURE CALCULATION LINE CONDITIONS",
            "PoroMultiphaseScatraOxyPartPressCalcCond",
            "PoroMultiphaseScatra Oxygen Partial Pressure Calculation line condition",
            DRT::Condition::PoroMultiphaseScatraOxyPartPressCalcCond, true, DRT::Condition::Line));
    Teuchos::RCP<ConditionDefinition> oxypartpresssurf = Teuchos::rcp(new ConditionDefinition(
        "DESIGN OXYGEN PARTIAL PRESSURE CALCULATION SURF CONDITIONS",
        "PoroMultiphaseScatraOxyPartPressCalcCond",
        "PoroMultiphaseScatra Oxygen Partial Pressure Calculation surface condition",
        DRT::Condition::PoroMultiphaseScatraOxyPartPressCalcCond, true, DRT::Condition::Surface));
    Teuchos::RCP<ConditionDefinition> oxypartpressvol = Teuchos::rcp(new ConditionDefinition(
        "DESIGN OXYGEN PARTIAL PRESSURE CALCULATION VOL CONDITIONS",
        "PoroMultiphaseScatraOxyPartPressCalcCond",
        "PoroMultiphaseScatra Oxygen Partial Pressure Calculation volume condition",
        DRT::Condition::PoroMultiphaseScatraOxyPartPressCalcCond, true, DRT::Condition::Volume));

    // equip condition definitions with input file line components
    std::vector<Teuchos::RCP<::INPUT::LineComponent>> oxypartpresscomponents;

    {
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("SCALARID")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::IntComponent("SCALARID")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("n")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::RealComponent("n")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("Pb50")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::RealComponent("Pb50")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("CaO2_max")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::RealComponent("CaO2_max")));
      oxypartpresscomponents.push_back(
          Teuchos::rcp(new ::INPUT::SeparatorComponent("alpha_bl_eff")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::RealComponent("alpha_bl_eff")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("rho_oxy")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::RealComponent("rho_oxy")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::SeparatorComponent("rho_bl")));
      oxypartpresscomponents.push_back(Teuchos::rcp(new ::INPUT::RealComponent("rho_bl")));
    }

    // insert input file line components into condition definitions
    for (unsigned i = 0; i < oxypartpresscomponents.size(); ++i)
    {
      oxypartpressline->AddComponent(oxypartpresscomponents[i]);
      oxypartpresssurf->AddComponent(oxypartpresscomponents[i]);
      oxypartpressvol->AddComponent(oxypartpresscomponents[i]);
    }

    // insert condition definitions into global list of valid condition definitions
    condlist.push_back(oxypartpressline);
    condlist.push_back(oxypartpresssurf);
    condlist.push_back(oxypartpressvol);
  }
  return;
}


#endif  // INPAR_INPAR_POROMULTIPHASE_SCATRA_CPP
