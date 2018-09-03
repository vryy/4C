/*----------------------------------------------------------------------*/
/*!
\file inpar_acou.cpp
\level 2
<pre>
\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/

/*----------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_acou.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::ACOU::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& acousticdyn = list->sublist(
      "ACOUSTIC DYNAMIC", false, "control parameters for acoustic or photoacoustic problems\n");

  // general time stepping settings
  DoubleParameter("TIMESTEP", 0.01, "Time increment dt", &acousticdyn);
  IntParameter("NUMSTEP", 100, "Total number of time steps", &acousticdyn);
  DoubleParameter("MAXTIME", 1.0, "Total simulation time", &acousticdyn);
  DoubleParameter("COURANTNUMBER", 0.3, "Total simulation time", &acousticdyn);

  // additional parameters
  IntParameter("CALCERRORFUNCNO", -1, "Function for Error Calculation", &acousticdyn);
  IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &acousticdyn);
  IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &acousticdyn);
  IntParameter(
      "LINEAR_SOLVER", -1, "Number of linear solver used for acoustical problem", &acousticdyn);
  IntParameter("STARTFUNCNO", -1, "Function for Initial Starting Field", &acousticdyn);
  IntParameter("SOURCETERMFUNCNO", -1, "Function for source term in volume", &acousticdyn);
  BoolParameter(
      "DOUBLEORFLOAT", "Yes", "Yes, if evaluation with double, no if with float", &acousticdyn);
  BoolParameter(
      "ALLELESEQUAL", "No", "Yes, if all elements have same shape and material", &acousticdyn);

  // distinguish viscous and lossless flows
  setStringToIntegralParameter<int>("PHYSICAL_TYPE", "lossless", "fluid properties",
      tuple<std::string>("lossless", "solid"), tuple<int>(acou_lossless, acou_solid), &acousticdyn);

  // for viscous flows, one can specify if the displacement gradient or the stresses are outputted
  BoolParameter(
      "WRITESTRESS", "Yes", "Output of stresses instead of displacement gradient", &acousticdyn);

  // photoacoustics
  BoolParameter(
      "PHOTOACOU", "No", "Coupling with Scatra for Diffusive Light Transport", &acousticdyn);
  BoolParameter("MESHCONFORM", "No", "Conformity of scatra and acoustical mesh", &acousticdyn);
  IntParameter("MINPATELEGID", -1, "minimum element gid of overlapping acou element", &acousticdyn);
  IntParameter("MAXPATELEGID", -1, "maximum element gid of overlapping acou element", &acousticdyn);


  // local postprocessing and p-adaptivity
  BoolParameter(
      "ERRORMAPS", "No", "Output of error maps obtained by local postprocessing", &acousticdyn);
  BoolParameter("P_ADAPTIVITY", "No", "p-adaptivity in time integration", &acousticdyn);
  DoubleParameter("P_ADAPT_TOL", 1.0e-15, "Error tolerance for p-adaptivity", &acousticdyn);

  // ader specific
  BoolParameter(
      "SPECTRAL_EVALUATION", "yes", "fast evaluation of Taylor-Cauchy-Kovalewski", &acousticdyn);
  BoolParameter("USE_ADER_POST", "no", "need superconvergent pressure results?", &acousticdyn);

  // ader lts specific
  IntParameter("MAX_NUM_CLUSTERS", 10, "Allowed number of clusters", &acousticdyn);
  IntParameter(
      "MAX_DIFF_CLUSTERS", 3, "Allowed difference of timestep between clusters", &acousticdyn);


  // time integration
  setStringToIntegralParameter<int>("TIMEINT", "impleuler", "Type of time integration scheme",
      tuple<std::string>("impleuler", "expleuler", "classrk4", "lsrk45reg2", "lsrk33reg2",
          "lsrk45reg3", "ssprk", "ader", "aderlts", "adertritet"),
      tuple<int>(acou_impleuler, acou_expleuler, acou_classrk4, acou_lsrk45reg2, acou_lsrk33reg2,
          acou_lsrk45reg3, acou_ssprk, acou_ader, acou_ader_lts, acou_ader_tritet),
      &acousticdyn);

  // output of "measured" pressure curves
  BoolParameter(
      "WRITEMONITOR", "No", "Write a monitor file for Pressure Monitor Condition", &acousticdyn);

  // inverse analysis
  setStringToIntegralParameter<int>("INV_ANALYSIS", "none",
      "Types of inverse analysis and on/off switch",
      tuple<std::string>("none", "patopti", "patoptisplit", "patoptisplitacousplit",
          "patoptisplitacouident", "patoptisplitacouidentsmart", "patreduction"),
      tuple<int>(pat_none, pat_opti, pat_optisplit, pat_optisplitacousplit, pat_optisplitacouident,
          pat_optisplitacouidentsmart, pat_reduction),
      &acousticdyn);

  // PML
  StringParameter("PML_DEFINITION_FILE", "none.txt",
      "Filename of file containing the pml definition", &acousticdyn);

  // sublist PA IMAGE RECONSTRUCTION
  Teuchos::ParameterList& acou_inv = acousticdyn.sublist("PA IMAGE RECONSTRUCTION", false, "");

  setStringToIntegralParameter<int>("OPTIMIZATION", "LBFGS", "types of optimization algorithm",
      tuple<std::string>("GradientDescent", "LBFGS"), tuple<int>(inv_gd, inv_lbfgs), &acou_inv);
  setStringToIntegralParameter<int>("REGULATYPE", "none", "types of regularization",
      tuple<std::string>("none", "tikh", "tvd", "tikhtvd"),
      tuple<int>(pat_regula_none, pat_regula_tikh, pat_regula_tvd, pat_regula_tikhtvd), &acou_inv);
  setStringToIntegralParameter<int>("PATCHTYPE", "none", "types of patch building",
      tuple<std::string>("none", "patchself", "patchreacgrad", "patchreacvals", "patchmixed"),
      tuple<int>(
          pat_patch_none, pat_patch_self, pat_patch_reacgrad, pat_patch_reacvals, pat_patch_mixed),
      &acou_inv);
  setStringToIntegralParameter<int>("TOMOGRAPHTYPE", "unknown", "type of tomograph",
      tuple<std::string>("circle", "unknown"), tuple<int>(pat_circle, pat_unknown), &acou_inv);

  BoolParameter("OVERWRITEOUTPUT", "Yes", "overwrite output (recommended)", &acou_inv);
  BoolParameter("ACOUOUTPUT", "No", "write BACI output for the acoustic problem", &acou_inv);
  StringParameter("MONITORFILE", "none.monitor",
      "Filename of file containing measured pressure values", &acou_inv);
  BoolParameter("FDCHECK", "No", "Finite difference check", &acou_inv);
  DoubleParameter(
      "INV_TOL", 1e-16, "Tolerance for objective function of inverse pat analysis", &acou_inv);
  IntParameter("INV_MAX_RUN", 10, "Maximal run number for inverse pat analysis", &acou_inv);
  IntParameter("INV_LS_MAX_RUN", 10, "Maximal run number for line search in inverse pat analysis",
      &acou_inv);
  DoubleParameter(
      "EPSILON", -1.0, "tolerated distance in which measured curve=nod curve", &acou_inv);
  StringParameter("ACOUPARAMLIST", "none",
      "list of std::string of acoustical parameters to be optimized", &acou_inv);
  StringParameter("OPTIPARAMLIST", "none",
      "list of std::string of optical parameters to be optimized", &acou_inv);
  StringParameter("SEGMENTATIONMATS", "none.material",
      "Filename of file containing table of materials", &acou_inv);
  IntParameter("SEQUENZE", -1, "sequential optimization for penalty and the rest", &acou_inv);
  BoolParameter("TIMEREVERSAL", "No", "Initial reaction guess with time reversal", &acou_inv);
  BoolParameter("SAMPLEOBJECTIVE", "No",
      "Sample objective function for several parameters (need be implemented)", &acou_inv);
  DoubleParameter("TIKHWEIGHT_MUA", 1.0, "tikhonow weight", &acou_inv);
  DoubleParameter("TVDWEIGHT_MUA", 1.0, "tvd weight", &acou_inv);
  DoubleParameter("TIKHWEIGHT_D", 1.0, "tikhonow weight", &acou_inv);
  DoubleParameter("TVDWEIGHT_D", 1.0, "tvd weight", &acou_inv);
  DoubleParameter("TIKHWEIGHT_C", 1.0, "tikhonow weight", &acou_inv);
  DoubleParameter("TVDWEIGHT_C", 1.0, "tvd weight", &acou_inv);
  DoubleParameter("TIKHWEIGHT_RHO", 1.0, "tikhonow weight", &acou_inv);
  DoubleParameter("TVDWEIGHT_RHO", 1.0, "tvd weight", &acou_inv);
  DoubleParameter("TVDEPS_MUA", 0.01, "tvd eps", &acou_inv);
  DoubleParameter("TVDEPS_D", 0.01, "tvd eps", &acou_inv);
  DoubleParameter("TVDEPS_C", 0.01, "tvd eps", &acou_inv);
  DoubleParameter("TVDEPS_RHO", 0.01, "tvd eps", &acou_inv);
  DoubleParameter(
      "IMPULSERESPONSE_DT", 0.0, "time step with which impulse response is recorded", &acou_inv);
  StringParameter("IMPULSERESPONSE", "none.impresp",
      "Filename of file containing the impulse response", &acou_inv);
  BoolParameter("ACOUIDENT_AVG", "No",
      "should acoustical properties set according to their two closest?", &acou_inv);
  DoubleParameter("REDUCTIONCUTTIME", 0.0,
      "the new monitor file shall only have values until this time", &acou_inv);
}


/// set specific acoustic conditions
void INPAR::ACOU::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // absorbing boundary condition for acoustic problems
  // line
  Teuchos::RCP<ConditionDefinition> absorbing_line =
      Teuchos::rcp(new ConditionDefinition("DESIGN ABSORBING LINE CONDITIONS", "Absorbing",
          "Absorbing line for acoustics", DRT::Condition::Absorb, true, DRT::Condition::Line));
  condlist.push_back(absorbing_line);

  // surface
  Teuchos::RCP<ConditionDefinition> absorbing_surface = Teuchos::rcp(new ConditionDefinition(
      "DESIGN ABSORBING SURF CONDITIONS", "Absorbing", "Absorbing surface for acoustics",
      DRT::Condition::Absorb, true, DRT::Condition::Surface));
  condlist.push_back(absorbing_surface);

  /*--------------------------------------------------------------------*/
  // monitor condition for acoustic problems
  // line
  Teuchos::RCP<ConditionDefinition> pressmon_line =
      Teuchos::rcp(new ConditionDefinition("DESIGN PRESSURE MONITOR LINE CONDITIONS",
          "PressureMonitor", "Pressure monitor line for acoustics", DRT::Condition::PressureMonitor,
          true, DRT::Condition::Line));
  condlist.push_back(pressmon_line);

  // surface
  Teuchos::RCP<ConditionDefinition> pressmon_surface =
      Teuchos::rcp(new ConditionDefinition("DESIGN PRESSURE MONITOR SURF CONDITIONS",
          "PressureMonitor", "Pressure monitor surface for acoustics",
          DRT::Condition::PressureMonitor, true, DRT::Condition::Surface));
  condlist.push_back(pressmon_surface);
}
