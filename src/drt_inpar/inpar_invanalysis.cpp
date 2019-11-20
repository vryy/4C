/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for inverse analysis

\level 2
\maintainer Harald Willmann

*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_invanalysis.H"
#include "inpar_statinvanalysis.H"
#include "inpar_structure.H"
#include "../drt_lib/drt_conditiondefinition.H"

#include "../drt_lib/drt_globalproblem_enums.H"
#include "inpar_problemtype.H"


void INPAR::INVANA::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::Array<std::string> yesnotuple =
      tuple<std::string>("Yes", "No", "yes", "no", "YES", "NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true, false, true, false, true, false);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& iap = list->sublist("INVERSE ANALYSIS", false, "");

  // Inverse Analysis
  setStringToIntegralParameter<int>("INV_ANALYSIS", "none",
      "types of inverse analysis and on/off switch", tuple<std::string>("none", "gen"),
      tuple<int>(INPAR::STR::inv_none, INPAR::STR::inv_generalized), &iap);

  // Special types of inv_analysis
  setStringToIntegralParameter<int>("SPECIAL_INV_ANALYSIS", "none",
      "types of special inverse analysis", tuple<std::string>("none", "coupled", "multiple"),
      tuple<int>(INPAR::STR::spec_inv_none, INPAR::STR::spec_inv_coup, INPAR::STR::spec_inv_mult),
      &iap);

  // tolerance for inv_analysis
  DoubleParameter("INV_ANA_TOL", 1.0, "tolerance for inverse analysis", &iap);
  IntParameter("INV_ANA_MAX_RUN", 100, "max iterations for inverse analysis", &iap);

  // perturbation parameters
  DoubleParameter("INV_ALPHA", 1.0e-3,
      "Absolute perturbation parameters for finite differencing of gradient", &iap);
  DoubleParameter("INV_BETA", 1.0e-3,
      "Relative perturbation parameters (to be scaled with the parameter) for finite differencing "
      "of gradient",
      &iap);

  // initial regularization parameter
  DoubleParameter("INV_INITREG", 1.0, "initial regularization parameter", &iap);

  // strategy to update regularization parameter
  setStringToIntegralParameter<int>("UPDATE_REG", "RES",
      "Update strategy for regularization parameter ",
      tuple<std::string>("RES", "res", "GRAD", "grad"),
      tuple<int>(INPAR::STR::reg_update_res, INPAR::STR::reg_update_res,
          INPAR::STR::reg_update_grad, INPAR::STR::reg_update_grad),
      &iap);

  IntParameter("EXP_ID_MULT_INV_ANA", -1,
      "ID of experiment when running multiple files in one inverse analysis. ID ranges from 0 to "
      "#experiments-1",
      &iap);

  StringParameter(
      "MONITORFILE", "none.monitor", "filename of file containing measured displacements", &iap);

  setNumericStringParameter("INV_LIST", "-1", "IDs of materials that have to be fitted", &iap);

  setStringToIntegralParameter<int>(
      "NEW_FILES", "yes", "new result files for each run", yesnotuple, yesnovalue, &iap);
  setStringToIntegralParameter<int>("PARAM_BOUNDS", "no",
      "Reset parameters if optstep predicts negative values", yesnotuple, yesnovalue, &iap);

  BoolParameter("PATCHES", "No", "Do you want to use smoothed patches?", &iap);
  StringParameter("DEFINEPATCHES", "MaterialNumber",
      "define how the patches are defined: MaterialNumber or Uniform", &iap);
  IntParameter("NUMPATCHES", 0, "number of patches", &iap);
  setNumericStringParameter(
      "INV_LIST_PATCHES", "-1", "IDs of materials that are included in the patches", &iap);
  IntParameter("SMOOTHINGSTEPSPATCHES", 1, "number of smoothing steps that are performed", &iap);
  setNumericStringParameter("STARTVALUESFORPATCHES", "1.0",
      "startvalues for the patches, only needed for Uniform Patches", &iap);

  {
    Teuchos::Array<std::string> name;
    Teuchos::Array<int> label;

    // fill the arrays
    {
      std::map<std::string, ProblemType> map = INPAR::PROBLEMTYPE::StringToProblemTypeMap();
      std::map<std::string, ProblemType>::const_iterator i;
      for (i = map.begin(); i != map.end(); ++i)
      {
        name.push_back(i->first);
        label.push_back(i->second);
      }
    }

    setStringToIntegralParameter<int>("FORWARD_PROBLEMTYP", "Structure",
        "defines which forward problem type is used for inverse analysis", name, label, &iap);

    // Type of measurement written in monitor file
    setStringToIntegralParameter<int>("MEASUREMENT_TYPE", "dofs",
        "Type of measurement written in monitor file", tuple<std::string>("dofs", "points"),
        tuple<int>(INPAR::STR::meas_dofs, INPAR::STR::meas_points), &iap);
  }

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& statinvp = list->sublist("STAT INVERSE ANALYSIS", false, "");

  // Statistical Inverse Analysis switch
  setStringToIntegralParameter<int>("STAT_INV_ANALYSIS", "none",
      "types of statistical inverse analysis and on/off switch",
      tuple<std::string>(
          "none", "MonteCarloSMC", "MonteCarloMH", "LBFGS", "BruteForce", "ParticlePrediction"),
      tuple<int>(stat_inv_none, stat_inv_smc, stat_inv_mh, stat_inv_lbfgs, stat_inv_bruteforce,
          stat_inv_prediction),
      &statinvp);

  // initial scaling for the LBFGS algorithm
  BoolParameter("LBFGSINITSCAL", "yes", "want initial scaling for the LBFGS?", &statinvp);

  // step to restart from
  IntParameter("FPRESTART", 0, "forward problem restart", &statinvp);

  StringParameter("FPOUTPUTFILENAME", "none",
      "controlfilename (without .control) which to use as forward problem output and restartfrom",
      &statinvp);

  // write restart info every so often
  IntParameter("RESTARTEVRY", 1, "write restart information every x-th step", &statinvp);

  // decide which parametrization of material parameters to use
  setStringToIntegralParameter<int>("PARAMETRIZATION", "none",
      "how to parametrize the parameter field",
      tuple<std::string>("none", "patchwise", "tvsvd", "elementwise", "uniform"),
      tuple<int>(stat_inv_mp_none, stat_inv_mp_patchwise, stat_inv_mp_tvsvd,
          stat_inv_mp_elementwise, stat_inv_mp_uniform),
      &statinvp);

  // number of levels/eigenvectors for the basis reduction
  IntParameter("NUM_REDUCT_LEVELS", 4,
      "number of levels for the basis reduction (patch-levels, or eigenvectors)", &statinvp);

  // anasazi's number of eigenvectors for the eigensolver in INVANA::MatParManagerTVSVD
  IntParameter("TVSVD_ANASAZI_NEV", 10,
      "anasazi's number of eigenvectors for the eigensolver in INVANA::MatParManagerTVSVD",
      &statinvp);

  // anasazi's number of blocks for the eigensolver in INVANA::MatParManagerTVSVD
  IntParameter("TVSVD_ANASAZI_NBLOCKS", 4,
      "anasazi's number of blocks for the eigensolver in INVANA::MatParManagerTVSVD", &statinvp);

  // anasazi's blocksize for the eigensolver in INVANA::MatParManagerTVSVD
  IntParameter("TVSVD_ANASAZI_BSIZE", 10,
      "anasazi's blocksize for the eigensolver in INVANA::MatParManagerTVSVD", &statinvp);

  // decide which weights to use for the graph of the elementwise parametrization
  setStringToIntegralParameter<int>("GRAPHWEIGHTS", "area",
      "weights for the elementwise graph creation", tuple<std::string>("area", "unity"),
      tuple<int>(stat_inv_graph_area, stat_inv_graph_unity), &statinvp);

  // want some regularization
  setStringToIntegralParameter<int>("REGULARIZATION", "none",
      "want regularization? ('tikhonov', 'totalvariation', 'none')",
      tuple<std::string>("none", "tikhonov", "totalvariation"),
      tuple<int>(stat_inv_reg_none, stat_inv_reg_tikhonov, stat_inv_reg_totalvariation), &statinvp);

  // objective function
  setStringToIntegralParameter<int>("OBJECTIVEFUNCT", "none",
      "choose type of objective function ('displacements', 'surfcurr')",
      tuple<std::string>("none", "displacements", "surfcurr"),
      tuple<int>(stat_inv_obj_none, stat_inv_obj_disp, stat_inv_obj_surfcurr), &statinvp);

  // scaling of objective function
  BoolParameter("OBJECTIVEFUNCTSCAL", "No", "want scaling of objective function?", &statinvp);

  // monitorfile to provide measurements
  StringParameter("MONITORFILE", "none.monitor",
      "filename of file containing measured displacements", &statinvp);

  // target discretization for surface currents
  StringParameter(
      "TARGETDISCRETIZATION", "none.dat", "datfile containing target discretization", &statinvp);

  // list of parameters for the respective material
  StringParameter("PARAMLIST", "none",
      "list of std::string of parameters to be optimized, order as in INV_LIST e.g. 1 YOUNG BETA",
      &statinvp);

  // number of optimization steps
  IntParameter("MAXITER", 100, "max iterations for inverse analysis", &statinvp);

  // number of optimization steps before using
  // parameter continuation in the forward problem
  IntParameter(
      "ITERTOPC", 10, "iterations before parameter continuation in the forward problem", &statinvp);

  // give prestressing method to be used for the adjoint formulation
  setStringToIntegralParameter<int>("PRESTRESS", "none", "prestressing takes values none mulf id",
      tuple<std::string>("none", "None", "NONE", "mulf", "Mulf", "MULF", "id", "Id", "ID"),
      tuple<int>(INPAR::STR::prestress_none, INPAR::STR::prestress_none, INPAR::STR::prestress_none,
          INPAR::STR::prestress_mulf, INPAR::STR::prestress_mulf, INPAR::STR::prestress_mulf,
          INPAR::STR::prestress_id, INPAR::STR::prestress_id, INPAR::STR::prestress_id),
      &statinvp);

  // stepsize for deterministic gradient based schemes
  DoubleParameter("STEPSIZE", 1.0, "stepsize for the gradient descent scheme", &statinvp);

  // linspace for the BruteForce optimizer
  setNumericStringParameter("BFLINSPACE", "0.0 1.0 100",
      "linear space (a,b,n) of n linearly spaced samples in the interval [a,b] ", &statinvp);

  // convergence criterion tolerance
  DoubleParameter("CONVTOL", 1.0e-06,
      "stop optimizaiton iterations for convergence criterion below this value", &statinvp);

  // weight of the Tikhonov regularization
  DoubleParameter("REG_WEIGHT", 1.0, "weight of the regularization functional", &statinvp);

  // regularization of the totalvariation functional
  DoubleParameter("TVD_EPS", 0.01, "differentiation epsilon for total variation", &statinvp);

  // number of optimization steps
  IntParameter("SIZESTORAGE", 20,
      "number of vectors to keep in storage; defaults to 20 (lbfgs usage only)", &statinvp);

  // number of SMC particles
  IntParameter(
      "NUM_PARTICLES", 1, "number of particles for the sequential monte carlo.", &statinvp);

  // meta parametrization of material parameters
  setStringToIntegralParameter<int>("METAPARAMS", "none",
      "choose type of metaparametrization (none/quad/exp/arctan)",
      tuple<std::string>("none", "quad", "exp", "arctan"),
      tuple<int>(stat_inv_meta_none, stat_inv_meta_quad, stat_inv_meta_exp, stat_inv_meta_arctan),
      &statinvp);


  // scale of the kernel functions used in surface currents
  DoubleParameter("KERNELSCALE", -1.0, "scale of the kernel function", &statinvp);

  // estimation of the variance of the measurement noise
  DoubleParameter("MEASVARESTIM", 1.0, "variance of the measurement noise", &statinvp);

  // add synthetic noise to the measurements (e.g. in case of synthetic data)
  BoolParameter("SYNTHNOISE", "No", "want noise on your synthetic measurements?", &statinvp);

  // seed used for synthetic noise geenration
  IntParameter("SYNTHNOISESEED", 1, "seed to be used for synthetic noise generation", &statinvp);

  // initial scale of the markov kernel used by the SMC algorithm
  DoubleParameter(
      "MC_INIT_SCALE", 0.1, "inital scaling factor for the markov kernel in the SMC", &statinvp);

  // number of kernel applications in smc rejuvenation step
  IntParameter(
      "SMC_KERNEL_ITER", 1, "number of kernel applications in smc rejuvenation", &statinvp);

  // scale the covariance matrix for monte carlo algorithms using it
  DoubleParameter(
      "MAP_PRIOR_SCALE", 1.0, "scaling for the prior covariance in the LogLikePrior", &statinvp);

  // file to read MAP approximation (as initial guess) from
  StringParameter("MAP_RESTARTFILE", "none",
      "control file to read the maximum a posterior approximation (as initial guess) from",
      &statinvp);

  // step from which to read the MAP (as initial guess) approximation
  IntParameter("MAP_RESTART", 0,
      "step to read the maximum a posterior approximation (as initial guess) from", &statinvp);

  // file to read MAP approximation from
  StringParameter("MAP_REDUCT_RESTARTFILE", "none",
      "control file to read the maximum a posterior approximation (to create a reduced basis) from",
      &statinvp);

  // step from which to read the MAP approximation
  IntParameter("MAP_REDUCT_RESTART", 0,
      "step to read the maximum a posterior approximation (to create a reduced basis) from",
      &statinvp);

  // target effective sample size reduction per time step
  DoubleParameter(
      "SMC_ESS_REDUCTION", 0.05, "targeted effective sample size reduction per step", &statinvp);

  // iterations used to adapt the acceptance ratio
  IntParameter("MH_ACCADAPT_ITER", 0, "iterations used to adapt the acceptance ratio", &statinvp);

  // adapt the acceptance ratio every x iterations
  IntParameter("MH_ACCADAPT_EVRY", 0, "adapt the acceptance ratio every x iterations", &statinvp);

  // use only every thin-th sample for the statistic
  IntParameter("MH_THIN", 0, "use only every thin-th sample for the statistic", &statinvp);

  // use samples in the statistic only after burnin
  IntParameter("MH_BURNIN", 0, "use samples in the statistic only after burnin", &statinvp);

  // decide how the initialize the optimization
  setStringToIntegralParameter<int>("INIT_TYPE", "dat", "how to initialize the optimization",
      tuple<std::string>("dat", "map"), tuple<int>(stat_inv_init_dat, stat_inv_init_map),
      &statinvp);
}



void INPAR::INVANA::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist)
{
  using namespace DRT::INPUT;

  /*--------------------------------------------------------------------*/
  // inverse analysis fitted surface

  std::vector<Teuchos::RCP<ConditionComponent>> invanacomponents;
  invanacomponents.push_back(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  Teuchos::RCP<ConditionDefinition> surfinvana = Teuchos::rcp(new ConditionDefinition(
      "DESIGN SURFACE INV ANALYSIS", "SurfInvAna", "Inverse Analysis Surface",
      DRT::Condition::InvAnaSurface, true, DRT::Condition::Surface));

  surfinvana->AddComponent(Teuchos::rcp(new IntConditionComponent("ConditionID")));

  condlist.push_back(surfinvana);

  /*--------------------------------------------------------------------*/
  // Surface current evaluation condition

  std::vector<Teuchos::RCP<ConditionComponent>> surfcurrcomponents;
  surfcurrcomponents.push_back(Teuchos::rcp(new IntConditionComponent("matching id")));

  Teuchos::RCP<ConditionDefinition> surfcurrcond = Teuchos::rcp(
      new ConditionDefinition("DESIGN SURFACE CURRENT EVALUATION CONDITION", "SurfaceCurrent",
          "Surface current", DRT::Condition::SurfaceCurrent, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < surfcurrcomponents.size(); ++i)
  {
    surfcurrcond->AddComponent(surfcurrcomponents[i]);
  }

  condlist.push_back(surfcurrcond);

  /*--------------------------------------------------------------------*/
  // Uncertain surface condition

  std::vector<Teuchos::RCP<ConditionComponent>> uncertsurfcomponents;
  uncertsurfcomponents.push_back(Teuchos::rcp(new IntConditionComponent("matching id")));

  Teuchos::RCP<ConditionDefinition> uncertsurfcond =
      Teuchos::rcp(new ConditionDefinition("DESIGN UNCERTAIN SURFACE CONDITION", "UncertainSurface",
          "Uncertain surface", DRT::Condition::UncertainSurface, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < uncertsurfcomponents.size(); ++i)
  {
    uncertsurfcond->AddComponent(uncertsurfcomponents[i]);
  }

  condlist.push_back(uncertsurfcond);
}
