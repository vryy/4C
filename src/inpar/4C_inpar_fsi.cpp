/*----------------------------------------------------------------------------*/
/*! \file

\level 1


\brief Input parameters for fluid structure interaction

*/
/*----------------------------------------------------------------------------*/

#include "4C_inpar_fsi.hpp"

#include "4C_fem_condition_definition.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
void Inpar::FSI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  Teuchos::ParameterList& fsidyn = list->sublist("FSI DYNAMIC", false,
      "Fluid Structure Interaction\n"
      "FSI solver with various coupling methods");

  Teuchos::Tuple<std::string, 25> name;
  Teuchos::Tuple<int, 25> label;

  name[0] = "basic_sequ_stagg";
  label[0] = fsi_basic_sequ_stagg;
  name[1] = "iter_stagg_fixed_rel_param";
  label[1] = fsi_iter_stagg_fixed_rel_param;
  name[2] = "iter_stagg_AITKEN_rel_param";
  label[2] = fsi_iter_stagg_AITKEN_rel_param;
  name[3] = "iter_stagg_steep_desc";
  label[3] = fsi_iter_stagg_steep_desc;
  name[4] = "iter_stagg_NLCG";
  label[4] = fsi_iter_stagg_NLCG;
  name[5] = "iter_stagg_MFNK_FD";
  label[5] = fsi_iter_stagg_MFNK_FD;
  name[6] = "iter_stagg_MFNK_FSI";
  label[6] = fsi_iter_stagg_MFNK_FSI;
  name[7] = "iter_stagg_MPE";
  label[7] = fsi_iter_stagg_MPE;
  name[8] = "iter_stagg_RRE";
  label[8] = fsi_iter_stagg_RRE;
  name[9] = "iter_monolithicfluidsplit";
  label[9] = fsi_iter_monolithicfluidsplit;
  name[10] = "iter_monolithicstructuresplit";
  label[10] = fsi_iter_monolithicstructuresplit;
  name[11] = "iter_lung_monolithicstructuresplit";
  label[11] = fsi_iter_lung_monolithicstructuresplit;
  name[12] = "iter_lung_monolithicfluidsplit";
  label[12] = fsi_iter_lung_monolithicfluidsplit;
  name[13] = "iter_xfem_monolithic";
  label[13] = fsi_iter_xfem_monolithic;
  name[14] = "iter_constr_monolithicfluidsplit";
  label[14] = fsi_iter_constr_monolithicfluidsplit;
  name[15] = "iter_constr_monolithicstructuresplit";
  label[15] = fsi_iter_constr_monolithicstructuresplit;
  name[16] = "iter_mortar_monolithicstructuresplit";
  label[16] = fsi_iter_mortar_monolithicstructuresplit;
  name[17] = "iter_mortar_monolithicfluidsplit";
  label[17] = fsi_iter_mortar_monolithicfluidsplit;
  name[18] = "iter_fluidfluid_monolithicstructuresplit";
  label[18] = fsi_iter_fluidfluid_monolithicstructuresplit;
  name[19] = "iter_fluidfluid_monolithicfluidsplit";
  label[19] = fsi_iter_fluidfluid_monolithicfluidsplit;
  name[20] = "iter_fluidfluid_monolithicstructuresplit_nonox";
  label[20] = fsi_iter_fluidfluid_monolithicstructuresplit_nonox;
  name[21] = "iter_fluidfluid_monolithicfluidsplit_nonox";
  label[21] = fsi_iter_fluidfluid_monolithicfluidsplit_nonox;
  name[22] = "iter_sliding_monolithicfluidsplit";
  label[22] = fsi_iter_sliding_monolithicfluidsplit;
  name[23] = "iter_sliding_monolithicstructuresplit";
  label[23] = fsi_iter_sliding_monolithicstructuresplit;
  name[24] = "iter_mortar_monolithicfluidsplit_saddlepoint";
  label[24] = fsi_iter_mortar_monolithicfluidsplit_saddlepoint;


  setStringToIntegralParameter<int>("COUPALGO", "iter_stagg_AITKEN_rel_param",
      "Iteration Scheme over the fields", name, label, &fsidyn);

  setStringToIntegralParameter<int>("DEBUGOUTPUT", "No",
      "Output of unconverged interface values during FSI iteration.\n"
      "There will be a new control file for each time step.\n"
      "This might be helpful to understand the coupling iteration.",
      tuple<std::string>(
          "No", "Yes", "no", "yes", "NO", "YES", "Interface", "Preconditioner", "All"),
      tuple<int>(0, 1, 0, 1, 0, 1, 1, 2, 256), &fsidyn);

  Core::UTILS::BoolParameter("MATCHGRID_FLUIDALE", "Yes", "is matching grid (fluid-ale)", &fsidyn);

  Core::UTILS::BoolParameter(
      "MATCHGRID_STRUCTALE", "Yes", "is matching grid (structure-ale)", &fsidyn);

  Core::UTILS::BoolParameter("MATCHALL", "Yes",
      "is matching grid (fluid-ale) and is full fluid-ale (without euler part)", &fsidyn);

  Core::UTILS::DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &fsidyn);
  Core::UTILS::IntParameter("NUMSTEP", 200, "Total number of Timesteps", &fsidyn);

  Core::UTILS::IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &fsidyn);

  Core::UTILS::BoolParameter("RESTART_FROM_PART_FSI", "No",
      "restart from partitioned fsi (e.g. from prestress calculations) instead of monolithic fsi",
      &fsidyn);

  Core::UTILS::BoolParameter("SECONDORDER", "No",
      "Second order displacement-velocity conversion at the interface.", &fsidyn);

  setStringToIntegralParameter<int>("SLIDEALEPROJ", "None",
      "Projection method to use for sliding FSI.",
      tuple<std::string>("None", "Curr", "Ref", "RotZ", "RotZSphere"),
      tuple<int>(Inpar::FSI::ALEprojection_none, Inpar::FSI::ALEprojection_curr,
          Inpar::FSI::ALEprojection_ref, Inpar::FSI::ALEprojection_rot_z,
          Inpar::FSI::ALEprojection_rot_zsphere),
      &fsidyn);

  Core::UTILS::DoubleParameter("TIMESTEP", 0.1, "Time increment dt", &fsidyn);

  Core::UTILS::IntParameter("RESULTSEVRY", 1, "Increment for writing solution", &fsidyn);

  setStringToIntegralParameter<int>("VERBOSITY", "full", "Verbosity of the FSI problem.",
      tuple<std::string>("full", "medium", "low", "subproblem"),
      tuple<int>(Inpar::FSI::verbosity_full, Inpar::FSI::verbosity_medium,
          Inpar::FSI::verbosity_low, Inpar::FSI::verbosity_subproblem),
      &fsidyn);

  /*----------------------------------------------------------------------*/
  /* parameters for time step size adaptivity in fsi dynamics */
  Teuchos::ParameterList& fsiadapt = fsidyn.sublist("TIMEADAPTIVITY", false, "");

  Core::UTILS::IntParameter("ADAPTSTEPMAX", 5,
      "Maximum number of repetitions of one time step for adapting/reducing the time step size "
      "(>0)",
      &fsiadapt);

  setStringToIntegralParameter<int>("AUXINTEGRATORFLUID", "AB2",
      "Method for error estimation in the fluid field",
      tuple<std::string>("None", "ExplicitEuler", "AB2"),
      tuple<int>(Inpar::FSI::timada_fld_none, Inpar::FSI::timada_fld_expleuler,
          Inpar::FSI::timada_fld_adamsbashforth2),
      &fsiadapt);

  Core::UTILS::StringParameter("AVERAGINGDT", "0.3 0.7",
      "Averaging of time step sizes in case of increasing time step size.\n"
      "Parameters are ordered from most recent weight to the most historic one.\n"
      "Number of parameters determines the number of previous time steps that are involved\n"
      "in the averaging procedure.",
      &fsiadapt);


  setStringToIntegralParameter<int>("DIVERCONT", "stop",
      "What to do if nonlinear solver does not converge?",
      tuple<std::string>("stop", "continue", "halve_step", "revert_dt"),
      tuple<int>(Inpar::FSI::divcont_stop, Inpar::FSI::divcont_continue,
          Inpar::FSI::divcont_halve_step, Inpar::FSI::divcont_revert_dt),
      &fsiadapt);

  Core::UTILS::DoubleParameter(
      "DTMAX", 0.1, "Limit maximally permitted time step size (>0)", &fsiadapt);
  Core::UTILS::DoubleParameter(
      "DTMIN", 1.0e-4, "Limit minimally allowed time step size (>0)", &fsiadapt);

  Core::UTILS::DoubleParameter(
      "LOCERRTOLFLUID", 1.0e-3, "Tolerance for the norm of local velocity error", &fsiadapt);


  Core::UTILS::IntParameter("NUMINCREASESTEPS", 0,
      "Number of consecutive steps that want to increase time step size before\n"
      "actually increasing it. Set 0 to deactivate this feature.",
      &fsiadapt);

  Core::UTILS::DoubleParameter("SAFETYFACTOR", 0.9,
      "This is a safety factor to scale theoretical optimal step size, \n"
      "should be lower than 1 and must be larger than 0",
      &fsiadapt);

  Core::UTILS::DoubleParameter("SIZERATIOMAX", 2.0,
      "Limit maximally permitted change of\n"
      "time step size compared to previous size (>0).",
      &fsiadapt);
  Core::UTILS::DoubleParameter("SIZERATIOMIN", 0.5,
      "Limit minimally permitted change of\n"
      "time step size compared to previous size (>0).",
      &fsiadapt);

  Core::UTILS::BoolParameter(
      "TIMEADAPTON", "No", "Activate or deactivate time step size adaptivity", &fsiadapt);

  /*--------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------*/
  /* parameters for monolithic FSI solvers */
  Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER", false, "");

  Core::UTILS::DoubleParameter("ADAPTIVEDIST", 0.0,
      "Required distance for adaptive convergence check in Newton-type FSI.\n"
      "This is the improvement we want to achieve in the linear extrapolation of the\n"
      "adaptive convergence check. Set to zero to avoid the adaptive check altogether.",
      &fsimono);

  Core::UTILS::DoubleParameter("BASETOL", 1e-3,
      "Basic tolerance for adaptive convergence check in monolithic FSI.\n"
      "This tolerance will be used for the linear solve of the FSI block system.\n"
      "The linear convergence test will always use the relative residual norm (AZ_r0).\n"
      "Not to be confused with the Newton tolerance (CONVTOL) that applies\n"
      "to the nonlinear convergence test using a absolute residual norm.",
      &fsimono);

  Core::UTILS::DoubleParameter("CONVTOL", 1e-6,
      "Nonlinear tolerance for lung/constraint/fluid-fluid FSI",
      &fsimono);  // ToDo remove

  Core::UTILS::BoolParameter("ENERGYFILE", "No",
      "Write artificial interface energy due to temporal discretization to file", &fsimono);

  Core::UTILS::BoolParameter(
      "FSIAMGANALYZE", "No", "run analysis on fsiamg multigrid scheme", &fsimono);

  setStringToIntegralParameter<int>("HYBRID_AS_TYPE", "ILU",
      "Type of inverse in additive part of hybrid preconditioner.",
      tuple<std::string>("ILU", "Amesos_LU"),
      tuple<int>(Inpar::FSI::hybrid_as_type_ILU, Inpar::FSI::hybrid_as_type_Amesos_LU), &fsimono);

  Core::UTILS::IntParameter(
      "HYBRID_FILL_LEVEL", 0, "Level of fill of the hybrid ILU preconditioner.", &fsimono);

  Core::UTILS::BoolParameter("HYBRIDFULL", "yes",
      "Apply Additive Schwarz Preconditioner on all procs (yes)\n"
      "or on interface procs only (no)",
      &fsimono);

  Core::UTILS::DoubleParameter("HYBRID_IMBALANCE_TOL", 1.1,
      "Relative imbalance of sudomain size for repartitioning in hybrid preconditioner", &fsimono);

  Core::UTILS::BoolParameter("INFNORMSCALING", "Yes", "Scale Blocks with row infnorm?", &fsimono);

  setStringToIntegralParameter<int>("INNERPREC", "PreconditionedKrylov",
      "Inner preconditioner used in a hybrid Schwarz setting.",
      tuple<std::string>("PreconditionedKrylov"), tuple<int>(Inpar::FSI::PreconditionedKrylov),
      &fsimono);

  Core::UTILS::IntParameter(
      "ITEMAX", 100, "Maximum allowed number of nonlinear iterations", &fsimono);

  Core::UTILS::IntParameter("KRYLOV_ITEMAX", 1000, "Max Iterations for linear solver.", &fsimono);

  Core::UTILS::IntParameter("KRYLOV_SIZE", 50, "Size of Krylov Subspace.", &fsimono);

  setStringToIntegralParameter<int>("LINEARBLOCKSOLVER", "PreconditionedKrylov",
      "Linear block preconditioner for block system in monolithic FSI.",
      tuple<std::string>("PreconditionedKrylov", "HybridSchwarz", "LinalgSolver"),
      tuple<int>(
          Inpar::FSI::PreconditionedKrylov, Inpar::FSI::HybridSchwarz, Inpar::FSI::LinalgSolver),
      &fsimono);

  Core::UTILS::IntParameter("LINEAR_SOLVER", -1,
      "Number of SOLVER block describing the linear solver and preconditioner", &fsimono);

  // Iteration parameters for convergence check of newton loop
  // for implementations without NOX
  setStringToIntegralParameter<int>("NORM_INC", "Rel",
      "type of norm for primary variables convergence check",
      tuple<std::string>("Abs", "Rel", "Mix"),
      tuple<int>(Inpar::FSI::convnorm_abs, Inpar::FSI::convnorm_rel, Inpar::FSI::convnorm_mix),
      &fsimono);

  // for implementations without NOX
  setStringToIntegralParameter<int>("NORM_RESF", "Rel",
      "type of norm for residual convergence check", tuple<std::string>("Abs", "Rel", "Mix"),
      tuple<int>(Inpar::FSI::convnorm_abs, Inpar::FSI::convnorm_rel, Inpar::FSI::convnorm_mix),
      &fsimono);

  // for implementations without NOX
  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC", "And",
      "binary operator to combine primary variables and residual force values",
      tuple<std::string>("And"), tuple<int>(Inpar::FSI::bop_and), &fsimono);

  Core::UTILS::IntParameter("PRECONDREUSE", 0,
      "Number of iterations in one time step reusing the preconditioner before rebuilding it",
      &fsimono);

  Core::UTILS::BoolParameter("REBUILDPRECEVERYSTEP", "Yes",
      "Enforce rebuilding the preconditioner at the beginning of every time step", &fsimono);

  setStringToIntegralParameter<int>("REDISTRIBUTE", "off", "Redistribute domain decomposition.",
      tuple<std::string>("off", "structure", "fluid", "both", "monolithic"),
      tuple<int>(Inpar::FSI::Redistribute_off, Inpar::FSI::Redistribute_structure,
          Inpar::FSI::Redistribute_fluid, Inpar::FSI::Redistribute_both,
          Inpar::FSI::Redistribute_monolithic),
      &fsimono);

  Core::UTILS::DoubleParameter(
      "REDIST_WEIGHT1", 1.0, "Weight for redistribution, general domain.", &fsimono);
  Core::UTILS::DoubleParameter(
      "REDIST_WEIGHT2", 50000.0, "Weight for redistribution, interface domain.", &fsimono);
  Core::UTILS::DoubleParameter(
      "REDIST_SECONDWEIGHT1", -1.0, "Weight for 2nd redistribution, general domain.", &fsimono);
  Core::UTILS::DoubleParameter(
      "REDIST_SECONDWEIGHT2", -1.0, "Weight for 2nd redistribution, interface domain.", &fsimono);

  Core::UTILS::BoolParameter("SHAPEDERIVATIVES", "No",
      "Include linearization with respect to mesh movement in Navier Stokes equation.", &fsimono);

  Core::UTILS::BoolParameter(
      "SYMMETRICPRECOND", "No", "Symmetric block GS preconditioner or ordinary GS", &fsimono);

  // monolithic preconditioner parameter

  Core::UTILS::StringParameter("ALEPCOMEGA", "1.0 1.0 1.0 1.0",
      "Relaxation factor for Richardson iteration on ale block in MFSI block preconditioner\n"
      "FSIAMG: each number belongs to a level\n"
      "PreconditiondKrylov: only first number is used for finest level",
      &fsimono);
  Core::UTILS::StringParameter("ALEPCITER", "1 1 1 1",
      "Number of Richardson iterations on ale block in MFSI block preconditioner\n"
      "FSIAMG: each number belongs to a level\n"
      "PreconditiondKrylov: only first number is used for finest level",
      &fsimono);
  Core::UTILS::StringParameter("FLUIDPCOMEGA", "1.0 1.0 1.0 1.0",
      "Relaxation factor for Richardson iteration on fluid block in MFSI block preconditioner\n"
      "FSIAMG: each number belongs to a level\n"
      "PreconditiondKrylov: only first number is used for finest level",
      &fsimono);
  Core::UTILS::StringParameter("FLUIDPCITER", "1 1 1 1",
      "Number of Richardson iterations on fluid block in MFSI block preconditioner\n"
      "FSIAMG: each number belongs to a level\n"
      "PreconditiondKrylov: only first number is used for finest level",
      &fsimono);
  Core::UTILS::StringParameter("STRUCTPCOMEGA", "1.0 1.0 1.0 1.0",
      "Relaxation factor for Richardson iteration on structural block in MFSI block "
      "preconditioner\n"
      "FSIAMG: each number belongs to a level\n"
      "PreconditiondKrylov: only first number is used for finest level",
      &fsimono);
  Core::UTILS::StringParameter("STRUCTPCITER", "1 1 1 1",
      "Number of Richardson iterations on structural block in MFSI block preconditioner\n"
      "FSIAMG: each number belongs to a level\n"
      "PreconditiondKrylov: only first number is used for finest level",
      &fsimono);

  Core::UTILS::StringParameter("PCOMEGA", "1.0 1.0 1.0",
      "Relaxation factor for Richardson iteration on whole MFSI block preconditioner\n"
      "FSIAMG: each number belongs to a level\n"
      "PreconditiondKrylov: only first number is used for finest level",
      &fsimono);
  Core::UTILS::StringParameter("PCITER", "1 1 1",
      "Number of Richardson iterations on whole MFSI block preconditioner\n"
      "FSIAMG: each number belongs to a level\n"
      "PreconditiondKrylov: only first number is used for finest level",
      &fsimono);

  Core::UTILS::StringParameter(
      "BLOCKSMOOTHER", "BGS BGS BGS", "Type of block smoother, can be BGS or Schur", &fsimono);

  Core::UTILS::StringParameter(
      "SCHUROMEGA", "0.001 0.01 0.1", "Damping factor for Schur complement construction", &fsimono);

  // tolerances for convergence check of nonlinear solver in monolithic FSI
  // structure displacements
  Core::UTILS::DoubleParameter("TOL_DIS_RES_L2", 1e-6,
      "Absolute tolerance for structure displacement residual in L2-norm", &fsimono);
  Core::UTILS::DoubleParameter("TOL_DIS_RES_INF", 1e-6,
      "Absolute tolerance for structure displacement residual in Inf-norm", &fsimono);
  Core::UTILS::DoubleParameter("TOL_DIS_INC_L2", 1e-6,
      "Absolute tolerance for structure displacement increment in L2-norm", &fsimono);
  Core::UTILS::DoubleParameter("TOL_DIS_INC_INF", 1e-6,
      "Absolute tolerance for structure displacement increment in Inf-norm", &fsimono);
  // interface tolerances
  Core::UTILS::DoubleParameter(
      "TOL_FSI_RES_L2", 1e-6, "Absolute tolerance for interface residual in L2-norm", &fsimono);
  Core::UTILS::DoubleParameter(
      "TOL_FSI_RES_INF", 1e-6, "Absolute tolerance for interface residual in Inf-norm", &fsimono);
  Core::UTILS::DoubleParameter(
      "TOL_FSI_INC_L2", 1e-6, "Absolute tolerance for interface increment in L2-norm", &fsimono);
  Core::UTILS::DoubleParameter(
      "TOL_FSI_INC_INF", 1e-6, "Absolute tolerance for interface increment in Inf-norm", &fsimono);
  // fluid pressure
  Core::UTILS::DoubleParameter("TOL_PRE_RES_L2", 1e-6,
      "Absolute tolerance for fluid pressure residual in L2-norm", &fsimono);
  Core::UTILS::DoubleParameter("TOL_PRE_RES_INF", 1e-6,
      "Absolute tolerance for fluid pressure residual in Inf-norm", &fsimono);
  Core::UTILS::DoubleParameter("TOL_PRE_INC_L2", 1e-6,
      "Absolute tolerance for fluid pressure increment in L2-norm", &fsimono);
  Core::UTILS::DoubleParameter("TOL_PRE_INC_INF", 1e-6,
      "Absolute tolerance for fluid pressure increment in Inf-norm", &fsimono);
  // fluid velocities
  Core::UTILS::DoubleParameter("TOL_VEL_RES_L2", 1e-6,
      "Absolute tolerance for fluid velocity residual in L2-norm", &fsimono);
  Core::UTILS::DoubleParameter("TOL_VEL_RES_INF", 1e-6,
      "Absolute tolerance for fluid velocity residual in Inf-norm", &fsimono);
  Core::UTILS::DoubleParameter("TOL_VEL_INC_L2", 1e-6,
      "Absolute tolerance for fluid velocity increment in L2-norm", &fsimono);
  Core::UTILS::DoubleParameter("TOL_VEL_INC_INF", 1e-6,
      "Absolute tolerance for fluid velocity increment in Inf-norm", &fsimono);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned FSI solvers */
  Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER", false, "");

  Core::UTILS::DoubleParameter("BASETOL", 1e-3,
      "Basic tolerance for adaptive convergence check in monolithic FSI.\n"
      "This tolerance will be used for the linear solve of the FSI block system.\n"
      "The linear convergence test will always use the relative residual norm (AZ_r0).\n"
      "Not to be confused with the Newton tolerance (CONVTOL) that applies\n"
      "to the nonlinear convergence test using a absolute residual norm.",
      &fsipart);

  Core::UTILS::DoubleParameter("CONVTOL", 1e-6,
      "Tolerance for iteration over fields in case of partitioned scheme", &fsipart);


  setStringToIntegralParameter<int>("COUPMETHOD", "conforming",
      "Coupling Method Mortar (mtr) or conforming nodes at interface",
      tuple<std::string>(  // ToDO introduce enum
          "MTR", "Mtr", "mtr", "conforming", "immersed"),
      tuple<int>(0, 0, 0, 1, 2), &fsipart);

  setStringToIntegralParameter<int>("COUPVARIABLE", "Displacement",

      "Coupling variable at the interface", tuple<std::string>("Displacement", "Force", "Velocity"),
      tuple<int>(Inpar::FSI::CoupVarPart::disp, Inpar::FSI::CoupVarPart::force,
          Inpar::FSI::CoupVarPart::vel),
      &fsipart);

  Core::UTILS::BoolParameter("DIVPROJECTION", "no",
      "Project velocity into divergence-free subspace for partitioned fsi", &fsipart);

  Core::UTILS::IntParameter("ITEMAX", 100, "Maximum number of iterations over fields", &fsipart);

  Core::UTILS::DoubleParameter("MAXOMEGA", 0.0,
      "largest omega allowed for Aitken relaxation (0.0 means no constraint)", &fsipart);

  Core::UTILS::DoubleParameter(
      "MINOMEGA", -1.0, "smallest omega allowed for Aitken relaxation (default is -1.0)", &fsipart);


  setStringToIntegralParameter<int>("PARTITIONED", "DirichletNeumann",
      "Coupling strategies for partitioned FSI solvers.",
      tuple<std::string>(
          "DirichletNeumann", "DirichletNeumannSlideALE", "DirichletNeumannVolCoupl"),
      tuple<int>(Inpar::FSI::DirichletNeumann, Inpar::FSI::DirichletNeumannSlideale,
          Inpar::FSI::DirichletNeumannVolCoupl),
      &fsipart);

  setStringToIntegralParameter<int>("PREDICTOR", "d(n)", "Predictor for interface displacements",
      tuple<std::string>(
          "d(n)", "d(n)+dt*(1.5*v(n)-0.5*v(n-1))", "d(n)+dt*v(n)", "d(n)+dt*v(n)+0.5*dt^2*a(n)"),
      tuple<int>(1, 2, 3, 4), &fsipart);

  Core::UTILS::DoubleParameter(
      "RELAX", 1.0, "fixed relaxation parameter for partitioned FSI solvers", &fsipart);

  /* ----------------------------------------------------------------------- */
  Teuchos::ParameterList& constrfsi = fsidyn.sublist("CONSTRAINT", false, "");

  setStringToIntegralParameter<int>("PRECONDITIONER", "Simple", "preconditioner to use",
      tuple<std::string>("Simple", "Simplec"), tuple<int>(Inpar::FSI::Simple, Inpar::FSI::Simplec),
      &constrfsi);
  Core::UTILS::IntParameter("SIMPLEITER", 2, "Number of iterations for simple pc", &constrfsi);
  Core::UTILS::DoubleParameter("ALPHA", 0.8, "alpha parameter for simple pc", &constrfsi);
}

/*----------------------------------------------------------------------------*/
void Inpar::FSI::SetValidConditions(
    std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist)
{
  using namespace Input;

  std::vector<Teuchos::RCP<Input::LineComponent>> fsicomponents;

  fsicomponents.push_back(Teuchos::rcp(new Input::IntComponent("coupling id")));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linefsi =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN FSI COUPLING LINE CONDITIONS",
          "FSICoupling", "FSI Coupling", Core::Conditions::FSICoupling, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surffsi =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition("DESIGN FSI COUPLING SURF CONDITIONS",
          "FSICoupling", "FSI Coupling", Core::Conditions::FSICoupling, true,
          Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < fsicomponents.size(); ++i)
  {
    linefsi->add_component(fsicomponents[i]);
    surffsi->add_component(fsicomponents[i]);
  }

  condlist.push_back(linefsi);
  condlist.push_back(surffsi);

  /*--------------------------------------------------------------------*/
  // FSI without sliding

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linefsins = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN FSI COUPLING NO SLIDE LINE CONDITIONS",
          "FSICouplingNoSlide", "FSI Coupling No Slide", Core::Conditions::FSICouplingNoSlide, true,
          Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surffsins = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN FSI COUPLING NO SLIDE SURF CONDITIONS",
          "FSICouplingNoSlide", "FSI Coupling No Slide", Core::Conditions::FSICouplingNoSlide, true,
          Core::Conditions::geometry_type_surface));

  condlist.push_back(linefsins);
  condlist.push_back(surffsins);

  /*--------------------------------------------------------------------*/
  // FSI define centerdisp for sliding interfaces

  Teuchos::RCP<Core::Conditions::ConditionDefinition> linefsicd = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN FSI COUPLING CENTER DISP LINE CONDITIONS",
          "FSICouplingCenterDisp", "FSI Coupling Center Disp",
          Core::Conditions::FSICouplingCenterDisp, true, Core::Conditions::geometry_type_line));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> surffsicd = Teuchos::rcp(
      new Core::Conditions::ConditionDefinition("DESIGN FSI COUPLING CENTER DISP SURF CONDITIONS",
          "FSICouplingCenterDisp", "FSI Coupling Center Disp",
          Core::Conditions::FSICouplingCenterDisp, true, Core::Conditions::geometry_type_surface));

  condlist.push_back(linefsicd);
  condlist.push_back(surffsicd);

  /*--------------------------------------------------------------------*/
  // Additional coupling of structure and ale fields (for lung fsi)

  std::vector<Teuchos::RCP<Input::LineComponent>> saccomponents;

  saccomponents.push_back(Teuchos::rcp(new Input::IntComponent("coupling id")));
  saccomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("field", "structure",
      Teuchos::tuple<std::string>("structure", "fluid"),
      Teuchos::tuple<std::string>("structure", "fluid"))));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfsac =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN STRUCTURE ALE COUPLING SURF CONDITIONS", "StructAleCoupling", "StructAleCoupling",
          Core::Conditions::StructAleCoupling, true, Core::Conditions::geometry_type_surface));

  for (unsigned i = 0; i < saccomponents.size(); ++i) surfsac->add_component(saccomponents[i]);

  condlist.push_back(surfsac);

  /*--------------------------------------------------------------------*/
  // Additional coupling of structure and fluid volumes (for lung fsi)

  std::vector<Teuchos::RCP<Input::LineComponent>> sfvcomponents;

  sfvcomponents.push_back(Teuchos::rcp(new Input::IntComponent("coupling id")));
  sfvcomponents.push_back(Teuchos::rcp(new Input::SelectionComponent("field", "structure",
      Teuchos::tuple<std::string>("structure", "fluid"),
      Teuchos::tuple<std::string>("structure", "fluid"))));

  Teuchos::RCP<Core::Conditions::ConditionDefinition> surfsfv =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN STRUCTURE FLUID VOLUME COUPLING SURF CONDITIONS", "StructFluidSurfCoupling",
          "StructFluidSurfCoupling", Core::Conditions::StructFluidSurfCoupling, true,
          Core::Conditions::geometry_type_surface));
  Teuchos::RCP<Core::Conditions::ConditionDefinition> volsfv =
      Teuchos::rcp(new Core::Conditions::ConditionDefinition(
          "DESIGN STRUCTURE FLUID VOLUME COUPLING VOL CONDITIONS", "StructFluidVolCoupling",
          "StructFluidVolCoupling", Core::Conditions::StructFluidVolCoupling, false,
          Core::Conditions::geometry_type_volume));

  for (unsigned i = 0; i < sfvcomponents.size(); ++i)
  {
    surfsfv->add_component(sfvcomponents[i]);
    volsfv->add_component(sfvcomponents[i]);
  }

  condlist.push_back(surfsfv);
  condlist.push_back(volsfv);
}

FOUR_C_NAMESPACE_CLOSE
