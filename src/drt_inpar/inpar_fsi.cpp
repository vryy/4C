/*----------------------------------------------------------------------------*/
/*!
\file inpar_fsi.cpp

\maintainer Matthias Mayr

\brief Input parameters for fluid structure interaction

*/
/*----------------------------------------------------------------------------*/

#include "drt_validparameters.H"
#include "inpar_fsi.H"
#include "../drt_lib/drt_conditiondefinition.H"

/*----------------------------------------------------------------------------*/
void INPAR::FSI::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::ParameterList& fsidyn = list->sublist("FSI DYNAMIC", false,
      "Fluid Structure Interaction\n"
          "FSI solver with various coupling methods");

  Teuchos::Tuple<std::string, 25> name;
  Teuchos::Tuple<int, 25> label;

  name[ 0] = "basic_sequ_stagg";                              label[ 0] = fsi_basic_sequ_stagg;
  name[ 1] = "iter_stagg_fixed_rel_param";                    label[ 1] = fsi_iter_stagg_fixed_rel_param;
  name[ 2] = "iter_stagg_AITKEN_rel_param";                   label[ 2] = fsi_iter_stagg_AITKEN_rel_param;
  name[ 3] = "iter_stagg_steep_desc";                         label[ 3] = fsi_iter_stagg_steep_desc;
  name[ 4] = "iter_stagg_NLCG";                               label[ 4] = fsi_iter_stagg_NLCG;
  name[ 5] = "iter_stagg_MFNK_FD";                            label[ 5] = fsi_iter_stagg_MFNK_FD;
  name[ 6] = "iter_stagg_MFNK_FSI";                           label[ 6] = fsi_iter_stagg_MFNK_FSI;
  name[ 7] = "iter_stagg_MPE";                                label[ 7] = fsi_iter_stagg_MPE;
  name[ 8] = "iter_stagg_RRE";                                label[ 8] = fsi_iter_stagg_RRE;
  name[ 9] = "iter_monolithicfluidsplit";                     label[ 9] = fsi_iter_monolithicfluidsplit;
  name[10] = "iter_monolithicstructuresplit";                 label[10] = fsi_iter_monolithicstructuresplit;
  name[11] = "iter_lung_monolithicstructuresplit";            label[11] = fsi_iter_lung_monolithicstructuresplit;
  name[12] = "iter_lung_monolithicfluidsplit";                label[12] = fsi_iter_lung_monolithicfluidsplit;
  name[13] = "iter_xfem_monolithic";                          label[13] = fsi_iter_xfem_monolithic;
  name[14] = "pseudo_structure";                              label[14] = fsi_pseudo_structureale;
  name[15] = "iter_constr_monolithicfluidsplit";              label[15] = fsi_iter_constr_monolithicfluidsplit;
  name[16] = "iter_constr_monolithicstructuresplit";          label[16] = fsi_iter_constr_monolithicstructuresplit;
  name[17] = "iter_mortar_monolithicstructuresplit";          label[17] = fsi_iter_mortar_monolithicstructuresplit;
  name[18] = "iter_mortar_monolithicfluidsplit";              label[18] = fsi_iter_mortar_monolithicfluidsplit;
  name[19] = "iter_fluidfluid_monolithicstructuresplit";      label[19] = fsi_iter_fluidfluid_monolithicstructuresplit;
  name[20] = "iter_fluidfluid_monolithicfluidsplit";          label[20] = fsi_iter_fluidfluid_monolithicfluidsplit;
  name[21] = "iter_fluidfluid_monolithicstructuresplit_nonox";label[21] = fsi_iter_fluidfluid_monolithicstructuresplit_nonox;
  name[22] = "iter_fluidfluid_monolithicfluidsplit_nonox";    label[22] = fsi_iter_fluidfluid_monolithicfluidsplit_nonox;
  name[23] = "iter_sliding_monolithicfluidsplit";             label[23] = fsi_iter_sliding_monolithicfluidsplit;
  name[24] = "iter_sliding_monolithicstructuresplit";         label[24] = fsi_iter_sliding_monolithicstructuresplit;


  setStringToIntegralParameter<int>("COUPALGO", "iter_stagg_AITKEN_rel_param",
      "Iteration Scheme over the fields", name, label, &fsidyn);

  setStringToIntegralParameter<int>("DEBUGOUTPUT","No",
      "Output of unconverged interface values during FSI iteration.\n"
      "There will be a new control file for each time step.\n"
      "This might be helpful to understand the coupling iteration.",
      tuple<std::string>(
        "No",
        "Yes",
        "no",
        "yes",
        "NO",
        "YES",
        "Interface",
        "Preconditioner",
        "All"
        ),
      tuple<int>(
        0,
        1,
        0,
        1,
        0,
        1,
        1,
        2,
        256
        ),
      &fsidyn);

  BoolParameter("MATCHGRID_FLUIDALE", "Yes", "is matching grid (fluid-ale)",
      &fsidyn);

  BoolParameter("MATCHALL", "Yes", "is matching grid (fluid-ale) and is full fluid-ale (without euler part)",
      &fsidyn);

  DoubleParameter("MAXTIME", 1000.0, "Total simulation time", &fsidyn);
  IntParameter("NUMSTEP", 200, "Total number of Timesteps", &fsidyn);

  IntParameter("RESTARTEVRY", 1, "Increment for writing restart", &fsidyn);

  BoolParameter("RESTART_FROM_PART_FSI", "No",
      "restart from partitioned fsi (e.g. from prestress calculations) instead of monolithic fsi",
        &fsidyn);

  BoolParameter("SECONDORDER", "No",
      "Second order displacement-velocity conversion at the interface.",
      &fsidyn);

  setStringToIntegralParameter<int>("SLIDEALEPROJ","None",
                                   "Projection method to use for sliding FSI.",
                                   tuple<std::string>(
                                       "None",
                                       "Curr",
                                       "Ref",
                                       "RotZ",
                                       "RotZSphere"),
                                   tuple<int>(
                                       INPAR::FSI::ALEprojection_none,
                                       INPAR::FSI::ALEprojection_curr,
                                       INPAR::FSI::ALEprojection_ref,
                                       INPAR::FSI::ALEprojection_rot_z,
                                       INPAR::FSI::ALEprojection_rot_zsphere),
                                   &fsidyn);

  DoubleParameter("TIMESTEP",0.1,"Time increment dt",&fsidyn);

  IntParameter("RESULTSEVRY",1,"Increment for writing solution",&fsidyn);

  setStringToIntegralParameter<int>("VERBOSITY","full",
                                 "Verbosity of the FSI problem.",
                                 tuple<std::string>(
                                     "full",
                                     "medium",
                                     "low",
                                     "subproblem"),
                                 tuple<int>(
                                     INPAR::FSI::verbosity_full,
                                     INPAR::FSI::verbosity_medium,
                                     INPAR::FSI::verbosity_low,
                                     INPAR::FSI::verbosity_subproblem),
                                 &fsidyn);

  /*----------------------------------------------------------------------*/
  /* parameters for time step size adaptivity in fsi dynamics */
  Teuchos::ParameterList& fsiadapt = fsidyn.sublist("TIMEADAPTIVITY",false,"");

  IntParameter("ADAPTSTEPMAX", 5,
      "Maximum number of repetitions of one time step for adapting/reducing the time step size (>0)",
      &fsiadapt);

  setStringToIntegralParameter<int>("AUXINTEGRATORFLUID", "AB2",
                                    "Method for error estimation in the fluid field",
                                    tuple<std::string>(
                                        "None",
                                        "ExplicitEuler",
                                        "AB2"),
                                    tuple<int>(
                                        INPAR::FSI::timada_fld_none,
                                        INPAR::FSI::timada_fld_expleuler,
                                        INPAR::FSI::timada_fld_adamsbashforth2),
                                    &fsiadapt);

  setNumericStringParameter("AVERAGINGDT", "0.3 0.7",
        "Averaging of time step sizes in case of increasing time step size.\n"
        "Parameters are ordered from most recent weight to the most historic one.\n"
        "Number of parameters determines the number of previous time steps that are involved\n"
        "in the averaging procedure.", &fsiadapt);


  setStringToIntegralParameter<int>("DIVERCONT", "stop",
        "What to do if nonlinear solver does not converge?",
        tuple<std::string>(
            "stop",
            "continue",
            "halve_step",
            "revert_dt"),
        tuple<int>(
            INPAR::FSI::divcont_stop,
            INPAR::FSI::divcont_continue,
            INPAR::FSI::divcont_halve_step,
            INPAR::FSI::divcont_revert_dt),
        &fsiadapt);

  DoubleParameter("DTMAX", 0.1, "Limit maximally permitted time step size (>0)",
      &fsiadapt);
  DoubleParameter("DTMIN", 1.0e-4,
      "Limit minimally allowed time step size (>0)", &fsiadapt);

  DoubleParameter("LOCERRTOLFLUID", 1.0e-3,
      "Tolerance for the norm of local velocity error", &fsiadapt);


  IntParameter("NUMINCREASESTEPS", 0,
      "Number of consecutive steps that want to increase time step size before\n"
      "actually increasing it. Set 0 to deactivate this feature.",
      &fsiadapt);

  DoubleParameter("SAFETYFACTOR", 0.9,
      "This is a safety factor to scale theoretical optimal step size, \n"
      "should be lower than 1 and must be larger than 0", &fsiadapt);

  DoubleParameter("SIZERATIOMAX", 2.0, "Limit maximally permitted change of\n"
      "time step size compared to previous size (>0).", &fsiadapt);
  DoubleParameter("SIZERATIOMIN", 0.5, "Limit minimally permitted change of\n"
      "time step size compared to previous size (>0).", &fsiadapt);

  BoolParameter("TIMEADAPTON", "No",
      "Activate or deactivate time step size adaptivity", &fsiadapt);

  /*--------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------*/
  /* parameters for monolithic FSI solvers */
  Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER",false,"");

  DoubleParameter("ADAPTIVEDIST", 0.0,
      "Required distance for adaptive convergence check in Newton-type FSI.\n"
          "This is the improvement we want to achieve in the linear extrapolation of the\n"
          "adaptive convergence check. Set to zero to avoid the adaptive check altogether.",
      &fsimono);

  DoubleParameter("BASETOL", 1e-3,
      "Basic tolerance for adaptive convergence check in monolithic FSI.\n"
          "This tolerance will be used for the linear solve of the FSI block system.\n"
          "The linear convergence test will always use the relative residual norm (AZ_r0).\n"
          "Not to be confused with the Newton tolerance (CONVTOL) that applies\n"
          "to the nonlinear convergence test using a absolute residual norm.",
      &fsimono);

  DoubleParameter("CONVTOL", 1e-6,
      "Nonlinear tolerance for lung/constraint/fluid-fluid FSI", &fsimono); // ToDo remove

  BoolParameter("ENERGYFILE", "No",
      "Write artificial interface energy due to temporal discretization to file",
      &fsimono);

  BoolParameter("FSIAMGANALYZE", "No",
      "run analysis on fsiamg multigrid scheme", &fsimono);

  setStringToIntegralParameter<int>(
        "HYBRID_AS_TYPE", "ILU",
        "Type of inverse in additive part of hybrid preconditioner.",
        tuple<std::string>(
            "ILU",
            "Amesos_LU"
            ),
        tuple<int>(
            INPAR::FSI::hybrid_as_type_ILU,
            INPAR::FSI::hybrid_as_type_Amesos_LU
            ),
        &fsimono);

  IntParameter("HYBRID_FILL_LEVEL", 0,
      "Level of fill of the hybrid ILU preconditioner.", &fsimono);

  BoolParameter("HYBRIDFULL", "yes",
      "Apply Additive Schwarz Preconditioner on all procs (yes)\n"
          "or on interface procs only (no)", &fsimono);

  BoolParameter("INFNORMSCALING", "Yes", "Scale Blocks with row infnorm?",
      &fsimono);

  setStringToIntegralParameter<int>(
      "INNERPREC", "PreconditionedKrylov",
      "Inner preconditioner used in a hybrid Schwarz setting.",
      tuple<std::string>(
          "PreconditionedKrylov",
          "FSIAMG"
          ),
      tuple<int>(
          INPAR::FSI::PreconditionedKrylov,
          INPAR::FSI::FSIAMG
          ),
      &fsimono);

  IntParameter("ITEMAX", 100, "Maximum allowed number of nonlinear iterations",
      &fsimono);

  IntParameter("KRYLOV_ITEMAX", 1000,
      "Max Iterations for linear solver.", &fsimono);

  IntParameter("KRYLOV_SIZE", 50,
      "Size of Krylov Subspace.", &fsimono);

  setStringToIntegralParameter<int>(
      "LINEARBLOCKSOLVER", "PreconditionedKrylov",
      "Linear block preconditioner for block system in monolithic FSI.",
      tuple<std::string>(
          "PreconditionedKrylov",
          "FSIAMG",
          "AMGnxn",
          "HybridSchwarz"
          ),
      tuple<int>(
          INPAR::FSI::PreconditionedKrylov,
          INPAR::FSI::FSIAMG,
          INPAR::FSI::AMGnxn,
          INPAR::FSI::HybridSchwarz
          ),
      &fsimono);

  // Iteration parameters for convergence check of newton loop
  // for implementations without NOX
  setStringToIntegralParameter<int>("NORM_INC", "Rel",
      "type of norm for primary variables convergence check",
      tuple<std::string>(
        "Abs",
        "Rel",
        "Mix"
        ),
      tuple<int>(
        INPAR::FSI::convnorm_abs,
        INPAR::FSI::convnorm_rel,
        INPAR::FSI::convnorm_mix
        ),
      &fsimono);

  // for implementations without NOX
  setStringToIntegralParameter<int>("NORM_RESF", "Rel",
      "type of norm for residual convergence check",
      tuple<std::string>(
        "Abs",
        "Rel",
        "Mix"
        ),
      tuple<int>(
        INPAR::FSI::convnorm_abs,
        INPAR::FSI::convnorm_rel,
        INPAR::FSI::convnorm_mix
        ),
      &fsimono);

  // for implementations without NOX
  setStringToIntegralParameter<int>("NORMCOMBI_RESFINC", "And",
      "binary operator to combine primary variables and residual force values",
      tuple<std::string>("And"), tuple<int>(INPAR::FSI::bop_and), &fsimono);

  IntParameter("PRECONDREUSE", 0,
      "Number of iterations in one time step reusing the preconditioner before rebuilding it",
      &fsimono);

  BoolParameter("REBUILDPRECEVERYSTEP", "Yes",
      "Enforce rebuilding the preconditioner at the beginning of every time step",
      &fsimono);

  setStringToIntegralParameter<int>(
      "REDISTRIBUTE", "off",
      "Redistribute domain decomposition.",
      tuple<std::string>(
          "off",
          "structure",
          "fluid",
          "both",
          "monolithic"),
      tuple<int>(
          INPAR::FSI::Redistribute_off,
          INPAR::FSI::Redistribute_structure,
          INPAR::FSI::Redistribute_fluid,
          INPAR::FSI::Redistribute_both,
          INPAR::FSI::Redistribute_monolithic
          ),
      &fsimono);

  DoubleParameter("REDIST_WEIGHT1",1.0,"Weight for redistribution, general domain.",&fsimono);
  DoubleParameter("REDIST_WEIGHT2",50000.0,"Weight for redistribution, interface domain.",&fsimono);
  DoubleParameter("REDIST_SECONDWEIGHT1",-1.0,"Weight for 2nd redistribution, general domain.",&fsimono);
  DoubleParameter("REDIST_SECONDWEIGHT2",-1.0,"Weight for 2nd redistribution, interface domain.",&fsimono);

  BoolParameter("SHAPEDERIVATIVES", "No",
      "Include linearization with respect to mesh movement in Navier Stokes equation.",
      &fsimono);

  BoolParameter("SYMMETRICPRECOND", "No",
      "Symmetric block GS preconditioner or ordinary GS", &fsimono);

  // monolithic preconditioner parameter

  setNumericStringParameter("ALEPCOMEGA","1.0 1.0 1.0 1.0",
                            "Relaxation factor for Richardson iteration on ale block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("ALEPCITER","1 1 1 1",
                            "Number of Richardson iterations on ale block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("FLUIDPCOMEGA","1.0 1.0 1.0 1.0",
                            "Relaxation factor for Richardson iteration on fluid block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("FLUIDPCITER","1 1 1 1",
                            "Number of Richardson iterations on fluid block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("STRUCTPCOMEGA","1.0 1.0 1.0 1.0",
                            "Relaxation factor for Richardson iteration on structural block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("STRUCTPCITER","1 1 1 1",
                            "Number of Richardson iterations on structural block in MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);

  setNumericStringParameter("PCOMEGA","1.0 1.0 1.0",
                            "Relaxation factor for Richardson iteration on whole MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);
  setNumericStringParameter("PCITER","1 1 1",
                            "Number of Richardson iterations on whole MFSI block preconditioner\n"
                            "FSIAMG: each number belongs to a level\n"
                            "PreconditiondKrylov: only first number is used for finest level",
                            &fsimono);

  StringParameter("BLOCKSMOOTHER", "BGS BGS BGS",
      "Type of block smoother, can be BGS or Schur", &fsimono);

  setNumericStringParameter("SCHUROMEGA", "0.001 0.01 0.1",
      "Damping factor for Schur complement construction", &fsimono);

  // tolerances for convergence check of nonlinear solver in monolithic FSI
  // structure displacements
  DoubleParameter("TOL_DIS_RES_L2",1e-6,"Absolute tolerance for structure displacement residual in L2-norm",&fsimono);
  DoubleParameter("TOL_DIS_RES_INF",1e-6,"Absolute tolerance for structure displacement residual in Inf-norm",&fsimono);
  DoubleParameter("TOL_DIS_INC_L2",1e-6,"Absolute tolerance for structure displacement increment in L2-norm",&fsimono);
  DoubleParameter("TOL_DIS_INC_INF",1e-6,"Absolute tolerance for structure displacement increment in Inf-norm",&fsimono);
  // interface tolerances
  DoubleParameter("TOL_FSI_RES_L2",1e-6,"Absolute tolerance for interface residual in L2-norm",&fsimono);
  DoubleParameter("TOL_FSI_RES_INF",1e-6,"Absolute tolerance for interface residual in Inf-norm",&fsimono);
  DoubleParameter("TOL_FSI_INC_L2",1e-6,"Absolute tolerance for interface increment in L2-norm",&fsimono);
  DoubleParameter("TOL_FSI_INC_INF",1e-6,"Absolute tolerance for interface increment in Inf-norm",&fsimono);
  // fluid pressure
  DoubleParameter("TOL_PRE_RES_L2",1e-6,"Absolute tolerance for fluid pressure residual in L2-norm",&fsimono);
  DoubleParameter("TOL_PRE_RES_INF",1e-6,"Absolute tolerance for fluid pressure residual in Inf-norm",&fsimono);
  DoubleParameter("TOL_PRE_INC_L2",1e-6,"Absolute tolerance for fluid pressure increment in L2-norm",&fsimono);
  DoubleParameter("TOL_PRE_INC_INF",1e-6,"Absolute tolerance for fluid pressure increment in Inf-norm",&fsimono);
  // fluid velocities
  DoubleParameter("TOL_VEL_RES_L2",1e-6,"Absolute tolerance for fluid velocity residual in L2-norm",&fsimono);
  DoubleParameter("TOL_VEL_RES_INF",1e-6,"Absolute tolerance for fluid velocity residual in Inf-norm",&fsimono);
  DoubleParameter("TOL_VEL_INC_L2",1e-6,"Absolute tolerance for fluid velocity increment in L2-norm",&fsimono);
  DoubleParameter("TOL_VEL_INC_INF",1e-6,"Absolute tolerance for fluid velocity increment in Inf-norm",&fsimono);

  /*----------------------------------------------------------------------*/
  /* parameters for partitioned FSI solvers */
  Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER",false,"");

  DoubleParameter("BASETOL", 1e-3,
      "Basic tolerance for adaptive convergence check in monolithic FSI.\n"
      "This tolerance will be used for the linear solve of the FSI block system.\n"
      "The linear convergence test will always use the relative residual norm (AZ_r0).\n"
      "Not to be confused with the Newton tolerance (CONVTOL) that applies\n"
      "to the nonlinear convergence test using a absolute residual norm.",
      &fsipart);

  DoubleParameter("CONVTOL", 1e-6,
      "Tolerance for iteration over fields in case of partitioned scheme",
      &fsipart);


  setStringToIntegralParameter<int>("COUPMETHOD","conforming",
      "Coupling Method Mortar (mtr) or conforming nodes at interface",
      tuple<std::string>( // ToDO introduce enum
          "MTR",
          "Mtr",
          "mtr",
          "conforming",
          "immersed"
      ),
      tuple<int>(0,0,0,1,2),
      &fsipart);

  setStringToIntegralParameter<int>("COUPVARIABLE","Displacement",
      "Coupling variable at the interface",
      tuple<std::string>("Displacement","Force"),
      tuple<int>(0,1),
      &fsipart);

  BoolParameter("DIVPROJECTION", "no",
      "Project velocity into divergence-free subspace for partitioned fsi",
      &fsipart);

  IntParameter("ITEMAX", 100, "Maximum number of iterations over fields",
        &fsipart);

  DoubleParameter("MAXOMEGA", 0.0,
      "largest omega allowed for Aitken relaxation (0.0 means no constraint)",
      &fsipart);

  setStringToIntegralParameter<int>("PARTITIONED","DirichletNeumann",
      "Coupling strategies for partitioned FSI solvers.",
      tuple<std::string>(
       "DirichletNeumann",
       "DirichletNeumannSlideALE"
       ),
      tuple<int>(
       INPAR::FSI::DirichletNeumann,
       INPAR::FSI::DirichletNeumannSlideale
       ),
      &fsipart);

  setStringToIntegralParameter<int>("PREDICTOR","d(n)",
      "Predictor for interface displacements",
      tuple<std::string>(
        "d(n)",
        "d(n)+dt*(1.5*v(n)-0.5*v(n-1))",
        "d(n)+dt*v(n)",
        "d(n)+dt*v(n)+0.5*dt^2*a(n)"
        ),
      tuple<int>(1,2,3,4),
      &fsipart);

  DoubleParameter("RELAX", 1.0,
      "fixed relaxation parameter for partitioned FSI solvers", &fsipart);

  /* ----------------------------------------------------------------------- */
  Teuchos::ParameterList& constrfsi = fsidyn.sublist("CONSTRAINT",false,"");

  setStringToIntegralParameter<int> ("PRECONDITIONER","Simple","preconditioner to use",
      tuple<std::string>("Simple","Simplec"),tuple<int>(INPAR::FSI::Simple,INPAR::FSI::Simplec),&constrfsi);
  IntParameter("SIMPLEITER",2,"Number of iterations for simple pc",&constrfsi);
  DoubleParameter("ALPHA",0.8,"alpha parameter for simple pc",&constrfsi);

}

/*----------------------------------------------------------------------------*/
void INPAR::FSI::SetValidConditions(
    std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  std::vector<Teuchos::RCP<ConditionComponent> > fsicomponents;

  fsicomponents.push_back(Teuchos::rcp(new IntConditionComponent("coupling id")));

  Teuchos::RCP<ConditionDefinition> linefsi = Teuchos::rcp(
      new ConditionDefinition("DESIGN FSI COUPLING LINE CONDITIONS",
          "FSICoupling", "FSI Coupling", DRT::Condition::FSICoupling, true,
          DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffsi = Teuchos::rcp(
      new ConditionDefinition("DESIGN FSI COUPLING SURF CONDITIONS",
          "FSICoupling", "FSI Coupling", DRT::Condition::FSICoupling, true,
          DRT::Condition::Surface));

  for (unsigned i = 0; i < fsicomponents.size(); ++i)
  {
    linefsi->AddComponent(fsicomponents[i]);
    surffsi->AddComponent(fsicomponents[i]);
  }

  condlist.push_back(linefsi);
  condlist.push_back(surffsi);

  /*--------------------------------------------------------------------*/
  // FSI without sliding

  Teuchos::RCP<ConditionDefinition> linefsins = Teuchos::rcp(
      new ConditionDefinition("DESIGN FSI COUPLING NO SLIDE LINE CONDITIONS",
          "FSICouplingNoSlide", "FSI Coupling No Slide",
          DRT::Condition::FSICouplingNoSlide, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffsins = Teuchos::rcp(
      new ConditionDefinition("DESIGN FSI COUPLING NO SLIDE SURF CONDITIONS",
          "FSICouplingNoSlide", "FSI Coupling No Slide",
          DRT::Condition::FSICouplingNoSlide, true, DRT::Condition::Surface));

  condlist.push_back(linefsins);
  condlist.push_back(surffsins);

  /*--------------------------------------------------------------------*/
  // FSI define centerdisp for sliding interfaces

  Teuchos::RCP<ConditionDefinition> linefsicd = Teuchos::rcp(
      new ConditionDefinition("DESIGN FSI COUPLING CENTER DISP LINE CONDITIONS",
          "FSICouplingCenterDisp", "FSI Coupling Center Disp",
          DRT::Condition::FSICouplingCenterDisp, true, DRT::Condition::Line));
  Teuchos::RCP<ConditionDefinition> surffsicd = Teuchos::rcp(
      new ConditionDefinition("DESIGN FSI COUPLING CENTER DISP SURF CONDITIONS",
          "FSICouplingCenterDisp", "FSI Coupling Center Disp",
          DRT::Condition::FSICouplingCenterDisp, true,
          DRT::Condition::Surface));

  condlist.push_back(linefsicd);
  condlist.push_back(surffsicd);

  /*--------------------------------------------------------------------*/
  // Additional coupling of structure and ale fields (for lung fsi)

  std::vector<Teuchos::RCP<ConditionComponent> > saccomponents;

  saccomponents.push_back(
      Teuchos::rcp(new IntConditionComponent("coupling id")));
  saccomponents.push_back(
      Teuchos::rcp(
          new StringConditionComponent("field", "structure",
              Teuchos::tuple<std::string>("structure", "fluid"),
              Teuchos::tuple<std::string>("structure", "fluid"))));

  Teuchos::RCP<ConditionDefinition> surfsac = Teuchos::rcp(
      new ConditionDefinition("DESIGN STRUCTURE ALE COUPLING SURF CONDITIONS",
          "StructAleCoupling", "StructAleCoupling",
          DRT::Condition::StructAleCoupling, true, DRT::Condition::Surface));

  for (unsigned i = 0; i < saccomponents.size(); ++i)
    surfsac->AddComponent(saccomponents[i]);

  condlist.push_back(surfsac);

  /*--------------------------------------------------------------------*/
  // Additional coupling of structure and fluid volumes (for lung fsi)

  std::vector<Teuchos::RCP<ConditionComponent> > sfvcomponents;

  sfvcomponents.push_back(
      Teuchos::rcp(new IntConditionComponent("coupling id")));
  sfvcomponents.push_back(
      Teuchos::rcp(
          new StringConditionComponent("field", "structure",
              Teuchos::tuple<std::string>("structure", "fluid"),
              Teuchos::tuple<std::string>("structure", "fluid"))));

  Teuchos::RCP<ConditionDefinition> surfsfv = Teuchos::rcp(
      new ConditionDefinition(
          "DESIGN STRUCTURE FLUID VOLUME COUPLING SURF CONDITIONS",
          "StructFluidSurfCoupling", "StructFluidSurfCoupling",
          DRT::Condition::StructFluidSurfCoupling, true,
          DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> volsfv = Teuchos::rcp(
      new ConditionDefinition(
          "DESIGN STRUCTURE FLUID VOLUME COUPLING VOL CONDITIONS",
          "StructFluidVolCoupling", "StructFluidVolCoupling",
          DRT::Condition::StructFluidVolCoupling, false,
          DRT::Condition::Volume));

  for (unsigned i = 0; i < sfvcomponents.size(); ++i)
  {
    surfsfv->AddComponent(sfvcomponents[i]);
    volsfv->AddComponent(sfvcomponents[i]);
  }

  condlist.push_back(surfsfv);
  condlist.push_back(volsfv);

}

