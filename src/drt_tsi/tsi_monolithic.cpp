/*----------------------------------------------------------------------*/
/*!

\brief  Basis of all monolithic TSI algorithms that perform a coupling between
        the linear momentum equation and the heat conduction equation

\level 2

\maintainer Christoph Meier

*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | headers                                                   dano 11/10 |
 *----------------------------------------------------------------------*/
#include "tsi_monolithic.H"
#include "tsi_defines.H"
#include "tsi_utils.H"
#include "../drt_thermo/thrtimint.H"

#include "../drt_thermo/thermo_ele_action.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_locsys.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_solver.H"

// contact
#include "../drt_contact/contact_abstract_strategy.H"
#include "../drt_contact/contact_tsi_lagrange_strategy.H"
#include "../drt_contact/contact_nitsche_strategy_tsi.H"
#include "../drt_contact/contact_interface.H"
#include "../drt_contact/contact_tsi_interface.H"
#include "../drt_contact/contact_node.H"
#include "../drt_contact/meshtying_contact_bridge.H"
#include "../drt_mortar/mortar_manager_base.H"

// for coupling of nonmatching meshes
#include "../drt_adapter/adapter_coupling_volmortar.H"

#include "../drt_structure_new/str_timint_base.H"
#include "../drt_structure_new/str_model_evaluator_contact.H"
#include "../drt_structure_new/str_model_evaluator_structure.H"
#include "../drt_structure_new/str_model_evaluator_data.H"
#include "../drt_lib/drt_elements_paramsminimal.H"

#include "../drt_mortar/mortar_multifield_coupling.H"

//! Note: The order of calling the two BaseAlgorithm-constructors is
//! important here! In here control file entries are written. And these entries
//! define the order in which the filters handle the Discretizations, which in
//! turn defines the dof number ordering of the Discretizations.


/*----------------------------------------------------------------------*
 | destructor (public)                                       dano 11/10 |
 *----------------------------------------------------------------------*/
TSI::Monolithic::~Monolithic() {}


/*----------------------------------------------------------------------*
 | monolithic                                                dano 11/10 |
 *----------------------------------------------------------------------*/
TSI::Monolithic::Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& sdynparams)
    : Algorithm(comm),
      solveradapttol_(DRT::INPUT::IntegralValue<int>(
                          ((DRT::Problem::Instance()->TSIDynamicParams()).sublist("MONOLITHIC")),
                          "ADAPTCONV") == 1),
      solveradaptolbetter_(((DRT::Problem::Instance()->TSIDynamicParams()).sublist("MONOLITHIC"))
                               .get<double>("ADAPTCONV_BETTER")),
      printiter_(true),      // ADD INPUT PARAMETER
      printerrfile_(false),  // ADD INPUT PARAMETER FOR 'true'
      errfile_(NULL),
      zeros_(Teuchos::null),
      strmethodname_(DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP")),
      tsidyn_(DRT::Problem::Instance()->TSIDynamicParams()),
      tsidynmono_((DRT::Problem::Instance()->TSIDynamicParams()).sublist("MONOLITHIC")),
      blockrowdofmap_(Teuchos::null),
      systemmatrix_(Teuchos::null),
      k_st_(Teuchos::null),
      k_ts_(Teuchos::null),
      merge_tsi_blockmatrix_(
          DRT::INPUT::IntegralValue<bool>(tsidynmono_, "MERGE_TSI_BLOCK_MATRIX")),
      soltech_(DRT::INPUT::IntegralValue<INPAR::TSI::NlnSolTech>(tsidynmono_, "NLNSOL")),
      iternorm_(DRT::INPUT::IntegralValue<INPAR::TSI::VectorNorm>(tsidynmono_, "ITERNORM")),
      iter_(0),
      sdyn_(sdynparams),
      timernewton_(comm),
      dtsolve_(0.),
      ptcdt_(tsidynmono_.get<double>("PTCDT")),
      dti_(1.0 / ptcdt_),
      ls_strategy_(
          DRT::INPUT::IntegralValue<INPAR::TSI::LineSearch>(tsidynmono_, "TSI_LINE_SEARCH")),
      vel_(Teuchos::null)
{
  FixTimeIntegrationParams();

  // another setup of structural time integration with the correct initial temperature is required,
  // so get the temperature
  if (ThermoField()->Tempnp() == Teuchos::null) dserror("this is NULL");

  if (matchinggrid_)
    StructureField()->Discretization()->SetState(1, "temperature", ThermoField()->Tempnp());
  else
    StructureField()->Discretization()->SetState(
        1, "temperature", volcoupl_->ApplyVectorMapping12(ThermoField()->Tempnp()));

  // setup structural time integrator with initial temperature
  structure_->Setup();
  StructureField()->Discretization()->ClearState(true);

  errfile_ = DRT::Problem::Instance()->ErrorFile()->Handle();
  if (errfile_) printerrfile_ = true;

  blockrowdofmap_ = Teuchos::rcp(new LINALG::MultiMapExtractor);

  // initialise internal varible with new velocities V_{n+1} at t_{n+1}
  vel_ = LINALG::CreateVector(*(StructureField()->DofRowMap(0)), true);

  // --------------------------------- TSI solver: create a linear solver

  // get iterative solver
  if (merge_tsi_blockmatrix_ == false) CreateLinearSolver();
  // get direct solver, e.g. UMFPACK
  else  // (merge_tsi_blockmatrix_ == true)
  {
#ifndef TFSI
    if (Comm().MyPID() == 0) std::cout << "Merged TSI block matrix is used!\n" << std::endl;
#endif

    // get solver parameter list of linear TSI solver
    const int linsolvernumber = tsidynmono_.get<int>("LINEAR_SOLVER");
    const Teuchos::ParameterList& tsisolverparams =
        DRT::Problem::Instance()->SolverParams(linsolvernumber);

    Teuchos::RCP<Teuchos::ParameterList> solverparams = Teuchos::rcp(new Teuchos::ParameterList);
    *solverparams = tsisolverparams;

    solver_ = Teuchos::rcp(
        new LINALG::Solver(*solverparams, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  }  // end BlockMatrixMerge

  // StructureField: check whether we have locsys BCs, i.e. inclined structural
  //  Dirichlet BC
  {
    std::vector<DRT::Condition*> locsysconditions(0);
    (StructureField()->Discretization())->GetCondition("Locsys", locsysconditions);

    // if there are inclined structural Dirichlet BC, get the structural LocSysManager
    if (locsysconditions.size())
    {
      locsysman_ = StructureField()->LocsysManager();
    }
    else
      locsysman_ = Teuchos::null;
  }

#ifndef TFSI
  if ((DRT::INPUT::IntegralValue<bool>(tsidynmono_, "CALC_NECKING_TSI_VALUES") == true) and
      (Comm().MyPID() == 0))
    std::cout
        << "CAUTION: calculation ONLY valid for necking of a cylindrical body!"
        << "\n Due to symmetry only 1/8 of the cylinder is simulated, i.e r/l = 6.413mm/53.334mm."
        << "\n The body is located between x: [0,6.413mm], y: [0,6.413mm], z: "
           "[-13.3335mm,13.3335mm]\n"
        << std::endl;
#endif  // TFSI

  // structural and thermal contact
  GetContactStrategy();

}  // Monolithic()


/*----------------------------------------------------------------------*
 | read restart information for given time step (public)     dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ReadRestart(int step)
{
  ThermoField()->ReadRestart(step);
  StructureField()->ReadRestart(step);

  // StructureField()->ReadRestart destroyed the old object and created
  // a new one, so we update the pointers
  GetContactStrategy();

  // pass the current coupling variables to the respective field
  // second ReadRestart needed due to the coupling variables
  ApplyStructCouplingState(StructureField()->Dispnp(), StructureField()->Velnp());
  ThermoField()->ReadRestart(step);
  ThermoField()->Discretization()->ClearState(true);

  ApplyThermoCouplingState(ThermoField()->Tempnp());
  StructureField()->ReadRestart(step);
  StructureField()->Discretization()->ClearState(true);

  // StructureField()->ReadRestart destroyed the old object and created
  // a new one, so we update the pointers
  GetContactStrategy();

  SetTimeStep(ThermoField()->TimeOld(), step);

  // Material pointers to other field were deleted during ReadRestart().
  // They need to be reset.
  if (matchinggrid_)
    TSI::UTILS::SetMaterialPointersMatchingGrid(
        StructureField()->Discretization(), ThermoField()->Discretization());
  else
  {
    Teuchos::RCP<TSI::UTILS::TSIMaterialStrategy> strategy =
        Teuchos::rcp(new TSI::UTILS::TSIMaterialStrategy());
    volcoupl_->AssignMaterials(
        StructureField()->Discretization(), ThermoField()->Discretization(), strategy);
  }

  Teuchos::ParameterList p;
  //! pointer to the model evaluator data container
  Teuchos::RCP<DRT::ELEMENTS::ParamsMinimal> EvalData =
      Teuchos::rcp(new DRT::ELEMENTS::ParamsMinimal());
  EvalData->SetActionType(DRT::ELEMENTS::struct_calc_reset_istep);
  p.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", EvalData);
  p.set<std::string>("action", "calc_struct_reset_istep");
  StructureField()->Discretization()->Evaluate(p);


  return;
}  // ReadRestart()


/*----------------------------------------------------------------------*
 | prepare time step (public)                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrepareTimeStep()
{
  // we may have changed the ghosting when redistributing contact
  // so we make sure all maps in the system are up to date
  if (contact_strategy_nitsche_ != Teuchos::null) SetupSystem();

  // counter and print header
  // increment time and step counter
  IncrementTimeAndStep();
  PrintHeader();

  // pass the current coupling variables to the respective fields
  ApplyStructCouplingState(StructureField()->Dispnp(), StructureField()->Velnp());
  ApplyThermoCouplingState(ThermoField()->Tempnp());

  // call the predictor
  StructureField()->PrepareTimeStep();
  ThermoField()->PrepareTimeStep();

}  // PrepareTimeStep()


/*----------------------------------------------------------------------*
 | create linear solver                                   wiesner 07/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::CreateLinearSolver()
{
  // get the solver number used for linear TSI solver
  const int linsolvernumber = tsidynmono_.get<int>("LINEAR_SOLVER");
  // check if the TSI solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for monolithic TSI. Please set LINEAR_SOLVER in TSI DYNAMIC to a "
        "valid number!");

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (slinsolvernumber == (-1))
    dserror(
        "no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL "
        "DYNAMIC to a valid number!");

  // get parameter list of thermal dynamics
  const Teuchos::ParameterList& tdyn = DRT::Problem::Instance()->ThermalDynamicParams();
  // use solver blocks for temperature (thermal field)
  // get the solver number used for thermal solver
  const int tlinsolvernumber = tdyn.get<int>("LINEAR_SOLVER");
  // check if the TSI solver has a valid solver number
  if (tlinsolvernumber == (-1))
    dserror(
        "no linear solver defined for thermal field. Please set LINEAR_SOLVER in THERMAL DYNAMIC "
        "to a valid number!");

  // get solver parameter list of linear TSI solver
  const Teuchos::ParameterList& tsisolverparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);

  const int solvertype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(tsisolverparams, "SOLVER");

  if ((solvertype != INPAR::SOLVER::aztec_msr) and (solvertype != INPAR::SOLVER::belos))
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " Note: the BGS2x2 preconditioner now " << std::endl;
    std::cout << " uses the structural solver and thermal solver blocks" << std::endl;
    std::cout << " for building the internal inverses" << std::endl;
    std::cout << " Remove the old BGS PRECONDITIONER BLOCK entries " << std::endl;
    std::cout << " in the dat files!" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    dserror("aztec solver expected");
  }
  const int azprectype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(tsisolverparams, "AZPREC");

  // plausibility check
  switch (azprectype)
  {
    case INPAR::SOLVER::azprec_BGS2x2:
      break;
    case INPAR::SOLVER::azprec_BGSnxn:
    case INPAR::SOLVER::azprec_TekoSIMPLE:
    {
#ifdef HAVE_TEKO
      // check if structural solver and thermal solver are Stratimikos based (Teko expects
      // stratimikos)
      int solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(
          DRT::Problem::Instance()->SolverParams(slinsolvernumber), "SOLVER");
      if ((solvertype != INPAR::SOLVER::stratimikos_amesos) and
          (solvertype != INPAR::SOLVER::stratimikos_aztec) and
          (solvertype != INPAR::SOLVER::stratimikos_belos))
        dserror("Teko expects a STRATIMIKOS solver object in STRUCTURE SOLVER");

      solvertype = DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(
          DRT::Problem::Instance()->SolverParams(tlinsolvernumber), "SOLVER");
      if ((solvertype != INPAR::SOLVER::stratimikos_amesos) and
          (solvertype != INPAR::SOLVER::stratimikos_aztec) and
          (solvertype != INPAR::SOLVER::stratimikos_belos))
        dserror("Teko expects a STRATIMIKOS solver object in thermal solver %3d", tlinsolvernumber);
#else
      dserror(
          "Teko preconditioners only available with HAVE_TEKO flag for TRILINOS_DEV (>Q1/2011)");
#endif
      break;
    }
    case INPAR::SOLVER::azprec_MueLuAMG_sym:
    case INPAR::SOLVER::azprec_AMGnxn:
    case INPAR::SOLVER::azprec_CheapSIMPLE:
    {
      // no plausibility checks here
      // if you forget to declare an xml file you will get an error message anyway
    }
    break;
    default:
      dserror(
          "Block Gauss-Seidel BGS2x2 preconditioner expected. Alternatively you can define your "
          "own AMG block preconditioner (using an xml file). This is experimental.");
      break;
  }


  // prepare linear solvers and preconditioners
  switch (azprectype)
  {
    case INPAR::SOLVER::azprec_BGS2x2:
    case INPAR::SOLVER::azprec_BGSnxn:
    case INPAR::SOLVER::azprec_TekoSIMPLE:
    case INPAR::SOLVER::azprec_AMGnxn:
    case INPAR::SOLVER::azprec_CheapSIMPLE:
    {
      // This should be the default case (well-tested and used)
      solver_ = Teuchos::rcp(new LINALG::Solver(tsisolverparams,
          // ggfs. explizit Comm von STR wie lungscatra
          Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

      // use solver blocks for structure and temperature (thermal field)
      const Teuchos::ParameterList& ssolverparams =
          DRT::Problem::Instance()->SolverParams(slinsolvernumber);
      const Teuchos::ParameterList& tsolverparams =
          DRT::Problem::Instance()->SolverParams(tlinsolvernumber);

      solver_->PutSolverParamsToSubParams("Inverse1", ssolverparams);
      solver_->PutSolverParamsToSubParams("Inverse2", tsolverparams);

      // prescribe rigid body modes
      StructureField()->Discretization()->ComputeNullSpaceIfNecessary(
          solver_->Params().sublist("Inverse1"));
      ThermoField()->Discretization()->ComputeNullSpaceIfNecessary(
          solver_->Params().sublist("Inverse2"));


      if (azprectype == INPAR::SOLVER::azprec_CheapSIMPLE)
      {
        // Tell to the LINALG::SOLVER::SimplePreconditioner that we use the general implementation
        solver_->Params().set<bool>("GENERAL", true);
      }

      break;
    }
    case INPAR::SOLVER::azprec_MueLuAMG_sym:
    {
      solver_ = Teuchos::rcp(new LINALG::Solver(tsisolverparams,
          // ggfs. explizit Comm von STR wie lungscatra
          Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

      // use solver blocks for structure and temperature (thermal field)
      const Teuchos::ParameterList& ssolverparams =
          DRT::Problem::Instance()->SolverParams(slinsolvernumber);
      const Teuchos::ParameterList& tsolverparams =
          DRT::Problem::Instance()->SolverParams(tlinsolvernumber);

      // This is not very elegant:
      // first read in solver parameters. These have to contain ML parameters such that...
      solver_->PutSolverParamsToSubParams("Inverse1", ssolverparams);
      solver_->PutSolverParamsToSubParams("Inverse2", tsolverparams);

      // ... BACI calculates the null space vectors. These are then stored in the sublists
      //     Inverse1 and Inverse2 from where they...
      StructureField()->Discretization()->ComputeNullSpaceIfNecessary(
          solver_->Params().sublist("Inverse1"));
      ThermoField()->Discretization()->ComputeNullSpaceIfNecessary(
          solver_->Params().sublist("Inverse2"));

      // ... are copied from here to ...
      const Teuchos::ParameterList& inv1source =
          solver_->Params().sublist("Inverse1").sublist("ML Parameters");
      const Teuchos::ParameterList& inv2source =
          solver_->Params().sublist("Inverse2").sublist("ML Parameters");

      // ... here. The "MueLu Parameters" sublists "Inverse1" and "Inverse2" only contain the basic
      //     information about the corresponding null space vectors, which are actually copied ...
      Teuchos::ParameterList& inv1 =
          solver_->Params().sublist("MueLu Parameters").sublist("Inverse1");
      Teuchos::ParameterList& inv2 =
          solver_->Params().sublist("MueLu Parameters").sublist("Inverse2");

      // ... here.
      inv1.set<int>("PDE equations", inv1source.get<int>("PDE equations"));
      inv2.set<int>("PDE equations", inv2source.get<int>("PDE equations"));
      inv1.set<int>("null space: dimension", inv1source.get<int>("null space: dimension"));
      inv2.set<int>("null space: dimension", inv2source.get<int>("null space: dimension"));
      inv1.set<double*>("null space: vectors", inv1source.get<double*>("null space: vectors"));
      inv2.set<double*>("null space: vectors", inv2source.get<double*>("null space: vectors"));
      inv1.set<Teuchos::RCP<std::vector<double>>>(
          "nullspace", inv1source.get<Teuchos::RCP<std::vector<double>>>("nullspace"));
      inv2.set<Teuchos::RCP<std::vector<double>>>(
          "nullspace", inv2source.get<Teuchos::RCP<std::vector<double>>>("nullspace"));

      solver_->Params().sublist("MueLu Parameters").set("TSI", true);
      break;
    }
    default:
      dserror("Block Gauss-Seidel BGS2x2 preconditioner expected");
      break;
  }

}  // CreateLinearSolver()


/*----------------------------------------------------------------------*
 | non-linear solve, i.e. (multiple) corrector (public)      dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::Solve()
{
  // choose solution technique according to input file
  switch (soltech_)
  {
    // Newton-Raphson iteration
    case INPAR::TSI::soltech_newtonfull:
      NewtonFull();
      break;
    // Pseudo-transient continuation
    case INPAR::TSI::soltech_ptc:
      PTC();
      break;
    // catch problems
    default:
      dserror("Solution technique \"%s\" is not implemented",
          INPAR::TSI::NlnSolTechString(soltech_).c_str());
      break;
  }  // end switch (soltechnique_)

  return;
}  // Solve()


/*----------------------------------------------------------------------*
 | time loop of the monolithic system                        dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::TimeLoop()
{
  // time loop
  while (NotFinished())
  {
    // counter and print header
    // predict solution of both field (call the adapter)
    PrepareTimeStep();

    // integrate time step
    Solve();

    // calculate stresses, strains, energies
    PrepareOutput();

    // update all single field solvers
    Update();

    // write output to screen and files
    Output();

#ifdef TSIMONOLITHASOUTPUT
    printf("Ende Timeloop ThermoField()->Tempnp[0] %12.8f\n", (*ThermoField()->Tempnp())[0]);
    printf("Ende Timeloop ThermoField()->Tempn[0] %12.8f\n", (*ThermoField()->Tempn())[0]);

    printf("Ende Timeloop disp %12.8f\n", (*StructureField()->Dispn())[0]);
    std::cout << "dispn\n" << *(StructureField()->Dispn()) << std::endl;
#endif  // TSIMONOLITHASOUTPUT

  }  // NotFinished
}  // TimeLoop()


/*----------------------------------------------------------------------*
 | solution with full Newton-Raphson iteration               dano 10/10 |
 | in tsi_algorithm: NewtonFull()                                       |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::NewtonFull()
{
#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0) std::cout << "TSI::Monolithic::NewtonFull()" << std::endl;
#endif  // TFSI
#endif  // TSI_DEBUG

  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic TSI tangent matrix

  // initialise equilibrium loop
  iter_ = 0;

  // incremental solution vector with length of all TSI dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);
  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);

  // compute residual forces #rhs_ and tangent #systemmatrix_
  // whose components are globally oriented
  // build linear system stiffness matrix and rhs/force residual for each
  // field, here e.g. for structure field: field want the iteration increment
  // 1.) Update(iterinc_),
  // 2.) EvaluateForceStiffResidual(),
  // 3.) PrepareSystemForNewtonSolve() --> if (locsysman_!=null) k_ss is rotated
  Evaluate(iterinc_);

  // create the linear system
  // \f$J(x_i) \Delta x_i = - R(x_i)\f$
  // create the systemmatrix
  SetupSystemMatrix();

  // check whether we have a sanely filled tangent matrix
  if (not systemmatrix_->Filled()) dserror("Effective tangent matrix must be filled here");

  // create full monolithic rhs vector
  // make negative residual not necessary: rhs_ is already negative
  // (STR/THR)-RHS is put negative in PrepareSystemForNewtonSolve()
  SetupRHS();

  // do the thermo contact modifications all at once
  if (contact_strategy_lagrange_ != Teuchos::null)
    contact_strategy_lagrange_->Evaluate(
        SystemMatrix(), rhs_, coupST_, StructureField()->Dispnp(), ThermoField()->Tempnp());
  ApplyDBC();

  // initialize with predictor values
  normrhsiter0_ = normrhs_ = CalculateVectorNorm(iternorm_, rhs_);
  normstrrhsiter0_ = normstrrhs_ = last_iter_res_.first =
      CalculateVectorNorm(iternormstr_, StructureField()->RHS());
  normthrrhsiter0_ = normthrrhs_ = last_iter_res_.second =
      CalculateVectorNorm(iternormthr_, ThermoField()->RHS());
  ls_step_length_ = 1.;
  normdisi_ = normtempi_ = norminc_ = 0.;
  normstrrhsiter0_ = normstrrhs_;
  norminciter0_ = normdisiiter0_ = normtempiiter0_ = -1.;

  // print stuff
  PrintNewtonIter();

  //------------------------------------------------------ iteration loop

  // equilibrium iteration loop (loop over k)
  while (((not Converged()) and (iter_ <= itermax_)) or (iter_ <= itermin_))
  {
    // increment equilibrium loop index
    ++iter_;
    ls_step_length_ = 1.;

    // reset timer
    timernewton_.ResetStartTime();

    // *********** time measurement ***********
    double dtcpu = timernewton_.WallTime();
    // *********** time measurement ***********
    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
    LinearSolve();
    // *********** time measurement ***********
    dtsolve_ = timernewton_.WallTime() - dtcpu;
    // *********** time measurement ***********

    // recover LM in the case of contact
    if (contact_strategy_lagrange_ != Teuchos::null) RecoverStructThermLM();

    // reset solver tolerance
    solver_->ResetTolerance();

    // compute residual forces #rhs_ and tangent #systemmatrix_
    // whose components are globally oriented
    // build linear system stiffness matrix and rhs/force residual for each
    // field, here e.g. for structure field: field want the iteration increment
    // 1.) Update(iterinc_),
    // 2.) EvaluateForceStiffResidual(),
    // 3.) PrepareSystemForNewtonSolve() --> if (locsysman_!=null) k_ss is rotated
    Evaluate(iterinc_);

    // create the linear system
    // \f$J(x_i) \Delta x_i = - R(x_i)\f$
    // create the systemmatrix
    SetupSystemMatrix();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled()) dserror("Effective tangent matrix must be filled here");

    // create full monolithic rhs vector
    // make negative residual not necessary: rhs_ is already negative
    // (STR/THR)-RHS is put negative in PrepareSystemForNewtonSolve()
    SetupRHS();

    // do the thermo contact modifications all at once
    if (contact_strategy_lagrange_ != Teuchos::null)
    {
      // *********** time measurement ***********
      double dtcpu = timernewton_.WallTime();
      // *********** time measurement ***********

      contact_strategy_lagrange_->Evaluate(
          SystemMatrix(), rhs_, coupST_, StructureField()->Dispnp(), ThermoField()->Tempnp());

      // *********** time measurement ***********
      dtcmt_ = timernewton_.WallTime() - dtcpu;
      // *********** time measurement ***********
    }
    ApplyDBC();

    // do line search
    switch (ls_strategy_)
    {
      case INPAR::TSI::LS_none:
        break;
      case INPAR::TSI::LS_structure:
      case INPAR::TSI::LS_thermo:
      case INPAR::TSI::LS_and:
      case INPAR::TSI::LS_or:
      {
        normstrrhs_ = CalculateVectorNorm(iternormstr_, StructureField()->RHS());
        normthrrhs_ = CalculateVectorNorm(iternormthr_, ThermoField()->RHS());
        iterinc_->Scale(-1.);

        while (ls_step_length_ > 1.e-8 && !LSadmissible())
        {
          iterinc_->Scale(.5);
          ls_step_length_ *= .5;
          Evaluate(iterinc_);
          normstrrhs_ = CalculateVectorNorm(iternormstr_, StructureField()->RHS());
          normthrrhs_ = CalculateVectorNorm(iternormthr_, ThermoField()->RHS());
        }

        last_iter_res_.first = CalculateVectorNorm(iternormstr_, StructureField()->RHS());
        last_iter_res_.second = CalculateVectorNorm(iternormthr_, ThermoField()->RHS());

        if (ls_step_length_ < 1.)
        {
          SetupSystemMatrix();
          SetupRHS();
          ApplyDBC();
        }
      }
      break;
      default:
        dserror("unknown line search strategy");
    }  // end line search

    // --------------------------------------------- build residual norms
    // include all stuff here related with convergence test
    normrhs_ = CalculateVectorNorm(iternorm_, rhs_);
    // vector of displacement and temperature residual
    Teuchos::RCP<Epetra_Vector> strrhs;
    Teuchos::RCP<Epetra_Vector> thrrhs;
    // extract field vectors
    ExtractFieldVectors(rhs_, strrhs, thrrhs);
    normstrrhs_ = CalculateVectorNorm(iternormstr_, strrhs);
    normthrrhs_ = CalculateVectorNorm(iternormthr_, thrrhs);

    // --------------------------------- build residual incremental norms
    // vector of displacement and temperature increments
    Teuchos::RCP<Epetra_Vector> sx;
    Teuchos::RCP<Epetra_Vector> tx;
    // extract field vectors
    ExtractFieldVectors(iterinc_, sx, tx);
    norminc_ = CalculateVectorNorm(iternorm_, iterinc_);
    normdisi_ = CalculateVectorNorm(iternormstr_, sx);
    normtempi_ = CalculateVectorNorm(iternormthr_, tx);

    // in case of 'Mix'-convergence criterion: save the norm of the 1st
    // iteration in (norm . iter0_)
    if (iter_ == 1)
    {
      // save initial incremental norms
      norminciter0_ = norminc_;
      normdisiiter0_ = normdisi_;
      normtempiiter0_ = normtempi_;

      // set the minimum of iter0_ and tolrhs_, because we want to prevent the
      // case of a zero characteristic initial norm
      if (normrhsiter0_ == 0.0) normrhsiter0_ = tolrhs_;
      if (normstrrhsiter0_ == 0.0) normstrrhsiter0_ = tolstrrhs_;
      if (normthrrhsiter0_ == 0.0) normthrrhsiter0_ = tolthrrhs_;
      if (norminciter0_ == 0.0) norminciter0_ = tolinc_;
      if (normdisiiter0_ == 0.0) normdisiiter0_ = toldisi_;
      if (normtempiiter0_ == 0.0) normtempiiter0_ = toltempi_;
    }

    // print stuff
    PrintNewtonIter();

  }  // end equilibrium loop

  // ----------------------------------------------------- iteration loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ((Converged()) and (Comm().MyPID() == 0))
  {
    PrintNewtonConv();
  }
  else if (iter_ >= itermax_)
  {
    if (DRT::INPUT::IntegralValue<INPAR::STR::DivContAct>(sdyn_, "DIVERCONT") ==
        INPAR::STR::divcont_continue)
      ;  // do nothing
    else
      dserror("Newton unconverged in %d iterations", iter_);
  }
  // for validation with literature calculate nodal TSI values
  if ((DRT::INPUT::IntegralValue<bool>(tsidynmono_, "CALC_NECKING_TSI_VALUES")) == true)
    CalculateNeckingTSIResults();

}  // NewtonFull()


/*----------------------------------------------------------------------*
 | solution with pseudo-transient continuation               dano 06/14 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PTC()
{
  // do a PTC iteration here
  // implementation is based on the work of Gee, Kelley, Lehouq (2009):
  // "Pseudo-transient continuation for nonlinear transient elasticity"

#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0) std::cout << "TSI::Monolithic::PTC()" << std::endl;
#endif  // TFSI
#endif  // TSI_DEBUG

  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic TSI tangent matrix

  // initialise equilibrium loop
  iter_ = 0;

  // incremental solution vector with length of all TSI dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);
  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);

  // compute residual forces #rhs_ and tangent #systemmatrix_
  // whose components are globally oriented
  // build linear system stiffness matrix and rhs/force residual for each
  // field, here e.g. for structure field: field want the iteration increment
  // 1.) Update(iterinc_),
  // 2.) EvaluateForceStiffResidual(),
  // 3.) PrepareSystemForNewtonSolve() --> if (locsysman_!=null) k_ss is rotated
  Evaluate(iterinc_);

  // create the linear system
  // \f$J(x_i) \Delta x_i = - R(x_i)\f$
  // create the systemmatrix
  SetupSystemMatrix();

  // check whether we have a sanely filled tangent matrix
  if (not systemmatrix_->Filled()) dserror("Effective tangent matrix must be filled here");

  // create full monolithic rhs vector
  // make negative residual not necessary: rhs_ is already negative
  // (STR/THR)-RHS is put negative in PrepareSystemForNewtonSolve()
  SetupRHS();

  ApplyDBC();

  // initialize with predictor values
  normrhsiter0_ = normrhs_ = CalculateVectorNorm(iternorm_, rhs_);
  normrhsiter0_ = normstrrhs_ = last_iter_res_.first =
      CalculateVectorNorm(iternormstr_, StructureField()->RHS());
  normthrrhsiter0_ = normthrrhs_ = last_iter_res_.second =
      CalculateVectorNorm(iternormthr_, ThermoField()->RHS());
  ls_step_length_ = 1.;
  normdisi_ = normtempi_ = norminc_ = 0.;
  normstrrhsiter0_ = normstrrhs_;
  normrhsiter0_ = normstrrhsiter0_ = normthrrhsiter0_ = norminciter0_ = normdisiiter0_ =
      normtempiiter0_ = -1.;

  // ----------------------------------------------- special stuff of PTC

  // compute the PTC parameters
  double ptcdt = ptcdt_;
  // norm of residual of old iteration step
  double nc = 0.0;
  // as recommended by Michael, apply PTC to the whole TSI system, i.e. use TSI
  // residual here
  // here the new convergence test stuff has to be included

  normrhs_ = CalculateVectorNorm(iternorm_, rhs_);
  rhs_->NormInf(&nc);
  // define the pseudo time step delta^{-1}
  double dti = 1 / ptcdt;

  // print stuff
  PrintNewtonIter();

  //------------------------------------------------------ iteration loop

  // equilibrium iteration loop (loop over k)
  while (((not Converged()) and (iter_ <= itermax_)) or (iter_ <= itermin_))
  {
    // increment equilibrium loop index
    ++iter_;

    // reset timer
    timernewton_.ResetStartTime();

    // ---------- modify diagonal blocks of systemmatrix according to PTC

    // modify structural diagonal block k_ss
    {
      Teuchos::RCP<Epetra_Vector> tmp_SS =
          LINALG::CreateVector(StructureField()->SystemMatrix()->RowMap(), false);
      tmp_SS->PutScalar(dti);
      Teuchos::RCP<Epetra_Vector> diag_SS =
          LINALG::CreateVector(StructureField()->SystemMatrix()->RowMap(), false);
      StructureField()->SystemMatrix()->ExtractDiagonalCopy(*diag_SS);

      diag_SS->Update(1.0, *tmp_SS, 1.0);

      StructureField()->SystemMatrix()->ReplaceDiagonalValues(*diag_SS);
    }
    // modify thermal diagonal block k_tt
    {
      Teuchos::RCP<Epetra_Vector> tmp_tt =
          LINALG::CreateVector(ThermoField()->SystemMatrix()->RowMap(), false);
      tmp_tt->PutScalar(dti);
      Teuchos::RCP<Epetra_Vector> diag_tt =
          LINALG::CreateVector(ThermoField()->SystemMatrix()->RowMap(), false);
      ThermoField()->SystemMatrix()->ExtractDiagonalCopy(*diag_tt);
      diag_tt->Update(1.0, *tmp_tt, 1.0);
      ThermoField()->SystemMatrix()->ReplaceDiagonalValues(*diag_tt);
    }

    // *********** time measurement ***********
    double dtcpu = timernewton_.WallTime();
    // *********** time measurement ***********
    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
    LinearSolve();
    // *********** time measurement ***********
    dtsolve_ = timernewton_.WallTime() - dtcpu;
    // *********** time measurement ***********

    // reset solver tolerance
    solver_->ResetTolerance();

    // compute residual forces #rhs_ and tangent #systemmatrix_
    // whose components are globally oriented
    // build linear TSI tangent matrix and rhs/force residual for each field,
    // here e.g. for structure field: STR field wants the iteration increment
    // 1.) Update(iterinc_),
    // 2.) EvaluateForceStiffResidual(),
    // 3.) PrepareSystemForNewtonSolve() --> if (locsysman_ != null) k_ss is rotated
    Evaluate(iterinc_);

    // create the linear system including PTC-modified systemmatrices k_ss and k_tt
    // \f$J(x_i) \Delta x_i = - R(x_i)\f$
    // create the systemmatrix
    SetupSystemMatrix();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled()) dserror("Effective tangent matrix must be filled here");

    // create full monolithic rhs vector
    // make negative residual not necessary: rhs_ is already negative
    // (STR/THR)-RHS is put negative in PrepareSystemForNewtonSolve()
    SetupRHS();

    // apply Dirichlet boundary conditions on System matrix and RHS
    ApplyDBC();

    // --------------------------------------------- build residual norms
    // include all stuff here related with convergence test
    normrhs_ = CalculateVectorNorm(iternorm_, rhs_);
    // vector of displacement and temperature residual
    Teuchos::RCP<Epetra_Vector> strrhs;
    Teuchos::RCP<Epetra_Vector> thrrhs;
    // extract field vectors
    ExtractFieldVectors(rhs_, strrhs, thrrhs);
    normstrrhs_ = CalculateVectorNorm(iternormstr_, strrhs);
    normthrrhs_ = CalculateVectorNorm(iternormthr_, thrrhs);

    // --------------------------------- build residual incremental norms
    // vector of displacement and temperature increments
    Teuchos::RCP<Epetra_Vector> sx;
    Teuchos::RCP<Epetra_Vector> tx;
    // extract field vectors
    ExtractFieldVectors(iterinc_, sx, tx);
    norminc_ = CalculateVectorNorm(iternorm_, iterinc_);
    normdisi_ = CalculateVectorNorm(iternormstr_, sx);
    normtempi_ = CalculateVectorNorm(iternormthr_, tx);

    // in case of 'Mix'-convergence criterion: save the norm of the 1st
    // iteration in norm*iter0_
    if (iter_ == 1)
    {
      // set residuals
      normrhsiter0_ = normrhs_;
      normstrrhsiter0_ = normstrrhs_;
      normthrrhsiter0_ = normthrrhs_;
      // set increments
      norminciter0_ = norminc_;
      normdisiiter0_ = normdisi_;
      normtempiiter0_ = normtempi_;

      // we set the minimum of iter0_ and tolrhs_, because
      // we want to prevent the case of a zero characteristic initial norm
      if (normrhsiter0_ == -1.) normrhsiter0_ = tolrhs_;
      if (normstrrhsiter0_ == -1.) normstrrhsiter0_ = tolstrrhs_;
      if (normthrrhsiter0_ == -1.) normthrrhsiter0_ = tolthrrhs_;
      if (norminciter0_ == -1.) norminciter0_ = tolinc_;
      if (normdisiiter0_ == -1.) normdisiiter0_ = toldisi_;
      if (normtempiiter0_ == -1.) normtempiiter0_ = toltempi_;
    }

    // print stuff
    PrintNewtonIter();

    // save old pseudo-time step in dti_
    dti_ = dti;

    // update ptc
    {
      double np = 0.0;
      rhs_->NormInf(&np);
      dti *= (np / nc);
      dti = std::max(dti, 0.0);
      nc = np;
    }

  }  // end equilibrium loop

  // ----------------------------------------------------- iteration loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ((Converged()) and (Comm().MyPID() == 0))
  {
    PrintNewtonConv();
  }
  else if (iter_ >= itermax_)
    dserror("PTC unconverged in %d iterations", iter_);

}  // PTC()


/*----------------------------------------------------------------------*
 | evaluate the single fields                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::Evaluate(Teuchos::RCP<Epetra_Vector> x)
{
#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0) std::cout << "\n TSI::Monolithic::Evaluate()" << std::endl;
#endif  // TFSI
#endif  // TSI_DEBUG

  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::Evaluate");

  // displacement and temperature incremental vector
  Teuchos::RCP<Epetra_Vector> sx;
  Teuchos::RCP<Epetra_Vector> tx;

  // if an increment vector exists
  if (x != Teuchos::null)
  {
    // extract displacement sx and temperature tx incremental vector of global
    // unknown incremental vector x
    ExtractFieldVectors(x, sx, tx);

#ifdef TSIMONOLITHASOUTPUT
    std::cout << "Recent thermal increment DT_n+1^i\n" << *(tx) << std::endl;
    std::cout << "Recent structural increment Dd_n+1^i\n" << *(sx) << std::endl;

    std::cout << "Until here only old solution of Newton step. No update applied\n"
              << *(ThermoField()->Tempnp()) << std::endl;
#endif  // TSIMONOLITHASOUTPUT
  }

  // else (x == Teuchos::null): initialise the system
#ifdef TSIMONOLITHASOUTPUT
  std::cout << "Tempnp vor UpdateNewton\n" << *(ThermoField()->Tempnp()) << std::endl;
  printf(
      "Tempnp vor UpdateNewton ThermoField()->Tempnp[0] %12.8f\n", (*ThermoField()->Tempnp())[0]);
#endif  // TSIMONOLITHASOUTPUT

  // Newton update of the thermo field
  // update temperature before passed to the structural field
  //   UpdateIterIncrementally(tx),
  ThermoField()->UpdateNewton(tx);

#ifdef TSIMONOLITHASOUTPUT
  std::cout << "Tempnp nach UpdateNewton\n" << *(ThermoField()->Tempnp()) << std::endl;
  printf(
      "Tempnp nach UpdateNewton ThermoField()->Tempnp[0] %12.8f\n", (*ThermoField()->Tempnp())[0]);
#endif  // TSIMONOLITHASOUTPUT

  // call all elements and assemble rhs and matrices
  /// structural field

  // structure Evaluate (builds tangent, residual and applies DBC)
  Epetra_Time timerstructure(Comm());

#ifndef MonTSIwithoutTHR
  // apply current temperature to structure
  ApplyThermoCouplingState(ThermoField()->Tempnp(), tx);
#endif

#ifdef TSIPARALLEL
  std::cout << Comm().MyPID() << " nach ApplyTemp!!" << std::endl;
#endif  // TSIPARALLEL

#ifdef TSIMONOLITHASOUTPUT
  Teuchos::RCP<Epetra_Vector> tempera =
      Teuchos::rcp(new Epetra_Vector(ThermoField()->Tempn()->Map(), true));

  if (ThermoField()->Tempnp() != Teuchos::null) tempera->Update(1.0, *ThermoField()->Tempnp(), 0.0);

  StructureField()->Discretization()->SetState(1, "temperature", tempera);
  StructureField()->Discretization()->SetState(1, "temperature", ThermoField()->Tempn());
#endif  // TSIMONOLITHASOUTPUT

  // Monolithic TSI accesses the linearised structure problem:
  //   UpdaterIterIncrementally(sx),
  //   EvaluateForceStiffResidual()
  //   PrepareSystemForNewtonSolve()
  //     blank residual DOFs that are on Dirichlet BC
  //     in case of local coordinate systems rotate the residual forth and back
  //     Be AWARE: ApplyDirichlettoSystem has to be called with rotated stiff_!
  if (iter_ == 0)
    StructureField()->Evaluate();
  else
    StructureField()->Evaluate(sx);
  StructureField()->Discretization()->ClearState(true);

#ifdef TSI_DEBUG
#ifndef TFSI
  std::cout << "  structure time for calling Evaluate: " << timerstructure.ElapsedTime() << "\n";
#endif  // TFSI
#endif  // TSI_DEBUG

#ifdef TSIMONOLITHASOUTPUT
  std::cout << "STR rhs_" << *StructureField()->RHS() << std::endl;
#endif  // TSIMONOLITHASOUTPUT

  /// thermal field

  // thermo Evaluate
  // (builds tangent, residual and applies DBC and recent coupling values)
  Epetra_Time timerthermo(Comm());

  // apply current displacements and velocities to the thermo field
  if (strmethodname_ == INPAR::STR::dyna_statics)
  {
    // calculate velocity V_n+1^k = (D_n+1^k-D_n)/Dt()
    vel_ = CalcVelocity(StructureField()->Dispnp());
  }
  // else: use velnp
  else
    vel_ = StructureField()->Velnp();

#ifndef MonTSIwithoutSTR
  // pass the structural values to the thermo field
  ApplyStructCouplingState(StructureField()->Dispnp(), vel_);
#endif

#ifdef TSIMONOLITHASOUTPUT
  std::cout << "d_n+1 inserted in THR field\n" << *(StructureField()->Dispnp()) << std::endl;
  std::cout << "v_n+1\n" << *vel_ << std::endl;
#endif  // TSIMONOLITHASOUTPUT

  // monolithic TSI accesses the linearised thermo problem
  //   EvaluateRhsTangResidual() and
  //   PrepareSystemForNewtonSolve()
  ThermoField()->Evaluate();
  ThermoField()->Discretization()->ClearState(true);
#ifdef TSI_DEBUG
#ifndef TFSI
  std::cout << "  thermo time for calling Evaluate: " << timerthermo.ElapsedTime() << "\n";
#endif  // TFSI
#endif  // TSI_DEBUG

}  // Evaluate()


/*----------------------------------------------------------------------*
 | extract field vectors for calling Evaluate() of the       dano 11/10 |
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ExtractFieldVectors(
    Teuchos::RCP<Epetra_Vector> x, Teuchos::RCP<Epetra_Vector>& sx, Teuchos::RCP<Epetra_Vector>& tx)
{
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::ExtractFieldVectors");

  // process structure unknowns of the first field
  sx = Extractor()->ExtractVector(x, 0);

  // process thermo unknowns of the second field
  tx = Extractor()->ExtractVector(x, 1);
}  // ExtractFieldVectors()


/*----------------------------------------------------------------------*
 | full monolithic dof row map                               dano 05/12 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> TSI::Monolithic::DofRowMap() const
{
  return Extractor()->FullMap();
}  // DofRowMap()


/*----------------------------------------------------------------------*
 | setup system (called in tsi_dyn)                          dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupSystem()
{
#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0) std::cout << " TSI::Monolithic::SetupSystem()" << std::endl;
#endif  // TFSI
#endif  // TSI_DEBUG

  // set parameters that remain the same in the whole calculation
  SetDefaultParameters();

#ifdef TSIPARALLEL
  std::cout << Comm().MyPID() << " :PID" << std::endl;
  std::cout << "structure dofmap" << std::endl;
  std::cout << *StructureField()->DofRowMap(0) << std::endl;
  std::cout << "thermo dofmap" << std::endl;
  std::cout << *StructureField()->DofRowMap(1) << std::endl;
#endif  // TSIPARALLEL

  SetDofRowMaps();

  /*----------------------------------------------------------------------*/
  // initialise TSI-systemmatrix_
  systemmatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *Extractor(), *Extractor(), 81, false, true));

  // create empty matrix
  k_st_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *(StructureField()->Discretization()->DofRowMap(0)), 81, false, true));

  // create empty matrix
  k_ts_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*(ThermoField()->Discretization()->DofRowMap(0)), 81, false, true));

}  // SetupSystem()


/*----------------------------------------------------------------------*
 | put the single maps to one full TSI map together          dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetDofRowMaps()
{
  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

  // use its own DofRowMap, that is the 0th map of the discretization
  vecSpaces.push_back(StructureField()->DofRowMap(0));
  vecSpaces.push_back(ThermoField()->DofRowMap(0));

  if (vecSpaces[0]->NumGlobalElements() == 0) dserror("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements() == 0) dserror("No temperature equation. Panic.");

  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(vecSpaces);

  // full TSI-blockmap
  Extractor()->Setup(*fullmap, vecSpaces);
}  // SetDofRowMaps()


/*----------------------------------------------------------------------*
 | setup system matrix of TSI                                dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupSystemMatrix()
{
#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0) std::cout << " TSI::Monolithic::SetupSystemMatrix()" << std::endl;
#endif  // TFSI
#endif  // TSI_DEBUG
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::SetupSystemMatrix");

  /*----------------------------------------------------------------------*/
  // pure structural part k_ss (3nx3n)

  // build pure structural block k_ss
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<LINALG::SparseMatrix> k_ss = StructureField()->SystemMatrix();

  // assign structure part to the TSI matrix
  systemmatrix_->Assign(0, 0, LINALG::View, *k_ss);

  /*----------------------------------------------------------------------*/
  // structural block k_st (3nxn)
  // build mechanical-thermal block

  k_st_->Reset();
  // call the element and calculate the matrix block
#ifndef MonTSIwithoutTHR
  ApplyStrCouplMatrix(k_st_);
#endif  // MonTSIwithoutTHR

  k_st_->Complete(*(StructureField()->Discretization()->DofRowMap(1)),
      *(StructureField()->Discretization()->DofRowMap(0)));

  if (!matchinggrid_) k_st_ = volcoupl_->ApplyMatrixMapping12(k_st_);

  k_st_->UnComplete();

  // assign thermo part to the TSI matrix
  systemmatrix_->Assign(0, 1, LINALG::View, *(k_st_));

  /*----------------------------------------------------------------------*/
  // pure thermo part k_tt (nxn)

  // build pure thermal block k_tt
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix systemmatrix_
  Teuchos::RCP<LINALG::SparseMatrix> k_tt = ThermoField()->SystemMatrix();

  // assign thermo part to the TSI matrix
  systemmatrix_->Assign(1, 1, LINALG::View, *(k_tt));

  /*----------------------------------------------------------------------*/
  // thermo part k_ts (nx3n)
  // build thermal-mechanical block

  k_ts_->Reset();

  // call the element and calculate the matrix block
#if ((!defined(MonTSIwithoutSTR)) and (!defined(COUPLEINITTEMPERATURE)))
  ApplyThrCouplMatrix(k_ts_);
  ApplyThrCouplMatrix_ConvBC(k_ts_);
#endif

  k_ts_->Complete(*(ThermoField()->Discretization()->DofRowMap(1)),
      *(ThermoField()->Discretization()->DofRowMap(0)));

  if (!matchinggrid_) k_ts_ = volcoupl_->ApplyMatrixMapping21(k_ts_);

  systemmatrix_->Assign(1, 0, LINALG::View, *k_ts_);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  systemmatrix_->Complete();

  // apply mortar coupling
  if (mortar_coupling_ != Teuchos::null) mortar_coupling_->CondenseMatrix(systemmatrix_);

}  // SetupSystemMatrix()


/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                                   dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupRHS()
{
#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0) std::cout << " TSI::Monolithic::SetupRHS()" << std::endl;
#endif  // TFSI
#endif  // TSI_DEBUG

  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::SetupRHS");

  // create full monolithic rhs vector
  rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  // get the structural rhs
  Teuchos::RCP<Epetra_Vector> str_rhs = Teuchos::rcp(new Epetra_Vector(*StructureField()->RHS()));
  if (DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(
          DRT::Problem::Instance()->StructuralDynamicParams(), "INT_STRATEGY") ==
      INPAR::STR::int_standard)
    str_rhs->Scale(-1.);

  // insert vectors to tsi rhs
  Extractor()->InsertVector(*str_rhs, 0, *rhs_);
  Extractor()->InsertVector(*ThermoField()->RHS(), 1, *rhs_);

  // apply mortar coupling
  if (mortar_coupling_ != Teuchos::null) mortar_coupling_->CondenseRhs(rhs_);

}  // SetupRHS()


/*----------------------------------------------------------------------*
 | solve linear TSI system                                   dano 04/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::LinearSolve()
{
#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0) std::cout << " TSI::Monolithic::LinearSolve()" << std::endl;
#endif  // TFSI
#endif  // TSI_DEBUG

  // Solve for inc_ = [disi_,tempi_]
  // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,T]
  // \f$x_{i+1} = x_i + \Delta x_i\f$
  if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normrhs_;
    double wanted = tolrhs_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }

  // Dirichlet boundary conditions are already applied to TSI system, i.e. TSI
  // system is prepared for solve, i.e. TSI systemmatrix, TSI rhs, TSI inc
  // --> in PrepareSystemForNewtonSolve(): done for rhs and diagonal blocks
  // --> in SetupSystemMatrix() done for off-diagonal blocks k_st, k_ts

  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  // default: use block matrix
  if (merge_tsi_blockmatrix_ == false)
  {
#ifdef TSI_DEBUG
#ifndef TFSI
    if (Comm().MyPID() == 0)
    {
      std::cout << " DBC applied to TSI system on proc" << Comm().MyPID() << std::endl;
    }
#endif  // TFSI
#endif  // TSI_DEBUG
    // Infnormscaling: scale system before solving
    ScaleSystem(*systemmatrix_, *rhs_);

    // solve the problem, work is done here!
    solver_->Solve(systemmatrix_->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);

    // Infnormscaling: unscale system after solving
    UnscaleSolution(*systemmatrix_, *iterinc_, *rhs_);

  }  // use block matrix

  else  // (merge_tsi_blockmatrix_ == true)
  {
    // merge blockmatrix to SparseMatrix and solve
    Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

    // standard solver call
    solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
  }  // MergeBlockMatrix

#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0)
  {
    std::cout << " Solved" << std::endl;
  }
#endif  // TFSI
#endif  // TSI_DEBUG

  // apply mortar coupling
  if (mortar_coupling_ != Teuchos::null) mortar_coupling_->RecoverIncr(iterinc_);

}  // LinearSolve()


/*----------------------------------------------------------------------*
 | initial guess of the displacements/temperatures           dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("TSI::Monolithic::InitialGuess");

  // InitalGuess() is called of the single fields and results are put in TSI
  // increment vector ig
  SetupVector(*ig,
      // returns residual displacements \f$\Delta D_{n+1}^{<k>}\f$ - disi_
      StructureField()->InitialGuess(),
      // returns residual temperatures or iterative thermal increment - tempi_
      ThermoField()->InitialGuess());

}  // InitialGuess()


/*----------------------------------------------------------------------*
 | setup vector of the structure and thermo field            dano 11/10 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetupVector(
    Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> tv)
{
  // extract dofs of the two fields
  // and put the structural/thermal field vector into the global vector f
  // noticing the block number
  Extractor()->InsertVector(*sv, 0, f);
  Extractor()->InsertVector(*tv, 1, f);

}  // SetupVector()


/*----------------------------------------------------------------------*
 | check convergence of Newton iteration (public)            dano 11/10 |
 *----------------------------------------------------------------------*/
bool TSI::Monolithic::Converged()
{
  // check for single norms
  bool convrhs = false;
  bool convinc = false;
  bool convstrrhs = false;
  bool convdisp = false;
  bool convthrrhs = false;
  bool convtemp = false;

  // ----------------------------------------------------------- TSI test
  // residual TSI forces
  switch (normtyperhs_)
  {
    case INPAR::TSI::convnorm_abs:
      convrhs = normrhs_ < tolrhs_;
      break;
    case INPAR::TSI::convnorm_rel:
      convrhs = normrhs_ < std::max(tolrhs_ * normrhsiter0_, 1.0e-15);
      break;
    case INPAR::TSI::convnorm_mix:
      convrhs = ((normrhs_ < tolrhs_) and (normrhs_ < std::max(normrhsiter0_ * tolrhs_, 1.0e-15)));
      break;
    default:
      dserror("Cannot check for convergence of residual forces!");
      break;
  }

  // residual TSI increments
  switch (normtypeinc_)
  {
    case INPAR::TSI::convnorm_abs:
      convinc = norminc_ < tolinc_;
      break;
    case INPAR::TSI::convnorm_rel:
      convinc = norminc_ < std::max(norminciter0_ * tolinc_, 1e-15);
      break;
    case INPAR::TSI::convnorm_mix:
      convinc = norminc_ < std::max(norminciter0_ * tolinc_, 1e-15);
      break;
    default:
      dserror("Cannot check for convergence of increments!");
      break;
  }  // switch (normtypeinc_)

  // -------------------------------------------------- single field test
  // ---------------------------------------------------------- structure
  // structural residual forces
  switch (normtypestrrhs_)
  {
    case INPAR::STR::convnorm_abs:
      convstrrhs = normstrrhs_ < tolstrrhs_;
      break;
    case INPAR::STR::convnorm_rel:
      convstrrhs = normstrrhs_ < std::max(normstrrhsiter0_ * tolstrrhs_, 1e-15);
      break;
    case INPAR::STR::convnorm_mix:
      convstrrhs = ((normstrrhs_ < tolstrrhs_) or
                    (normstrrhs_ < std::max(normstrrhsiter0_ * tolstrrhs_, 1e-15)));
      break;
    default:
      dserror("Cannot check for convergence of residual forces!");
      break;
  }  // switch (normtypestrrhs_)

  // residual displacements
  switch (normtypedisi_)
  {
    case INPAR::STR::convnorm_abs:
      convdisp = normdisi_ < toldisi_;
      break;
    case INPAR::STR::convnorm_rel:
      convdisp = normdisi_ < std::max(normdisiiter0_ * toldisi_, 1e-15);
      break;
    case INPAR::STR::convnorm_mix:
      convdisp =
          ((normdisi_ < toldisi_) or (normdisi_ < std::max(normdisiiter0_ * toldisi_, 1e-15)));
      break;
    default:
      dserror("Cannot check for convergence of displacements!");
      break;
  }  // switch (normtypedisi_)

  // ------------------------------------------------------------- thermo
  // thermal residual forces
  switch (normtypethrrhs_)
  {
    case INPAR::THR::convnorm_abs:
      convthrrhs = normthrrhs_ < tolthrrhs_;
      break;
    case INPAR::THR::convnorm_rel:
      convthrrhs = normthrrhs_ < normthrrhsiter0_ * tolthrrhs_;
      break;
    case INPAR::THR::convnorm_mix:
      convthrrhs = ((normthrrhs_ < tolthrrhs_) or (normthrrhs_ < normthrrhsiter0_ * tolthrrhs_));
      break;
    default:
      dserror("Cannot check for convergence of residual forces!");
      break;
  }  // switch (normtypethrrhs_)

  // residual temperatures
  switch (normtypetempi_)
  {
    case INPAR::THR::convnorm_abs:
      convtemp = normtempi_ < toltempi_;
      break;
    case INPAR::THR::convnorm_rel:
      convtemp = normtempi_ < normtempiiter0_ * toltempi_;
      break;
    case INPAR::THR::convnorm_mix:
      convtemp = ((normtempi_ < toltempi_) or (normtempi_ < normtempiiter0_ * toltempi_));
      break;
    default:
      dserror("Cannot check for convergence of temperatures!");
      break;
  }  // switch (normtypetempi_)

  // -------------------------------------------------------- convergence
  // combine increment-like and force-like residuals, combine TSI and single
  // field values
  bool conv = false;
  if (combincrhs_ == INPAR::TSI::bop_and)
    conv = convinc and convrhs;
  else if (combincrhs_ == INPAR::TSI::bop_or)
    conv = convinc or convrhs;
  else if (combincrhs_ == INPAR::TSI::bop_coupl_and_singl)
    conv = convinc and convrhs and convdisp and convstrrhs and convtemp and convthrrhs;
  else if (combincrhs_ == INPAR::TSI::bop_coupl_or_singl)
    conv = (convinc and convrhs) or (convdisp and convstrrhs and convtemp and convthrrhs);
  else if (combincrhs_ == INPAR::TSI::bop_and_singl)
    conv = convdisp and convstrrhs and convtemp and convthrrhs;
  else if (combincrhs_ == INPAR::TSI::bop_or_singl)
    conv = (convdisp or convstrrhs or convtemp or convthrrhs);
  else
    dserror("Something went terribly wrong with binary operator!");

  // convergence of active contact set
  if (contact_strategy_lagrange_ != Teuchos::null)
  {
    conv = conv && (contact_strategy_lagrange_->mech_contact_res_ <
                       contact_strategy_lagrange_->Params().get<double>("TOLCONTCONSTR"));
    conv = conv && (contact_strategy_lagrange_->mech_contact_incr_ <
                       contact_strategy_lagrange_->Params().get<double>("TOLLAGR"));
    conv = conv && contact_strategy_lagrange_->ActiveSetSemiSmoothConverged();
  }

  // return things
  return conv;

}  // Converged()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file   dano 11/10 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrintNewtonIter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ((Comm().MyPID() == 0) and PrintScreenEvry() and (Step() % PrintScreenEvry() == 0) and
      printiter_)
  {
    if (iter_ == 0) PrintNewtonIterHeader(stdout);
    PrintNewtonIterText(stdout);
  }

  // print to error file
  if (printerrfile_ and printiter_)
  {
    if (iter_ == 1) PrintNewtonIterHeader(errfile_);
    PrintNewtonIterText(errfile_);
  }

  // see you
  return;
}  // PrintNewtonIter()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file   dano 11/10 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrintNewtonIterHeader(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(5) << "# iter";

  // line search
  if (ls_strategy_) oss << std::setw(11) << " ls_step";

  // ---------------------------------------------------------------- TSI
  // different style due relative or absolute error checking
  // displacement
  switch (normtyperhs_)
  {
    case INPAR::TSI::convnorm_abs:
      oss << std::setw(15) << "abs-res-norm";
      break;
    case INPAR::TSI::convnorm_rel:
      oss << std::setw(15) << "rel-res-norm";
      break;
    case INPAR::TSI::convnorm_mix:
      oss << std::setw(15) << "mix-res-norm";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::TSI::convnorm_abs:
      oss << std::setw(15) << "abs-inc-norm";
      break;
    case INPAR::TSI::convnorm_rel:
      oss << std::setw(15) << "rel-inc-norm";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  // -------------------------------------------------- single field test
  // ---------------------------------------------------------- structure
  switch (normtypestrrhs_)
  {
    case INPAR::STR::convnorm_rel:
      oss << std::setw(18) << "rel-str-res-norm";
      break;
    case INPAR::STR::convnorm_abs:
      oss << std::setw(18) << "abs-str-res-norm";
      break;
    case INPAR::STR::convnorm_mix:
      oss << std::setw(18) << "mix-str-res-norm";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }  // switch (normtypestrrhs_)

  switch (normtypedisi_)
  {
    case INPAR::STR::convnorm_rel:
      oss << std::setw(16) << "rel-dis-norm";
      break;
    case INPAR::STR::convnorm_abs:
      oss << std::setw(16) << "abs-dis-norm";
      break;
    case INPAR::STR::convnorm_mix:
      oss << std::setw(16) << "mix-dis-norm";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }  // switch (normtypedisi_)

  // ------------------------------------------------------------- thermo
  switch (normtypethrrhs_)
  {
    case INPAR::THR::convnorm_rel:
      oss << std::setw(18) << "rel-thr-res-norm";
      break;
    case INPAR::THR::convnorm_abs:
      oss << std::setw(18) << "abs-thr-res-norm";
      break;
    case INPAR::THR::convnorm_mix:
      oss << std::setw(18) << "mix-thr-res-norm";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }  // switch (normtypethrrhs_)

  switch (normtypetempi_)
  {
    case INPAR::THR::convnorm_rel:
      oss << std::setw(16) << "rel-temp-norm";
      break;
    case INPAR::THR::convnorm_abs:
      oss << std::setw(16) << "abs-temp-norm";
      break;
    case INPAR::THR::convnorm_mix:
      oss << std::setw(16) << "mix-temp-norm";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }  // switch (normtypetempi_)

  if (contact_strategy_lagrange_ != Teuchos::null)
  {
    oss << std::setw(16) << "sLMres-norm";
    oss << std::setw(16) << "sLMinc-norm";
    oss << std::setw(16) << "tLMinc-norm";
  }


  if (soltech_ == INPAR::TSI::soltech_ptc)
  {
    oss << std::setw(16) << "        PTC-dti";
  }

  // add solution time
  oss << std::setw(12) << "ts";
  oss << std::setw(12) << "wct";

  // add contact set information
  if (contact_strategy_lagrange_ != Teuchos::null)
  {
    oss << std::setw(12) << "tc";
    oss << std::setw(11) << "#active";
    if (contact_strategy_lagrange_->Friction()) oss << std::setw(10) << "#slip";
  }

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == NULL) dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // PrintNewtonIterHeader()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen                  dano 11/10 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrintNewtonIterText(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(6) << iter_;

  // line search step
  if (ls_strategy_)
    oss << std::setw(11) << std::setprecision(3) << std::scientific << ls_step_length_;

  // different style due relative or absolute error checking

  // ----------------------------------------------- test coupled problem
  switch (normtyperhs_)
  {
    case INPAR::TSI::convnorm_abs:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhs_;
      break;
    case INPAR::TSI::convnorm_rel:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhs_ / normrhsiter0_;
      break;
    case INPAR::TSI::convnorm_mix:
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << std::min(normrhs_, normrhs_ / normrhsiter0_);
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::TSI::convnorm_abs:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << norminc_;
      break;
    case INPAR::TSI::convnorm_rel:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << norminc_ / norminciter0_;
      break;
    case INPAR::TSI::convnorm_mix:
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << std::min(norminc_, norminc_ / norminciter0_);
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }  // switch (normtypeinc_)

  // ------------------------------------------------- test single fields
  // ---------------------------------------------------------- structure
  // different style due relative or absolute error checking
  // displacement
  switch (normtypestrrhs_)
  {
    case INPAR::STR::convnorm_abs:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normstrrhs_;
      break;
    case INPAR::STR::convnorm_rel:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << normstrrhs_ / normstrrhsiter0_;
      break;
    case INPAR::STR::convnorm_mix:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << std::min(normstrrhs_, normstrrhs_ / normstrrhsiter0_);
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }  // switch (normtypestrrhs_)

  switch (normtypedisi_)
  {
    case INPAR::STR::convnorm_abs:
      oss << std::setw(16) << std::setprecision(5) << std::scientific << normdisi_;
      break;
    case INPAR::STR::convnorm_rel:
      oss << std::setw(16) << std::setprecision(5) << std::scientific << normdisi_ / normdisiiter0_;
      break;
    case INPAR::STR::convnorm_mix:
      oss << std::setw(16) << std::setprecision(5) << std::scientific
          << std::min(normdisi_, normdisi_ / normdisiiter0_);
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }  // switch (normtypedisi_)

  // ------------------------------------------------------------- thermo
  switch (normtypethrrhs_)
  {
    case INPAR::THR::convnorm_abs:
      oss << std::setw(18) << std::setprecision(5) << std::scientific << normthrrhs_;
      break;
    case INPAR::THR::convnorm_rel:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << normthrrhs_ / normthrrhsiter0_;
      break;
    case INPAR::THR::convnorm_mix:
      oss << std::setw(18) << std::setprecision(5) << std::scientific
          << std::min(normthrrhs_, normthrrhs_ / normthrrhsiter0_);
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }  // switch (normtypethrrhs_)

  switch (normtypetempi_)
  {
    case INPAR::THR::convnorm_abs:
      oss << std::setw(16) << std::setprecision(5) << std::scientific << normtempi_;
      break;
    case INPAR::THR::convnorm_rel:
      oss << std::setw(16) << std::setprecision(5) << std::scientific
          << normtempi_ / normtempiiter0_;
      break;
    case INPAR::THR::convnorm_mix:
      oss << std::setw(16) << std::setprecision(5) << std::scientific
          << std::min(normtempi_, normtempi_ / normtempiiter0_);
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }  // switch (normtypetempi_)

  if (contact_strategy_lagrange_ != Teuchos::null)
  {
    oss << std::setw(16) << std::setprecision(5) << std::scientific
        << contact_strategy_lagrange_->mech_contact_res_;
    oss << std::setw(16) << std::setprecision(5) << std::scientific
        << contact_strategy_lagrange_->mech_contact_incr_;
    oss << std::setw(16) << std::setprecision(5) << std::scientific
        << contact_strategy_lagrange_->thr_contact_incr_;
  }

  if (soltech_ == INPAR::TSI::soltech_ptc)
  {
    oss << std::setw(16) << std::setprecision(5) << std::scientific << dti_;
  }

  // add solution time of to print to screen
  oss << std::setw(12) << std::setprecision(2) << std::scientific << dtsolve_;
  oss << std::setw(12) << std::setprecision(2) << std::scientific << timernewton_.ElapsedTime();

  // add contact information
  if (contact_strategy_lagrange_ != Teuchos::null)
  {
    oss << std::setw(12) << std::setprecision(2) << std::scientific << dtcmt_;
    oss << std::setw(11) << contact_strategy_lagrange_->NumberOfActiveNodes();
    if (contact_strategy_lagrange_->Friction())
      oss << std::setw(10) << contact_strategy_lagrange_->NumberOfSlipNodes();
  }

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == NULL) dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // PrintNewtonIterText


/*----------------------------------------------------------------------*
 | print statistics of converged NRI                         dano 11/10 |
 | orignially by bborn 08/09                                            |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrintNewtonConv()
{
  // somebody did the door
  return;
}  // PrintNewtonConv()


/*----------------------------------------------------------------------*
 | evaluate mechanical-thermal system matrix at state        dano 03/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyStrCouplMatrix(
    Teuchos::RCP<LINALG::SparseMatrix> k_st  //!< off-diagonal tangent matrix term
)
{
#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0) std::cout << " TSI::Monolithic::ApplyStrCouplMatrix()" << std::endl;
#endif  // TFSI
#endif  // TSI_DEBUG

  // create the parameters for the discretization
  Teuchos::ParameterList sparams;

  //! pointer to the model evaluator data container
  Teuchos::RCP<DRT::ELEMENTS::ParamsMinimal> EvalData =
      Teuchos::rcp(new DRT::ELEMENTS::ParamsMinimal());

  // set parameters needed for element evalutation
  EvalData->SetActionType(DRT::ELEMENTS::struct_calc_stifftemp);
  EvalData->SetTotalTime(Time());
  EvalData->SetDeltaTime(Dt());

  const std::string action = "calc_struct_stifftemp";
  sparams.set("action", action);
  // other parameters that might be needed by the elements
  sparams.set("delta time", Dt());
  sparams.set("total time", Time());

  StructureField()->Discretization()->ClearState(true);
  StructureField()->Discretization()->SetState(0, "displacement", StructureField()->Dispnp());

  // in case of temperature-dependent material parameters, here E(T), T_{n+1} is required in STR
  sparams.set<int>(
      "young_temp", (DRT::INPUT::IntegralValue<int>(sdyn_, "YOUNG_IS_TEMP_DEPENDENT")));

  ApplyThermoCouplingState(ThermoField()->Tempnp());

  // build specific assemble strategy for mechanical-thermal system matrix
  // from the point of view of StructureField:
  // structdofset = 0, thermdofset = 1
  DRT::AssembleStrategy structuralstrategy(0,  // structdofset for row
      1,                                       // thermdofset for column
      k_st,                                    // build mechanical-thermal matrix
      Teuchos::null,                           // no other matrix or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);


  sparams.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", EvalData);
  StructureField()->Discretization()->Evaluate(sparams, structuralstrategy);
  StructureField()->Discretization()->ClearState(true);

  // add nitsche contact integral
  if (contact_strategy_nitsche_ != Teuchos::null)
    k_st->Add(
        *contact_strategy_nitsche_->GetMatrixBlockPtr(DRT::UTILS::block_displ_temp), false, 1., 1.);

  // TODO 2013-11-11 move scaling to the so3_thermo element
  // --> consistent with thermo element and clearer, more consistent

  // for consistent linearisation scale k_st with time factor
  // major switch to different time integrators
  switch (strmethodname_)
  {
    case INPAR::STR::dyna_statics:
    {
      // continue
      break;
    }
    case INPAR::STR::dyna_onesteptheta:
    {
      double theta = sdyn_.sublist("ONESTEPTHETA").get<double>("THETA");
      // K_Teffdyn(T_n+1^i) = theta * k_st
      k_st->Scale(theta);
      break;
    }
    case INPAR::STR::dyna_genalpha:
    {
      double alphaf = sdyn_.sublist("GENALPHA").get<double>("ALPHA_F");
      // K_Teffdyn(T_n+1) = (1-alphaf_) . kst
      // Lin(dT_n+1-alphaf_/ dT_n+1) = (1-alphaf_)
      k_st->Scale(1.0 - alphaf);
      break;
    }
    default:
      dserror("Don't know what to do...");
      break;
  }  // end of switch(strmethodname_)

}  // ApplyStrCouplMatrix()


/*----------------------------------------------------------------------*
 | evaluate thermal-mechanical system matrix at state        dano 03/11 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyThrCouplMatrix(
    Teuchos::RCP<LINALG::SparseMatrix> k_ts  //!< off-diagonal tangent matrix term
)
{
#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0) std::cout << " TSI::Monolithic::ApplyThrCouplMatrix()" << std::endl;
#endif  // TFSI
#endif  // TSI_DEBUG

  // create the parameters for the discretization
  Teuchos::ParameterList tparams;
  // action for elements
  const THR::Action action = THR::calc_thermo_coupltang;
  tparams.set<int>("action", action);
  // other parameters that might be needed by the elements
  tparams.set("delta time", Dt());
  tparams.set("total time", Time());
  tparams.set<int>(
      "young_temp", (DRT::INPUT::IntegralValue<int>(sdyn_, "YOUNG_IS_TEMP_DEPENDENT")));

  // create specific time integrator
  const Teuchos::ParameterList& tdyn = DRT::Problem::Instance()->ThermalDynamicParams();
  tparams.set<int>(
      "time integrator", DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn, "DYNAMICTYP"));
  tparams.set<int>("structural time integrator", strmethodname_);
  double timefac = -1.;
  switch (DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn, "DYNAMICTYP"))
  {
    // static analysis
    case INPAR::THR::dyna_statics:
    {
      timefac = 1.;
      // continue
      break;
    }
    // dynamic analysis
    case INPAR::THR::dyna_onesteptheta:
    {
      // K_Td = theta . k_Td^e
      double theta = tdyn.sublist("ONESTEPTHETA").get<double>("THETA");
      tparams.set("theta", theta);
      timefac = theta;
      break;
    }
    case INPAR::THR::dyna_genalpha:
    {
      double alphaf = tdyn.sublist("GENALPHA").get<double>("ALPHA_F");
      tparams.set("alphaf", alphaf);

      timefac = alphaf;
      break;
    }
    case INPAR::THR::dyna_undefined:
    default:
    {
      dserror("Don't know what to do...");
      break;
    }
  }  // switch (THR::DynamicType)

  switch (strmethodname_)
  {
    case INPAR::STR::dyna_statics:
    {
      // continue
      break;
    }
    case INPAR::STR::dyna_onesteptheta:
    {
      // put the structural theta value to the thermal parameter list
      double str_theta = sdyn_.sublist("ONESTEPTHETA").get<double>("THETA");
      tparams.set("str_theta", str_theta);
      break;
    }
    case INPAR::STR::dyna_genalpha:
    {
      // put the structural theta value to the thermal parameter list
      double str_beta = sdyn_.sublist("GENALPHA").get<double>("BETA");
      double str_gamma = sdyn_.sublist("GENALPHA").get<double>("GAMMA");

      tparams.set("str_beta", str_beta);
      tparams.set("str_gamma", str_gamma);
      break;
    }
    default:
      dserror("Don't know what to do...");
      break;
  }


  ThermoField()->Discretization()->ClearState(true);
  // set the variables that are needed by the elements
  ThermoField()->Discretization()->SetState(0, "temperature", ThermoField()->Tempnp());

  ApplyStructCouplingState(StructureField()->Dispnp(), vel_);

  // build specific assemble strategy for the thermal-mechanical system matrix
  // from the point of view of ThermoField:
  // thermdofset = 0, structdofset = 1
  DRT::AssembleStrategy thermostrategy(0,  // thermdofset for row
      1,                                   // structdofset for column
      k_ts,                                // thermal-mechanical matrix
      Teuchos::null,                       // no other matrix or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate the thermal-mechanical system matrix on the thermal element
  ThermoField()->Discretization()->Evaluate(tparams, thermostrategy);
  ThermoField()->Discretization()->ClearState(true);

  // add nitsche contact integral
  if (contact_strategy_nitsche_ != Teuchos::null)
    k_ts->Add(*contact_strategy_nitsche_->GetMatrixBlockPtr(DRT::UTILS::block_temp_displ), false,
        timefac, 1.);
}  // ApplyThrCouplMatrix()


/*----------------------------------------------------------------------*
 | evaluate thermal-mechanical system matrix at state        dano 12/12 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyThrCouplMatrix_ConvBC(
    Teuchos::RCP<LINALG::SparseMatrix> k_ts  //!< off-diagonal tangent matrix term
)
{
#ifdef TSI_DEBUG
#ifndef TFSI
  if (Comm().MyPID() == 0)
    std::cout << " TSI::Monolithic::ApplyThrCouplMatrix_ConvBC()" << std::endl;
#endif  // TFSI
#endif  // TSI_DEBUG

  std::vector<DRT::Condition*> cond;
  std::string condstring("ThermoConvections");
  ThermoField()->Discretization()->GetCondition(condstring, cond);
  if (cond.size() > 0)
  {
    // create the parameters for the discretization
    Teuchos::ParameterList tparams;
    // action for elements
    const THR::BoundaryAction action = THR::calc_thermo_fextconvection_coupltang;
    tparams.set<int>("action", action);
    // other parameters that might be needed by the elements
    tparams.set("delta time", Dt());
    tparams.set("total time", Time());
    // create specific time integrator
    const Teuchos::ParameterList& tdyn = DRT::Problem::Instance()->ThermalDynamicParams();
    tparams.set<int>(
        "time integrator", DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn, "DYNAMICTYP"));
    tparams.set<int>("structural time integrator", strmethodname_);
    switch (DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(tdyn, "DYNAMICTYP"))
    {
      // static analysis
      case INPAR::THR::dyna_statics:
      {
        break;
      }
      // dynamic analysis
      case INPAR::THR::dyna_onesteptheta:
      {
        // K_Td = theta . k_Td^e
        double theta = tdyn.sublist("ONESTEPTHETA").get<double>("THETA");
        tparams.set("theta", theta);
        // put the structural theta value to the thermal parameter list
        double str_theta = sdyn_.sublist("ONESTEPTHETA").get<double>("THETA");
        tparams.set("str_theta", str_theta);
        break;
      }
      case INPAR::THR::dyna_genalpha:
      {
        // K_Td = alphaf . k_Td^e
        double alphaf = tdyn.sublist("GENALPHA").get<double>("ALPHA_F");
        tparams.set("alphaf", alphaf);

        // put the structural theta value to the thermal parameter list
        double str_beta = sdyn_.sublist("GENALPHA").get<double>("BETA");
        double str_gamma = sdyn_.sublist("GENALPHA").get<double>("GAMMA");
        tparams.set("str_beta", str_beta);
        tparams.set("str_gamma", str_gamma);
        break;
      }
      case INPAR::THR::dyna_undefined:
      default:
      {
        dserror("Don't know what to do...");
        break;
      }
    }  // end(switch)
    // clear all states set in discretization
    ThermoField()->Discretization()->ClearState(true);
    // set the variables that are needed by the elements
    ThermoField()->Discretization()->SetState(0, "temperature", ThermoField()->Tempnp());
    ApplyStructCouplingState(StructureField()->Dispnp(), vel_);

    // build specific assemble strategy for the thermal-mechanical system matrix
    // from the point of view of ThermoField:
    // thermdofset = 0, structdofset = 1
    DRT::AssembleStrategy thermostrategy(0,  // thermdofset for row
        1,                                   // structdofset for column
        k_ts,                                // thermal-mechanical matrix
        Teuchos::null,                       // no other matrix or vectors
        Teuchos::null, Teuchos::null, Teuchos::null);

    // evaluate the thermal-mechanical system matrix on the thermal element
    ThermoField()->Discretization()->EvaluateCondition(tparams, thermostrategy, condstring);
    // clear all states set in discretization
    ThermoField()->Discretization()->ClearState(true);
  }  // cond.size()>0

}  // ApplyThrCouplMatrix()


/*----------------------------------------------------------------------*
 | map containing the dofs with Dirichlet BC                 dano 03/11 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> TSI::Monolithic::CombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map> scondmap = StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> tcondmap = ThermoField()->GetDBCMapExtractor()->CondMap();
  Teuchos::RCP<Epetra_Map> condmap = LINALG::MergeMap(scondmap, tcondmap, false);
  return condmap;

}  // CombinedDBCMap()



/*----------------------------------------------------------------------*
 | recover structural and thermal Lagrange multipliers from  seitz 11/15|
 | displacements and temperature                                        |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::RecoverStructThermLM()
{
  // only in the case of contact
  if (contact_strategy_lagrange_ == Teuchos::null) return;

  // split the increment
  Teuchos::RCP<Epetra_Vector> sx;
  Teuchos::RCP<Epetra_Vector> tx;

  // extract field vectors
  ExtractFieldVectors(iterinc_, sx, tx);

  contact_strategy_lagrange_->RecoverCoupled(sx, tx, coupST_);

  return;
}


/*----------------------------------------------------------------------*
 | scale system, i.e. apply infnorm scaling to linear        dano 02/13 |
 | block system before solving system                                   |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  // should we scale the system?
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(tsidynmono_, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if ((A->LeftScale(*srowsum_)) or (A->RightScale(*scolsum_)) or
        (mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_)) or
        (mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_)))
      dserror("structure scaling failed");

    A = mat.Matrix(1, 1).EpetraMatrix();
    trowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    tcolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*trowsum_);
    A->InvColSums(*tcolsum_);
    if ((A->LeftScale(*trowsum_)) or (A->RightScale(*tcolsum_)) or
        (mat.Matrix(1, 0).EpetraMatrix()->LeftScale(*trowsum_)) or
        (mat.Matrix(0, 1).EpetraMatrix()->RightScale(*tcolsum_)))
      dserror("thermo scaling failed");

    Teuchos::RCP<Epetra_Vector> sx = Extractor()->ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> tx = Extractor()->ExtractVector(b, 1);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0)) dserror("structure scaling failed");
    if (tx->Multiply(1.0, *trowsum_, *tx, 0.0)) dserror("thermo scaling failed");

    Extractor()->InsertVector(*sx, 0, b);
    Extractor()->InsertVector(*tx, 1, b);
  }
}  // ScaleSystem


/*----------------------------------------------------------------------*
 | unscale solution after solving the linear system          dano 02/13 |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::UnscaleSolution(
    LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(tsidynmono_, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = Extractor()->ExtractVector(x, 0);
    Teuchos::RCP<Epetra_Vector> ty = Extractor()->ExtractVector(x, 1);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0)) dserror("structure scaling failed");
    if (ty->Multiply(1.0, *tcolsum_, *ty, 0.0)) dserror("thermo scaling failed");

    Extractor()->InsertVector(*sy, 0, x);
    Extractor()->InsertVector(*ty, 1, x);

    Teuchos::RCP<Epetra_Vector> sx = Extractor()->ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> tx = Extractor()->ExtractVector(b, 1);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0)) dserror("structure scaling failed");
    if (tx->ReciprocalMultiply(1.0, *trowsum_, *tx, 0.0)) dserror("thermo scaling failed");

    Extractor()->InsertVector(*sx, 0, b);
    Extractor()->InsertVector(*tx, 1, b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if ((A->LeftScale(*srowsum_)) or (A->RightScale(*scolsum_)) or
        (mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_)) or
        (mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_)))
      dserror("structure scaling failed");

    A = mat.Matrix(1, 1).EpetraMatrix();
    trowsum_->Reciprocal(*trowsum_);
    tcolsum_->Reciprocal(*tcolsum_);
    if ((A->LeftScale(*trowsum_)) or (A->RightScale(*tcolsum_)) or
        (mat.Matrix(1, 0).EpetraMatrix()->LeftScale(*trowsum_)) or
        (mat.Matrix(0, 1).EpetraMatrix()->RightScale(*tcolsum_)))
      dserror("thermo scaling failed");

  }  // if (scaling_infnorm)

}  // UnscaleSolution()


/*----------------------------------------------------------------------*
 | calculate vector norm                                     dano 04/13 |
 *----------------------------------------------------------------------*/
double TSI::Monolithic::CalculateVectorNorm(
    const enum INPAR::TSI::VectorNorm norm, const Teuchos::RCP<const Epetra_Vector> vect)
{
  // L1 norm
  // norm = sum_0^i vect[i]
  if (norm == INPAR::TSI::norm_l1)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm;
  }
  // L2/Euclidian norm
  // norm = sqrt{sum_0^i vect[i]^2 }
  else if (norm == INPAR::TSI::norm_l2)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm;
  }
  // RMS norm
  // norm = sqrt{sum_0^i vect[i]^2 }/ sqrt{length_vect}
  else if (norm == INPAR::TSI::norm_rms)
  {
    double vectnorm;
    vect->Norm2(&vectnorm);
    return vectnorm / sqrt((double)vect->GlobalLength());
  }
  // infinity/maximum norm
  // norm = max( vect[i] )
  else if (norm == INPAR::TSI::norm_inf)
  {
    double vectnorm;
    vect->NormInf(&vectnorm);
    return vectnorm;
  }
  // norm = sum_0^i vect[i]/length_vect
  else if (norm == INPAR::TSI::norm_l1_scaled)
  {
    double vectnorm;
    vect->Norm1(&vectnorm);
    return vectnorm / ((double)vect->GlobalLength());
  }
  else
  {
    dserror("Cannot handle vector norm");
    return 0;
  }
}  // CalculateVectorNorm()


/*----------------------------------------------------------------------*
 | set parameters for TSI remaining constant over whole      dano 04/13 |
 | simulation                                                           |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::SetDefaultParameters()
{
  // time parameters
  // call the TSI parameter list
  const Teuchos::ParameterList& tdyn = DRT::Problem::Instance()->ThermalDynamicParams();

  // get the parameters for the Newton iteration
  itermax_ = tsidyn_.get<int>("ITEMAX");
  itermin_ = tsidyn_.get<int>("ITEMIN");

  // what kind of norm do we wanna test for coupled TSI problem
  normtypeinc_ = DRT::INPUT::IntegralValue<INPAR::TSI::ConvNorm>(tsidyn_, "NORM_INC");
  normtyperhs_ = DRT::INPUT::IntegralValue<INPAR::TSI::ConvNorm>(tsidynmono_, "NORM_RESF");
  // what kind of norm do we wanna test for the single fields
  normtypedisi_ = DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(sdyn_, "NORM_DISP");
  normtypestrrhs_ = DRT::INPUT::IntegralValue<INPAR::STR::ConvNorm>(sdyn_, "NORM_RESF");
  enum INPAR::STR::VectorNorm striternorm =
      DRT::INPUT::IntegralValue<INPAR::STR::VectorNorm>(sdyn_, "ITERNORM");
  normtypetempi_ = DRT::INPUT::IntegralValue<INPAR::THR::ConvNorm>(tdyn, "NORM_TEMP");
  normtypethrrhs_ = DRT::INPUT::IntegralValue<INPAR::THR::ConvNorm>(tdyn, "NORM_RESF");
  enum INPAR::THR::VectorNorm thriternorm =
      DRT::INPUT::IntegralValue<INPAR::THR::VectorNorm>(tdyn, "ITERNORM");
  // in total when do we reach a converged state for complete problem
  combincrhs_ = DRT::INPUT::IntegralValue<INPAR::TSI::BinaryOp>(tsidynmono_, "NORMCOMBI_RESFINC");

#ifndef TFSI
  switch (combincrhs_)
  {
    case INPAR::TSI::bop_and:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of TSI:\n res, inc with 'AND'." << std::endl;
      break;
    }
    case INPAR::TSI::bop_or:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of TSI:\n res, inc with 'OR'." << std::endl;
      break;
    }
    case INPAR::TSI::bop_coupl_and_singl:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of TSI:\n res, inc, str-res, thr-res, dis, temp with 'AND'."
                  << std::endl;
      break;
    }
    case INPAR::TSI::bop_coupl_or_singl:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of TSI:\n (res, inc) or (str-res, thr-res, dis, temp)."
                  << std::endl;
      break;
    }
    case INPAR::TSI::bop_and_singl:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of TSI:\n str-res, thr-res, dis, temp with 'AND'."
                  << std::endl;
      break;
    }
    case INPAR::TSI::bop_or_singl:
    {
      if (Comm().MyPID() == 0)
        std::cout << "Convergence test of TSI:\n str-res, thr-res, dis, temp with 'OR'."
                  << std::endl;
      break;
    }
    default:
    {
      dserror("Something went terribly wrong with binary operator!");
      break;
    }
  }  // switch (combincrhs_)
#endif

  // convert the single field norms to be used within TSI
  // what norm is used for structure
  switch (striternorm)
  {
    case INPAR::STR::norm_l1:
      iternormstr_ = INPAR::TSI::norm_l1;
      break;
    case INPAR::STR::norm_l2:
      iternormstr_ = INPAR::TSI::norm_l2;
      break;
    case INPAR::STR::norm_rms:
      iternormstr_ = INPAR::TSI::norm_rms;
      break;
    case INPAR::STR::norm_inf:
      iternormstr_ = INPAR::TSI::norm_inf;
      break;
    case INPAR::STR::norm_vague:
    default:
      dserror("STR norm is not determined");
      break;
  }  // switch (striternorm)

  // what norm is used for thermo
  switch (thriternorm)
  {
    case INPAR::THR::norm_l1:
      iternormthr_ = INPAR::TSI::norm_l1;
      break;
    case INPAR::THR::norm_l2:
      iternormthr_ = INPAR::TSI::norm_l2;
      break;
    case INPAR::THR::norm_rms:
      iternormthr_ = INPAR::TSI::norm_rms;
      break;
    case INPAR::THR::norm_inf:
      iternormthr_ = INPAR::TSI::norm_inf;
      break;
    case INPAR::THR::norm_vague:
    default:
    {
      dserror("THR norm is not determined.");
      break;
    }
  }  // switch (thriternorm)

  // if scaled L1-norm is wished to be used
  if ((iternorm_ == INPAR::TSI::norm_l1_scaled) and
      ((combincrhs_ == INPAR::TSI::bop_coupl_and_singl) or
          (combincrhs_ == INPAR::TSI::bop_coupl_or_singl)))
  {
    iternormstr_ = INPAR::TSI::norm_l1_scaled;
    iternormthr_ = INPAR::TSI::norm_l1_scaled;
  }

  // test the TSI-residual and the TSI-increment
  tolinc_ = tsidynmono_.get<double>("TOLINC");
  tolrhs_ = tsidynmono_.get<double>("CONVTOL");

  // get the single field tolerances from this field itselves
  toldisi_ = sdyn_.get<double>("TOLDISP");
  tolstrrhs_ = sdyn_.get<double>("TOLRES");
  toltempi_ = tdyn.get<double>("TOLTEMP");
  tolthrrhs_ = tdyn.get<double>("TOLRES");

  // initialise norms for coupled TSI problem
  normrhs_ = 0.0;
  normrhsiter0_ = 0.0;
  norminc_ = 0.0;
  norminciter0_ = 0.0;

  // initialise norms for single field tests
  normdisi_ = 0.0;
  normstrrhs_ = 0.0;
  normstrrhsiter0_ = 0.0;
  normtempi_ = 0.0;
  normthrrhs_ = 0.0;
  normthrrhsiter0_ = 0.0;

  return;

}  // SetDefaultParameter()


/*----------------------------------------------------------------------*
 | calculate nodal TSI results for evaluation                dano 09/13 |
 | used for thermoplasticity and necking to calculate displacements,    |
 | temperatures, and reaction forces at different points                |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::CalculateNeckingTSIResults()
{
  // be aware, parallel bug arises using baci-debug

  // --------------------------------------------- initialise/ define constants

  // initialise initial temperature required to calculate temperature increase
  const double inittemp = 293.0;

  // necking point A
  const double necking_x = 6.413;
  const double necking_y = 0.0;
  const double necking_z = -13.3335;
  // --> point A is extracted via restriction of coordinates

  // top point B
  const double top_x = 6.413;
  const double top_y = 0.0;
  const double top_z = 13.3335;
  // --> point B is extracted via DBC

  //---------------------------------------------------------------------------
  // -------------- get the nodes with STRUCTURAL Dirichlet boundary conditions
  //---------------------------------------------------------------------------

  // initialise a vector containing all structural DBC
  std::vector<DRT::Condition*> dbc(0);
  StructureField()->Discretization()->GetCondition("Dirichlet", dbc);

  // initialise a vector contatining all DBC in a special direction (here: in z)
  std::vector<int> one_dof_in_dbc(1);
  one_dof_in_dbc.at(0) = -1;

  // local list of found structural DOF IDs that have a DBC
  // in this special case DBC are applied at the nodes which are interested for
  // evaluation
  std::vector<int> sdata(0);

  // loop over Dirichlet boundary conditions
  for (int i = 0; i < (int)dbc.size(); ++i)
  {
    // get nodes which have DBCs
    const std::vector<int>* nodeids_withdbc = dbc[i]->Nodes();
    if (!nodeids_withdbc) dserror("Condition does not have Node Ids");

    // loop over DBC nodes
    for (int k = 0; k < (int)(*nodeids_withdbc).size(); ++k)
    {
      int gid = (*nodeids_withdbc)[k];
      // do only nodes which are in my discretisation
      if (StructureField()->Discretization()->NodeRowMap()->MyGID(gid) == false) continue;

      // -------------------- evaluation in special direction, here z-direction
      // get node with global id gid
      DRT::Node* node = StructureField()->Discretization()->gNode(gid);
      if (!node) dserror("Cannot find node with gid %", gid);
      // check coordinates in z-direction, i.e. third value of X()
      double zcoord = node->X()[2];
      // possible push-back
      bool this_is_new_gid = true;
      // get the z-displacement DOFS located at the top surface (z=13.3335mm)
      if (abs(zcoord - top_z) < 1.0e-8)  // change here value for different geometries
      {
        for (unsigned j = 0; j < sdata.size(); j++)
        {
          if (sdata.at(j) == StructureField()->Discretization()->Dof(0, node, 2))
            this_is_new_gid = false;
        }
        if (this_is_new_gid) sdata.push_back(StructureField()->Discretization()->Dof(0, node, 2));
        one_dof_in_dbc.at(0) = StructureField()->Discretization()->Dof(0, node, 2);
      }  // top surface
    }    // loop over DBC nodes
  }      // loop over all STRUCTURAL DBC conditions

  // map containing all z-displacement DOFs which have a DBC
  Teuchos::RCP<Epetra_Map> newdofmap = Teuchos::rcp(new Epetra_Map(
      -1, (int)sdata.size(), &sdata[0], 0, StructureField()->Discretization()->Comm()));

  //---------------------------------------------------------------------------
  // ------------------------------------ initialse STRUCTURAL output variables
  //---------------------------------------------------------------------------

  // ---------------------------------------------------------------- top force
  // nodal reaction force at outer edge for whole support area
  Teuchos::RCP<Epetra_Vector> tension =
      Teuchos::rcp(new Epetra_Vector(*newdofmap,  // map containing
                                                  // all DOFs at top surf with DBC
          false));
  // copy the structural reaction force to tension
  LINALG::Export(*(StructureField()->Freact()), *tension);
  double top_force_local = 0.0;  // local force
  for (int i = 0; i < tension->MyLength(); i++) top_force_local -= (*tension)[i];

  // complete force pointing in axial direction
  double top_force_global = 0.0;

  // sum all nodal forces (top_force_local) in one global vector (top_force_global)
  StructureField()->Discretization()->Comm().SumAll(&top_force_local, &top_force_global, 1);

  // --------------------------------------------- reaction force of whole body
  // due to symmetry only 1/8 is simulated, i.e. only 1/4 of the surface is considered
  // R = force_top * 4
  double top_reaction_force = 4 * top_force_global;

  // --------------------------------------------------------- top displacement
  // top displacement, i.e. displacement at outer edge in axial-direction
  // (here: in z-direction)

  // get the DOFs corresponding to z-displacements and located at the top surface
  std::vector<double> top_disp_local(1);
  top_disp_local.at(0) = 0.0;

  std::vector<int> one_dof_in_dbc_global(1);
  one_dof_in_dbc_global.at(0) = -1;

  StructureField()->Discretization()->Comm().MaxAll(
      &one_dof_in_dbc.at(0), &one_dof_in_dbc_global.at(0), 1);

  // extract axial displacements (here z-displacements) of top surface
  if (StructureField()->Discretization()->DofRowMap()->MyGID(one_dof_in_dbc_global.at(0)))
  {
    DRT::UTILS::ExtractMyValues(
        *(StructureField()->Dispnp()), top_disp_local, one_dof_in_dbc_global);
  }

  // initialse the top displacement
  double top_disp_global = 0.0;
  // sum all nodal displacements (top_disp_local) in one global vector (top_disp_global)
  StructureField()->Discretization()->Comm().SumAll(&top_disp_local.at(0), &top_disp_global, 1);

  // ------------------------------------------------ necking radius at point A
  // necking, i.e. radial displacements in centre plane (here: xy-plane)

  // get the necking DOFs
  // i.e. displacement-DOFs of the xy-plane in the middle of the body
  std::vector<int> necking_radius_dof(1);
  necking_radius_dof.at(0) = -1.;
  for (int k = 0; k < (int)StructureField()->Discretization()->NodeRowMap()->NumMyElements(); k++)
  {
    DRT::Node* node = StructureField()->Discretization()->lRowNode(k);
    // change here value for different geometries
    if ((abs(node->X()[0] - necking_x) < 1.e-8)      // x-direction
        and (abs(node->X()[1] - necking_y) < 1.e-8)  // y-direction
        and (abs(node->X()[2] - necking_z) < 1.e-8)  // z-direction
    )
    {
      // we choose point A (6.413mm / 0mm / -13.3335mm)
      necking_radius_dof.at(0) = StructureField()->Discretization()->Dof(0, node, 0);
      break;  // we only look for one specific node, if we have found it: stop

    }  // end point A(6.413/0/-13.3335)
  }    // sum all nodes
  std::vector<double> necking_radius(1);
  necking_radius.at(0) = 0.0;
  if (necking_radius_dof.at(0) != -1)
  {
    DRT::UTILS::ExtractMyValues(*(StructureField()->Dispnp()), necking_radius, necking_radius_dof);
  }

  // sum necking deformations in the global variable necking_radius_global
  double necking_radius_global = 0.0;
  StructureField()->Discretization()->Comm().SumAll(
      &necking_radius.at(0), &necking_radius_global, 1);

  //---------------------------------------------------------------------------
  // -------------------------------------------------- initialise TEMPERATURES
  //---------------------------------------------------------------------------

  // ------------------------------------------- necking temperature at point A

  // necking temperature at point A, i.e. temperature at the outer side in the middle
  // A (6.413mm / 0.0mm / 13.3335mm)
  std::vector<int> neck_temperature_dof(1);
  neck_temperature_dof.at(0) = -1.0;
  for (int k = 0; k < (int)ThermoField()->Discretization()->NodeRowMap()->NumMyElements(); k++)
  {
    DRT::Node* node = ThermoField()->Discretization()->lRowNode(k);
    // change here value for different geometries
    if ((abs(node->X()[0] - necking_x) < 1.e-8)      // x-direction
        and (abs(node->X()[1] - necking_y) < 1.e-8)  // y-direction
        and (abs(node->X()[2] - necking_z) < 1.e-8)  // z-direction
    )
    {
      neck_temperature_dof.at(0) = ThermoField()->Discretization()->Dof(0, node, 0);
      break;  // we only look for one specific node, if we have found it: stop
    }
  }
  std::vector<double> temperature(1);
  temperature.at(0) = 0.0;
  if (neck_temperature_dof.at(0) != -1)
  {
    DRT::UTILS::ExtractMyValues(*(ThermoField()->Tempnp()),  // global (i)
        temperature,                                         // local (o)
        neck_temperature_dof                                 // global ids to be extracted
    );
  }
  // sum necking temperatures in the variable temperature_global
  double necking_temperature_global = 0.0;
  ThermoField()->Discretization()->Comm().SumAll(
      &temperature.at(0), &necking_temperature_global, 1);

  // -----------------------------------------temperatures at top, i.e. point B

  // necking temperature at point B, i.e. temperature at the outer side at top
  // B (6.413mm / 0.0mm / -13.3335mm)
  std::vector<int> top_temperature_dof(1);
  top_temperature_dof.at(0) = -1.0;
  for (int k = 0; k < (int)ThermoField()->Discretization()->NodeRowMap()->NumMyElements(); k++)
  {
    DRT::Node* node = ThermoField()->Discretization()->lRowNode(k);
    // change here value for different geometries
    if ((abs(node->X()[0] - top_x) < 1.e-8)      // x-direction
        and (abs(node->X()[1] - top_y) < 1.e-8)  // y-direction
        and (abs(node->X()[2] - top_z) < 1.e-8)  // z-direction
    )
    {
      top_temperature_dof.at(0) = ThermoField()->Discretization()->Dof(0, node, 0);
      break;  // we only look for one specific node, if we have found it: stop
    }
  }  // loop over thermal nodes

  // extract top-temperatures of top surface out of global temperature vector
  std::vector<double> top_temperature_local(1);
  top_temperature_local.at(0) = 0.0;
  if (top_temperature_dof.at(0) != -1.)
  {
    DRT::UTILS::ExtractMyValues(*(ThermoField()->Tempnp()),  // global vector (i)
        top_temperature_local,                               // local, i.e. at specific position (o)
        top_temperature_dof                                  // global ids to be extracted
    );
  }

  // sum top-temperatures in the variable top_temperature_global
  double top_temperature_global = 0.0;
  ThermoField()->Discretization()->Comm().SumAll(
      &top_temperature_local.at(0), &top_temperature_global, 1);

  // -------------------------------------------------- print results to screen
  std::cout.precision(7);
  std::cout << std::scientific;
  std::cout << std::fixed;
  if (ThermoField()->Discretization()->Comm().MyPID() == 0)
  {
    std::cout << "OUTPUT:\ttop-disp \ttop-Freact \tneck-disp \tneck-tempi \ttop-tempi \ttop-force\n"
              << "\t" << top_disp_global << "\t" << top_reaction_force << "\t"
              << necking_radius_global << "\t" << (necking_temperature_global - inittemp) << "\t"
              << (top_temperature_global - inittemp) << "\t" << top_force_global << std::endl;
  }

}  // CalculateNeckingTSIResults()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void TSI::Monolithic::PrepareOutput()
{
  // set temperatures on structure field for evaluating stresses
  ApplyThermoCouplingState(ThermoField()->Tempnp());
  // prepare output (i.e. calculate stresses, strains, energies)
  StructureField()->PrepareOutput();

  // reset states
  StructureField()->Discretization()->ClearState(true);
}
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyThermoCouplingState(
    Teuchos::RCP<const Epetra_Vector> temp, Teuchos::RCP<const Epetra_Vector> temp_res)
{
  TSI::Algorithm::ApplyThermoCouplingState(temp, temp_res);

  // set new temperatures to contact
  if (contact_strategy_lagrange_ != Teuchos::null || contact_strategy_nitsche_ != Teuchos::null)
  {
    if (contact_strategy_lagrange_ != Teuchos::null)
      contact_strategy_lagrange_->SetState(
          MORTAR::state_temperature, *coupST_()->SlaveToMaster(ThermoField()->Tempnp()));
    if (contact_strategy_nitsche_ != Teuchos::null)
      contact_strategy_nitsche_->SetState(MORTAR::state_temperature, *ThermoField()->Tempnp());
  }
}  // ApplyThermoCouplingState()


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyStructCouplingState(
    Teuchos::RCP<const Epetra_Vector> disp, Teuchos::RCP<const Epetra_Vector> vel)
{
  if (matchinggrid_)
  {
    if (disp != Teuchos::null) ThermoField()->Discretization()->SetState(1, "displacement", disp);
    if (vel != Teuchos::null) ThermoField()->Discretization()->SetState(1, "velocity", vel);
  }
  else
  {
    if (disp != Teuchos::null)
      ThermoField()->Discretization()->SetState(
          1, "displacement", volcoupl_->ApplyVectorMapping21(disp));
    if (vel != Teuchos::null)
      ThermoField()->Discretization()->SetState(
          1, "velocity", volcoupl_->ApplyVectorMapping21(vel));
  }
}  // ApplyStructCouplingState()



/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::FixTimeIntegrationParams()
{
  if (DRT::INPUT::IntegralValue<INPAR::THR::DynamicType>(
          DRT::Problem::Instance()->ThermalDynamicParams(), "DYNAMICTYP") ==
      INPAR::THR::dyna_genalpha)
  {
    Teuchos::ParameterList& ga = DRT::Problem::Instance()
                                     ->getNonconstParameterList()
                                     ->sublist("THERMAL DYNAMIC")
                                     .sublist("GENALPHA");
    double rhoinf = ga.get<double>("RHO_INF");

    if (rhoinf != -1.)
    {
      if ((rhoinf < 0.0) or (rhoinf > 1.0)) dserror("rho_inf out of range [0.0,1.0]");
      double alpham = 0.5 * (3.0 - rhoinf) / (rhoinf + 1.0);
      double alphaf = 1.0 / (rhoinf + 1.0);
      double gamma = 0.5 + alpham - alphaf;
      ga.set<double>("GAMMA", gamma);
      ga.set<double>("ALPHA_F", alphaf);
      ga.set<double>("ALPHA_M", alpham);
    }
  }

  if (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(
          DRT::Problem::Instance()->StructuralDynamicParams(), "DYNAMICTYP") ==
      INPAR::STR::dyna_genalpha)
  {
    Teuchos::ParameterList& ga = DRT::Problem::Instance()
                                     ->getNonconstParameterList()
                                     ->sublist("STRUCTURAL DYNAMIC")
                                     .sublist("GENALPHA");
    double rhoinf = ga.get<double>("RHO_INF");

    if (rhoinf != -1.)
    {
      double alpham = (2.0 * rhoinf - 1.0) / (rhoinf + 1.0);
      double alphaf = rhoinf / (rhoinf + 1.0);
      double beta = 0.25 * (1.0 - alpham + alphaf) * (1.0 - alpham + alphaf);
      double gamma = 0.5 - alpham + alphaf;

      ga.set<double>("BETA", beta);
      ga.set<double>("GAMMA", gamma);
      ga.set<double>("ALPHA_F", alphaf);
      ga.set<double>("ALPHA_M", alpham);
    }
  }
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void TSI::Monolithic::ApplyDBC()
{
  Teuchos::RCP<LINALG::SparseMatrix> k_ss =
      Teuchos::rcp(new LINALG::SparseMatrix(systemmatrix_->Matrix(0, 0).EpetraMatrix(),
          LINALG::Copy, true, false, LINALG::SparseMatrix::CRS_MATRIX));
  Teuchos::RCP<LINALG::SparseMatrix> k_st =
      Teuchos::rcp(new LINALG::SparseMatrix(systemmatrix_->Matrix(0, 1)));
  Teuchos::RCP<LINALG::SparseMatrix> k_ts =
      Teuchos::rcp(new LINALG::SparseMatrix(systemmatrix_->Matrix(1, 0)));
  Teuchos::RCP<LINALG::SparseMatrix> k_tt =
      Teuchos::rcp(new LINALG::SparseMatrix(systemmatrix_->Matrix(1, 1)));
  if (locsysman_ != Teuchos::null)
  {
    {
      locsysman_->RotateGlobalToLocal(k_ss);
      k_ss->ApplyDirichletWithTrafo(
          locsysman_->Trafo(), *StructureField()->GetDBCMapExtractor()->CondMap(), true);
      locsysman_->RotateLocalToGlobal(k_ss);
    }
    {
      locsysman_->RotateGlobalToLocal(k_st);
      k_st->ApplyDirichletWithTrafo(
          locsysman_->Trafo(), *StructureField()->GetDBCMapExtractor()->CondMap(), false);
      locsysman_->RotateLocalToGlobal(k_st);
    }
  }
  else
  {
    k_ss->ApplyDirichlet(*StructureField()->GetDBCMapExtractor()->CondMap(), true);
    k_st->ApplyDirichlet(*StructureField()->GetDBCMapExtractor()->CondMap(), false);
  }
  k_ts->ApplyDirichlet(*ThermoField()->GetDBCMapExtractor()->CondMap(), false);
  k_tt->ApplyDirichlet(*ThermoField()->GetDBCMapExtractor()->CondMap(), true);


  systemmatrix_->UnComplete();
  systemmatrix_->Assign(0, 0, LINALG::View, *k_ss);
  systemmatrix_->Assign(0, 1, LINALG::View, *k_st);
  systemmatrix_->Assign(1, 0, LINALG::View, *k_ts);
  systemmatrix_->Assign(1, 1, LINALG::View, *k_tt);
  systemmatrix_->Complete();


  if (locsysman_ != Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> s_rhs, t_rhs;
    ExtractFieldVectors(rhs_, s_rhs, t_rhs);
    locsysman_->RotateGlobalToLocal(s_rhs);
    LINALG::ApplyDirichlettoSystem(
        s_rhs, zeros_, *StructureField()->GetDBCMapExtractor()->CondMap());
    locsysman_->RotateLocalToGlobal(s_rhs);

    LINALG::ApplyDirichlettoSystem(t_rhs, zeros_, *ThermoField()->GetDBCMapExtractor()->CondMap());

    Extractor()->InsertVector(*s_rhs, 0, *rhs_);
    Extractor()->InsertVector(*t_rhs, 1, *rhs_);
  }
  else
  {
    LINALG::ApplyDirichlettoSystem(
        rhs_, zeros_, *StructureField()->GetDBCMapExtractor()->CondMap());
    LINALG::ApplyDirichlettoSystem(rhs_, zeros_, *ThermoField()->GetDBCMapExtractor()->CondMap());
  }
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
bool TSI::Monolithic::LSadmissible()
{
  switch (ls_strategy_)
  {
    case INPAR::TSI::LS_structure:
      return normstrrhs_ < last_iter_res_.first;
    case INPAR::TSI::LS_thermo:
      return normthrrhs_ < last_iter_res_.second;
    case INPAR::TSI::LS_or:
      return (normstrrhs_ < last_iter_res_.first || normthrrhs_ < last_iter_res_.second);
    case INPAR::TSI::LS_and:
      return (normstrrhs_ < last_iter_res_.first && normthrrhs_ < last_iter_res_.second);
    case INPAR::TSI::LS_none:
    default:
      dserror("you should not be here");
      return false;
  }
}
