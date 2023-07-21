/*----------------------------------------------------------------------*/
/*! \file

\brief  Basis of all monolithic poroelasticity algorithms

\level 2


*/
#include "baci_poroelast_monolithic.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
// needed for PrintNewton
#include <sstream>

#include "baci_poroelast_defines.H"

// adapters
#include "baci_coupling_adapter_volmortar.H"
#include "baci_adapter_fld_base_algorithm.H"
#include "baci_adapter_fld_poro.H"
#include "baci_adapter_str_fpsiwrapper.H"

// contact
#include "baci_contact_poro_lagrange_strategy.H"
#include "baci_contact_meshtying_poro_lagrange_strategy.H"
#include "baci_contact_meshtying_contact_bridge.H"
#include "baci_contact_nitsche_strategy_poro.H"

#include "baci_fluid_ele_action.H"
#include "baci_fluid_utils_mapextractor.H"

// include this header for coupling stiffness terms
#include "baci_lib_assemblestrategy.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_condition_utils.H"

#include "baci_io_control.H"
#include "baci_inpar_solver.H"

#include "baci_structure_aux.H"

#include "baci_mortar_manager_base.H"

#include "baci_linalg_equilibrate.H"
#include "baci_linalg_utils_sparse_algebra_assemble.H"
#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_linalg_utils_sparse_algebra_manipulation.H"
#include "baci_linear_solver_method_linalg.H"

#include "baci_lib_elements_paramsminimal.H"



POROELAST::Monolithic::Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams,
    Teuchos::RCP<CORE::LINALG::MapExtractor> porosity_splitter)
    : PoroBase(comm, timeparams, porosity_splitter),
      printscreen_(true),   // ADD INPUT PARAMETER
      printiter_(true),     // ADD INPUT PARAMETER
      printerrfile_(true),  // ADD INPUT PARAMETER FOR 'true'
      errfile_(DRT::Problem::Instance()->ErrorFile()->Handle()),
      zeros_(Teuchos::null),
      blockrowdofmap_(Teuchos::null),
      normtypeinc_(INPAR::POROELAST::convnorm_undefined),
      normtypefres_(INPAR::POROELAST::convnorm_undefined),
      combincfres_(INPAR::POROELAST::bop_undefined),
      vectornormfres_(INPAR::POROELAST::norm_undefined),
      vectornorminc_(INPAR::POROELAST::norm_undefined),
      tolinc_(0.0),
      tolfres_(0.0),
      tolinc_struct_(0.0),
      tolfres_struct_(0.0),
      tolinc_velocity_(0.0),
      tolfres_velocity_(0.0),
      tolinc_pressure_(0.0),
      tolfres_pressure_(0.0),
      tolinc_porosity_(0.0),
      tolfres_porosity_(0.0),
      itermax_(0),
      itermin_(0),
      normrhs_(0.0),
      norminc_(0.0),
      normrhsfluidvel_(0.0),
      normincfluidvel_(0.0),
      normrhsfluidpres_(0.0),
      normincfluidpres_(0.0),
      normrhsfluid_(0.0),
      normincfluid_(0.0),
      normrhsstruct_(0.0),
      normincstruct_(0.0),
      normrhsporo_(0.0),
      normincporo_(0.0),
      timer_(Teuchos::rcp(new Teuchos::Time("", false))),
      iter_(-1),
      iterinc_(Teuchos::null),
      directsolve_(true),
      del_(Teuchos::null),
      delhist_(Teuchos::null),
      mu_(0.0),
      equilibration_(Teuchos::null),
      equilibration_method_(CORE::LINALG::EquilibrationMethod::none)
{
  const Teuchos::ParameterList& sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();

  // some solver paramaters are red form the structure dynamic list (this is not the best way to do
  // it ...)
  solveradapttol_ = (DRT::INPUT::IntegralValue<int>(sdynparams, "ADAPTCONV") == 1);
  solveradaptolbetter_ = (sdynparams.get<double>("ADAPTCONV_BETTER"));

  const Teuchos::ParameterList& poroparams = DRT::Problem::Instance()->PoroelastDynamicParams();
  equilibration_method_ =
      Teuchos::getIntegralValue<CORE::LINALG::EquilibrationMethod>(poroparams, "EQUILIBRATION");

  strmethodname_ = DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdynparams, "DYNAMICTYP");
  no_penetration_ = false;
  nit_contact_ = false;
  // if inpar is set to nopenetration for contact!!! to be done!
  // clean up as soon as old time integration is unused!
  if (oldstructimint_)
  {
    if (StructureField()->MeshtyingContactBridge() != Teuchos::null)
    {
      if (StructureField()->MeshtyingContactBridge()->HaveContact())
      {
        auto* poro_lagrange_strat = dynamic_cast<CONTACT::PoroLagrangeStrategy*>(
            &StructureField()->MeshtyingContactBridge()->ContactManager()->GetStrategy());
        if (poro_lagrange_strat)
          no_penetration_ = poro_lagrange_strat->HasPoroNoPenetration();
        else
        {
          auto* co_nitsche_strategy = dynamic_cast<CONTACT::CoNitscheStrategyPoro*>(
              &StructureField()->MeshtyingContactBridge()->ContactManager()->GetStrategy());
          if (co_nitsche_strategy)
          {
            nit_contact_ = true;
            no_penetration_ = co_nitsche_strategy->HasPoroNoPenetration();
          }
        }
      }
    }
  }
  blockrowdofmap_ = Teuchos::rcp(new CORE::LINALG::MultiMapExtractor);

  // contact no penetration constraint not yet works for non-matching structure and fluid
  // discretizations
  if (no_penetration_ and (not matchinggrid_))
  {
    dserror(
        "The contact no penetration constraint does not yet work for non-matching "
        "discretizations!");
  }
}

void POROELAST::Monolithic::DoTimeStep()
{
  // counter and print header
  // predict solution of both field (call the adapter)
  PrepareTimeStep();

  // Newton-Raphson iteration
  Solve();

  // calculate stresses, strains, energies
  constexpr bool force_prepare = false;
  PrepareOutput(force_prepare);

  // update all single field solvers
  Update();

  // write output to screen and files
  Output();
}

void POROELAST::Monolithic::Solve()
{
  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic tangent matrix

  //---------------------------------------- initialise equilibrium loop and norms
  SetupNewton();

  //---------------------------------------------- iteration loop
  // equilibrium iteration loop (loop over k)
  while (((not Converged()) and (iter_ <= itermax_)) or (iter_ <= itermin_))
  {
    timer_->start();
    Teuchos::Time timer("eval", true);
    timer.start();
    // compute residual forces #rhs_ and tangent #tang_
    // whose components are globally oriented
    // build linear system stiffness matrix and rhs/force residual for each
    // field, here e.g. for structure field: field want the iteration increment
    // 1.) Update(iterinc_),
    // 2.) EvaluateForceStiffResidual(),
    // 3.) PrepareSystemForNewtonSolve()
    Evaluate(iterinc_, iter_ == 1);
    // std::cout << "  time for Evaluate : " <<  timer.totalElapsedTime(true) << "\n";
    // timer.reset();

    // if (iter_>1 and Step()>2 )
    // PoroFDCheck();

    // build norms
    BuildConvergenceNorms();
    if ((not Converged()) or combincfres_ == INPAR::POROELAST::bop_or)
    {
      // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
      // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
      LinearSolve();
      // std::cout << "  time for Evaluate LinearSolve: " << timer.totalElapsedTime(true) << "\n";
      // timer.reset();

      // reset solver tolerance
      solver_->ResetTolerance();

      // rebuild norms
      BuildConvergenceNorms();
    }

    // print stuff
    PrintNewtonIter();

    // Aitken();

    // Recover Lagrangean Multiplier in Newton Iteration (for contact & contact no penetration!)
    // adjust LM Recovery for the meshtying case
    RecoverLagrangeMultiplierAfterNewtonStep(iterinc_);

    // increment equilibrium loop index
    iter_ += 1;
  }  // end equilibrium loop

  //---------------------------------------------- output of number of iterations
  //  {
  //    std::ostringstream s;
  //    s << std::right << std::setw(16) << std::scientific << Time()
  //      << std::right << std::setw(10) << std::scientific << Step()
  //      << std::right << std::setw(10) << std::scientific << iter_-1;
  //
  //    std::ofstream f;
  //    const std::string fname = DRT::Problem::Instance()->OutputControlFile()->FileName()
  //                            + "_numiter.txt";
  //
  //    if (Step() <= 1)
  //      f.open(fname.c_str(),std::fstream::trunc); //f << header.str() << std::endl;
  //    else
  //      f.open(fname.c_str(),std::fstream::ate | std::fstream::app);
  //
  //    f << s.str() << "\n";
  //    f.close();
  //  }

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ((Converged()) and (Comm().MyPID() == 0))
  {
    PrintNewtonConv();
  }
  else if (iter_ >= itermax_)
  {
    dserror("Newton unconverged in %d iterations", iter_);
  }
}

void POROELAST::Monolithic::UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::UpdateStateIncrementally");

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> sx = Teuchos::null;
  Teuchos::RCP<const Epetra_Vector> fx = Teuchos::null;

  // do nothing in the case no increment vector is given
  if (x == Teuchos::null)
  {
    return;
  }
  else
  {
    // extract displacement sx and fluid fx incremental vector of global
    // unknown incremental vector x (different for splits)
    ExtractFieldVectors(x, sx, fx);

    // update poro iterinc
    if (is_part_of_multifield_problem_) UpdatePoroIterinc(x);
  }

  UpdateStateIncrementally(sx, fx);
}

void POROELAST::Monolithic::UpdateStateIncrementally(
    Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx)
{
  // Newton update of the fluid field
  // update velocities and pressures before passed to the structural field
  //  UpdateIterIncrementally(fx),
  FluidField()->UpdateNewton(fx);

  // call all elements and assemble rhs and matrices
  // structural field

  // structure Evaluate (builds tangent, residual and applies DBC)

  // apply current velocity and pressures to structure
  SetFluidSolution();

  // apply current velocity of fluid to ContactMangager if contact problem
  if (no_penetration_)
    SetPoroContactStates();  // ATM svel is set in structure evaluate as the vel of the structure is
                             // evaluated there ...

  // Monolithic Poroelasticity accesses the linearised structure problem:
  //   UpdaterIterIncrementally(sx),
  StructureField()->UpdateStateIncrementally(sx);
  // fluid field

  // fluid Evaluate
  // (builds tangent, residual and applies DBC and recent coupling values)

  // set structure displacements onto the fluid
  SetStructSolution();
}

void POROELAST::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x, bool firstiter)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::Evaluate");

  // Evaluate Stiffness and RHS of Structural & Fluid field
  EvaluateFields(x);

  // fill off diagonal blocks and fill global system matrix
  SetupSystemMatrix();

  // check whether we have a sanely filled tangent matrix
  if (not systemmatrix_->Filled())
  {
    dserror("Effective tangent matrix must be filled here");
  }

  // create full monolithic rhs vector
  SetupRHS(firstiter);

  // Modify System for Contact or Meshtying!
  EvalPoroMortar();
}

void POROELAST::Monolithic::Evaluate(
    Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx, bool firstiter)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::Evaluate");

  // Evaluate Stiffness and RHS of Structural & Fluid field
  EvaluateFields(sx, fx);

  // fill off diagonal blocks and fill global system matrix
  SetupSystemMatrix();

  // check whether we have a sanely filled tangent matrix
  if (not systemmatrix_->Filled())
  {
    dserror("Effective tangent matrix must be filled here");
  }

  // create full monolithic rhs vector
  SetupRHS(firstiter);

  // Modify System for Contact or Meshtying!
  EvalPoroMortar();
}

void POROELAST::Monolithic::EvaluateFields(
    Teuchos::RCP<const Epetra_Vector> sx, Teuchos::RCP<const Epetra_Vector> fx)
{
  // Update State incrementally
  UpdateStateIncrementally(sx, fx);

  // Monolithic Poroelasticity accesses the linearised structure problem:
  //   EvaluateForceStiffResidual()
  //   PrepareSystemForNewtonSolve()
  StructureField()->Evaluate(Teuchos::null);

  // monolithic Poroelasticity accesses the linearised fluid problem
  //   EvaluateRhsTangResidual() and
  //   PrepareSystemForNewtonSolve()
  FluidField()->Evaluate(Teuchos::null);
}

void POROELAST::Monolithic::EvaluateFields(Teuchos::RCP<const Epetra_Vector> x)
{
  // Update State incrementally
  UpdateStateIncrementally(x);

  // Monolithic Poroelasticity accesses the linearised structure problem:
  //   EvaluateForceStiffResidual()
  //   PrepareSystemForNewtonSolve()
  StructureField()->Evaluate(Teuchos::null);

  // monolithic Poroelasticity accesses the linearised fluid problem
  //   EvaluateRhsTangResidual() and
  //   PrepareSystemForNewtonSolve()
  FluidField()->Evaluate(Teuchos::null);
}

void POROELAST::Monolithic::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::ExtractFieldVectors");

  // process structure unknowns of the first field
  sx = Extractor()->ExtractVector(x, 0);

  // process fluid unknowns of the second field
  fx = Extractor()->ExtractVector(x, 1);
}

void POROELAST::Monolithic::SetupSystem()
{
  {
    // -------------------------------------------------------------create combined map
    std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

    // Note:
    // when using constraints applied via Lagrange-Multipliers there is a
    // difference between StructureField()->DofRowMap() and StructureField()->DofRowMap(0).
    // StructureField()->DofRowMap(0) returns the DofRowMap
    // known to the discretization (without lagrange multipliers)
    // while StructureField()->DofRowMap() returns the DofRowMap known to
    // the constraint manager (with lagrange multipliers)
    // In the constrained case we want the "whole" RowDofMap,
    // otherwise both calls are equivalent

    vecSpaces.push_back(StructureField()->DofRowMap());
    // use its own DofRowMap, that is the 0th map of the discretization
    vecSpaces.push_back(FluidField()->DofRowMap(0));

    if (vecSpaces[0]->NumGlobalElements() == 0) dserror("No structure equation. Panic.");
    if (vecSpaces[1]->NumGlobalElements() == 0) dserror("No fluid equation. Panic.");

    // full Poroelasticity-map
    fullmap_ = CORE::LINALG::MultiMapExtractor::MergeMaps(vecSpaces);
    // full Poroelasticity-blockmap
    blockrowdofmap_->Setup(*fullmap_, vecSpaces);
  }
  // -------------------------------------------------------------

  //-----------------------------------build map of global dofs with DBC
  BuildCombinedDBCMap();
  // -------------------------------------------------------------

  // initialize Poroelasticity-systemmatrix_
  systemmatrix_ =
      Teuchos::rcp(new CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          *Extractor(), *Extractor(), 81, false, true));

  k_sf_ = Teuchos::rcp(
      new CORE::LINALG::SparseMatrix(*(StructureField()->DofRowMap()), 81, true, true));
  k_fs_ =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*(FluidField()->Discretization()->DofRowMap(0)),
          //*(FluidField()->DofRowMap()),
          81, true, true));

  nopen_handle_->Setup(DofRowMap(), (FluidField()->Discretization()->DofRowMap(0)));

  SetupEquilibration();
}  // SetupSystem()

void POROELAST::Monolithic::SetupEquilibration()
{
  // instantiate appropriate equilibration class
  auto equilibration_method =
      std::vector<CORE::LINALG::EquilibrationMethod>(1, equilibration_method_);
  equilibration_ = CORE::LINALG::BuildEquilibration(
      CORE::LINALG::MatrixType::block_field, equilibration_method, fullmap_);
}

void POROELAST::Monolithic::SetupSystemMatrix(CORE::LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::SetupSystemMatrix");

  // pure structural part k_ss (3nx3n)

  // build pure structural block k_ss
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<CORE::LINALG::SparseMatrix> k_ss = StructureField()->SystemMatrix();

  if (k_ss == Teuchos::null) dserror("structure system matrix null pointer!");

  /*----------------------------------------------------------------------*/
  // structural part k_sf (3nxn)
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<CORE::LINALG::SparseMatrix> k_sf = StructFluidCouplingMatrix();

  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_sf);

  if (!matchinggrid_)
  {
    k_sf->Complete(*(StructureField()->Discretization()->DofRowMap(1)),
        *(StructureField()->Discretization()->DofRowMap(0)));

    k_sf = volcoupl_->ApplyMatrixMapping12(k_sf);
  }

  /*----------------------------------------------------------------------*/
  // pure fluid part k_ff ( (3n+1)x(3n+1) )

  // build pure fluid block k_ff
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<CORE::LINALG::SparseMatrix> k_ff = FluidField()->SystemMatrix();

  if (k_ff == Teuchos::null) dserror("fuid system matrix null pointer!");

  if (nopen_handle_->HasCond())
  {
    // Evaluate poroelasticity specific conditions
    EvaluateCondition(k_ff);
  }

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (3n+1)x3n )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<CORE::LINALG::SparseMatrix> k_fs = FluidStructCouplingMatrix();

  // call the element and calculate the matrix block
  ApplyFluidCouplMatrix(k_fs);

  if (!matchinggrid_)
  {
    k_fs->Complete(*(FluidField()->Discretization()->DofRowMap(1)),
        *(FluidField()->Discretization()->DofRowMap(0)));

    k_fs = volcoupl_->ApplyMatrixMapping21(k_fs);
  }

  // assign structure part to the Poroelasticity matrix
  mat.Assign(0, 0, CORE::LINALG::View, *k_ss);
  // assign coupling part to the Poroelasticity matrix
  mat.Assign(0, 1, CORE::LINALG::View, *k_sf);
  // assign fluid part to the poroelasticity matrix
  mat.Assign(1, 1, CORE::LINALG::View, *k_ff);
  // assign coupling part to the Poroelasticity matrix
  mat.Assign(1, 0, CORE::LINALG::View, *k_fs);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();
}

void POROELAST::Monolithic::SetupRHS(bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::SetupRHS");

  // create full monolithic rhs vector
  if (rhs_ == Teuchos::null) rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  // fill the Poroelasticity rhs vector rhs_ with the single field rhss
  SetupVector(*rhs_, StructureField()->RHS(), FluidField()->RHS());

  // add rhs terms due to no penetration condition
  nopen_handle_->ApplyCondRHS(iterinc_, rhs_);
}

void POROELAST::Monolithic::PrepareTimeStep() { PoroBase::PrepareTimeStep(); }

void POROELAST::Monolithic::LinearSolve()
{
  if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normrhs_;
    double wanted = tolfres_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  // equilibrate global system of equations if necessary
  equilibration_->EquilibrateSystem(systemmatrix_, rhs_, blockrowdofmap_);

  if (directsolve_)
  {
    // merge blockmatrix to SparseMatrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

    // apply dirichlet boundary conditions
    CORE::LINALG::ApplyDirichletToSystem(*sparse, *iterinc_, *rhs_, *zeros_, *CombinedDBCMap());

    // standard solver call
    solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
  }
  else  // use bgs2x2_operator
  {
    // apply dirichlet boundary conditions
    CORE::LINALG::ApplyDirichletToSystem(
        *systemmatrix_, *iterinc_, *rhs_, *zeros_, *CombinedDBCMap());

    // standard solver call
    solver_->Solve(systemmatrix_->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
  }

  equilibration_->UnequilibrateIncrement(iterinc_);
}

void POROELAST::Monolithic::CreateLinearSolver()
{
  // get dynamic section
  const Teuchos::ParameterList& porodyn = DRT::Problem::Instance()->PoroelastDynamicParams();

  // get the linear solver number
  const int linsolvernumber = porodyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
  {
    dserror(
        "no linear solver defined for monolithic Poroelasticity. Please set LINEAR_SOLVER in "
        "POROELASTICITY DYNAMIC to a valid number!");
  }

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (slinsolvernumber == (-1))
  {
    dserror(
        "no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL "
        "DYNAMIC to a valid number!");
  }

  // get parameter list of fluid dynamics
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  // use solver blocks for fluid
  // get the solver number used for fluid solver
  const int flinsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  // check if the fluid solver has a valid solver number
  if (flinsolvernumber == (-1))
  {
    dserror(
        "no linear solver defined for fluid field. Please set LINEAR_SOLVER in FLUID DYNAMIC to a "
        "valid number!");
  }

  // get solver parameter list of linear Poroelasticity solver
  const Teuchos::ParameterList& porosolverparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);

  const auto solvertype =
      Teuchos::getIntegralValue<INPAR::SOLVER::SolverType>(porosolverparams, "SOLVER");

  if (solvertype != INPAR::SOLVER::SolverType::belos)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " Note: the BGS2x2 preconditioner now " << std::endl;
    std::cout << " uses the structural solver and fluid solver blocks" << std::endl;
    std::cout << " for building the internal inverses" << std::endl;
    std::cout << " Remove the old BGS PRECONDITIONER BLOCK entries " << std::endl;
    std::cout << " in the dat files!" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    dserror("Iterative solver expected");
  }
  const auto azprectype =
      Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(porosolverparams, "AZPREC");

  // plausibility check
  switch (azprectype)
  {
    case INPAR::SOLVER::PreconditionerType::block_gauss_seidel_2x2:
      break;
    case INPAR::SOLVER::PreconditionerType::multigrid_nxn:
    {
      // no plausibility checks here
      // if you forget to declare an xml file you will get an error message anyway
    }
    break;
    default:
      dserror("Block Gauss-Seidel BGS2x2 preconditioner expected");
      break;
  }

  solver_ = Teuchos::rcp(new CORE::LINALG::Solver(
      porosolverparams, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

  // use solver blocks for structure and fluid
  const Teuchos::ParameterList& ssolverparams =
      DRT::Problem::Instance()->SolverParams(slinsolvernumber);
  const Teuchos::ParameterList& fsolverparams =
      DRT::Problem::Instance()->SolverParams(flinsolvernumber);

  solver_->PutSolverParamsToSubParams("Inverse1", ssolverparams);
  solver_->PutSolverParamsToSubParams("Inverse2", fsolverparams);

  // prescribe rigid body modes
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse1"));
  FluidField()->Discretization()->ComputeNullSpaceIfNecessary(
      solver_->Params().sublist("Inverse2"));
}

void POROELAST::Monolithic::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::InitialGuess");

  // InitalGuess() is called of the single fields and results are put in
  // increment vector ig
  SetupVector(*ig,
      // returns residual displacements \f$\Delta D_{n+1}^{<k>}\f$ - disi_
      StructureField()->InitialGuess(),
      // returns residual velocities or iterative fluid increment - incvel_
      FluidField()->InitialGuess());
}

void POROELAST::Monolithic::SetupVector(
    Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv)
{
  // extract dofs of the two fields
  // and put the structural/fluid field vector into the global vector f
  // noticing the block number

  Extractor()->InsertVector(*sv, 0, f);
  if (not oldstructimint_) f.Scale(-1);
  Extractor()->InsertVector(*fv, 1, f);
}

bool POROELAST::Monolithic::Converged()
{
  // check for single norms
  bool convinc = false;
  bool convfres = false;

  // residual increments
  switch (normtypeinc_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      convinc = norminc_ < tolinc_;
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      convinc = (normincstruct_ < tolinc_struct_ and normincfluidvel_ < tolinc_velocity_ and
                 normincfluidpres_ < tolinc_pressure_ and normincporo_ < tolinc_porosity_);
      break;
    default:
      dserror("Cannot check for convergence of residual values!");
      break;
  }

  // residual forces
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      convfres = normrhs_ < tolfres_;
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      convfres = (normrhsstruct_ < tolfres_struct_ and normrhsfluidvel_ < tolfres_velocity_ and
                  normrhsfluidpres_ < tolfres_pressure_ and normrhsporo_ < tolfres_porosity_);
      break;
    default:
      dserror("Cannot check for convergence of residual forces!");
      break;
  }

  // combine increments and forces
  bool conv = false;
  if (combincfres_ == INPAR::POROELAST::bop_and)
    conv = convinc and convfres;
  else if (combincfres_ == INPAR::POROELAST::bop_or)
    conv = convinc or convfres;
  else
    dserror("Something went terribly wrong with binary operator!");

  if (conv)
  {
    // clean up as soon as old time integration is unused!
    if (oldstructimint_)
    {
      if (StructureField()->MeshtyingContactBridge() != Teuchos::null)
      // assume dual mortar lagmult contact!!! ... for general case store member into algo!
      {
        if (StructureField()->MeshtyingContactBridge()->HaveContact())
        {
          conv = StructureField()
                     ->MeshtyingContactBridge()
                     ->GetStrategy()
                     .ActiveSetSemiSmoothConverged();
        }
      }
    }
  }

  // return things
  return conv;
}

void POROELAST::Monolithic::PrintNewtonIter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ((Comm().MyPID() == 0) and printscreen_ and (Step() % printscreen_ == 0) and printiter_)
  {
    if (iter_ == 1) PrintNewtonIterHeader(stdout);
    PrintNewtonIterText(stdout);
  }

  // print to error file
  if (printerrfile_ and printiter_)
  {
    if (iter_ == 1) PrintNewtonIterHeader(errfile_);
    PrintNewtonIterText(errfile_);
  }
}

void POROELAST::Monolithic::PrintNewtonIterHeader(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  PrintNewtonIterHeaderStream(oss);

  // add solution time
  oss << std::setw(14) << "wct";
  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == nullptr) dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);
}

void POROELAST::Monolithic::PrintNewtonIterHeaderStream(std::ostringstream& oss)
{
  oss << "------------------------------------------------------------" << std::endl;
  oss << "                   Newton-Raphson Scheme                    " << std::endl;
  oss << "                NormRES " << VectorNormString(vectornormfres_);
  oss << "     NormINC " << VectorNormString(vectornorminc_) << "                    " << std::endl;
  oss << "------------------------------------------------------------" << std::endl;

  // enter converged state etc
  oss << "numiter";

  // different style due relative or absolute error checking

  // --------------------------------------------------------global system test
  // residual forces
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      oss << std::setw(15) << "abs-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_ << ")";
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      oss << std::setw(15) << "abs-s-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_struct_ << ")";
      if (porosity_dof_)
        oss << std::setw(15) << "abs-poro-res"
            << "(" << std::setw(5) << std::setprecision(2) << tolfres_porosity_ << ")";
      oss << std::setw(15) << "abs-fvel-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_velocity_ << ")";
      oss << std::setw(15) << "abs-fpres-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_pressure_ << ")";
      break;
    default:
      dserror("Unknown or undefined convergence form for residual.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      oss << std::setw(15) << "abs-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_ << ")";
      break;
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      oss << std::setw(15) << "abs-s-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_struct_ << ")";
      if (porosity_dof_)
        oss << std::setw(15) << "abs-poro-inc"
            << "(" << std::setw(5) << std::setprecision(2) << tolinc_porosity_ << ")";
      oss << std::setw(15) << "abs-fvel-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_velocity_ << ")";
      oss << std::setw(15) << "abs-fpres-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_pressure_ << ")";
      break;
    default:
      dserror("Unknown or undefined convergence form for increment.");
      break;
  }
}

void POROELAST::Monolithic::PrintNewtonIterText(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  PrintNewtonIterTextStream(oss);

  // add solution time
  oss << std::setw(14) << std::setprecision(2) << std::scientific << timer_->totalElapsedTime(true);
  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == nullptr) dserror("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);
}

void POROELAST::Monolithic::PrintNewtonIterTextStream(std::ostringstream& oss)
{
  // enter converged state etc
  oss << std::setw(7) << iter_;

  // different style due relative or absolute error checking

  // --------------------------------------------------------global system test
  // residual forces
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhs_;
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      break;
    default:
      dserror("Unknown or undefined convergence form for global residual.");
      break;
  }
  // increments
  switch (normtypeinc_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      oss << std::setw(22) << std::setprecision(5) << std::scientific << norminc_;
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      break;
    default:
      dserror("Unknown or undefined convergence form for global increment.");
      break;
  }

  // --------------------------------------------------------single field test
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs_singlefields:
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsstruct_;
      if (porosity_dof_)
        oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsporo_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsfluidvel_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsfluidpres_;
      break;
    case INPAR::POROELAST::convnorm_abs_global:
      break;
    default:
      dserror("Unknown or undefined convergence form for single field residual.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::POROELAST::convnorm_abs_singlefields:
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincstruct_;
      if (porosity_dof_)
        oss << std::setw(22) << std::setprecision(5) << std::scientific << normincporo_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincfluidvel_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincfluidpres_;
      break;
    case INPAR::POROELAST::convnorm_abs_global:
      break;
    default:
      dserror("Unknown or undefined convergence form for single field increment.");
      break;
  }
}

void POROELAST::Monolithic::PrintNewtonConv() {}

void POROELAST::Monolithic::ApplyStrCouplMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> k_sf)
{
  k_sf->Zero();

  // create the parameters for the discretization
  Teuchos::ParameterList sparams;

  if (oldstructimint_)
  {
    const std::string action = "struct_poro_calc_fluidcoupling";
    sparams.set("action", action);
    // other parameters that might be needed by the elements
    sparams.set("delta time", Dt());
    sparams.set("total time", Time());
  }
  else
  {
    //! pointer to the model evaluator data container
    Teuchos::RCP<DRT::ELEMENTS::ParamsMinimal> params =
        Teuchos::rcp(new DRT::ELEMENTS::ParamsMinimal());

    // set parameters needed for element evalutation
    params->SetActionType(DRT::ELEMENTS::struct_poro_calc_fluidcoupling);
    params->SetTotalTime(Time());
    params->SetDeltaTime(Dt());

    sparams.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", params);
  }

  StructureField()->Discretization()->ClearState();
  StructureField()->Discretization()->SetState(0, "displacement", StructureField()->Dispnp());
  StructureField()->Discretization()->SetState(0, "velocity", StructureField()->Velnp());

  // build specific assemble strategy for mechanical-fluid system matrix
  // from the point of view of StructureField:
  // structdofset = 0, fluiddofset = 1
  DRT::AssembleStrategy structuralstrategy(0,  // structdofset for row
      1,                                       // fluiddofset for column
      k_sf,                                    // mechanical-fluid coupling matrix
      Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate the mechanical-fluid system matrix on the structural element
  StructureField()->Discretization()->EvaluateCondition(
      sparams, structuralstrategy, "PoroCoupling");
  // StructureField()->Discretization()->Evaluate( sparams, structuralstrategy);
  StructureField()->Discretization()->ClearState();

  // scale with time integration factor
  k_sf->Scale(1.0 - StructureField()->TimIntParam());
}

void POROELAST::Monolithic::ApplyFluidCouplMatrix(Teuchos::RCP<CORE::LINALG::SparseOperator> k_fs)
{
  k_fs->Zero();

  // create the parameters for the discretization
  Teuchos::ParameterList fparams;
  // action for elements
  fparams.set<int>("action", FLD::calc_porousflow_fluid_coupling);
  // physical type
  fparams.set<int>("Physical Type", FluidField()->PhysicalType());
  // other parameters that might be needed by the elements
  fparams.set("delta time", Dt());
  fparams.set("total time", Time());

  FluidField()->Discretization()->ClearState();

  // set general vector values needed by elements
  FluidField()->Discretization()->SetState(0, "dispnp", FluidField()->Dispnp());
  FluidField()->Discretization()->SetState(0, "dispn", FluidField()->Dispn());
  FluidField()->Discretization()->SetState(0, "gridv", FluidField()->GridVel());
  FluidField()->Discretization()->SetState(0, "gridvn", FluidField()->GridVeln());
  FluidField()->Discretization()->SetState(0, "veln", FluidField()->Veln());
  FluidField()->Discretization()->SetState(0, "accnp", FluidField()->Accnp());
  FluidField()->Discretization()->SetState(0, "accam", FluidField()->Accam());
  FluidField()->Discretization()->SetState(0, "accn", FluidField()->Accn());

  FluidField()->Discretization()->SetState(0, "scaaf", FluidField()->Scaaf());

  FluidField()->Discretization()->SetState(0, "hist", FluidField()->Hist());

  // set scheme-specific element parameters and vector values
  if (FluidField()->TimIntScheme() == INPAR::FLUID::timeint_npgenalpha or
      FluidField()->TimIntScheme() == INPAR::FLUID::timeint_npgenalpha)
    FluidField()->Discretization()->SetState(0, "velaf", FluidField()->Velaf());
  else
    FluidField()->Discretization()->SetState(0, "velaf", FluidField()->Velnp());

  FluidField()->Discretization()->SetState(0, "velnp", FluidField()->Velnp());

  // build specific assemble strategy for the fluid-mechanical system matrix
  // from the point of view of FluidField:
  // fluiddofset = 0, structdofset = 1
  DRT::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
      1,                                  // structdofset for column
      k_fs,                               // fluid-mechanical matrix
      Teuchos::null,                      // no other matrix or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate the fluid-mechanical system matrix on the fluid element
  FluidField()->Discretization()->EvaluateCondition(fparams, fluidstrategy, "PoroCoupling");
  // FluidField()->Discretization()->Evaluate( fparams, fluidstrategy );

  // evaluate coupling terms from partial integration of continuity equation
  if (part_int_cond_)
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_boundary);
    params.set("total time", Time());
    params.set("delta time", Dt());
    params.set<POROELAST::coupltype>("coupling", POROELAST::fluidstructure);
    params.set("timescale", FluidField()->ResidualScaling());
    params.set<int>("Physical Type", FluidField()->PhysicalType());

    FluidField()->Discretization()->ClearState();
    FluidField()->Discretization()->SetState(0, "dispnp", FluidField()->Dispnp());
    FluidField()->Discretization()->SetState(0, "gridv", FluidField()->GridVel());
    FluidField()->Discretization()->SetState(0, "velnp", FluidField()->Velnp());
    FluidField()->Discretization()->SetState(0, "scaaf", FluidField()->Scaaf());

    FluidField()->Discretization()->EvaluateCondition(params, fluidstrategy, "PoroPartInt");

    FluidField()->Discretization()->ClearState();
  }
  if (pres_int_cond_)
  {
    // create the parameters for the discretization
    Teuchos::ParameterList params;
    // action for elements
    params.set<int>("action", FLD::poro_prescoupl);
    params.set<POROELAST::coupltype>("coupling", POROELAST::fluidstructure);
    params.set<int>("Physical Type", FluidField()->PhysicalType());

    FluidField()->Discretization()->ClearState();
    FluidField()->Discretization()->SetState(0, "dispnp", FluidField()->Dispnp());
    FluidField()->Discretization()->SetState(0, "velnp", FluidField()->Velnp());

    FluidField()->Discretization()->EvaluateCondition(params, fluidstrategy, "PoroPresInt");

    FluidField()->Discretization()->ClearState();
  }

  // apply normal flux condition on coupling part
  if (nopen_handle_->HasCond())
  {
    k_fs->Complete(
        StructureField()->SystemMatrix()->RangeMap(), FluidField()->SystemMatrix()->RangeMap());
    EvaluateCondition(k_fs, POROELAST::fluidstructure);
  }

  FluidField()->Discretization()->ClearState();
}

[[maybe_unused]] void POROELAST::Monolithic::PoroFDCheck()
{
  std::cout << "\n******************finite difference check***************" << std::endl;

  int dof_struct = (DofRowMapStructure()->NumGlobalElements());
  int dof_fluid = (DofRowMapFluid()->NumGlobalElements());

  std::cout << "structure field has " << dof_struct << " DOFs" << std::endl;
  std::cout << "fluid field has " << dof_fluid << " DOFs" << std::endl;

  Teuchos::RCP<Epetra_Vector> iterinc = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> abs_iterinc = Teuchos::null;
  iterinc = CORE::LINALG::CreateVector(*DofRowMap(), true);
  abs_iterinc = CORE::LINALG::CreateVector(*DofRowMap(), true);

  const int dofs = iterinc->GlobalLength();
  std::cout << "in total " << dofs << " DOFs" << std::endl;
  const double delta = 1e-8;

  iterinc->PutScalar(0.0);

  iterinc->ReplaceGlobalValue(0, 0, delta);

  abs_iterinc->Update(1.0, *iterinc_, 0.0);

  Teuchos::RCP<Epetra_CrsMatrix> stiff_approx = Teuchos::null;
  stiff_approx = CORE::LINALG::CreateMatrix(*DofRowMap(), 81);

  Teuchos::RCP<Epetra_Vector> rhs_old = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));
  rhs_old->Update(1.0, *rhs_, 0.0);
  Teuchos::RCP<Epetra_Vector> rhs_copy = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  Teuchos::RCP<CORE::LINALG::SparseMatrix> sparse = systemmatrix_->Merge();
  Teuchos::RCP<CORE::LINALG::SparseMatrix> sparse_copy =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(sparse->EpetraMatrix(), CORE::LINALG::Copy));

  bool output = false;
  if (output)
  {
    std::cout << "iterinc_" << std::endl << *iterinc_ << std::endl;
    std::cout << "iterinc" << std::endl << *iterinc << std::endl;
    std::cout << "meshdisp: " << std::endl << *(FluidField()->Dispnp());
    std::cout << "disp: " << std::endl << *(StructureField()->Dispnp());
    std::cout << "fluid vel" << std::endl << *(FluidField()->Velnp());
    std::cout << "fluid acc" << std::endl << *(FluidField()->Accnp());
    std::cout << "gridvel fluid" << std::endl << *(FluidField()->GridVel());
    std::cout << "gridvel struct" << std::endl << *(StructureField()->Velnp());
  }

  const int row_number = -1;
  const int column_number = -1;
  for (int i = 0; i < dofs; ++i)
  {
    if (CombinedDBCMap()->MyGID(i))
    {
      iterinc->ReplaceGlobalValue(i, 0, 0.0);
    }
    abs_iterinc->Update(1.0, *iterinc, 1.0);

    if (i == column_number)
      std::cout << "\n******************" << column_number + 1 << ". Spalte!!***************"
                << std::endl;

    Evaluate(iterinc, iter_ == 1);

    rhs_copy->Update(1.0, *rhs_, 0.0);

    iterinc_->PutScalar(0.0);  // Useful? depends on solver and more
    CORE::LINALG::ApplyDirichletToSystem(
        *sparse_copy, *iterinc_, *rhs_copy, *zeros_, *CombinedDBCMap());


    if (i == column_number)
    {
      std::cout << "rhs_: " << (*rhs_copy)[row_number] << std::endl;
      std::cout << "rhs_old: " << (*rhs_old)[row_number] << std::endl;
    }

    rhs_copy->Update(-1.0, *rhs_old, 1.0);
    rhs_copy->Scale(-1.0 / delta);

    int* index = &i;
    for (int j = 0; j < dofs; ++j)
    {
      double value = (*rhs_copy)[j];
      stiff_approx->InsertGlobalValues(j, 1, &value, index);

      if ((j == row_number) and (i == column_number))
      {
        std::cout << "\n******************" << row_number + 1 << ". Row!!***************"
                  << std::endl;
        std::cout << "iterinc_" << std::endl << *iterinc_ << std::endl;
        std::cout << "iterinc" << std::endl << *iterinc << std::endl;
        std::cout << "meshdisp: " << std::endl << *(FluidField()->Dispnp());
        std::cout << "disp: " << std::endl << *(StructureField()->Dispnp());
        std::cout << "fluid vel" << std::endl << *(FluidField()->Velnp());
        std::cout << "fluid acc" << std::endl << *(FluidField()->Accnp());
        std::cout << "gridvel fluid" << std::endl << *(FluidField()->GridVel());
        std::cout << "gridvel struct" << std::endl << *(StructureField()->Velnp());

        std::cout << "stiff_apprx(" << row_number << "," << column_number
                  << "): " << (*rhs_copy)[row_number] << std::endl;

        std::cout << "value(" << row_number << "," << column_number << "): " << value << std::endl;
        std::cout << "\n******************" << row_number + 1 << ". Row End!!***************"
                  << std::endl;
      }
    }

    if (not CombinedDBCMap()->MyGID(i)) iterinc->ReplaceGlobalValue(i, 0, -delta);

    iterinc->ReplaceGlobalValue(i - 1, 0, 0.0);

    if (i != dofs - 1) iterinc->ReplaceGlobalValue(i + 1, 0, delta);

    if (i == column_number)
      std::cout << "\n******************" << column_number + 1 << ". Column End!!***************"
                << std::endl;
  }

  Evaluate(iterinc, iter_ == 1);

  stiff_approx->FillComplete();

  Teuchos::RCP<CORE::LINALG::SparseMatrix> stiff_approx_sparse = Teuchos::null;
  stiff_approx_sparse =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(stiff_approx, CORE::LINALG::Copy));

  stiff_approx_sparse->Add(*sparse_copy, false, -1.0, 1.0);

  Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = sparse_copy->EpetraMatrix();

  Teuchos::RCP<Epetra_CrsMatrix> error_crs = stiff_approx_sparse->EpetraMatrix();

  error_crs->FillComplete();
  sparse_crs->FillComplete();

  bool success = true;
  double error_max = 0.0;
  double abs_error_max = 0.0;
  for (int i = 0; i < dofs; ++i)
  {
    if (not CombinedDBCMap()->MyGID(i))
    {
      for (int j = 0; j < dofs; ++j)
      {
        if (not CombinedDBCMap()->MyGID(j))
        {
          double stiff_approx_ij = 0.0;
          double sparse_ij = 0.0;
          double error_ij = 0.0;

          {
            // get error_crs entry ij
            int errornumentries;
            int errorlength = error_crs->NumGlobalEntries(i);
            std::vector<double> errorvalues(errorlength);
            std::vector<int> errorindices(errorlength);
            // int errorextractionstatus =
            error_crs->ExtractGlobalRowCopy(
                i, errorlength, errornumentries, errorvalues.data(), errorindices.data());
            for (int k = 0; k < errorlength; ++k)
            {
              if (errorindices[k] == j)
              {
                error_ij = errorvalues[k];
                break;
              }
              else
                error_ij = 0.0;
            }
          }

          // get sparse_ij entry ij
          {
            int sparsenumentries;
            int sparselength = sparse_crs->NumGlobalEntries(i);
            std::vector<double> sparsevalues(sparselength);
            std::vector<int> sparseindices(sparselength);
            // int sparseextractionstatus =
            sparse_crs->ExtractGlobalRowCopy(
                i, sparselength, sparsenumentries, sparsevalues.data(), sparseindices.data());
            for (int k = 0; k < sparselength; ++k)
            {
              if (sparseindices[k] == j)
              {
                sparse_ij = sparsevalues[k];
                break;
              }
              else
                sparse_ij = 0.0;
            }
          }

          // get stiff_approx entry ij
          {
            int approxnumentries;
            int approxlength = stiff_approx->NumGlobalEntries(i);
            std::vector<double> approxvalues(approxlength);
            std::vector<int> approxindices(approxlength);
            // int approxextractionstatus =
            stiff_approx->ExtractGlobalRowCopy(
                i, approxlength, approxnumentries, approxvalues.data(), approxindices.data());
            for (int k = 0; k < approxlength; ++k)
            {
              if (approxindices[k] == j)
              {
                stiff_approx_ij = approxvalues[k];
                break;
              }
              else
                stiff_approx_ij = 0.0;
            }
          }

          double error = 0.0;
          if (abs(stiff_approx_ij) > 1e-5)
            error = error_ij / (stiff_approx_ij);
          else if (abs(sparse_ij) > 1e-5)
            error = error_ij / (sparse_ij);

          if (abs(error) > abs(error_max)) error_max = abs(error);
          if (abs(error_ij) > abs(abs_error_max)) abs_error_max = abs(error_ij);

          if ((abs(error) > 1e-4))
          {
            if ((abs(error_ij) > 1e-5))
            //  if( (sparse_ij>1e-1) or (stiff_approx_ij>1e-1) )
            {
              std::cout << "finite difference check failed entry (" << i << "," << j
                        << ")! stiff: " << sparse_ij << ", approx: " << stiff_approx_ij
                        << " ,abs. error: " << error_ij << " , rel. error: " << error << std::endl;

              success = false;
            }
          }
        }
      }
    }
  }

  if (success)
  {
    std::cout << "finite difference check successful, max. rel. error: " << error_max
              << "  (max. abs. error: " << abs_error_max << ")" << std::endl;
    std::cout << "******************finite difference check done***************\n\n" << std::endl;
  }
  else
    dserror("PoroFDCheck failed in step: %d, iter: %d", Step(), iter_);
}

void POROELAST::Monolithic::EvaluateCondition(
    Teuchos::RCP<CORE::LINALG::SparseOperator> Sysmat, POROELAST::coupltype coupltype)
{
  nopen_handle_->Clear(coupltype);
  Teuchos::RCP<CORE::LINALG::SparseMatrix> ConstraintMatrix =
      nopen_handle_->ConstraintMatrix(coupltype);
  Teuchos::RCP<CORE::LINALG::SparseMatrix> StructVelConstraintMatrix =
      nopen_handle_->StructVelConstraintMatrix(coupltype);

  // evaluate condition on elements and assemble matrices
  FluidField()->EvaluateNoPenetrationCond(nopen_handle_->RHS(), ConstraintMatrix,
      StructVelConstraintMatrix, nopen_handle_->CondVector(), nopen_handle_->CondIDs(), coupltype);

  if (coupltype == POROELAST::fluidfluid)  // fluid fluid part
  {
    ConstraintMatrix->Complete();
    nopen_handle_->BuidNoPenetrationMap(
        FluidField()->Discretization()->Comm(), FluidField()->DofRowMap());
  }
  else  // fluid structure part
  {
    // double timescale = FluidField()->TimeScaling();
    double timescale = FluidField()->ResidualScaling();
    StructVelConstraintMatrix->Scale(timescale);
    StructVelConstraintMatrix->Complete(
        StructureField()->SystemMatrix()->RangeMap(), FluidField()->SystemMatrix()->RangeMap());
    ConstraintMatrix->Add(*StructVelConstraintMatrix, false, 1.0, 1.0);
    ConstraintMatrix->Complete(
        StructureField()->SystemMatrix()->RangeMap(), FluidField()->SystemMatrix()->RangeMap());
  }

  const Teuchos::RCP<const Epetra_Map>& nopenetrationmap = nopen_handle_->Extractor()->Map(1);
  const Teuchos::RCP<const Epetra_Map>& othermap = nopen_handle_->Extractor()->Map(0);
  ConstraintMatrix->ApplyDirichlet(*othermap, false);
  Sysmat->ApplyDirichlet(*nopenetrationmap, false);
  Sysmat->UnComplete();
  Sysmat->Add(*ConstraintMatrix, false, 1.0, 1.0);
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> POROELAST::Monolithic::StructFluidCouplingMatrix()
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> sparse =
      Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(k_sf_);
  if (sparse == Teuchos::null) dserror("cast to CORE::LINALG::SparseMatrix failed!");

  return sparse;
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> POROELAST::Monolithic::FluidStructCouplingMatrix()
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> sparse =
      Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(k_fs_);
  if (sparse == Teuchos::null) dserror("cast to CORE::LINALG::SparseMatrix failed!");

  return sparse;
}

Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>
POROELAST::Monolithic::StructFluidCouplingBlockMatrix()
{
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> blocksparse =
      Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(k_sf_);
  if (blocksparse == Teuchos::null) dserror("cast to CORE::LINALG::BlockSparseMatrixBase failed!");

  return blocksparse;
}


Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase>
POROELAST::Monolithic::FluidStructCouplingBlockMatrix()
{
  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> blocksparse =
      Teuchos::rcp_dynamic_cast<CORE::LINALG::BlockSparseMatrixBase>(k_fs_);
  if (blocksparse == Teuchos::null) dserror("cast to CORE::LINALG::BlockSparseMatrixBase failed!");

  return blocksparse;
}

void POROELAST::Monolithic::SetupNewton()
{
  // initialise equilibrium loop and norms
  iter_ = 1;
  normrhs_ = 0.0;
  norminc_ = 0.0;
  normrhsfluid_ = 0.0;
  normincfluid_ = 0.0;
  normrhsfluidvel_ = 0.0;
  normincfluidvel_ = 0.0;
  normrhsfluidpres_ = 0.0;
  normincfluidpres_ = 0.0;
  normrhsstruct_ = 0.0;
  normincstruct_ = 0.0;
  normrhsporo_ = 0.0;
  normincporo_ = 0.0;

  // incremental solution vector with length of all dofs
  if (iterinc_ == Teuchos::null)
    iterinc_ = CORE::LINALG::CreateVector(*DofRowMap(), true);
  else
    iterinc_->PutScalar(0.0);

  // a zero vector of full length
  if (zeros_ == Teuchos::null)
    zeros_ = CORE::LINALG::CreateVector(*DofRowMap(), true);
  else
    zeros_->PutScalar(0.0);

  // AitkenReset();
}

void POROELAST::Monolithic::BuildConvergenceNorms()
{
  //------------------------------------------------------------ build residual force norms
  normrhs_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_);
  Teuchos::RCP<const Epetra_Vector> rhs_s;
  Teuchos::RCP<const Epetra_Vector> rhs_f;
  Teuchos::RCP<const Epetra_Vector> rhs_fvel;
  Teuchos::RCP<const Epetra_Vector> rhs_fpres;

  // process structure unknowns of the first field
  rhs_s = Extractor()->ExtractVector(rhs_, 0);
  // process fluid unknowns of the second field
  rhs_f = Extractor()->ExtractVector(rhs_, 1);
  rhs_fvel = FluidField()->ExtractVelocityPart(rhs_f);
  rhs_fpres = FluidField()->ExtractPressurePart(rhs_f);

  if (porosity_dof_)
  {
    Teuchos::RCP<const Epetra_Vector> rhs_poro = porosity_splitter_->ExtractCondVector(rhs_s);
    Teuchos::RCP<const Epetra_Vector> rhs_sdisp = porosity_splitter_->ExtractOtherVector(rhs_s);

    normrhsstruct_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_sdisp);
    normrhsporo_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_poro);
  }
  else
    normrhsstruct_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_s);

  normrhsfluid_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_f);
  normrhsfluidvel_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_fvel);
  normrhsfluidpres_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_fpres);


  //------------------------------------------------------------- build residual increment norms
  iterinc_->Norm2(&norminc_);

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> interincs;
  Teuchos::RCP<const Epetra_Vector> interincf;
  Teuchos::RCP<const Epetra_Vector> interincfvel;
  Teuchos::RCP<const Epetra_Vector> interincfpres;
  // process structure unknowns of the first field
  interincs = Extractor()->ExtractVector(iterinc_, 0);
  // process fluid unknowns of the second field
  interincf = Extractor()->ExtractVector(iterinc_, 1);
  interincfvel = FluidField()->ExtractVelocityPart(interincf);
  interincfpres = FluidField()->ExtractPressurePart(interincf);

  if (porosity_dof_)
  {
    Teuchos::RCP<const Epetra_Vector> interincporo =
        porosity_splitter_->ExtractCondVector(interincs);
    Teuchos::RCP<const Epetra_Vector> interincsdisp =
        porosity_splitter_->ExtractOtherVector(interincs);

    normincstruct_ = UTILS::CalculateVectorNorm(vectornorminc_, interincsdisp);
    normincporo_ = UTILS::CalculateVectorNorm(vectornorminc_, interincporo);
  }
  else
    normincstruct_ = UTILS::CalculateVectorNorm(vectornorminc_, interincs);

  normincfluid_ = UTILS::CalculateVectorNorm(vectornorminc_, interincf);
  normincfluidvel_ = UTILS::CalculateVectorNorm(vectornorminc_, interincfvel);
  normincfluidpres_ = UTILS::CalculateVectorNorm(vectornorminc_, interincfpres);
}

bool POROELAST::Monolithic::SetupSolver()
{
  //  solver
  // create a linear solver
  // get dynamic section of poroelasticity
  const Teuchos::ParameterList& poroelastdyn = DRT::Problem::Instance()->PoroelastDynamicParams();
  // get the solver number used for linear poroelasticity solver
  const int linsolvernumber = poroelastdyn.get<int>("LINEAR_SOLVER");
  // check if the poroelasticity solver has a valid solver number
  if (linsolvernumber == (-1))
  {
    dserror(
        "no linear solver defined for poroelasticity. Please set LINEAR_SOLVER in POROELASTICITY "
        "DYNAMIC to a valid number!");
  }
  const Teuchos::ParameterList& solverparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const auto solvertype =
      Teuchos::getIntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");

  directsolve_ = (solvertype == INPAR::SOLVER::SolverType::umfpack or
                  solvertype == INPAR::SOLVER::SolverType::superlu);

  if (directsolve_)
  {
    solver_ = Teuchos::rcp(new CORE::LINALG::Solver(
        solverparams, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  }
  else
    // create a linear solver
    CreateLinearSolver();

  // Get the parameters for the Newton iteration
  itermax_ = poroelastdyn.get<int>("ITEMAX");
  itermin_ = poroelastdyn.get<int>("ITEMIN");
  normtypeinc_ = DRT::INPUT::IntegralValue<INPAR::POROELAST::ConvNorm>(poroelastdyn, "NORM_INC");
  normtypefres_ = DRT::INPUT::IntegralValue<INPAR::POROELAST::ConvNorm>(poroelastdyn, "NORM_RESF");
  combincfres_ =
      DRT::INPUT::IntegralValue<INPAR::POROELAST::BinaryOp>(poroelastdyn, "NORMCOMBI_RESFINC");
  vectornormfres_ =
      DRT::INPUT::IntegralValue<INPAR::POROELAST::VectorNorm>(poroelastdyn, "VECTORNORM_RESF");
  vectornorminc_ =
      DRT::INPUT::IntegralValue<INPAR::POROELAST::VectorNorm>(poroelastdyn, "VECTORNORM_INC");

  tolinc_ = poroelastdyn.get<double>("TOLINC_GLOBAL");
  tolfres_ = poroelastdyn.get<double>("TOLRES_GLOBAL");

  tolinc_struct_ = poroelastdyn.get<double>("TOLINC_DISP");
  tolinc_velocity_ = poroelastdyn.get<double>("TOLINC_VEL");
  tolinc_pressure_ = poroelastdyn.get<double>("TOLINC_PRES");
  tolinc_porosity_ = poroelastdyn.get<double>("TOLINC_PORO");

  tolfres_struct_ = poroelastdyn.get<double>("TOLRES_DISP");
  tolfres_velocity_ = poroelastdyn.get<double>("TOLRES_VEL");
  tolfres_pressure_ = poroelastdyn.get<double>("TOLRES_PRES");
  tolfres_porosity_ = poroelastdyn.get<double>("TOLRES_PORO");

  return true;
}

Teuchos::RCP<CORE::LINALG::SparseMatrix> POROELAST::Monolithic::SystemMatrix()
{
  return systemmatrix_->Merge();
}

void POROELAST::Monolithic::IncrementPoroIter() { iter_ += 1; }

void POROELAST::Monolithic::UpdatePoroIterinc(Teuchos::RCP<const Epetra_Vector> poroinc)
{
  iterinc_->PutScalar(0.0);
  iterinc_->Update(1.0, *poroinc, 0.0);
}

void POROELAST::Monolithic::ClearPoroIterinc() { iterinc_->PutScalar(0.0); }

void POROELAST::Monolithic::Aitken()
{
  // initialise increment vector with solution of last iteration (i)
  // update del_ with current residual vector
  // difference of last two solutions
  if (del_ == Teuchos::null)  // first iteration, itnum==1
  {
    del_ = CORE::LINALG::CreateVector(*DofRowMap(), true);
    delhist_ = CORE::LINALG::CreateVector(*DofRowMap(), true);
    del_->PutScalar(1.0e20);
    delhist_->PutScalar(0.0);
  }

  // calculate difference of current (i+1) and old (i) residual vector
  // delhist = ( r^{i+1}_{n+1} - r^i_{n+1} )
  // update history vector old increment r^i_{n+1}
  delhist_->Update(1.0, *del_, 0.0);         // r^i_{n+1}
  delhist_->Update(1.0, *iterinc_, (-1.0));  // update r^{i+1}_{n+1}


  // del_ = r^{i+1}_{n+1} = T^{i+1}_{n+1} - T^{i}_{n+1}
  del_->Update(1.0, *iterinc_, 0.0);
  // den = |r^{i+1} - r^{i}|^2
  double den = 0.0;
  delhist_->Norm2(&den);
  // calculate dot product
  // dot = delhist_ . del_ = ( r^{i+1}_{n+1} - r^i_{n+1} )^T . r^{i+1}_{n+1}
  double top = 0.0;
  delhist_->Dot(*del_, &top);

  // mu_: Aikten factor in Mok's version
  // mu_: relaxation parameter in Irons & Tuck
  // mu^{i+1} = mu^i + (mu^i -1) . (r^{i+1} - r^i)^T . (-r^{i+1}) / |r^{i+1} - r^{i}|^2
  // top = ( r^{i+1} - r^i )^T . r^{i+1} --> use -top
  // Uli's implementation: mu_ = mu_ + (mu_ - 1.0) * top / (den*den). with '-' included in top
  mu_ = mu_ + (mu_ - 1) * (-top) / (den * den);

  iterinc_->Scale(1.0 - mu_);
}

[[maybe_unused]] void POROELAST::Monolithic::AitkenReset()
{
  if (del_ == Teuchos::null)  // first iteration, itnum==1
  {
    del_ = CORE::LINALG::CreateVector(*DofRowMap(), true);
    delhist_ = CORE::LINALG::CreateVector(*DofRowMap(), true);
  }
  del_->PutScalar(1.0e20);
  delhist_->PutScalar(0.0);
  mu_ = 0.0;
}

const Epetra_Map& POROELAST::Monolithic::FluidRangeMap()
{
  return FluidField()->SystemMatrix()->RangeMap();
}

const Epetra_Map& POROELAST::Monolithic::FluidDomainMap()
{
  return FluidField()->SystemMatrix()->DomainMap();
}

const Epetra_Map& POROELAST::Monolithic::StructureDomainMap()
{
  return StructureField()->DomainMap();
}

Teuchos::RCP<const Epetra_Map> POROELAST::Monolithic::DofRowMap()
{
  return blockrowdofmap_->FullMap();
}

Teuchos::RCP<const Epetra_Map> POROELAST::Monolithic::DofRowMapStructure()
{
  return blockrowdofmap_->Map(0);
}

Teuchos::RCP<const Epetra_Map> POROELAST::Monolithic::DofRowMapFluid()
{
  return blockrowdofmap_->Map(1);
}

void POROELAST::Monolithic::RecoverLagrangeMultiplierAfterNewtonStep(
    Teuchos::RCP<const Epetra_Vector> x)
{
  // clean up as soon as old time integration is unused!
  if (oldstructimint_)
  {
    if (StructureField()->MeshtyingContactBridge() != Teuchos::null)
    {
      if (StructureField()->MeshtyingContactBridge()->HaveContact() && !nit_contact_)
      {
        // Recover structural contact lagrange multiplier !!! For Poro & FPSI Problems this is
        // deactivated in the Structure!!!

        CONTACT::PoroLagrangeStrategy& costrategy = static_cast<CONTACT::PoroLagrangeStrategy&>(
            StructureField()->MeshtyingContactBridge()->ContactManager()->GetStrategy());

        // displacement and fluid velocity & pressure incremental vector
        Teuchos::RCP<const Epetra_Vector> sx;
        Teuchos::RCP<const Epetra_Vector> fx;
        ExtractFieldVectors(x, sx, fx);

        // RecoverStructuralLM
        Teuchos::RCP<Epetra_Vector> tmpsx = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*sx));
        Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*fx));

        costrategy.RecoverCoupled(tmpsx, tmpfx);
        if (no_penetration_) costrategy.RecoverPoroNoPen(tmpsx, tmpfx);
      }
      else if (StructureField()->MeshtyingContactBridge()->HaveMeshtying())
      {  // if meshtying  //h.Willmann

        CONTACT::PoroMtLagrangeStrategy& costrategy = static_cast<CONTACT::PoroMtLagrangeStrategy&>(
            StructureField()->MeshtyingContactBridge()->MtManager()->GetStrategy());

        // displacement and fluid velocity & pressure incremental vector
        Teuchos::RCP<const Epetra_Vector> sx;
        Teuchos::RCP<const Epetra_Vector> fx;
        ExtractFieldVectors(x, sx, fx);

        Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp<Epetra_Vector>(new Epetra_Vector(*fx));

        // Recover part of LM stemming from offdiagonal coupling matrix
        costrategy.RecoverCouplingMatrixPartofLMP(tmpfx);
      }
    }
  }
}

void POROELAST::Monolithic::SetPoroContactStates()
{
  // clean up as soon as old time integration is unused!
  if (oldstructimint_)
  {
    if (StructureField()->MeshtyingContactBridge() != Teuchos::null)
    {
      if (StructureField()->MeshtyingContactBridge()->HaveContact())
      {
        if (!nit_contact_)  // Lagmult contact
        {
          CONTACT::PoroLagrangeStrategy& costrategy = static_cast<CONTACT::PoroLagrangeStrategy&>(
              StructureField()->MeshtyingContactBridge()->ContactManager()->GetStrategy());
          Teuchos::RCP<Epetra_Vector> fvel = Teuchos::rcp(
              new Epetra_Vector(*FluidField()->ExtractVelocityPart(FluidField()->Velnp())));
          fvel = FluidStructureCoupling().SlaveToMaster(fvel);
          costrategy.SetState(MORTAR::state_fvelocity, *fvel);

          // To get pressure dofs into first structural component!!! - any idea for nice
          // implementation?
          Teuchos::RCP<const Epetra_Vector> fpres =
              FluidField()->ExtractPressurePart(FluidField()->Velnp());
          Teuchos::RCP<Epetra_Vector> modfpres =
              Teuchos::rcp(new Epetra_Vector(*FluidField()->VelocityRowMap(), true));

          int* mygids = fpres->Map().MyGlobalElements();
          double* val = fpres->Values();
          const int ndim = DRT::Problem::Instance()->NDim();
          for (int i = 0; i < fpres->MyLength(); ++i)
          {
            int gid = mygids[i] - ndim;
            modfpres->ReplaceGlobalValues(1, &val[i], &gid);
          }

          modfpres = FluidStructureCoupling().SlaveToMaster(modfpres);
          costrategy.SetState(MORTAR::state_fpressure, *modfpres);

          Teuchos::RCP<Epetra_Vector> dis =
              Teuchos::rcp(new Epetra_Vector(*StructureField()->Dispnp()));
          costrategy.SetParentState("displacement", dis,
              StructureField()->Discretization());  // add displacements of the parent element!!!
        }
        else
        {
          CONTACT::CoNitscheStrategyPoro& costrategy =
              dynamic_cast<CONTACT::CoNitscheStrategyPoro&>(
                  StructureField()->MeshtyingContactBridge()->ContactManager()->GetStrategy());
          costrategy.SetParentState(MORTAR::state_fvelocity, *FluidField()->Velnp());
        }
      }
    }
  }
}

void POROELAST::Monolithic::EvalPoroMortar()
{
  // clean up as soon as old time integration is unused!
  if (oldstructimint_)
  {
    if (StructureField()->MeshtyingContactBridge() != Teuchos::null)
    {
      if (StructureField()->MeshtyingContactBridge()->HaveContact())
      {
        if (!nit_contact_)  // Lagmult contact
        {
          CONTACT::PoroLagrangeStrategy& costrategy = static_cast<CONTACT::PoroLagrangeStrategy&>(
              StructureField()->MeshtyingContactBridge()->ContactManager()->GetStrategy());
          //---Modifiy coupling matrix k_sf

          // Get matrix block!
          Teuchos::RCP<CORE::LINALG::SparseOperator> k_ss =
              Teuchos::rcp<CORE::LINALG::SparseMatrix>(
                  new CORE::LINALG::SparseMatrix(systemmatrix_->Matrix(0, 0)));
          Teuchos::RCP<CORE::LINALG::SparseOperator> k_sf =
              Teuchos::rcp<CORE::LINALG::SparseMatrix>(
                  new CORE::LINALG::SparseMatrix(systemmatrix_->Matrix(0, 1)));
          Teuchos::RCP<Epetra_Vector> rhs_s = Extractor()->ExtractVector(rhs_, 0);

          // Evaluate Poro Contact Condensation for K_ss, K_sf
          costrategy.ApplyForceStiffCmtCoupled(
              StructureField()->WriteAccessDispnp(), k_ss, k_sf, rhs_s, Step(), iter_, false);

          // Assign modified matrixes & vectors
          systemmatrix_->Assign(0, 0, CORE::LINALG::Copy,
              *Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(k_ss));
          systemmatrix_->Assign(0, 1, CORE::LINALG::Copy,
              *Teuchos::rcp_dynamic_cast<CORE::LINALG::SparseMatrix>(k_sf));
          Extractor()->InsertVector(rhs_s, 0, rhs_);

          //---Modify fluid matrix, coupling matrix k_fs and rhs of fluid
          if (no_penetration_)
          {
            //---Initialize Poro Contact
            costrategy.PoroInitialize(FluidStructureCoupling(),
                FluidField()->DofRowMap());  // true stands for the no_penetration condition !!!
            // Get matrix blocks & rhs vector!
            Teuchos::RCP<CORE::LINALG::SparseMatrix> f = Teuchos::rcp<CORE::LINALG::SparseMatrix>(
                new CORE::LINALG::SparseMatrix(systemmatrix_->Matrix(1, 1)));
            Teuchos::RCP<CORE::LINALG::SparseMatrix> k_fs =
                Teuchos::rcp<CORE::LINALG::SparseMatrix>(
                    new CORE::LINALG::SparseMatrix(systemmatrix_->Matrix(1, 0)));

            Teuchos::RCP<Epetra_Vector> frhs = Extractor()->ExtractVector(rhs_, 1);

            // Evaluate Poro No Penetration Contact Condensation
            costrategy.EvaluatePoroNoPenContact(k_fs, f, frhs);

            // Assign modified matrixes & vectors
            systemmatrix_->Assign(1, 1, CORE::LINALG::Copy, *f);
            systemmatrix_->Assign(1, 0, CORE::LINALG::Copy, *k_fs);

            Extractor()->InsertVector(*frhs, 1, *rhs_);
          }
        }
        else  // add nitsche contributions to sysmat and rhs
        {
          CONTACT::CoNitscheStrategyPoro* pstrat = dynamic_cast<CONTACT::CoNitscheStrategyPoro*>(
              &StructureField()->MeshtyingContactBridge()->ContactManager()->GetStrategy());

          systemmatrix_->UnComplete();
          systemmatrix_->Matrix(0, 1).Add(
              *pstrat->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::displ_porofluid), false, 1.0,
              1.0);
          systemmatrix_->Matrix(1, 1).Add(
              *pstrat->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::porofluid_porofluid), false, 1.0,
              1.0);
          systemmatrix_->Matrix(1, 0).Add(
              *pstrat->GetMatrixBlockPtr(DRT::UTILS::MatBlockType::porofluid_displ), false, 1.0,
              1.0);
          systemmatrix_->Complete();

          Extractor()->AddVector(
              *pstrat->GetRhsBlockPtr(DRT::UTILS::VecBlockType::porofluid), 1, *rhs_);
        }
      }
      else if (StructureField()->MeshtyingContactBridge()->HaveMeshtying())
      {  // if meshtying  //h.Willmann

        CONTACT::PoroMtLagrangeStrategy& costrategy = static_cast<CONTACT::PoroMtLagrangeStrategy&>(
            StructureField()->MeshtyingContactBridge()->MtManager()->GetStrategy());


        //---Modifiy coupling matrix k_sf

        // Get matrix block!
        Teuchos::RCP<CORE::LINALG::SparseMatrix> k_sf = Teuchos::rcp<CORE::LINALG::SparseMatrix>(
            new CORE::LINALG::SparseMatrix(systemmatrix_->Matrix(0, 1)));

        // initialize poro meshtying
        costrategy.InitializePoroMt(k_sf);

        // Evaluate Poro Meshtying Condensation for K_sf
        costrategy.EvaluateMeshtyingPoroOffDiag(k_sf);

        // Assign modified matrix
        systemmatrix_->Assign(0, 1, CORE::LINALG::Copy, *k_sf);
      }
    }
  }
}

void POROELAST::Monolithic::BuildCombinedDBCMap()
{
  const Teuchos::RCP<const Epetra_Map> scondmap = StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> fcondmap = FluidField()->GetDBCMapExtractor()->CondMap();
  combinedDBCMap_ = CORE::LINALG::MergeMap(scondmap, fcondmap, false);
}

void POROELAST::Monolithic::ReadRestart(const int step)
{
  // call base class
  POROELAST::PoroBase::ReadRestart(step);

  // set states for porous contact
  // apply current velocity of fluid to ContactMangager if contact problem
  if (no_penetration_) SetPoroContactStates();
}
