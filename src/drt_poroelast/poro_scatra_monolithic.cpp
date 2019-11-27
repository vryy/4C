/*----------------------------------------------------------------------*/
/*! \file

 \brief monolithic coupling algorithm for scalar transport within porous medium

\level 2

\maintainer Johannes Kremheller
            kremheller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249

 *----------------------------------------------------------------------*/

#include "poro_scatra_monolithic.H"

#include <Teuchos_TimeMonitor.hpp>

#include "poro_base.H"
#include "poroelast_utils.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"

#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../drt_io/io_control.H"
#include "../drt_inpar/inpar_poroelast.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_assemblestrategy.H"

#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_densematrix_manipulation.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_blocksparsematrix.H"



/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
POROELAST::PoroScatraMono::PoroScatraMono(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraBase(comm, timeparams),
      printscreen_(true),   // ADD INPUT PARAMETER
      printiter_(true),     // ADD INPUT PARAMETER
      printerrfile_(true),  // ADD INPUT PARAMETER FOR 'true'
      errfile_(DRT::Problem::Instance()->ErrorFile()->Handle()),
      timer_(comm),
      iterinc_(Teuchos::null),
      zeros_(Teuchos::null),
      systemmatrix_(Teuchos::null),
      rhs_(Teuchos::null),
      blockrowdofmap_(Teuchos::null),
      directsolve_(true)
{
  const Teuchos::ParameterList& sdynparams = DRT::Problem::Instance()->StructuralDynamicParams();

  // some solver paramaters are red form the structure dynamic list (this is not the best way to do
  // it ...)
  solveradapttol_ = (DRT::INPUT::IntegralValue<int>(sdynparams, "ADAPTCONV") == 1);
  solveradaptolbetter_ = (sdynparams.get<double>("ADAPTCONV_BETTER"));

  blockrowdofmap_ = Teuchos::rcp(new LINALG::MultiMapExtractor);
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::Timeloop()
{
  while (NotFinished())
  {
    DoTimeStep();
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::DoTimeStep()
{
  // counter and print header
  // predict solution of both field (call the adapter)
  PrepareTimeStep();

  // Newton-Raphson iteration
  Solve();

  // calculate stresses, strains, energies
  PrepareOutput();

  // update all single field solvers
  Update();

  // write output to screen and files
  Output();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::ReadRestart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    SetScatraSolution();
    SetPoroSolution();

    PoroField()->ReadRestart(restart);
    ScaTraField()->ReadRestart(restart);

    // in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if (PoroField()->HasSubmeshes())
      ReplaceDofSets(StructureField()->Discretization(), FluidField()->Discretization(),
          ScaTraField()->Discretization());

    // the variables need to be set on other field
    SetScatraSolution();
    SetPoroSolution();

    // second restart needed due to two way coupling.
    ScaTraField()->ReadRestart(restart);
    PoroField()->ReadRestart(restart);

    // in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if (PoroField()->HasSubmeshes())
      ReplaceDofSets(StructureField()->Discretization(), FluidField()->Discretization(),
          ScaTraField()->Discretization());

    SetTimeStep(PoroField()->Time(), restart);

    // Material pointers to other field were deleted during ReadRestart().
    // They need to be reset.
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(
        PoroField()->StructureField()->Discretization(), ScaTraField()->Discretization());
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(
        PoroField()->FluidField()->Discretization(), ScaTraField()->Discretization());
  }
}

/*----------------------------------------------------------------------*/
// prepare time step
/*----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::PrepareTimeStep(bool printheader)
{
  // the global control routine has its own time_ and step_ variables, as well as the single fields
  // keep them in sinc!
  IncrementTimeAndStep();

  if (printheader) PrintHeader();

  SetPoroSolution();
  ScaTraField()->PrepareTimeStep();
  // set structure-based scalar transport values
  SetScatraSolution();

  PoroField()->PrepareTimeStep();
  SetPoroSolution();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::PrepareOutput() { PoroField()->PrepareOutput(); }

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::Output()
{
  PoroField()->Output();
  ScaTraField()->Output();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::Update()
{
  PoroField()->Update();

  ScaTraField()->Update();
  ScaTraField()->EvaluateErrorComparedToAnalyticalSol();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::Solve()
{
  // we do a Newton-Raphson iteration here.
  // the specific time integration has set the following
  // --> On #rhs_ is the positive force residuum
  // --> On #systemmatrix_ is the effective dynamic tangent matrix

  //-------------------------------------- initialize variables needed in newton loop
  iter_ = 1;
  normrhs_ = 0.0;
  norminc_ = 0.0;

  // incremental solution vector with length of all dofs
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_->PutScalar(0.0);

  // a zero vector of full length
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_->PutScalar(0.0);

  //---------------------------------------------- iteration loop

  // equilibrium iteration loop (loop over k)
  while (((not Converged()) and (iter_ <= itermax_)) or (iter_ <= itermin_))
  {
    timer_.ResetStartTime();
    Epetra_Time timer(Comm());

    // compute residual forces #rhs_ and tangent #tang_
    // whose components are globally oriented
    // build linear system stiffness matrix and rhs/force residual for each
    // field, here e.g. for structure field: field want the iteration increment
    // 1.) Update(iterinc_),
    // 2.) EvaluateForceStiffResidual(),
    // 3.) PrepareSystemForNewtonSolve()
    Evaluate(iterinc_);
    // cout << "  time for Evaluate diagonal blocks: " << timer.ElapsedTime() << "\n";
    // timer.ResetStartTime();

    // check whether we have a sanely filled tangent matrix
    if (not systemmatrix_->Filled())
    {
      dserror("Effective tangent matrix must be filled here");
    }

    // create full monolithic rhs vector
    SetupRHS(iter_ == 1);
    // cout << "  time for Evaluate SetupRHS: " << timer.ElapsedTime() << "\n";
    // timer.ResetStartTime();

    // FDCheck();

    // build norms
    BuildConvergenceNorms();

    if ((not Converged()) or combincfres_ == INPAR::POROELAST::bop_or)
    {
      // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
      // is done in PrepareSystemForNewtonSolve() within Evaluate(iterinc_)
      LinearSolve();
      // cout << "  time for Evaluate LinearSolve: " << timer.ElapsedTime() << "\n";
      // timer.ResetStartTime();

      // reset solver tolerance
      solver_->ResetTolerance();

      // build norms
      BuildConvergenceNorms();
    }

    // print stuff
    PrintNewtonIter();

    // increment equilibrium loop index
    iter_ += 1;
  }  // end equilibrium loop

  //---------------------------------------------- iteration loop

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

}  // Solve()

/*----------------------------------------------------------------------*
 | evaluate the single fields                              vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("PoroScatraMono::Monolithic::Evaluate");

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> porostructinc;
  Teuchos::RCP<const Epetra_Vector> porofluidinc;
  Teuchos::RCP<const Epetra_Vector> scatrainc;

  // if an increment vector exists
  if (x != Teuchos::null)
  {
    // process structure unknowns of the first field
    porostructinc = Extractor()->ExtractVector(x, 0);
    porofluidinc = Extractor()->ExtractVector(x, 1);

    // process fluid unknowns of the second field
    scatrainc = Extractor()->ExtractVector(x, 2);
  }

  // Newton update of the fluid field
  // update velocities and pressures before passed to the structural field
  //  UpdateIterIncrementally(fx),
  ScaTraField()->UpdateIter(scatrainc);

  // call all elements and assemble rhs and matrices
  /// poro field

  // structure Evaluate (builds tangent, residual and applies DBC)
  Epetra_Time timerporo(Comm());

  // apply current velocity and pressures to structure
  SetScatraSolution();

  // access poro problem to build poro-poro block
  PoroField()->Evaluate(porostructinc, porofluidinc, iter_ == 1);

  /// scatra field

  // Scatra Evaluate
  // (builds tangent, residual and applies DBC and recent coupling values)
  SetPoroSolution();

  // access ScaTra problem to build scatra-scatra block
  ScaTraField()->PrepareLinearSolve();

  // fill off diagonal blocks and build monolithic system matrix
  SetupSystemMatrix();
}

/*----------------------------------------------------------------------*
 | setup system (called in poro_dyn.cpp)                 vuong 01/12    |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::SetupSystem()
{
  // setup the poro subsystem first
  PoroField()->SetupSystem();

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

  {
    // vecSpaces.push_back(PoroField()->DofRowMap());
    vecSpaces.push_back(PoroField()->DofRowMapStructure());
    vecSpaces.push_back(PoroField()->DofRowMapFluid());
    const Epetra_Map* dofrowmapscatra = (ScaTraField()->Discretization())->DofRowMap(0);
    vecSpaces.push_back(Teuchos::rcp(dofrowmapscatra, false));
  }

  if (vecSpaces[0]->NumGlobalElements() == 0) dserror("No poro structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements() == 0) dserror("No poro fluid equation. Panic.");
  if (vecSpaces[2]->NumGlobalElements() == 0) dserror("No scatra equation. Panic.");

  // build dof row map of monolithic system
  SetDofRowMaps(vecSpaces);

  // build dbc map of monolithic system
  {
    const Teuchos::RCP<const Epetra_Map> porocondmap = PoroField()->CombinedDBCMap();
    const Teuchos::RCP<const Epetra_Map> scatracondmap = ScaTraField()->DirichMaps()->CondMap();
    Teuchos::RCP<const Epetra_Map> dbcmap = LINALG::MergeMap(porocondmap, scatracondmap, false);

    // Finally, create the global FSI Dirichlet map extractor
    dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor(*DofRowMap(), dbcmap, true));
    if (dbcmaps_ == Teuchos::null) dserror("Creation of Dirichlet map extractor failed.");
  }

  // initialize Poroscatra-systemmatrix_
  systemmatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *Extractor(), *Extractor(), 81, false, true));

  {
    std::vector<Teuchos::RCP<const Epetra_Map>> scatravecSpaces;
    const Epetra_Map* dofrowmapscatra = (ScaTraField()->Discretization())->DofRowMap(0);
    scatravecSpaces.push_back(Teuchos::rcp(dofrowmapscatra, false));
    scatrarowdofmap_.Setup(*dofrowmapscatra, scatravecSpaces);
  }

  k_pss_ =
      Teuchos::rcp(new LINALG::SparseMatrix(*(PoroField()->DofRowMapStructure()), 81, true, true));
  k_pfs_ = Teuchos::rcp(new LINALG::SparseMatrix(*(PoroField()->DofRowMapFluid()),
      //*(FluidField()->DofRowMap()),
      81, true, true));

  k_sps_ = Teuchos::rcp(
      new LINALG::SparseMatrix(*(ScaTraField()->Discretization()->DofRowMap()), 81, true, true));
  k_spf_ = Teuchos::rcp(new LINALG::SparseMatrix(*(ScaTraField()->Discretization()->DofRowMap()),
      //*(FluidField()->DofRowMap()),
      81, true, true));
}

/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                                 vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::SetupRHS(bool firstcall)
{
  // create full monolithic rhs vector
  if (rhs_ == Teuchos::null) rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  // fill the Poroelasticity rhs vector rhs_ with the single field rhss
  SetupVector(*rhs_, PoroField()->RHS(), ScaTraField()->Residual());
}

/*----------------------------------------------------------------------*
 | setup vector of the structure and fluid field            vuong 01/12|
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::SetupVector(
    Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> pv, Teuchos::RCP<const Epetra_Vector> sv)
{
  // extract dofs of the two fields
  // and put the poro/scatra field vector into the global vector f
  // noticing the block number

  //  Teuchos::RCP<const Epetra_Vector> psx;
  //  Teuchos::RCP<const Epetra_Vector> pfx;

  Extractor()->InsertVector(*(PoroField()->Extractor()->ExtractVector(pv, 0)), 0, f);
  Extractor()->InsertVector(*(PoroField()->Extractor()->ExtractVector(pv, 1)), 1, f);
  Extractor()->InsertVector(*sv, 2, f);
}

/*----------------------------------------------------------------------*
 | setup system matrix of poroelasticity                   vuong 01/12  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::SetupSystemMatrix()
{
  // set loma block matrix to zero
  systemmatrix_->Zero();

  //----------------------------------------------------------------------
  // 1st diagonal block (upper left): poro weighting - poro solution
  //----------------------------------------------------------------------
  // get matrix block
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mat_pp = PoroField()->BlockSystemMatrix();

  // uncomplete matrix block (appears to be required in certain cases)
  mat_pp->UnComplete();

  // assign matrix block
  systemmatrix_->Assign(0, 0, LINALG::View, mat_pp->Matrix(0, 0));
  systemmatrix_->Assign(0, 1, LINALG::View, mat_pp->Matrix(0, 1));
  systemmatrix_->Assign(1, 0, LINALG::View, mat_pp->Matrix(1, 0));
  systemmatrix_->Assign(1, 1, LINALG::View, mat_pp->Matrix(1, 1));

  //----------------------------------------------------------------------
  // 2nd diagonal block (lower right): scatra weighting - scatra solution
  //----------------------------------------------------------------------
  // get matrix block
  Teuchos::RCP<LINALG::SparseMatrix> mat_ss = ScaTraField()->SystemMatrix();

  // uncomplete matrix block (appears to be required in certain cases)
  mat_ss->UnComplete();

  // assign matrix block
  systemmatrix_->Assign(2, 2, LINALG::View, *mat_ss);

  // complete scatra block matrix
  systemmatrix_->Complete();

  //----------------------------------------------------------------------
  // 1st off-diagonal block (upper right): poro weighting - scatra solution
  //----------------------------------------------------------------------

  // evaluate off-diagonal matrix block in fluid
  EvaluateODBlockMatPoro();

  // k_ps_->Complete(mat_pp->DomainMap(),mat_ss->RangeMap());
  //  k_ps_->Complete();
  //
  //  Teuchos::RCP<LINALG::SparseMatrix> k_ps_sparse = k_ps_->Merge();

  // uncomplete matrix block (appears to be required in certain cases)
  // k_ps_sparse->UnComplete();
  k_pss_->UnComplete();
  k_pfs_->UnComplete();

  // assign matrix block
  systemmatrix_->Assign(0, 2, LINALG::View, *(k_pss_));
  systemmatrix_->Assign(1, 2, LINALG::View, *(k_pfs_));

  //----------------------------------------------------------------------
  // 2nd off-diagonal block (lower left): scatra weighting - poro solution
  //----------------------------------------------------------------------

  // evaluate off-diagonal matrix block in scatra
  EvaluateODBlockMatScatra();

  //  k_sp_->Complete();
  //
  //  Teuchos::RCP<LINALG::SparseMatrix> k_sp_sparse = k_sp_->Merge();

  // uncomplete matrix block (appears to be required in certain cases)
  // k_sp_sparse->UnComplete();
  k_sps_->UnComplete();
  k_spf_->UnComplete();

  // assign matrix block
  systemmatrix_->Assign(2, 0, LINALG::View, *(k_sps_));
  systemmatrix_->Assign(2, 1, LINALG::View, *(k_spf_));

  // complete block matrix
  systemmatrix_->Complete();
}

/*----------------------------------------------------------------------*
 | Solve linear Poroelasticity system                      vuong 01/12   |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::LinearSolve()
{
  // Solve for inc_ = [disi_,tempi_]
  // Solve K_Teffdyn . IncX = -R  ===>  IncX_{n+1} with X=[d,T]
  // \f$x_{i+1} = x_i + \Delta x_i\f$
  if (solveradapttol_ and (iter_ > 1))
  {
    double worst = normrhs_;
    double wanted = tolfres_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }
  // apply Dirichlet BCs to system of equations
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  if (directsolve_)
  {
    // merge blockmatrix to SparseMatrix and solve
    Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

    LINALG::ApplyDirichlettoSystem(
        sparse, iterinc_, rhs_, Teuchos::null, zeros_, *CombinedDBCMap());
    //  if ( Comm().MyPID()==0 ) { cout << " DBC applied to system" << endl; }

    // standard solver call
    solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
    //  if ( Comm().MyPID()==0 ) { cout << " Solved" << endl; }
  }
  else
  {
    // in case of inclined boundary conditions
    // rotate systemmatrix_ using GetLocSysTrafo()!=Teuchos::null
    LINALG::ApplyDirichlettoSystem(
        systemmatrix_, iterinc_, rhs_, Teuchos::null, zeros_, *CombinedDBCMap());

    solver_->Solve(systemmatrix_->EpetraOperator(), iterinc_, rhs_, true, iter_ == 1);
  }
}

/*----------------------------------------------------------------------*
 | setup solver for monolithic system                    vuong 01/12     |
 *----------------------------------------------------------------------*/
bool POROELAST::PoroScatraMono::SetupSolver()
{
  //  solver
  // create a linear solver
  // get dynamic section of poroelasticity
  const Teuchos::ParameterList& poroscatradyn = DRT::Problem::Instance()->PoroScatraControlParams();
  // get the solver number used for linear poroelasticity solver
  const int linsolvernumber = poroscatradyn.get<int>("LINEAR_SOLVER");
  // check if the poroelasticity solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for scalar transport in porous media. Please set LINEAR_SOLVER "
        "in POROSCATRA CONTROL to a valid number!");
  const Teuchos::ParameterList& solverparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const int solvertype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");

  directsolve_ = (solvertype == INPAR::SOLVER::umfpack);

  if (directsolve_)
  {
    solver_ = Teuchos::rcp(
        new LINALG::Solver(solverparams, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  }
  else
    // create a linear solver
    // CreateLinearSolver();
    dserror("no implicit solver supported yet!");

  // Get the parameters for the Newton iteration
  itermax_ = poroscatradyn.get<int>("ITEMAX");
  itermin_ = poroscatradyn.get<int>("ITEMIN");
  normtypeinc_ = DRT::INPUT::IntegralValue<INPAR::POROELAST::ConvNorm>(poroscatradyn, "NORM_INC");
  normtypeinc_ = DRT::INPUT::IntegralValue<INPAR::POROELAST::ConvNorm>(poroscatradyn, "NORM_INC");
  normtypefres_ = DRT::INPUT::IntegralValue<INPAR::POROELAST::ConvNorm>(poroscatradyn, "NORM_RESF");
  combincfres_ =
      DRT::INPUT::IntegralValue<INPAR::POROELAST::BinaryOp>(poroscatradyn, "NORMCOMBI_RESFINC");
  vectornormfres_ =
      DRT::INPUT::IntegralValue<INPAR::POROELAST::VectorNorm>(poroscatradyn, "VECTORNORM_RESF");
  vectornorminc_ =
      DRT::INPUT::IntegralValue<INPAR::POROELAST::VectorNorm>(poroscatradyn, "VECTORNORM_INC");

  tolinc_ = poroscatradyn.get<double>("TOLINC_GLOBAL");
  tolfres_ = poroscatradyn.get<double>("TOLRES_GLOBAL");

  tolinc_struct_ = poroscatradyn.get<double>("TOLINC_DISP");
  tolinc_velocity_ = poroscatradyn.get<double>("TOLINC_VEL");
  tolinc_pressure_ = poroscatradyn.get<double>("TOLINC_PRES");
  tolinc_scalar_ = poroscatradyn.get<double>("TOLINC_SCALAR");
  // tolinc_porosity_= poroscatradyn.get<double> ("TOLINC_PORO");

  tolfres_struct_ = poroscatradyn.get<double>("TOLRES_DISP");
  tolfres_velocity_ = poroscatradyn.get<double>("TOLRES_VEL");
  tolfres_pressure_ = poroscatradyn.get<double>("TOLRES_PRES");
  tolfres_scalar_ = poroscatradyn.get<double>("TOLINC_SCALAR");
  // tolfres_porosity_= poroscatradyn.get<double> ("TOLRES_PORO");

  return true;
}

/*----------------------------------------------------------------------*
 | check convergence of Newton iteration (public)         vuong 01/12   |
 *----------------------------------------------------------------------*/
bool POROELAST::PoroScatraMono::Converged()
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
                 normincfluidpres_ < tolinc_pressure_ and normincscalar_ < tolinc_scalar_
          //                  normincporo_      < tolinc_porosity_
      );
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
                  normrhsfluidpres_ < tolfres_pressure_ and normrhsscalar_ < tolfres_scalar_
          //                   normrhsporo_      < tolfres_porosity_
      );
      break;
    default:
      dserror("Cannot check for convergence of residual forces!");
      break;
  }

  // combine increments and forces
  bool conv = false;
  switch (combincfres_)
  {
    case INPAR::POROELAST::bop_and:
      conv = convinc and convfres;
      break;
    case INPAR::POROELAST::bop_or:
      conv = convinc or convfres;
      break;
    default:
      dserror("Something went terribly wrong with binary operator!");
      break;
  }

  // return things
  return conv;
}  // Converged()

/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file    vuong 08/13   |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::PrintNewtonIter()
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

  return;
}  // PrintNewtonIter()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen and error file    vuong 08/13  |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::PrintNewtonIterHeader(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

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
      //    case INPAR::POROELAST::convnorm_rel_global:
      //      oss << std::setw(18) << "rel-res";
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      //  case INPAR::POROELAST::convnorm_rel_singlefields:
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      oss << std::setw(15) << "abs-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_ << ")";
      break;
      //    case INPAR::POROELAST::convnorm_rel_global:
      //      oss << std::setw(18) << "rel-inc";
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      // case INPAR::POROELAST::convnorm_rel_singlefields:
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  // --------------------------------------------------------single field test
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      oss << std::setw(15) << "abs-s-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_struct_ << ")";
      //      if(porositydof_)
      //        oss <<std::setw(15)<< "abs-poro-res"<<"("<<std::setprecision(2)
      //        <<tolfres_porosity_<<")";
      oss << std::setw(15) << "abs-fvel-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_velocity_ << ")";
      oss << std::setw(15) << "abs-fpres-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_pressure_ << ")";
      oss << std::setw(15) << "abs-sca-res"
          << "(" << std::setw(5) << std::setprecision(2) << tolfres_scalar_ << ")";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::POROELAST::convnorm_abs_global:
      break;
    case INPAR::POROELAST::convnorm_abs_singlefields:
      oss << std::setw(15) << "abs-s-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_struct_ << ")";
      //      if(porositydof_)
      //        oss <<std::setw(15)<< "abs-poro-inc"<<"("<<std::setprecision(2)
      //        <<tolinc_porosity_<<")";
      oss << std::setw(15) << "abs-fvel-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_velocity_ << ")";
      oss << std::setw(15) << "abs-fpres-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_pressure_ << ")";
      oss << std::setw(15) << "abs-sca-inc"
          << "(" << std::setw(5) << std::setprecision(2) << tolinc_scalar_ << ")";
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  // add solution time
  oss << std::setw(14) << "wct";

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
 | print Newton-Raphson iteration to screen              vuong 08/13 |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::PrintNewtonIterText(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

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
      dserror("You should not turn up here.");
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
      dserror("You should not turn up here.");
      break;
  }

  // --------------------------------------------------------single field test
  switch (normtypefres_)
  {
    case INPAR::POROELAST::convnorm_abs_singlefields:
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsstruct_;
      //      if(porositydof_)
      //        oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsporo_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsfluidvel_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsfluidpres_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normrhsscalar_;
      break;
    case INPAR::POROELAST::convnorm_abs_global:
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case INPAR::POROELAST::convnorm_abs_singlefields:
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincstruct_;
      //      if(porositydof_)
      //        oss << std::setw(22) << std::setprecision(5) << std::scientific << normincporo_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincfluidvel_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincfluidpres_;
      oss << std::setw(22) << std::setprecision(5) << std::scientific << normincscalar_;
      break;
    case INPAR::POROELAST::convnorm_abs_global:
      break;
    default:
      dserror("You should not turn up here.");
      break;
  }

  // add solution time
  oss << std::setw(14) << std::setprecision(2) << std::scientific << timer_.ElapsedTime();

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
 | print statistics of converged NRI                      vuong 08/13    |
 | orignially by bborn 08/09                                            |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::PrintNewtonConv()
{
  // somebody did the door
  return;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::BuildConvergenceNorms()
{
  //------------------------------------------------------------ build residual force norms

  // global norm
  normrhs_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_);

  // split vectors
  Teuchos::RCP<const Epetra_Vector> rhs_s;
  Teuchos::RCP<const Epetra_Vector> rhs_f;
  Teuchos::RCP<const Epetra_Vector> rhs_fvel;
  Teuchos::RCP<const Epetra_Vector> rhs_fpres;
  Teuchos::RCP<const Epetra_Vector> rhs_scalar;

  // process structure unknowns of the first field
  rhs_s = Extractor()->ExtractVector(rhs_, 0);

  // process fluid unknowns of the second field
  rhs_f = Extractor()->ExtractVector(rhs_, 1);
  rhs_fvel = PoroField()->FluidField()->ExtractVelocityPart(rhs_f);
  rhs_fpres = PoroField()->FluidField()->ExtractPressurePart(rhs_f);

  // process scalar unknowns of the third field
  rhs_scalar = Extractor()->ExtractVector(rhs_, 2);

  //  if(porositydof_)
  //  {
  //    Teuchos::RCP<const Epetra_Vector> rhs_poro = porositysplitter_->ExtractCondVector(rhs_s);
  //    Teuchos::RCP<const Epetra_Vector> rhs_sdisp = porositysplitter_->ExtractOtherVector(rhs_s);
  //
  //    normrhsstruct_ = UTILS::CalculateVectorNorm(vectornormfres_,rhs_sdisp);
  //    normrhsporo_ = UTILS::CalculateVectorNorm(vectornormfres_,rhs_poro);
  //  }
  //  else
  normrhsstruct_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_s);

  normrhsfluid_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_f);
  normrhsfluidvel_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_fvel);
  normrhsfluidpres_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_fpres);

  normrhsscalar_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_scalar);


  //------------------------------------------------------------- build residual increment norms
  iterinc_->Norm2(&norminc_);

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> interincs;
  Teuchos::RCP<const Epetra_Vector> interincf;
  Teuchos::RCP<const Epetra_Vector> interincfvel;
  Teuchos::RCP<const Epetra_Vector> interincfpres;
  Teuchos::RCP<const Epetra_Vector> interincscalar;

  // process structure unknowns of the first field
  interincs = Extractor()->ExtractVector(iterinc_, 0);

  // process fluid unknowns of the second field
  interincf = Extractor()->ExtractVector(iterinc_, 1);
  interincfvel = PoroField()->FluidField()->ExtractVelocityPart(interincf);
  interincfpres = PoroField()->FluidField()->ExtractPressurePart(interincf);

  // process scalar unknowns of the third field
  interincscalar = Extractor()->ExtractVector(iterinc_, 2);

  //  if(porositydof_)
  //  {
  //    Teuchos::RCP<const Epetra_Vector> interincporo =
  //    porositysplitter_->ExtractCondVector(interincs); Teuchos::RCP<const Epetra_Vector>
  //    interincsdisp = porositysplitter_->ExtractOtherVector(interincs);
  //
  //    normincstruct_     = UTILS::CalculateVectorNorm(vectornorminc_,interincsdisp);
  //    normincporo_       = UTILS::CalculateVectorNorm(vectornorminc_,interincporo);
  //  }
  //  else
  normincstruct_ = UTILS::CalculateVectorNorm(vectornorminc_, interincs);

  normincfluid_ = UTILS::CalculateVectorNorm(vectornorminc_, interincf);
  normincfluidvel_ = UTILS::CalculateVectorNorm(vectornorminc_, interincfvel);
  normincfluidpres_ = UTILS::CalculateVectorNorm(vectornorminc_, interincfpres);

  normincscalar_ = UTILS::CalculateVectorNorm(vectornorminc_, interincscalar);

  return;
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROELAST::PoroScatraMono::DofRowMap() const
{
  return blockrowdofmap_->FullMap();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROELAST::PoroScatraMono::CombinedDBCMap() const
{
  return dbcmaps_->CondMap();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> POROELAST::PoroScatraMono::SystemMatrix()
{
  return systemmatrix_->Merge();
}

/*----------------------------------------------------------------------*
 | put the single maps to one full                                      |
 | Poroelasticity map together                              vuong 01/12 |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::SetDofRowMaps(
    const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);

  // full monolithic-blockmap
  blockrowdofmap_->Setup(*fullmap, maps);
}


/*----------------------------------------------------------------------*
 |  Evaluate off diagonal matrix in poro row                  vuong 08/13   |
 *---------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::EvaluateODBlockMatPoro()
{
  k_pfs_->Zero();

  // create the parameters for the discretization
  Teuchos::ParameterList fparams;
  // action for elements
  fparams.set<int>("action", FLD::calc_poroscatra_mono_odblock);
  // physical type
  fparams.set<int>("Physical Type", PoroField()->FluidField()->PhysicalType());

  // other parameters that might be needed by the elements
  fparams.set("delta time", Dt());
  fparams.set("total time", Time());

  const Teuchos::RCP<DRT::Discretization>& porofluiddis =
      PoroField()->FluidField()->Discretization();
  porofluiddis->ClearState();

  // set general vector values needed by elements
  porofluiddis->SetState(0, "dispnp", PoroField()->FluidField()->Dispnp());
  porofluiddis->SetState(0, "gridv", PoroField()->FluidField()->GridVel());
  porofluiddis->SetState(0, "veln", PoroField()->FluidField()->Veln());
  porofluiddis->SetState(0, "accnp", PoroField()->FluidField()->Accnp());
  porofluiddis->SetState(0, "hist", PoroField()->FluidField()->Hist());

  PoroField()->FluidField()->Discretization()->SetState(
      0, "scaaf", PoroField()->FluidField()->Scaaf());

  porofluiddis->SetState(0, "velaf", PoroField()->FluidField()->Velnp());
  porofluiddis->SetState(0, "velnp", PoroField()->FluidField()->Velnp());

  // build specific assemble strategy for the fluid-mechanical system matrix
  // from the point of view of FluidField:
  // fluiddofset = 0, structdofset = 1
  DRT::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
      2,                                  // scatradofset for column
      k_pfs_,                             // scatra-mechanical matrix
      Teuchos::null,                      // no other matrix or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate the fluid-mechanical system matrix on the fluid element
  PoroField()->FluidField()->Discretization()->EvaluateCondition(
      fparams, fluidstrategy, "PoroCoupling");
  // PoroField()->FluidField()->Discretization()->Evaluate( fparams, fluidstrategy );

  porofluiddis->ClearState(true);

  //************************************************************************************
  //************************************************************************************

  k_pss_->Zero();


  // create the parameters for the discretization
  Teuchos::ParameterList sparams;

  const std::string action = "struct_poro_calc_scatracoupling";
  sparams.set("action", action);
  // other parameters that might be needed by the elements
  sparams.set("delta time", Dt());
  sparams.set("total time", Time());

  PoroField()->StructureField()->Discretization()->ClearState();
  PoroField()->StructureField()->Discretization()->SetState(
      0, "displacement", PoroField()->StructureField()->Dispnp());
  PoroField()->StructureField()->Discretization()->SetState(
      0, "velocity", PoroField()->StructureField()->Velnp());

  // build specific assemble strategy for mechanical-fluid system matrix
  // from the point of view of StructureField:
  // structdofset = 0, fluiddofset = 1
  DRT::AssembleStrategy structuralstrategy(0,  // structdofset for row
      2,                                       // scatradofset for column
      k_pss_,                                  // mechanical-scatra coupling matrix
      Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate the mechanical-fluid system matrix on the structural element
  PoroField()->StructureField()->Discretization()->EvaluateCondition(
      sparams, structuralstrategy, "PoroCoupling");
  // StructureField()->Discretization()->Evaluate( sparams, structuralstrategy);
  PoroField()->StructureField()->Discretization()->ClearState(true);

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate off diagonal matrix in scatra row                    |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::EvaluateODBlockMatScatra()
{
  // create the parameters for the discretization
  Teuchos::ParameterList sparams_struct;

  k_sps_->Zero();

  sparams_struct.set<int>("action", SCATRA::calc_scatra_mono_odblock_mesh);
  // other parameters that might be needed by the elements
  sparams_struct.set("delta time", Dt());
  sparams_struct.set("total time", Time());

  // provide element parameter list with numbers of dofsets associated with displacement and
  // velocity dofs on scatra discretization
  sparams_struct.set<int>("ndsdisp", ScaTraField()->NdsDisp());
  sparams_struct.set<int>("ndsvel", ScaTraField()->NdsVel());

  ScaTraField()->Discretization()->ClearState();
  ScaTraField()->Discretization()->SetState(0, "hist", ScaTraField()->Hist());
  ScaTraField()->Discretization()->SetState(0, "phinp", ScaTraField()->Phinp());

  // build specific assemble strategy for mechanical-fluid system matrix
  // from the point of view of StructureField:
  // structdofset = 0, fluiddofset = 1
  DRT::AssembleStrategy scatrastrategy_struct(0,  // scatradofset for row
      1,                                          // structuredofset for column
      k_sps_,                                     // scatra-structure coupling matrix
      Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate the mechanical-fluid system matrix on the structural element
  ScaTraField()->Discretization()->EvaluateCondition(
      sparams_struct, scatrastrategy_struct, "PoroCoupling");
  // StructureField()->Discretization()->Evaluate( sparams, structuralstrategy);

  ScaTraField()->Discretization()->ClearState();
  // dserror("stop");
  //************************************************************************************
  //************************************************************************************

  k_spf_->Zero();

  // create the parameters for the discretization
  Teuchos::ParameterList sparams_fluid;

  sparams_fluid.set<int>("action", SCATRA::calc_scatra_mono_odblock_fluid);
  // other parameters that might be needed by the elements
  sparams_fluid.set("delta time", Dt());
  sparams_fluid.set("total time", Time());

  // provide element parameter list with numbers of dofsets associated with displacement and
  // velocity dofs on scatra discretization
  sparams_fluid.set<int>("ndsdisp", ScaTraField()->NdsDisp());
  sparams_fluid.set<int>("ndsvel", ScaTraField()->NdsVel());

  ScaTraField()->Discretization()->ClearState();
  ScaTraField()->Discretization()->SetState(0, "hist", ScaTraField()->Hist());
  ScaTraField()->Discretization()->SetState(0, "phinp", ScaTraField()->Phinp());

  // build specific assemble strategy for mechanical-fluid system matrix
  // from the point of view of StructureField:
  // structdofset = 0, fluiddofset = 1
  DRT::AssembleStrategy scatrastrategy_fluid(0,  // scatradofset for row
      2,                                         // fluiddofset for column
      k_spf_,                                    // scatra-fluid coupling matrix
      Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  // evaluate the mechanical-fluid system matrix on the structural element
  ScaTraField()->Discretization()->EvaluateCondition(
      sparams_fluid, scatrastrategy_fluid, "PoroCoupling");
  // StructureField()->Discretization()->Evaluate( sparams, structuralstrategy);
  ScaTraField()->Discretization()->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 |  check tangent stiffness matrix vie finite differences      vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraMono::FDCheck()
{
  std::cout << "\n******************finite difference check***************" << std::endl;

  int dof_struct = (PoroField()->StructureField()->Discretization()->NumGlobalNodes()) * 3;
  int dof_fluid = (PoroField()->FluidField()->Discretization()->NumGlobalNodes()) * 4;
  int dof_scatra = (ScaTraField()->Discretization()->NumGlobalNodes());

  std::cout << "structure field has " << dof_struct << " DOFs" << std::endl;
  std::cout << "fluid field has " << dof_fluid << " DOFs" << std::endl;
  std::cout << "scatra field has " << dof_scatra << " DOFs" << std::endl;

  Teuchos::RCP<Epetra_Vector> iterinc = Teuchos::null;
  iterinc = LINALG::CreateVector(*DofRowMap(), true);

  const int dofs = iterinc->GlobalLength();
  std::cout << "in total " << dofs << " DOFs" << std::endl;
  const double delta = 1e-8;

  iterinc->PutScalar(0.0);

  iterinc->ReplaceGlobalValue(0, 0, delta);

  Teuchos::RCP<Epetra_CrsMatrix> stiff_approx = Teuchos::null;
  stiff_approx = LINALG::CreateMatrix(*DofRowMap(), 81);

  Teuchos::RCP<Epetra_Vector> rhs_old = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));
  rhs_old->Update(1.0, *rhs_, 0.0);
  Teuchos::RCP<Epetra_Vector> rhs_copy = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();
  Teuchos::RCP<LINALG::SparseMatrix> sparse_copy =
      Teuchos::rcp(new LINALG::SparseMatrix(sparse->EpetraMatrix(), LINALG::Copy));

  if (false)
  {
    std::cout << "iterinc_" << std::endl << *iterinc_ << std::endl;
    std::cout << "iterinc" << std::endl << *iterinc << std::endl;
    std::cout << "meshdisp: " << std::endl << *(PoroField()->FluidField()->Dispnp());
    std::cout << "disp: " << std::endl << *(PoroField()->StructureField()->Dispnp());
    std::cout << "fluid vel" << std::endl << *(PoroField()->FluidField()->Velnp());
    std::cout << "fluid acc" << std::endl << *(PoroField()->FluidField()->Accnp());
    std::cout << "gridvel fluid" << std::endl << *(PoroField()->FluidField()->GridVel());
    std::cout << "gridvel struct" << std::endl << *(PoroField()->StructureField()->Velnp());
  }

  const int zeilennr = -1;
  const int spaltenr = -1;
  for (int i = 0; i < dofs; ++i)
  {
    if (CombinedDBCMap()->MyGID(i))
    {
      iterinc->ReplaceGlobalValue(i, 0, 0.0);
    }

    if (i == spaltenr)
      std::cout << "\n******************" << spaltenr + 1 << ". Spalte!!***************"
                << std::endl;

    Evaluate(iterinc);
    SetupRHS();

    rhs_copy->Update(1.0, *rhs_, 0.0);

    iterinc_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(
        sparse_copy, iterinc_, rhs_copy, Teuchos::null, zeros_, *CombinedDBCMap());


    if (i == spaltenr)
    {
      std::cout << "rhs_: " << (*rhs_copy)[zeilennr] << std::endl;
      std::cout << "rhs_old: " << (*rhs_old)[zeilennr] << std::endl;
    }

    rhs_copy->Update(-1.0, *rhs_old, 1.0);
    rhs_copy->Scale(-1.0 / delta);

    int* index = &i;
    for (int j = 0; j < dofs; ++j)
    {
      double value = (*rhs_copy)[j];
      stiff_approx->InsertGlobalValues(j, 1, &value, index);

      if ((j == zeilennr) and (i == spaltenr))
      {
        std::cout << "\n******************" << zeilennr + 1 << ". Zeile!!***************"
                  << std::endl;
        std::cout << "iterinc_" << std::endl << *iterinc_ << std::endl;
        std::cout << "iterinc" << std::endl << *iterinc << std::endl;
        std::cout << "meshdisp: " << std::endl << *(PoroField()->FluidField()->Dispnp());
        std::cout << "meshdisp scatra: " << std::endl
                  << *(ScaTraField()->Discretization()->GetState(
                         ScaTraField()->NdsDisp(), "dispnp"));
        std::cout << "disp: " << std::endl << *(PoroField()->StructureField()->Dispnp());
        std::cout << "fluid vel" << std::endl << *(PoroField()->FluidField()->Velnp());
        std::cout << "scatra vel" << std::endl
                  << *(ScaTraField()->Discretization()->GetState(
                         ScaTraField()->NdsVel(), "velocity field"));
        std::cout << "fluid acc" << std::endl << *(PoroField()->FluidField()->Accnp());
        std::cout << "gridvel fluid" << std::endl << *(PoroField()->FluidField()->GridVel());
        std::cout << "gridvel struct" << std::endl << *(PoroField()->StructureField()->Velnp());

        std::cout << "stiff_apprx(" << zeilennr << "," << spaltenr << "): " << (*rhs_copy)[zeilennr]
                  << std::endl;

        std::cout << "value(" << zeilennr << "," << spaltenr << "): " << value << std::endl;
        std::cout << "\n******************" << zeilennr + 1 << ". Zeile Ende!!***************"
                  << std::endl;
      }
    }

    if (not CombinedDBCMap()->MyGID(i)) iterinc->ReplaceGlobalValue(i, 0, -delta);

    iterinc->ReplaceGlobalValue(i - 1, 0, 0.0);

    if (i != dofs - 1) iterinc->ReplaceGlobalValue(i + 1, 0, delta);

    if (i == spaltenr)
      std::cout << "\n******************" << spaltenr + 1 << ". Spalte Ende!!***************"
                << std::endl;
  }

  Evaluate(iterinc);
  SetupRHS();

  stiff_approx->FillComplete();

  Teuchos::RCP<LINALG::SparseMatrix> stiff_approx_sparse = Teuchos::null;
  stiff_approx_sparse = Teuchos::rcp(new LINALG::SparseMatrix(stiff_approx, LINALG::Copy));

  stiff_approx_sparse->Add(*sparse_copy, false, -1.0, 1.0);

  Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = sparse_copy->EpetraMatrix();

  Teuchos::RCP<Epetra_CrsMatrix> error_crs = stiff_approx_sparse->EpetraMatrix();

  error_crs->FillComplete();
  sparse_crs->FillComplete();

  bool success = true;
  double error_max_rel = 0.0;
  double error_max_abs = 0.0;
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
                i, errorlength, errornumentries, &errorvalues[0], &errorindices[0]);
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
                i, sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);
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
                i, approxlength, approxnumentries, &approxvalues[0], &approxindices[0]);
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

          if (abs(error) > abs(error_max_rel)) error_max_rel = abs(error);
          if (abs(error_ij) > abs(error_max_abs)) error_max_abs = abs(error_ij);

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
    std::cout << "finite difference check successful, max. rel. error: " << error_max_rel
              << " , max. abs. error: " << error_max_abs << std::endl;
    std::cout << "******************finite difference check done***************\n\n" << std::endl;
  }
  else
    dserror("PoroFDCheck failed");

  return;
}
