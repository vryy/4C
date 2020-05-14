/*----------------------------------------------------------------------*/
/*! \file
 \brief two-way coupled monolithic solution algorithm
        for porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#include "poromultiphase_monolithic_twoway.H"
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_adapter/ad_porofluidmultiphase_wrapper.H"
#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/ad_art_net.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/inpar_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_utils_sparse_algebra_print.H"
#include "../linalg/linalg_solver.H"
#include "../drt_lib/drt_locsys.H"


#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_elements_paramsminimal.H"

#include "poromultiphase_utils.H"


/*----------------------------------------------------------------------*
 | constructor                                          kremheller 03/17 |
 *----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::PoroMultiPhaseMonolithicTwoWay(
    const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : PoroMultiPhaseMonolithic(comm, globaltimeparams),
      directsolve_(true),
      ittolinc_(0.0),
      ittolres_(0.0),
      itmax_(0),
      itmin_(1),
      itnum_(0),
      solveradaptolbetter_(0.0),
      solveradapttol_(false),
      blockrowdofmap_(Teuchos::null),
      equilibration_(INPAR::POROMULTIPHASE::EquilibrationMethods::equilibration_none),
      invrowsums_(Teuchos::null),
      tolinc_(0.0),
      tolfres_(0.0),
      tolinc_struct_(0.0),
      tolfres_struct_(0.0),
      tolinc_fluid_(0.0),
      tolfres_fluid_(0.0),
      normrhs_(0.0),
      norminc_(0.0),
      normrhsfluid_(0.0),
      normincfluid_(0.0),
      normrhsstruct_(0.0),
      normincstruct_(0.0),
      vectornormfres_(INPAR::POROMULTIPHASE::norm_undefined),
      vectornorminc_(INPAR::POROMULTIPHASE::norm_undefined),
      timernewton_(comm),
      dtsolve_(0.0),
      dtele_(0.0),
      fdcheck_(INPAR::POROMULTIPHASE::FDCheck::fdcheck_none)
{
}

/*----------------------------------------------------------------------*
 | initialization                                       kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::Init(
    const Teuchos::ParameterList& globaltimeparams, const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& structparams, const Teuchos::ParameterList& fluidparams,
    const std::string& struct_disname, const std::string& fluid_disname, bool isale, int nds_disp,
    int nds_vel, int nds_solidpressure, int ndsporofluid_scatra,
    const std::map<int, std::set<int>>* nearbyelepairs)
{
  // call base class
  POROMULTIPHASE::PoroMultiPhaseMonolithic::Init(globaltimeparams, algoparams, structparams,
      fluidparams, struct_disname, fluid_disname, isale, nds_disp, nds_vel, nds_solidpressure,
      ndsporofluid_scatra, nearbyelepairs);

  // inform user that structure field will not be solved but displacements will just be set to zero
  if (not solve_structure_) PrintStructureDisabledInfo();

  // Get the parameters for the ConvergenceCheck
  itmax_ = algoparams.get<int>("ITEMAX");
  ittolres_ = algoparams.sublist("MONOLITHIC").get<double>("TOLRES_GLOBAL");
  ittolinc_ = algoparams.sublist("MONOLITHIC").get<double>("TOLINC_GLOBAL");

  blockrowdofmap_ = Teuchos::rcp(new LINALG::MultiMapExtractor);

  fdcheck_ = DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASE::FDCheck>(
      algoparams.sublist("MONOLITHIC"), "FDCHECK");

  equilibration_ = DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASE::EquilibrationMethods>(
      algoparams.sublist("MONOLITHIC"), "EQUILIBRATION");

  solveradaptolbetter_ = algoparams.sublist("MONOLITHIC").get<double>("ADAPTCONV_BETTER");
  solveradapttol_ =
      (DRT::INPUT::IntegralValue<int>(algoparams.sublist("MONOLITHIC"), "ADAPTCONV") == 1);
}

/*----------------------------------------------------------------------*
 | setup the system if necessary (called in poromultiphase_dyn.cpp)     |
 |                                                     kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::SetupSystem()
{
  // -------------------------------------------------------------create combined map
  SetupMaps();

  // check global map extractor
  blockrowdofmap_->CheckForValidMapExtractor();

  //-----------------------------------build map of global dofs with DBC
  BuildCombinedDBCMap();
  // -------------------------------------------------------------

  // initialize Poromultiphase-elasticity-systemmatrix_
  systemmatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *Extractor(), *Extractor(), 81, false, true));

  // Initialize rhs
  rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  k_sf_ = Teuchos::rcp(new LINALG::SparseMatrix(*(StructDofRowMap()), 81, true, true));
  k_fs_ = Teuchos::rcp(new LINALG::SparseMatrix(*(FluidDofRowMap()), 81, true, true));

  // perform initialization associated with equilibration of global system of equations
  switch (equilibration_)
  {
    case INPAR::POROMULTIPHASE::EquilibrationMethods::equilibration_none:
    {
      // do nothing
      break;
    }

    case INPAR::POROMULTIPHASE::EquilibrationMethods::equilibration_rows_full:
    case INPAR::POROMULTIPHASE::EquilibrationMethods::equilibration_rows_maindiag:
    {
      // initialize vector for row sums of global system matrix if necessary
      invrowsums_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), false));
      break;
    }

    default:
    {
      dserror("Equilibration method not yet implemented!");
      break;
    }
  }

  // StructureField: check whether we have locsys BCs, i.e. inclined structural
  //  Dirichlet BC

  std::vector<DRT::Condition*> locsysconditions(0);
  (StructureField()->Discretization())->GetCondition("Locsys", locsysconditions);

  // if there are inclined structural Dirichlet BC, get the structural LocSysManager
  if (locsysconditions.size()) locsysman_ = StructureField()->LocsysManager();

  return;
}

/*----------------------------------------------------------------------*
 | setup the map                                       kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::SetupMaps()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

  vecSpaces.push_back(StructDofRowMap());

  vecSpaces.push_back(FluidDofRowMap());

  if (vecSpaces[0]->NumGlobalElements() == 0) dserror("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements() == 0) dserror("No fluid equation. Panic.");

  // full Poromultiphase-elasticity-map
  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(vecSpaces);

  // full Poromultiphase-elasticity-blockmap
  blockrowdofmap_->Setup(*fullmap_, vecSpaces);

  return;
}

/*----------------------------------------------------------------------*
 | Monolithic Time Step                                kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::TimeStep()
{
  // Prepare stuff
  SetupNewton();
  PrintHeader();

  // Evaluate
  Evaluate(iterinc_);

  // Newton-Loop
  while ((not Converged() and itnum_ < itmax_) or (itnum_ < itmin_))
  {
    // increment number of iteration
    itnum_++;

    // Solve
    LinearSolve();
    solver_->ResetTolerance();

    // Build Convergence Norms
    BuildConvergenceNorms();

    if (not Converged())
    {
      // Evaluate
      Evaluate(iterinc_);

      // perform FD Check of monolithic system matrix
      if (fdcheck_ == INPAR::POROMULTIPHASE::fdcheck_global) PoroFDCheck();
    }
    else
    {
      // convergence check is based on residual(phi_i) < tol and phi_i+1 - phi_i < tol
      // in this function we update phi_i+1 as phi_i+1 = phi_i + iterinc for all fields
      // even though we have not evaluated the residual of phi_i+1 it will still be more exact than
      // the one at phi_i
      UpdateFieldsAfterConvergence();
    }

    // print output
    NewtonOutput();
  }

  // Error-Check
  NewtonErrorCheck();

  return;
}

/*-----------------------------------------------------------------------/
|  build the combined dbcmap                           kremheller 03/17  |
/-----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::BuildCombinedDBCMap()
{
  // get structure and fluid dbc maps
  const Teuchos::RCP<const Epetra_Map> scondmap = StructureField()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> fcondmap = FluidField()->GetDBCMapExtractor()->CondMap();
  // merge them
  combinedDBCMap_ = LINALG::MergeMap(scondmap, fcondmap, false);

  return;
}

/*----------------------------------------------------------------------*
 | Evaluate (build global Matrix and RHS)            kremheller 03/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::Evaluate");

  // reset timer
  timernewton_.ResetStartTime();
  // *********** time measurement ***********
  double dtcpu = timernewton_.WallTime();
  // *********** time measurement ***********

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  ExtractFieldVectors(x, sx, fx);

  Evaluate(sx, fx, itnum_ == 0);

  // *********** time measurement ***********
  dtele_ = timernewton_.WallTime() - dtcpu;
  // *********** time measurement ***********

  return;
}
/*----------------------------------------------------------------------*
 | Evaluate (build global Matrix and RHS, public --> allows access      |
 | from outside --> monolithic scatra-coupling)        kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::Evaluate(Teuchos::RCP<const Epetra_Vector> sx,
    Teuchos::RCP<const Epetra_Vector> fx, const bool firstcall)
{
  // (1) Update fluid Field and reconstruct pressures and saturations
  FluidField()->UpdateIter(fx);

  if (solve_structure_)
  {
    // (2) set fluid solution in structure field
    StructureField()->Discretization()->SetState(1, "porofluid", FluidField()->Phinp());

    // (3) evaluate structure
    if (firstcall)  // first call (iterinc_ = 0) --> sx = 0
      StructureField()->Evaluate();
    else  //(this call will also update displacements and velocities)
      StructureField()->Evaluate(sx);

    // (4) Set structure solution on fluid field
    SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());
  }
  else
  {
    // (4) Set structure solution on fluid field
    SetStructSolution(struct_zeros_, struct_zeros_);
    StructureField()->SystemMatrix()->Zero();
    StructureField()->SystemMatrix()->Complete(
        StructureField()->SystemMatrix()->RangeMap(), StructureField()->SystemMatrix()->RangeMap());
  }

  // (5) Evaluate the fluid
  FluidField()->Evaluate();

  // (6) Build the monolithic system matrix
  SetupSystemMatrix();

  // check whether we have a sanely filled tangent matrix
  if (not systemmatrix_->Filled())
  {
    dserror("Effective tangent matrix must be filled here");
  }

  // (7) Build the monolithic system vector
  SetupRHS();
}

/*----------------------------------------------------------------------*
 | setup system matrix of poromultiphase-elasticity   kremheller 03/17  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::SetupSystemMatrix(
    LINALG::BlockSparseMatrixBase& mat)
{
  // TEUCHOS_FUNC_TIME_MONITOR("POROELAST::Monolithic::SetupSystemMatrix");

  // pure structural part k_ss ((ndim*n_nodes)x(ndim*n_nodes))

  // build pure structural block k_ss
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<LINALG::SparseMatrix> k_ss = StructureField()->SystemMatrix();

  if (k_ss == Teuchos::null) dserror("structure system matrix null pointer!");

  // Copy from TSI
  if (locsysman_ != Teuchos::null)
  {
    // rotate k_ss to local coordinate system --> k_ss^{~}
    locsysman_->RotateGlobalToLocal(k_ss);
    // apply ApplyDirichletWithTrafo() on rotated block k_ss^{~}
    // --> if dof has an inclined DBC: blank the complete row, the '1.0' is set
    //     on diagonal of row, i.e. on diagonal of k_ss
    k_ss->ApplyDirichletWithTrafo(
        locsysman_->Trafo(), *StructureField()->GetDBCMapExtractor()->CondMap(), true);
  }  // end locsys
  // default: (locsysman_ == Teuchos::null), i.e. NO inclined Dirichlet BC
  else
    k_ss->ApplyDirichlet(*StructureField()->GetDBCMapExtractor()->CondMap(), true);

  /*----------------------------------------------------------------------*/
  // structural part k_sf ((ndim*n_nodes)x(n_phases*n_nodes))
  // build mechanical-fluid block

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_sf = StructFluidCouplingMatrix();

  // Epetra_Time timerstrcoupl(Comm());
  // call the element and calculate the matrix block
  ApplyStrCouplMatrix(k_sf);

  // Copy from TSI
  // apply dirichlet boundary conditions properly on matrix k_sf, i.e. blank row
  // if dof is a structural DBC
  // Normally, DBC should be applied on complete systemmatrix mat, but for
  // diagonal blocks (here k_ss, k_tt) DBC are ALREADY applied in
  // PrepareSystemForNewtonSolve() included in Evaluate(sx)
  //
  // to avoid double work, we only call ApplyDirichlet for the off-diagonal blocks,
  // here k_sf
  // k_sf is an off-diagonal block --> pass the bool diagonal==false
  // ApplyDirichlet*() expect filled matrix
  //
  // in case of inclined STR-DBC
  //   1.) transform the off-diagonal block k_sf to the local system --> k_st^{~}
  //   2.) apply ApplyDirichletWithTrafo() on rotated block k_sf^{~}
  //              --> blank the row, which has a DBC

  // to apply Multiply in LocSys, k_st has to be FillCompleted
  k_sf->Complete(
      FluidField()->SystemMatrix()->RangeMap(), StructureField()->SystemMatrix()->RangeMap());

  if (locsysman_ != Teuchos::null)
  {
    // rotate k_st to local coordinate system --> k_st^{~}
    locsysman_->RotateGlobalToLocal(k_sf);
    // apply ApplyDirichletWithTrafo() on rotated block k_st^{~}
    // --> if dof has an inclined DBC: blank the complete row, the '1.0' is set
    //     on diagonal of row, i.e. on diagonal of k_ss
    k_sf->ApplyDirichletWithTrafo(
        locsysman_->Trafo(), *StructureField()->GetDBCMapExtractor()->CondMap(), false);
  }  // end locsys
  // default: (locsysman_ == Teuchos::null), i.e. NO inclined Dirichlet BC
  else
    k_sf->ApplyDirichlet(*StructureField()->GetDBCMapExtractor()->CondMap(), false);

  /*----------------------------------------------------------------------*/
  // pure fluid part k_ff ( (n_phases*n_nodes)x(n_phases*n_nodes) )

  // build pure fluid block k_ff
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W
  // NOTE: DBC's have already been applied within Evaluate (PrepareSystemForNewtonSolve())
  Teuchos::RCP<LINALG::SparseMatrix> k_ff = FluidField()->SystemMatrix();

  if (k_ff == Teuchos::null) dserror("fluid system matrix null pointer!");

  /*----------------------------------------------------------------------*/
  // fluid part k_fs ( (n_phases*n_nodes)x(ndim*n_nodes) )
  // build fluid-mechanical block

  // create empty matrix
  Teuchos::RCP<LINALG::SparseMatrix> k_fs = FluidStructCouplingMatrix();

  // call the element and calculate the matrix block
  ApplyFluidCouplMatrix(k_fs);

  // apply DBC's also on off-diagonal fluid-structure coupling block
  k_fs->ApplyDirichlet(*FluidField()->GetDBCMapExtractor()->CondMap(), false);

  // uncomplete matrix block (appears to be required in certain cases (locsys+iterative solver))
  if (solve_structure_)
  {
    k_ss->UnComplete();
    k_sf->UnComplete();
  }
  k_fs->UnComplete();
  k_ff->UnComplete();

  // assign structure part to the Poroelasticity matrix
  mat.Assign(0, 0, LINALG::View, *k_ss);
  // assign coupling part to the Poroelasticity matrix
  mat.Assign(0, 1, LINALG::View, *k_sf);
  // assign fluid part to the poroelasticity matrix
  mat.Assign(1, 1, LINALG::View, *k_ff);
  // assign coupling part to the Poroelasticity matrix
  mat.Assign(1, 0, LINALG::View, *k_fs);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();

  // Debug: matlab output of system matrix
  bool matlab = false;
  if (matlab)
  {
    std::string filename = "../o/mymatrix.dat";
    LINALG::PrintBlockMatrixInMatlabFormat(filename, mat);
    dserror("exit");
  }

}  // SetupSystemMatrix

/*----------------------------------------------------------------------*
 | get fluid structure-coupling sparse matrix           kremheller 03/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix>
POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::FluidStructCouplingMatrix()
{
  Teuchos::RCP<LINALG::SparseMatrix> sparse =
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_fs_);
  if (sparse == Teuchos::null) dserror("cast to LINALG::SparseMatrix failed!");

  return sparse;
}  // FluidStructCouplingMatrix()

/*----------------------------------------------------------------------*
 | get structure fluid-coupling sparse matrix           kremheller 03/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix>
POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::StructFluidCouplingMatrix()
{
  Teuchos::RCP<LINALG::SparseMatrix> sparse =
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_sf_);
  if (sparse == Teuchos::null) dserror("cast to LINALG::SparseMatrix failed!");

  return sparse;
}  // FluidStructCouplingMatrix()

/*----------------------------------------------------------------------*
 | evaluate fluid-structural system matrix at state    kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::ApplyFluidCouplMatrix(
    Teuchos::RCP<LINALG::SparseOperator> k_fs  //!< off-diagonal tangent matrix term
)
{
  // reset
  k_fs->Zero();
  if (solve_structure_) FluidField()->AssembleFluidStructCouplingMat(k_fs);
  k_fs->Complete(
      StructureField()->SystemMatrix()->RangeMap(), FluidField()->SystemMatrix()->RangeMap());

  return;
}

/*-----------------------------------------------------------------------------*
 | update fields after convergence as phi_i+1=phi_i+iterinc   kremheller 07/17 |
 *-----------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::UpdateFieldsAfterConvergence()
{
  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  ExtractFieldVectors(iterinc_, sx, fx);

  UpdateFieldsAfterConvergence(sx, fx);

  return;
}

/*-----------------------------------------------------------------------------*
 | update fields after convergence as phi_i+1=phi_i+iterinc   kremheller 07/17 |
 *-----------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::UpdateFieldsAfterConvergence(
    Teuchos::RCP<const Epetra_Vector>& sx, Teuchos::RCP<const Epetra_Vector>& fx)
{
  // (1) Update fluid Field and reconstruct pressures and saturations
  FluidField()->UpdateIter(fx);
  FluidField()->ReconstructPressuresAndSaturations();
  FluidField()->ReconstructFlux();

  if (solve_structure_) StructureField()->Evaluate(sx);

  // (4) Set structure solution on fluid field
  SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());

  return;
}

/*----------------------------------------------------------------------*
 | evaluate structural-fluid system matrix at state    kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::ApplyStrCouplMatrix(
    Teuchos::RCP<LINALG::SparseOperator> k_sf  //!< off-diagonal tangent matrix term
)
{
  k_sf->Zero();

  if (solve_structure_)
  {
    // create the parameters for the discretization
    Teuchos::ParameterList sparams;

    //! pointer to the model evaluator data container
    Teuchos::RCP<DRT::ELEMENTS::ParamsMinimal> params =
        Teuchos::rcp(new DRT::ELEMENTS::ParamsMinimal());

    // set parameters needed for element evalutation
    params->SetActionType(DRT::ELEMENTS::struct_poro_calc_fluidcoupling);
    params->SetTotalTime(Time());
    params->SetDeltaTime(Dt());
    // std::cout << Dt() << std::endl;

    sparams.set<Teuchos::RCP<DRT::ELEMENTS::ParamsInterface>>("interface", params);

    StructureField()->Discretization()->ClearState();
    StructureField()->Discretization()->SetState(0, "displacement", StructureField()->Dispnp());
    StructureField()->Discretization()->SetState(0, "velocity", StructureField()->Velnp());
    StructureField()->Discretization()->SetState(1, "porofluid", FluidField()->Phinp());

    // build specific assemble strategy for mechanical-fluid system matrix
    // from the point of view of StructureField:
    // structdofset = 0, fluiddofset = 1
    DRT::AssembleStrategy structuralstrategy(0,  // structdofset for row
        1,                                       // fluiddofset for column
        k_sf,                                    // mechanical-fluid coupling matrix
        Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

    // evaluate the mechanical-fluid system matrix on the structural element
    StructureField()->Discretization()->Evaluate(sparams, structuralstrategy);

    StructureField()->Discretization()->ClearState();

    // scale with time integration factor
    k_sf->Scale(1.0 - StructureField()->TimIntParam());
  }

  return;
}
/*----------------------------------------------------------------------*
 | setup solver for monolithic problem                kremheller 03/17  |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::SetupSolver()
{
  //  solver
  // create a linear solver
  // get dynamic section of poroelasticity
  const Teuchos::ParameterList& poromultdyn =
      DRT::Problem::Instance()->PoroMultiPhaseDynamicParams();
  // get the solver number used for linear poroelasticity solver
  const int linsolvernumber = poromultdyn.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER");
  // check if the poroelasticity solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for poromultiphaseflow. Please set LINEAR_SOLVER in "
        "POROMULTIPHASE DYNAMIC to a valid number!");
  const Teuchos::ParameterList& solverparams =
      DRT::Problem::Instance()->SolverParams(linsolvernumber);
  const int solvertype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(solverparams, "SOLVER");

  directsolve_ = (solvertype == INPAR::SOLVER::umfpack or solvertype == INPAR::SOLVER::superlu or
                  solvertype == INPAR::SOLVER::amesos_klu_nonsym);

  if (directsolve_)
  {
    solver_ = Teuchos::rcp(
        new LINALG::Solver(solverparams, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));
  }
  else
    CreateLinearSolver(solverparams, solvertype);

  vectornormfres_ = DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASE::VectorNorm>(
      poromultdyn.sublist("MONOLITHIC"), "VECTORNORM_RESF");
  vectornorminc_ = DRT::INPUT::IntegralValue<INPAR::POROMULTIPHASE::VectorNorm>(
      poromultdyn.sublist("MONOLITHIC"), "VECTORNORM_INC");

  return true;
}

/*----------------------------------------------------------------------*
 | Create linear (iterative) solver                  kremheller 08/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::CreateLinearSolver(
    const Teuchos::ParameterList& solverparams, const int solvertype)
{
  if (solvertype != INPAR::SOLVER::aztec_msr && solvertype != INPAR::SOLVER::belos)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " Note: the BGS2x2 preconditioner now " << std::endl;
    std::cout << " uses the structural solver and fluid solver blocks" << std::endl;
    std::cout << " for building the internal inverses" << std::endl;
    std::cout << " Remove the old BGS PRECONDITIONER BLOCK entries " << std::endl;
    std::cout << " in the dat files!" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    dserror("aztec solver expected");
  }
  const int azprectype =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(solverparams, "AZPREC");

  // plausibility check
  switch (azprectype)
  {
    case INPAR::SOLVER::azprec_AMGnxn:
    {
      // no plausibility checks here
      // if you forget to declare an xml file you will get an error message anyway
    }
    break;
    default:
      dserror("AMGnxn preconditioner expected");
      break;
  }

  solver_ = Teuchos::rcp(
      new LINALG::Solver(solverparams, Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

  // build the null spaces of the single blocks
  BuildBlockNullSpaces(solver_);
}

/*-----------------------------------------------------------------------------------*
 | build null spaces associated with blocks of global system matrix kremheller 08/17 |
 *-----------------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::BuildBlockNullSpaces(
    Teuchos::RCP<LINALG::Solver>& solver)
{
  // equip smoother for structure matrix block with empty parameter sublists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams = solver->Params().sublist("Inverse1");
  blocksmootherparams.sublist("Aztec Parameters");
  blocksmootherparams.sublist("MueLu Parameters");

  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

  // equip smoother for fluid matrix block with empty parameter sublists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams2 = solver->Params().sublist("Inverse2");
  blocksmootherparams2.sublist("Aztec Parameters");
  blocksmootherparams2.sublist("MueLu Parameters");

  FluidField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams2);

  return;
}

/*----------------------------------------------------------------------*
 | Setup Newton-Raphson iteration                    kremheller 03/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::SetupNewton()
{
  // initialise equilibrium loop and norms
  itnum_ = 0;
  normrhs_ = 0.0;
  norminc_ = 0.0;
  normrhsfluid_ = 0.0;
  normincfluid_ = 0.0;
  normrhsstruct_ = 0.0;
  normincstruct_ = 0.0;
  tolinc_ = 0.0;
  tolfres_ = 0.0;
  tolinc_struct_ = 0.0;
  tolfres_struct_ = 0.0;
  tolinc_fluid_ = 0.0;
  tolfres_fluid_ = 0.0;

  // incremental solution vector with length of all dofs
  if (iterinc_ == Teuchos::null)
    iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  else
    iterinc_->PutScalar(0.0);

  // a zero vector of full length
  if (zeros_ == Teuchos::null)
    zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  else
    zeros_->PutScalar(0.0);

  // AitkenReset();

  return;
}

/*----------------------------------------------------------------------*
 | Print Header                                      kremheller 03/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::PrintHeader()
{
  if (Comm().MyPID() == 0)
  {
    if (!solve_structure_) PrintStructureDisabledInfo();
    std::cout << "+--------------------------------------------------------------------------------"
                 "----------+"
              << std::endl;
    std::cout << "| MONOLITHIC POROMULTIPHASE SOLVER                                               "
                 "          |"
              << std::endl;
    std::cout << "| STEP: " << std::setw(5) << std::setprecision(4) << std::scientific << Step()
              << "/" << std::setw(5) << std::setprecision(4) << std::scientific << NStep()
              << ", Time: " << std::setw(11) << std::setprecision(4) << std::scientific << Time()
              << "/" << std::setw(11) << std::setprecision(4) << std::scientific << MaxTime()
              << ", Dt: " << std::setw(11) << std::setprecision(4) << std::scientific << Dt()
              << "                        |" << std::endl;
  }
}

/*----------------------------------------------------------------------*
 | Build necessary norms                               kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::BuildConvergenceNorms()
{
  //------------------------------------------------------------ build residual force norms
  normrhs_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_);
  Teuchos::RCP<const Epetra_Vector> rhs_s;
  Teuchos::RCP<const Epetra_Vector> rhs_f;

  // get structure and fluid RHS
  ExtractFieldVectors(rhs_, rhs_s, rhs_f);

  // build also norms for fluid and structure
  normrhsstruct_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_s);
  normrhsfluid_ = UTILS::CalculateVectorNorm(vectornormfres_, rhs_f);

  //------------------------------------------------------------- build residual increment norms
  norminc_ = UTILS::CalculateVectorNorm(vectornorminc_, iterinc_);

  // displacement and fluid velocity & pressure incremental vector
  Teuchos::RCP<const Epetra_Vector> iterincs;
  Teuchos::RCP<const Epetra_Vector> iterincf;

  // get structure and fluid increment
  ExtractFieldVectors(iterinc_, iterincs, iterincf);

  // build also norms for fluid and structure
  normincstruct_ = UTILS::CalculateVectorNorm(vectornorminc_, iterincs);
  normincfluid_ = UTILS::CalculateVectorNorm(vectornorminc_, iterincf);

  Teuchos::RCP<Epetra_Vector> sol_vec = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));
  SetupVector(*sol_vec, StructureField()->Dispnp(), FluidField()->Phinp());

  double dispnorm = UTILS::CalculateVectorNorm(vectornorminc_, (StructureField()->Dispnp()));
  double fluidnorm = UTILS::CalculateVectorNorm(vectornorminc_, (FluidField()->Phinp()));
  double totalnorm = UTILS::CalculateVectorNorm(vectornorminc_, sol_vec);

  // take care of very small norms
  if (dispnorm < 1.0e-6) dispnorm = 1.0;
  if (fluidnorm < 1.0e-6) fluidnorm = 1.0;
  if (totalnorm < 1.0e-6) totalnorm = 1.0;

  // build relative increment norm
  normincstruct_ /= dispnorm;
  normincfluid_ /= fluidnorm;
  norminc_ /= totalnorm;


  return;
}

/*----------------------------------------------------------------------*
 | Newton Output (adapted form tsi)                    kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::NewtonOutput()
{
  // print the incremental based convergence check to the screen
  if (Comm().MyPID() == 0)
  {
    if (itnum_ == 1)
      printf(
          "+--------------+----------------+------------------+--------------------+---------------"
          "---+\n");
    printf(
        "|-  step/max  -|-  total-inc   -|-  fluid-inc     -|-  disp-inc        -|-  norm-rhs      "
        "-| (ts =%10.3E,",
        dtsolve_);
    printf("\n");
    printf(
        "|   %3d/%3d    |  %10.3E    |  %10.3E      |  %10.3E        |  %10.3E      |  te =%10.3E)",
        itnum_, itmax_, norminc_, normincfluid_, normincstruct_, normrhs_, dtele_);
    printf("\n");
    printf(
        "+--------------+----------------+------------------+--------------------+-----------------"
        "-+\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 | Error-Check and final output                        kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::NewtonErrorCheck()
{
  // build the maximum value of the residuals and increments
  const double maxinc = std::max(norminc_, std::max(normincfluid_, normincstruct_));
  const double maxres = std::max(normrhs_, std::max(normrhsfluid_, normrhsstruct_));

  // print the incremental based convergence check to the screen
  if (Converged())  // norminc_ < ittol_ && normrhs_ < ittol_ && normincfluid_ < ittol_ &&
                    // normincstruct_ < ittol_
  {
    if (Comm().MyPID() == 0)
    {
      printf(
          "|  Monolithic iteration loop converged after iteration %3d/%3d !                        "
          "   |\n",
          itnum_, itmax_);
      printf(
          "|  Quantity           [norm]:                 TOL                                       "
          "   |\n");
      printf(
          "|  Max. rel. increment [%3s]:  %10.3E  < %10.3E                                    |\n",
          VectorNormString(vectornorminc_).c_str(), maxinc, ittolinc_);
      printf(
          "|  Maximum    residual [%3s]:  %10.3E  < %10.3E                                    |\n",
          VectorNormString(vectornormfres_).c_str(), maxres, ittolres_);
      printf(
          "+--------------+----------------+------------------+--------------------+---------------"
          "---+\n");
      printf("\n");
    }
  }
  else
  {
    if ((Comm().MyPID() == 0))
    {
      printf(
          "|     >>>>>> not converged in %3d steps!                                                "
          "   |\n",
          itmax_);
      printf(
          "|  Max. rel. increment [%3s]:  %10.3E    %10.3E                                    |\n",
          VectorNormString(vectornorminc_).c_str(), maxinc, ittolinc_);
      printf(
          "|  Maximum    residual [%3s]:  %10.3E    %10.3E                                    |\n",
          VectorNormString(vectornormfres_).c_str(), maxres, ittolres_);
      printf(
          "+--------------+----------------+------------------+--------------------+---------------"
          "---+\n");
      printf("\n");
      printf("\n");
    }
    dserror("The monolithic solver did not converge in ITEMAX steps!");
  }


  return;
}

/*----------------------------------------------------------------------*
 | simple convergence check                            kremheller 03/17 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::Converged()
{
  return (norminc_ < ittolinc_ && normincfluid_ < ittolinc_ && normincstruct_ < ittolinc_ &&
          normrhs_ < ittolres_ && normrhsstruct_ < ittolres_ && normrhsfluid_ < ittolres_);
}

/*----------------------------------------------------------------------*
 | Solve linear Poromultiphase-elasticity system     kremheller 03/17   |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::LinearSolve()
{
  // reset timer
  timernewton_.ResetStartTime();
  // *********** time measurement ***********
  double dtcpu = timernewton_.WallTime();
  // *********** time measurement ***********

  if (solveradapttol_ and (itnum_ > 1))
  {
    double worst = std::max(normrhs_, norminc_);
    double wanted = tolfres_;
    solver_->AdaptTolerance(wanted, worst, solveradaptolbetter_);
  }
  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more

  // equilibrate global system of equations if necessary
  EquilibrateSystem();

  if (directsolve_)
  {
    // merge blockmatrix to SparseMatrix
    Teuchos::RCP<LINALG::SparseMatrix> sparse = systemmatrix_->Merge();

    // standard solver call
    // system is ready to solve since Dirichlet Boundary conditions have been applied in
    // SetupSystemMatrix or Evaluate
    solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, true, itnum_ == 1);
  }
  else
  {
    // standard solver call
    solver_->Solve(systemmatrix_->EpetraOperator(), iterinc_, rhs_, true, itnum_ == 1);
  }

  // *********** time measurement ***********
  dtsolve_ = timernewton_.WallTime() - dtcpu;
  // *********** time measurement ***********

  return;
}

/*----------------------------------------------------------------------*
 | equilibrate global system of equations if necessary kremheller 08/17 |
 | (adapated from sti_algorithm.cpp)                                    |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::EquilibrateSystem()
{
  switch (equilibration_)
  {
    case INPAR::POROMULTIPHASE::EquilibrationMethods::equilibration_none:
    {
      // do nothing
      break;
    }

    case INPAR::POROMULTIPHASE::EquilibrationMethods::equilibration_rows_full:
    case INPAR::POROMULTIPHASE::EquilibrationMethods::equilibration_rows_maindiag:
    {
      // perform row equilibration
      for (int i = 0; i < systemmatrix_->Rows(); ++i)
      {
        // initialize vector for inverse row sums
        const Teuchos::RCP<Epetra_Vector> invrowsums(
            Teuchos::rcp(new Epetra_Vector(systemmatrix_->Matrix(i, i).RowMap())));

        // compute inverse row sums of current main diagonal matrix block
        if (equilibration_ ==
            INPAR::POROMULTIPHASE::EquilibrationMethods::equilibration_rows_maindiag)
          ComputeInvRowSums(systemmatrix_->Matrix(i, i), invrowsums);

        // compute inverse row sums of current row block of global system matrix
        else
        {
          // loop over all column blocks of global system matrix
          for (int j = 0; j < systemmatrix_->Cols(); ++j)
          {
            // extract current block of global system matrix
            const LINALG::SparseMatrix& matrix = systemmatrix_->Matrix(i, j);

            // loop over all rows of current matrix block
            for (int irow = 0; irow < matrix.RowMap().NumMyElements(); ++irow)
            {
              // determine length of current matrix row
              const int length = matrix.EpetraMatrix()->NumMyEntries(irow);

              if (length > 0)
              {
                // extract current matrix row from matrix block
                int numentries(0);
                std::vector<double> values(length, 0.);
                if (matrix.EpetraMatrix()->ExtractMyRowCopy(irow, length, numentries, &values[0]))
                  dserror("Cannot extract matrix row with local ID %d from matrix block!", irow);

                // compute and store current row sum
                double rowsum(0.);
                for (int ientry = 0; ientry < numentries; ++ientry)
                  rowsum += std::abs(values[ientry]);
                (*invrowsums)[irow] += rowsum;
              }
            }
          }

          // invert row sums
          if (invrowsums->Reciprocal(*invrowsums)) dserror("Vector could not be inverted!");
        }

        // perform row equilibration of matrix blocks in current row block of global system matrix
        for (int j = 0; j < systemmatrix_->Cols(); ++j)
          EquilibrateMatrixRows(systemmatrix_->Matrix(i, j), invrowsums);

        // insert inverse row sums of current main diagonal matrix block into global vector
        Extractor()->InsertVector(invrowsums, i, invrowsums_);
      }

      // perform equilibration of global residual vector
      if (rhs_->Multiply(1., *invrowsums_, *rhs_, 0.))
        dserror("Equilibration of global residual vector failed!");

      break;
    }

    default:
    {
      dserror("Equilibration method not yet implemented!");
      break;
    }
  }

  return;
}

/*--------------------------------------------------------------------------------*
 | compute inverse sums of absolute values of matrix row entries kremheller 08/17 |
 | (copy from sti_algorithm.cpp)                                                  |
 *--------------------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::ComputeInvRowSums(
    const LINALG::SparseMatrix& matrix,  //!< matrix
    const Teuchos::RCP<Epetra_Vector>&
        invrowsums  //!< inverse sums of absolute values of row entries in matrix
    ) const
{
  // compute inverse row sums of matrix
  if (matrix.EpetraMatrix()->InvRowSums(*invrowsums))
    dserror("Inverse row sums of matrix could not be successfully computed!");

  return;
}  // POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::ComputeInvRowSums

/*----------------------------------------------------------------------*
 | equilibrate matrix rows                             kremheller 08/17 |
 | (copy from sti_algorithm.cpp)                                        |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::EquilibrateMatrixRows(
    LINALG::SparseMatrix& matrix,  //!< matrix
    const Teuchos::RCP<Epetra_Vector>&
        invrowsums  //!< sums of absolute values of row entries in matrix
    ) const
{
  if (matrix.LeftScale(*invrowsums)) dserror("Row equilibration of matrix failed!");

  return;
}  // POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::EquilibrateMatrixRows

/*----------------------------------------------------------------------*
 | get the dof row map                                 kremheller 03/17 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::DofRowMap()
{
  return blockrowdofmap_->FullMap();
}

/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                             kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::SetupRHS()
{
  // get structure part
  Teuchos::RCP<Epetra_Vector> str_rhs = SetupStructurePartofRHS();

  // fill the Poroelasticity rhs vector rhs_ with the single field rhss
  SetupVector(*rhs_, str_rhs, FluidField()->RHS());

}  // SetupRHS()

/*----------------------------------------------------------------------*
 | setup RHS (like fsimon)                             kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::SetupStructurePartofRHS()
{
  // Copy from TSI
  Teuchos::RCP<Epetra_Vector> str_rhs = struct_zeros_;
  if (solve_structure_) str_rhs = Teuchos::rcp(new Epetra_Vector(*StructureField()->RHS()));
  if (locsysman_ != Teuchos::null) locsysman_->RotateGlobalToLocal(str_rhs);

  return str_rhs;
}

/*----------------------------------------------------------------------*
 | setup vector of the structure and fluid field         kremheller 03/17|
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::SetupVector(
    Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> sv, Teuchos::RCP<const Epetra_Vector> fv)
{
  Extractor()->InsertVector(*sv, 0, f);

  f.Scale(-1);
  Extractor()->InsertVector(*fv, 1, f);
}
/*----------------------------------------------------------------------*
 | extract field vectors for calling Evaluate() of the  kremheller 03/17|
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::ExtractFieldVectors(
    Teuchos::RCP<const Epetra_Vector> x, Teuchos::RCP<const Epetra_Vector>& sx,
    Teuchos::RCP<const Epetra_Vector>& fx)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::ExtractFieldVectors");

  // process structure unknowns of the first field
  sx = Extractor()->ExtractVector(x, 0);

  // process fluid unknowns of the second field
  fx = Extractor()->ExtractVector(x, 1);
}

/*----------------------------------------------------------------------*
 | inform user that structure is not solved            kremheller 08/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::PrintStructureDisabledInfo()
{
  // print out Info
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n";
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 "++++++++++++++++++++++++++++++++\n";
    std::cout << " INFO:    STRUCTURE FIELD IS NOT SOLVED; MAKE SURE YOU HAVE CONSTRAINED ALL DOFS "
                 "IN YOUR STRUCTURE WITH A DBC\n";
    std::cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                 "++++++++++++++++++++++++++++++++\n";
  }
}

/*----------------------------------------------------------------------*
 |  check tangent stiffness matrix via finite differences     vuong 01/12 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWay::PoroFDCheck()
{
  std::cout << "\n******************finite difference check***************" << std::endl;

  int dof_struct = (StructureField()->DofRowMap()->NumGlobalElements());
  int dof_fluid = (FluidField()->DofRowMap()->NumGlobalElements());

  std::cout << "structure field has " << dof_struct << " DOFs" << std::endl;
  std::cout << "fluid field has " << dof_fluid << " DOFs" << std::endl;

  Teuchos::RCP<Epetra_Vector> iterinc = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> abs_iterinc = Teuchos::null;
  iterinc = LINALG::CreateVector(*DofRowMap(), true);
  abs_iterinc = LINALG::CreateVector(*DofRowMap(), true);

  const int dofs = iterinc->GlobalLength();
  std::cout << "in total " << dofs << " DOFs" << std::endl;
  const double delta = 1e-8;

  iterinc->PutScalar(0.0);

  iterinc->ReplaceGlobalValue(0, 0, delta);

  abs_iterinc->Update(1.0, *iterinc_, 0.0);

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
    // std::cout << "meshdisp: " << std::endl << *(FluidField()->Dispnp());
    std::cout << "disp: " << std::endl << *(StructureField()->Dispnp());
    // std::cout << "fluid vel" << std::endl << *(FluidField()->Velnp());
    // std::cout << "fluid acc" << std::endl << *(FluidField()->Accnp());
    // std::cout << "gridvel fluid" << std::endl << *(FluidField()->GridVel());
    std::cout << "gridvel struct" << std::endl << *(StructureField()->Velnp());
  }

  const int zeilennr = -1;
  const int spaltenr = -1;
  for (int i = 0; i < dofs; ++i)
  {
    if (CombinedDBCMap()->MyGID(i))
    {
      iterinc->ReplaceGlobalValue(i, 0, 0.0);
    }
    abs_iterinc->Update(1.0, *iterinc, 1.0);

    if (i == spaltenr)
      std::cout << "\n******************" << spaltenr + 1 << ". Spalte!!***************"
                << std::endl;

    Evaluate(iterinc);

    rhs_copy->Update(1.0, *rhs_, 0.0);

    iterinc_->PutScalar(0.0);  // Useful? depends on solver and more
    LINALG::ApplyDirichlettoSystem(
        sparse_copy, iterinc_, rhs_copy, Teuchos::null, zeros_, *CombinedDBCMap());
    Teuchos::RCP<Epetra_CrsMatrix> test_crs = sparse_copy->EpetraMatrix();
    int sparsenumentries;
    int sparselength = test_crs->NumGlobalEntries(i);
    std::vector<double> sparsevalues(sparselength);
    std::vector<int> sparseindices(sparselength);
    // int sparseextractionstatus =
    test_crs->ExtractGlobalRowCopy(
        i, sparselength, sparsenumentries, &sparsevalues[0], &sparseindices[0]);


    if (i == spaltenr)
    {
      std::cout << "rhs_: " << (*rhs_copy)[zeilennr] << std::endl;
      std::cout << "rhs_old: " << (*rhs_old)[zeilennr] << std::endl;
    }
    // rhs_copy = ( rhs_disturb - rhs_old ) . (-1)/delta with rhs_copy==rhs_disturb
    rhs_copy->Update(-1.0, *rhs_old, 1.0);
    rhs_copy->Scale(-1.0 / delta);

    if (i == spaltenr)
    {
      std::cout << "( rhs_disturb - rhs_old )               "
                << (*rhs_copy)[zeilennr] * (-1.0) * delta << std::endl;
      std::cout << "( rhs_disturb - rhs_old ) . (-1)/delta: " << (*rhs_copy)[zeilennr] << std::endl;
      // LINALG::PrintMatrixInMatlabFormat("../o/mymatrix.dat",*test_crs);
      // dserror("exit here");
      // std::cout << sparse_copy(4,35) << std::endl;
    }
    int* index = &i;
    for (int j = 0; j < dofs; ++j)
    {
      double value = (*rhs_copy)[j];
      stiff_approx->InsertGlobalValues(j, 1, &value, index);

      if ((j == zeilennr) and (i == spaltenr))
      {
        std::cout << "\n******************" << zeilennr + 1 << ". Zeile!!***************"
                  << std::endl;
        // std::cout << "iterinc_" << std::endl << *iterinc_ << std::endl;
        // std::cout << "iterinc" << std::endl << *iterinc << std::endl;
        // std::cout << "meshdisp: " << std::endl << *(FluidField()->Dispnp());
        std::cout << "disp: " << std::endl << *(StructureField()->Dispnp());
        // std::cout << "fluid vel" << std::endl << *(FluidField()->Velnp());
        // std::cout << "fluid acc" << std::endl << *(FluidField()->Accnp());
        // std::cout << "gridvel fluid" << std::endl << *(FluidField()->GridVel());
        std::cout << "gridvel struct" << std::endl << *(StructureField()->Velnp());

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

  stiff_approx->FillComplete();

  Teuchos::RCP<LINALG::SparseMatrix> stiff_approx_sparse = Teuchos::null;
  stiff_approx_sparse = Teuchos::rcp(new LINALG::SparseMatrix(stiff_approx, LINALG::Copy));

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

          if (abs(error) > abs(error_max)) error_max = abs(error);
          if (abs(error_ij) > abs(abs_error_max)) abs_error_max = abs(error_ij);

          if ((abs(error) > 1e-4))
          {
            if ((abs(error_ij) > 1e-5))
            {
              // if( (sparse_ij>1e-1) or (stiff_approx_ij>1e-1) )
              {
                std::cout << "finite difference check failed entry (" << i << "," << j
                          << ")! stiff: " << sparse_ij << ", approx: " << stiff_approx_ij
                          << " ,abs. error: " << error_ij << " , rel. error: " << error
                          << std::endl;

                success = false;
              }
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
    dserror("PoroFDCheck failed in step: %d, iter: %d", Step(), itnum_);

  return;
}

/*----------------------------------------------------------------------*
 | constructor                                         kremheller 05/18 |
 *----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::
    PoroMultiPhaseMonolithicTwoWayArteryCoupling(
        const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : PoroMultiPhaseMonolithicTwoWay(comm, globaltimeparams)
{
  blockrowdofmap_artporo_ = Teuchos::rcp(new LINALG::MultiMapExtractor);

  return;
}

/*----------------------------------------------------------------------*
 | setup the map                                       kremheller 04/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::SetupMaps()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;

  vecSpaces.push_back(StructDofRowMap());

  vecSpaces.push_back(FluidDofRowMap());

  vecSpaces.push_back(ArteryDofRowMap());

  if (vecSpaces[0]->NumGlobalElements() == 0) dserror("No structure equation. Panic.");
  if (vecSpaces[1]->NumGlobalElements() == 0) dserror("No fluid equation. Panic.");
  if (vecSpaces[2]->NumGlobalElements() == 0) dserror("No fluid equation. Panic.");

  // full Poromultiphase-elasticity-map
  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(vecSpaces);

  // full Poromultiphase-elasticity-blockmap
  blockrowdofmap_->Setup(*fullmap_, vecSpaces);

  // full map of artery and poromulti DOFs
  fullmap_artporo_ = LINALG::MultiMapExtractor::MergeMaps({vecSpaces[1], vecSpaces[2]});

  // full artery-poromulti-blockmap
  blockrowdofmap_artporo_->Setup(*fullmap_artporo_, {vecSpaces[1], vecSpaces[2]});

  return;
}

/*----------------------------------------------------------------------*
 | extract field vectors for calling Evaluate() of the  kremheller 04/18|
 | single fields                                                        |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::ExtractFieldVectors(
    Teuchos::RCP<const Epetra_Vector> x, Teuchos::RCP<const Epetra_Vector>& sx,
    Teuchos::RCP<const Epetra_Vector>& fx)
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::ExtractFieldVectors");

  // process structure unknowns of the first field
  sx = Extractor()->ExtractVector(x, 0);

  // process artery and porofluid unknowns
  Teuchos::RCP<const Epetra_Vector> porofluid = Extractor()->ExtractVector(x, 1);
  Teuchos::RCP<const Epetra_Vector> artery = Extractor()->ExtractVector(x, 2);

  Teuchos::RCP<Epetra_Vector> dummy = Teuchos::rcp(new Epetra_Vector(*fullmap_artporo_));

  blockrowdofmap_artporo_->InsertVector(porofluid, 0, dummy);
  blockrowdofmap_artporo_->InsertVector(artery, 1, dummy);

  fx = dummy;

  return;
}

/*----------------------------------------------------------------------*
 | setup system matrix of poromultiphase-elasticity with artery         |
 | coupling                                           kremheller 05/18  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::SetupSystemMatrix(
    LINALG::BlockSparseMatrixBase& mat)
{
  PoroMultiPhaseMonolithicTwoWay::SetupSystemMatrix(mat);

  // pure artery part
  mat.Assign(2, 2, LINALG::View, ArteryPorofluidSysmat()->Matrix(1, 1));
  // artery-porofluid part
  mat.Assign(2, 1, LINALG::View, ArteryPorofluidSysmat()->Matrix(1, 0));
  // porofluid-artery part
  mat.Assign(1, 2, LINALG::View, ArteryPorofluidSysmat()->Matrix(0, 1));

  // Debug: matlab output of system matrix
  bool matlab = false;
  if (matlab)
  {
    std::string filename = "../o/mymatrix.dat";
    LINALG::PrintBlockMatrixInMatlabFormat(filename, mat);
    dserror("exit");
  }

  return;
}

/*----------------------------------------------------------------------*
 | setup rhs of poromultiphase-elasticity with artery coupling          |
 |                                                    kremheller 05/18  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::SetupRHS()
{
  // get structure part
  Teuchos::RCP<Epetra_Vector> str_rhs = SetupStructurePartofRHS();

  // insert and scale
  Extractor()->InsertVector(*str_rhs, 0, *rhs_);
  rhs_->Scale(-1.0);

  // insert artery part and porofluid part
  Extractor()->InsertVector(
      *(blockrowdofmap_artporo_->ExtractVector(FluidField()->ArteryPorofluidRHS(), 0)), 1, *rhs_);
  Extractor()->InsertVector(
      *(blockrowdofmap_artporo_->ExtractVector(FluidField()->ArteryPorofluidRHS(), 1)), 2, *rhs_);

  return;
}

/*-----------------------------------------------------------------------/
|  build the combined dbcmap                           kremheller 05/18  |
/-----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::BuildCombinedDBCMap()
{
  PoroMultiPhaseMonolithicTwoWay::BuildCombinedDBCMap();

  const Teuchos::RCP<const Epetra_Map> artcondmap =
      FluidField()->ArtNetTimInt()->GetDBCMapExtractor()->CondMap();

  // merge them
  combinedDBCMap_ = LINALG::MergeMap(combinedDBCMap_, artcondmap, false);

  return;
}
/*----------------------------------------------------------------------------*
 | build null space for artery block of global system matrix kremheller 05/18 |
 *--------------------------------------------------------------------------- */
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::BuildArteryBlockNullSpace(
    Teuchos::RCP<LINALG::Solver>& solver, const int& arteryblocknum)
{
  // equip smoother for fluid matrix block with empty parameter sublists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams3 =
      solver->Params().sublist("Inverse" + std::to_string(arteryblocknum));
  blocksmootherparams3.sublist("Aztec Parameters");
  blocksmootherparams3.sublist("MueLu Parameters");

  // build null space of complete discretization
  FluidField()->ArtNetTimInt()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams3);
  // fix the null space if some DOFs are condensed out
  LINALG::Solver::FixMLNullspace("Artery",
      *(FluidField()->ArtNetTimInt()->Discretization()->DofRowMap(0)),
      *(FluidField()->ArteryDofRowMap()), blocksmootherparams3);

  return;
}

/*----------------------------------------------------------------------*
 | Create linear (iterative) solver                    kremheller 05/18 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhaseMonolithicTwoWayArteryCoupling::CreateLinearSolver(
    const Teuchos::ParameterList& solverparams, const int solvertype)
{
  PoroMultiPhaseMonolithicTwoWay::CreateLinearSolver(solverparams, solvertype);

  // build also the artery null space
  BuildArteryBlockNullSpace(solver_, 3);
}
