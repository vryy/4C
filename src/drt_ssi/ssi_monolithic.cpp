/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "ssi_monolithic.H"

#include "ssi_coupling.H"
#include "ssi_monolithic_assemble_strategy.H"
#include "ssi_monolithic_convcheck_strategies.H"
#include "ssi_monolithic_evaluate_OffDiag.H"
#include "ssi_resulttest.H"
#include "ssi_str_model_evaluator_monolithic.H"
#include "ssi_utils.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/ad_str_structure_new.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_inpar/inpar_scatra.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_scatra/scatra_timint_meshtying_strategy_s2i.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_equilibrate.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

#include <Epetra_Time.h>


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
SSI::SSI_Mono::SSI_Mono(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams)
    : SSI_Base(comm, globaltimeparams),
      dtele_(0.0),
      dtsolve_(0.0),
      equilibration_method_(Teuchos::getIntegralValue<LINALG::EquilibrationMethod>(
          globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION")),
      increment_(Teuchos::null),
      map_structure_(Teuchos::null),
      maps_scatra_(Teuchos::null),
      maps_sub_problems_(Teuchos::null),
      maps_systemmatrix_(Teuchos::null),
      matrixtype_(Teuchos::getIntegralValue<LINALG::MatrixType>(
          globaltimeparams.sublist("MONOLITHIC"), "MATRIXTYPE")),
      meshtying_strategy_s2i_(Teuchos::null),
      residual_(Teuchos::null),
      scatrastructureOffDiagcoupling_(Teuchos::null),
      solver_(Teuchos::rcp(
          new LINALG::Solver(DRT::Problem::Instance()->SolverParams(
                                 globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
              comm, DRT::Problem::Instance()->ErrorFile()->Handle()))),
      ssi_matrices_(Teuchos::null),
      strategy_assemble_(Teuchos::null),
      strategy_convcheck_(Teuchos::null),
      strategy_equilibration_(Teuchos::null),
      timer_(Teuchos::rcp(new Epetra_Time(comm)))

{
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::AssembleMatAndRHS()
{
  // needed to communicate to NOX state
  StructureField()->SetState(StructureField()->WriteAccessDispnp());

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());

  // pass scalar transport degrees of freedom to structural discretization
  SetScatraSolution(ScaTraField()->Phinp());

  // evaluate temperature from function and set to structural discretization
  EvaluateAndSetTemperatureField();

  // build system matrix and residual for structure field
  StructureField()->Evaluate();

  // build system matrix and residual for scalar transport field
  ScaTraField()->PrepareLinearSolve();

  // evaluate off-diagonal scatra-structure block (domain contributions) of global system matrix
  scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraStructureDomain(
      ssi_matrices_->ScaTraStructureDomain());

  // evaluate off-diagonal scatra-structure block (interface contributions) of global system matrix
  if (SSIInterfaceMeshtying())
    scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraStructureInterface(
        ssi_matrices_->ScaTraStructureInterface());

  // evaluate off-diagonal structure-scatra block (we only have domain contributions so far) of
  // global system matrix
  scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockStructureScatraDomain(
      ssi_matrices_->StructureScaTraDomain());

  // assemble scatra block into system matrix
  strategy_assemble_->AssembleScatraDomain(
      ssi_matrices_->SystemMatrix(), ScaTraField()->SystemMatrixOperator());

  // assemble scatra-strucutre block (domain contributions) into system matrix
  strategy_assemble_->AssembleScatraStructureDomain(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraStructureDomain());

  // assemble scatra-strucutre block (interface contributions) into system matrix
  if (SSIInterfaceMeshtying())
    strategy_assemble_->AssembleScatraStructureInterface(
        ssi_matrices_->SystemMatrix(), ssi_matrices_->ScaTraStructureInterface());

  // assemble structure-scatra block (domain contributions) into system matrix
  strategy_assemble_->AssembleStructureScatraDomain(
      ssi_matrices_->SystemMatrix(), ssi_matrices_->StructureScaTraDomain());

  // assemble strucutre block into system matrix
  strategy_assemble_->AssembleStructureDomain(
      ssi_matrices_->SystemMatrix(), StructureField()->SystemMatrix());

  // apply meshtying
  strategy_assemble_->ApplyMeshtyingSystemMatrix(ssi_matrices_->SystemMatrix());

  // finalize global system matrix
  ssi_matrices_->SystemMatrix()->Complete();

  // apply scatra Dirichlet
  ssi_matrices_->SystemMatrix()->ApplyDirichlet(*ScaTraField()->DirichMaps()->CondMap(), true);

  // apply structural Dirichlet conditions
  strategy_assemble_->ApplyStructuralDBCSystemMatrix(ssi_matrices_->SystemMatrix());

  // assemble monolithic RHS
  strategy_assemble_->AssembleRHS(residual_, ScaTraField()->Residual(), StructureField()->RHS());
}

/*-------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------*/
void SSI::SSI_Mono::BuildNullSpaces() const
{
  switch (ScaTraField()->MatrixType())
  {
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      // equip smoother for scatra matrix blocks with null space
      ScaTraField()->BuildBlockNullSpaces(solver_, 0);

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sub lists to trigger null space
      // computation
      Teuchos::ParameterList& blocksmootherparams = solver_->Params().sublist("Inverse1");
      blocksmootherparams.sublist("Aztec Parameters");
      blocksmootherparams.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // store number of matrix block associated with structural field as string
  std::stringstream iblockstr;
  iblockstr << MapsSystemMatrix()->NumMaps();

  // equip smoother for structural matrix block with empty parameter sub lists to trigger null space
  // computation
  Teuchos::ParameterList& blocksmootherparams =
      solver_->Params().sublist("Inverse" + iblockstr.str());
  blocksmootherparams.sublist("Aztec Parameters");
  blocksmootherparams.sublist("MueLu Parameters");

  // equip smoother for structural matrix block with null space associated with all degrees of
  // freedom on structural discretization
  StructureField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);
}  // SSI::SSI_Mono::BuildNullSpaces


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& SSI::SSI_Mono::DofRowMap() const
{
  return MapsSubProblems()->FullMap();
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::FDCheck()
{
  // initial screen output
  if (Comm().MyPID() == 0)
    std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR SSI SYSTEM MATRIX" << std::endl;

  // create global state vector
  Teuchos::RCP<Epetra_Vector> statenp(LINALG::CreateVector(*DofRowMap(), true));
  MapsSubProblems()->InsertVector(ScaTraField()->Phinp(), 0, statenp);
  MapsSubProblems()->InsertVector(StructureField()->Dispnp(), 1, statenp);

  // make a copy of global state vector to undo perturbations later
  Teuchos::RCP<Epetra_Vector> statenp_original = Teuchos::rcp(new Epetra_Vector(*statenp));

  // make a copy of system matrix as Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> sysmat_original = Teuchos::null;
  if (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(ssi_matrices_->SystemMatrix()) !=
      Teuchos::null)
  {
    sysmat_original =
        (new LINALG::SparseMatrix(
             *(Teuchos::rcp_static_cast<LINALG::SparseMatrix>(ssi_matrices_->SystemMatrix()))))
            ->EpetraMatrix();
  }
  else if (Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(
               ssi_matrices_->SystemMatrix()) != Teuchos::null)
  {
    sysmat_original =
        (new LINALG::SparseMatrix(*(
             Teuchos::rcp_static_cast<LINALG::BlockSparseMatrixBase>(ssi_matrices_->SystemMatrix())
                 ->Merge())))
            ->EpetraMatrix();
  }
  else
    dserror("Type of system matrix unknown!");
  sysmat_original->FillComplete();

  // make a copy of system right-hand side vector
  Teuchos::RCP<Epetra_Vector> rhs_original = Teuchos::rcp(new Epetra_Vector(*residual_));

  // copy and zero out system increment vector if necessary
  Teuchos::RCP<Epetra_Vector> increment_original(Teuchos::null);
  if (IterationCount() != 1)
  {
    increment_original = Teuchos::rcp(new Epetra_Vector(*increment_));
    increment_->PutScalar(0.0);
  }

  // initialize counter for system matrix entries with failing finite difference check
  int counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  for (int colgid = 0; colgid <= sysmat_original->ColMap().MaxAllGID(); ++colgid)
  {
    // continue loop if current column index is not a valid global column index
    int collid(sysmat_original->ColMap().LID(colgid));
    int maxcollid(-1);
    Comm().MaxAll(&collid, &maxcollid, 1);
    if (maxcollid < 0) continue;

    // continue loop if current column index is associated with slave side of structural
    // meshtying interface
    if (SSIInterfaceMeshtying())
    {
      collid = InterfaceCouplingAdapterStructure()->SlaveDofMap()->LID(colgid);
      Comm().MaxAll(&collid, &maxcollid, 1);
      if (maxcollid >= 0) continue;
    }

    // fill global state vector with original state variables
    statenp->Update(1., *statenp_original, 0.);

    // impose perturbation
    if (statenp->Map().MyGID(colgid))
      if (statenp->SumIntoGlobalValue(colgid, 0, ScaTraField()->FDCheckEps()))
        dserror("Perturbation could not be imposed on state vector for finite difference check!");
    if (SSIInterfaceMeshtying() and
        InterfaceCouplingAdapterStructure()->PermMasterDofMap()->MyGID(colgid))
    {
      if (statenp->SumIntoGlobalValue(
              InterfaceCouplingAdapterStructure()->SlaveDofMap()->GID(
                  InterfaceCouplingAdapterStructure()->PermMasterDofMap()->LID(colgid)),
              0, ScaTraField()->FDCheckEps()))
        dserror("Perturbation could not be imposed on state vector for finite difference check!");
    }
    ScaTraField()->Phinp()->Update(1., *MapsSubProblems()->ExtractVector(statenp, 0), 0.);
    StructureField()->SetState(MapsSubProblems()->ExtractVector(statenp, 1));

    // calculate element right-hand side vector for perturbed state
    if (IterationCount() != 1) UpdateIterStructure();
    AssembleMatAndRHS();

    // Now we compare the difference between the current entries in the system matrix
    // and their finite difference approximations according to
    // entries ?= (residual_perturbed - residual_original) / epsilon

    // Note that the residual_ vector actually denotes the right-hand side of the linear
    // system of equations, i.e., the negative system residual.
    // To account for errors due to numerical cancellation, we additionally consider
    // entries + residual_original / epsilon ?= residual_perturbed / epsilon

    // Note that we still need to evaluate the first comparison as well. For small entries in
    // the system matrix, the second comparison might yield good agreement in spite of the
    // entries being wrong!
    for (int rowlid = 0; rowlid < DofRowMap()->NumMyElements(); ++rowlid)
    {
      // get global index of current matrix row
      const int rowgid = sysmat_original->RowMap().GID(rowlid);
      if (rowgid < 0) dserror("Invalid global ID of matrix row!");

      // skip matrix rows associated with Dirichlet boundary conditions and slave side of
      // structural meshtying interface
      if (ScaTraField()->DirichMaps()->CondMap()->MyGID(rowgid) or
          StructureField()->GetDBCMapExtractor()->CondMap()->MyGID(rowgid) or
          (SSIInterfaceMeshtying() and
              InterfaceCouplingAdapterStructure()->SlaveDofMap()->MyGID(rowgid)))
        continue;

      // get relevant entry in current row of original system matrix
      double entry(0.);
      int length = sysmat_original->NumMyEntries(rowlid);
      int numentries;
      std::vector<double> values(length);
      std::vector<int> indices(length);
      sysmat_original->ExtractMyRowCopy(rowlid, length, numentries, &values[0], &indices[0]);
      for (int ientry = 0; ientry < length; ++ientry)
      {
        if (sysmat_original->ColMap().GID(indices[ientry]) == colgid)
        {
          entry = values[ientry];
          break;
        }
      }

      // finite difference suggestion (first divide by epsilon and then add for better
      // conditioning)
      const double fdval = -(*residual_)[rowlid] / ScaTraField()->FDCheckEps() +
                           (*rhs_original)[rowlid] / ScaTraField()->FDCheckEps();

      // confirm accuracy of first comparison
      if (abs(fdval) > 1.e-20 and abs(fdval) < 1.e-15)
      {
        // output warning
        std::cout
            << "WARNING: Finite difference check involves values very close to numerical zero!"
            << std::endl;

        // skip comparison if current entry is very small
        if (abs(entry) < 1.e-15) continue;
      }

      // absolute and relative errors in first comparison
      const double abserr1 = entry - fdval;
      if (abs(abserr1) > maxabserr) maxabserr = abs(abserr1);
      double relerr1(0.);
      if (abs(entry) > 1.e-17)
        relerr1 = abserr1 / abs(entry);
      else if (abs(fdval) > 1.e-17)
        relerr1 = abserr1 / abs(fdval);
      if (abs(relerr1) > maxrelerr) maxrelerr = abs(relerr1);

      // evaluate first comparison
      if (abs(relerr1) > ScaTraField()->FDCheckTol())
      {
        std::cout << "sysmat[" << rowgid << "," << colgid << "]:  " << entry << "   ";
        std::cout << "finite difference suggestion:  " << fdval << "   ";
        std::cout << "absolute error:  " << abserr1 << "   ";
        std::cout << "relative error:  " << relerr1 << std::endl;

        counter++;
      }

      // first comparison OK
      else
      {
        // left-hand side in second comparison
        const double left = entry - (*rhs_original)[rowlid] / ScaTraField()->FDCheckEps();

        // right-hand side in second comparison
        const double right = -(*residual_)[rowlid] / ScaTraField()->FDCheckEps();

        // confirm accuracy of second comparison
        if (abs(right) > 1.e-20 and abs(right) < 1.e-15)
        {
          // output warning
          std::cout << "WARNING: Finite difference check involves values very close to "
                       "numerical zero!"
                    << std::endl;

          // skip comparison if current left-hand side is very small
          if (abs(left) < 1.e-15) continue;
        }

        // absolute and relative errors in second comparison
        const double abserr2 = left - right;
        if (abs(abserr2) > maxabserr) maxabserr = abs(abserr2);
        double relerr2(0.);
        if (abs(left) > 1.e-17)
          relerr2 = abserr2 / abs(left);
        else if (abs(right) > 1.e-17)
          relerr2 = abserr2 / abs(right);
        if (abs(relerr2) > maxrelerr) maxrelerr = abs(relerr2);

        // evaluate second comparison
        if (abs(relerr2) > ScaTraField()->FDCheckTol())
        {
          std::cout << "sysmat[" << rowgid << "," << colgid << "]-rhs[" << rowgid
                    << "]/eps:  " << left << "   ";
          std::cout << "-rhs_perturbed[" << rowgid << "]/eps:  " << right << "   ";
          std::cout << "absolute error:  " << abserr2 << "   ";
          std::cout << "relative error:  " << relerr2 << std::endl;

          counter++;
        }
      }
    }
  }

  // communicate tracking variables
  int counterglobal(0);
  Comm().SumAll(&counter, &counterglobal, 1);
  double maxabserrglobal(0.);
  Comm().MaxAll(&maxabserr, &maxabserrglobal, 1);
  double maxrelerrglobal(0.);
  Comm().MaxAll(&maxrelerr, &maxrelerrglobal, 1);

  // final screen output
  if (Comm().MyPID() == 0)
  {
    if (counterglobal)
    {
      printf(
          "--> FAILED AS LISTED ABOVE WITH %d CRITICAL MATRIX ENTRIES IN TOTAL\n\n", counterglobal);
      printf(
          "--> FAILED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR "
          "%+12.5e\n\n",
          maxabserrglobal, maxrelerrglobal);
      dserror("Finite difference check failed for SSI system matrix!");
    }
    else
    {
      printf(
          "--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR "
          "%+12.5e\n\n",
          maxabserrglobal, maxrelerrglobal);
    }
  }

  // undo perturbations of state variables
  ScaTraField()->Phinp()->Update(1., *MapsSubProblems()->ExtractVector(statenp_original, 0), 0.0);
  StructureField()->SetState(MapsSubProblems()->ExtractVector(statenp_original, 1));

  // recompute system matrix and right-hand side vector based on original state variables
  if (IterationCount() != 1) UpdateIterStructure();
  AssembleMatAndRHS();

  // restore system increment vector if necessary
  if (increment_original != Teuchos::null) increment_ = increment_original;
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
int SSI::SSI_Mono::Init(const Epetra_Comm& comm, const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams, const Teuchos::ParameterList& structparams,
    const std::string struct_disname, const std::string scatra_disname, bool isAle)
{
  // check input parameters for scalar transport field
  if (DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatraparams, "VELOCITYFIELD") !=
      INPAR::SCATRA::velocity_Navier_Stokes)
    dserror("Invalid type of velocity field for scalar-structure interaction!");

  // initialize strategy for Newton-Raphson convergence check
  switch (
      DRT::INPUT::IntegralValue<INPAR::SSI::ScaTraTimIntType>(globaltimeparams, "SCATRATIMINTTYPE"))
  {
    case INPAR::SSI::scatratiminttype_elch:
    {
      strategy_convcheck_ =
          Teuchos::rcp(new SSI::SSI_Mono::ConvCheckStrategyElch(globaltimeparams));
      break;
    }

    case INPAR::SSI::scatratiminttype_standard:
    {
      strategy_convcheck_ = Teuchos::rcp(new SSI::SSI_Mono::ConvCheckStrategyStd(globaltimeparams));
      break;
    }

    default:
    {
      dserror("Type of scalar transport time integrator currently not supported!");
      break;
    }
  }

  // call base class routine
  return SSI_Base::Init(
      comm, globaltimeparams, scatraparams, structparams, struct_disname, scatra_disname, isAle);
}


/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Output()
{
  // output scalar transport field
  ScaTraField()->Output();

  // output structure field
  StructureField()->Output();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());

  // prepare time step for scalar transport field
  ScaTraField()->PrepareTimeStep();

  // if adaptive time stepping and different time step size: calculate time step in scatra
  // (PrepareTimeStep() of Scatra) and pass to structure
  if (ScaTraField()->TimeStepAdapted()) SetDtFromScaTraToStructure();

  // pass scalar transport degrees of freedom to structural discretization
  // has to be called AFTER ScaTraField()->PrepareTimeStep() to ensure
  // consistent scalar transport state vector with valid Dirichlet conditions
  SetScatraSolution(ScaTraField()->Phinp());

  // evaluate temperature from function and set to structural discretization
  EvaluateAndSetTemperatureField();

  // prepare time step for structural field
  StructureField()->PrepareTimeStep();

  // print time step information to screen
  ScaTraField()->PrintTimeStepInfo();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Setup()
{
  // call base class routine
  SSI_Base::Setup();

  // safety checks
  if (ScaTraField()->NumScal() != 1)
  {
    dserror(
        "Since the ssi_monolithic framework is only implemented for usage in combination with "
        "volume change laws 'MAT_InelasticDefgradLinScalarIso' or "
        "'MAT_InelasticDefgradLinScalarAniso' so far and these laws are implemented for only "
        "one transported scalar at the moment it is not reasonable to use them with more than one "
        "transported scalar. So you need to cope with it or change implementation! ;-)");
  }
  if (ScaTraField()->EquilibrationMethod() != LINALG::EquilibrationMethod::none)
  {
    dserror(
        "You are within the monolithic solid scatra interaction framework but activated a pure "
        "scatra equilibration method. Delete this from 'SCALAR TRANSPORT DYNAMIC' section and set "
        "it in 'SSI CONTROL/MONOLITHIC' instead.");
  }

  if (!ScaTraField()->IsIncremental())
    dserror("Must have incremental solution approach for monolithic scalar-structure interaction!");

  // set up scatra-scatra interface meshtying if set in input file
  if (SSIInterfaceMeshtying())
  {
    // extract meshtying strategy for scatra-scatra interface coupling on scatra discretization
    meshtying_strategy_s2i_ =
        Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(ScaTraField()->Strategy());

    // safety checks
    if (meshtying_strategy_s2i_ == Teuchos::null)
      dserror("Invalid scatra-scatra interface coupling strategy!");
    if (meshtying_strategy_s2i_->CouplingType() != INPAR::S2I::coupling_matching_nodes)
    {
      dserror(
          "Monolithic scalar-structure interaction only implemented for scatra-scatra "
          "interface coupling with matching interface nodes!");
    }
  }
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::SetupSystem()
{
  // merge slave and master side block maps for interface matrix for scatra
  Teuchos::RCP<Epetra_Map> interface_map_scatra(Teuchos::null);

  if (SSIInterfaceMeshtying())
  {
    // check whether slave-side degrees of freedom are Dirichlet-free
    std::vector<Teuchos::RCP<const Epetra_Map>> maps(2, Teuchos::null);
    maps[0] = InterfaceCouplingAdapterStructure()->SlaveDofMap();
    maps[1] = StructureField()->GetDBCMapExtractor()->CondMap();
    if (LINALG::MultiMapExtractor::IntersectMaps(maps)->NumGlobalElements() > 0)
      dserror("Must not apply Dirichlet conditions to slave-side structural displacements!");

    interface_map_scatra = LINALG::MultiMapExtractor::MergeMaps(
        {meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(),
            meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap()});
  }

  // initialize global map extractor
  maps_sub_problems_ = Teuchos::rcp(new LINALG::MapExtractor(
      *LINALG::MergeMap(*ScaTraField()->DofRowMap(), *StructureField()->DofRowMap(), false),
      StructureField()->DofRowMap(), ScaTraField()->DofRowMap()));

  // check global map extractor
  maps_sub_problems_->CheckForValidMapExtractor();

  // initialize global increment vector for Newton-Raphson iteration
  increment_ = LINALG::CreateVector(*DofRowMap(), true);

  // initialize global residual vector
  residual_ = LINALG::CreateVector(*DofRowMap(), true);

  // initialize map extractors associated with blocks of global system matrix
  switch (ScaTraField()->MatrixType())
  {
    // one single main-diagonal matrix block associated with scalar transport field
    case LINALG::MatrixType::sparse:
    {
      maps_systemmatrix_ = MapsSubProblems();
      break;
    }

    // several main-diagonal matrix blocks associated with scalar transport field
    case LINALG::MatrixType::block_condition:
    case LINALG::MatrixType::block_condition_dof:
    {
      // store an RCP to the block maps of the scatra field
      maps_scatra_ = Teuchos::rcpFromRef(ScaTraField()->BlockMaps());

      // safety check
      maps_scatra_->CheckForValidMapExtractor();

      // extract maps underlying main-diagonal matrix blocks associated with scalar transport  field
      const int maps_systemmatrix_scatra = maps_scatra_->NumMaps();
      std::vector<Teuchos::RCP<const Epetra_Map>> maps_systemmatrix(maps_systemmatrix_scatra + 1);
      for (int imap = 0; imap < maps_systemmatrix_scatra; ++imap)
        maps_systemmatrix[imap] = maps_scatra_->Map(imap);

      // extract map underlying single main-diagonal matrix block associated with structural field
      maps_systemmatrix[maps_systemmatrix_scatra] = StructureField()->DofRowMap();

      // initialize map extractor associated with blocks of global system matrix
      maps_systemmatrix_ =
          Teuchos::rcp(new LINALG::MultiMapExtractor(*DofRowMap(), maps_systemmatrix));

      // initialize map extractor associated with all degrees of freedom inside structural field
      map_structure_ = Teuchos::rcp(
          new LINALG::MultiMapExtractor(*StructureField()->Discretization()->DofRowMap(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(1, StructureField()->DofRowMap())));

      // safety check
      map_structure_->CheckForValidMapExtractor();

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // safety check
  maps_systemmatrix_->CheckForValidMapExtractor();

  // perform initializations associated with global system matrix
  switch (matrixtype_)
  {
    case LINALG::MatrixType::block_field:
    {
      // safety check
      if (!solver_->Params().isSublist("AMGnxn Parameters"))
        dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

      // feed AMGnxn block preconditioner with null space information for each block of global block
      // system matrix
      BuildNullSpaces();

      break;
    }

    case LINALG::MatrixType::sparse:
    {
      // safety check
      if (ScaTraField()->SystemMatrix() == Teuchos::null)
        dserror("Incompatible matrix type associated with scalar transport field!");
      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // initialize subblocks and system matrix
  ssi_matrices_ = Teuchos::rcp(new SSI::UTILS::SSIMatrices(*this, interface_map_scatra));

  // initialize strategy for assembly
  strategy_assemble_ = SSI::BuildAssembleStrategy(
      Teuchos::rcp(this, false), matrixtype_, ScaTraField()->MatrixType());

  // initialize object, that performs evaluations of OD coupling
  scatrastructureOffDiagcoupling_ =
      Teuchos::rcp(new SSI::ScatraStructureOffDiagCoupling(MapsStructure(),
          MapsSubProblems()->Map(0), MapsSubProblems()->Map(1), InterfaceCouplingAdapterStructure(),
          interface_map_scatra, meshtying_strategy_s2i_, ScaTraBaseAlgorithm(), StructureField()));

  // instantiate appropriate equilibration class
  strategy_equilibration_ =
      LINALG::BuildEquilibration(matrixtype_, equilibration_method_, MapsSubProblems()->FullMap());
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSI_Mono::SetupModelEvaluator() const
{
  // construct and register structural model evaluator if necessary
  if (DRT::INPUT::IntegralValue<INPAR::STR::StressType>(
          DRT::Problem::Instance()->IOParams(), "STRUCT_STRESS") != INPAR::STR::stress_none and
      SSIInterfaceMeshtying())
    StructureBaseAlgorithm()->RegisterModelEvaluator("Monolithic Coupling Model",
        Teuchos::rcp(new STR::MODELEVALUATOR::MonolithicSSI(Teuchos::rcp(this, false))));
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
void SSI::SSI_Mono::SolveLinearSystem()
{
  strategy_equilibration_->EquilibrateSystem(
      ssi_matrices_->SystemMatrix(), residual_, *MapsSystemMatrix());

  // solve global system of equations
  // Dirichlet boundary conditions have already been applied to global system of equations
  solver_->Solve(ssi_matrices_->SystemMatrix()->EpetraOperator(), increment_, residual_, true,
      IterationCount() == 1);

  strategy_equilibration_->UnequilibrateIncrement(increment_);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::NewtonLoop()
{
  // reset counter for Newton-Raphson iteration
  ResetIterationCount();

  // start Newton-Raphson iteration
  while (true)
  {
    // update iteration counter
    IncrementIterationCount();

    // reset timer
    timer_->ResetStartTime();

    // store time before evaluating elements and assembling global system of equations
    double time = timer_->WallTime();

    // assemble global system of equations
    AssembleMatAndRHS();

    // determine time needed for evaluating elements and assembling global system of
    // equations, and take maximum over all processors via communication
    double mydtele = timer_->WallTime() - time;
    Comm().MaxAll(&mydtele, &dtele_, 1);

    // safety check
    if (!ssi_matrices_->SystemMatrix()->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // perform finite difference check on time integrator level
    if ((ScaTraField()->FDCheckType() == INPAR::SCATRA::fdcheck_global) and (Step() > 1)) FDCheck();

    // check termination criterion for Newton-Raphson iteration
    if (strategy_convcheck_->ExitNewtonRaphson(*this)) break;

    // initialize global increment vector
    increment_->PutScalar(0.0);

    // store time before solving global system of equations
    time = timer_->WallTime();

    SolveLinearSystem();

    // determine time needed for solving global system of equations,
    // and take maximum over all processors via communication
    double mydtsolve = timer_->WallTime() - time;
    Comm().MaxAll(&mydtsolve, &dtsolve_, 1);

    // output performance statistics associated with linear solver into text file if
    // applicable
    if (DRT::INPUT::IntegralValue<bool>(
            *ScaTraField()->ScatraParameterList(), "OUTPUTLINSOLVERSTATS"))
      ScaTraField()->OutputLinSolverStats(
          *solver_, dtsolve_, Step(), IterationCount(), residual_->Map().NumGlobalElements());

    // update states for next Newton iteration
    UpdateIterScaTra();
    UpdateIterStructure();

  }  // Newton-Raphson iteration
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Timeloop()
{
  // output initial scalar transport solution to screen and files
  if (Step() == 0)
  {
    SetStructSolution(StructureField()->Dispnp(), StructureField()->Velnp());
    ScaTraField()->Output();
  }

  // time loop
  while (NotFinished() and ScaTraField()->NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // store time before calling nonlinear solver
    const double time = timer_->WallTime();

    // evaluate time step
    NewtonLoop();

    // determine time spent by nonlinear solver and take maximum over all processors via
    // communication
    double mydtnonlinsolve(timer_->WallTime() - time), dtnonlinsolve(0.);
    Comm().MaxAll(&mydtnonlinsolve, &dtnonlinsolve, 1);

    // output performance statistics associated with nonlinear solver into *.csv file if
    // applicable
    if (DRT::INPUT::IntegralValue<int>(
            *ScaTraField()->ScatraParameterList(), "OUTPUTNONLINSOLVERSTATS"))
      ScaTraField()->OutputNonlinSolverStats(IterationCount(), dtnonlinsolve, Step(), Comm());

    // prepare structure output
    StructureField()->PrepareOutput();

    // update scalar transport and structure fields
    Update();

    // output solution to screen and files
    Output();
  }
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSI_Mono::Update()
{
  // update scalar transport field
  ScaTraField()->Update();

  // update structure field
  StructureField()->Update();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSI_Mono::UpdateIterScaTra()
{
  // update scalar transport field
  ScaTraField()->UpdateIter(MapsSubProblems()->ExtractVector(increment_, 0));
  ScaTraField()->ComputeIntermediateValues();
}

/*--------------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------------*/
void SSI::SSI_Mono::UpdateIterStructure()
{
  // set up structural increment vector
  const Teuchos::RCP<Epetra_Vector> increment_structure =
      MapsSubProblems()->ExtractVector(increment_, 1);

  // consider structural meshtying. Copy master increments and displacements to slave side.
  if (SSIInterfaceMeshtying())
  {
    MapsStructure()->InsertVector(
        InterfaceCouplingAdapterStructure()->MasterToSlave(
            MapsStructure()->ExtractVector(StructureField()->Dispnp(), 2)),
        1, StructureField()->WriteAccessDispnp());
    StructureField()->SetState(StructureField()->WriteAccessDispnp());
    MapsStructure()->InsertVector(InterfaceCouplingAdapterStructure()->MasterToSlave(
                                      MapsStructure()->ExtractVector(increment_structure, 2)),
        1, increment_structure);
  }

  // update displacement of structure field
  StructureField()->UpdateStateIncrementally(increment_structure);
}
