/*--------------------------------------------------------------------------*/
/*! \file
\brief monolithic scalar-structure interaction

\level 2

\maintainer Christoph Schmidt

*/
/*--------------------------------------------------------------------------*/
#include "ssi_monolithic.H"
#include "ssi_coupling.H"
#include "ssi_monolithic_convcheck_strategies.H"
#include "ssi_resulttest.H"
#include "ssi_str_model_evaluator_monolithic.H"
#include "ssi_monolithic_assemble_strategy.H"
#include <Epetra_Time.h>

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

#include "../drt_scatra_ele/scatra_ele_action.H"

#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_matrixtransform.H"
#include "ssi_monolithic_evaluate_OffDiag.H"

/*--------------------------------------------------------------------------*
 | constructor                                                   fang 08/17 |
 *--------------------------------------------------------------------------*/
SSI::SSI_Mono::SSI_Mono(const Epetra_Comm& comm,    //!< communicator
    const Teuchos::ParameterList& globaltimeparams  //!< parameter list for time integration
    )
    : SSI_Base(comm, globaltimeparams),
      dtele_(0.0),
      dtsolve_(0.0),
      maps_scatra_(Teuchos::null),
      map_structure_(Teuchos::null),
      maps_(Teuchos::null),
      maps_systemmatrix_(Teuchos::null),
      equilibration_method_(Teuchos::getIntegralValue<INPAR::SSI::EquilibrationMethod>(
          globaltimeparams.sublist("MONOLITHIC"), "EQUILIBRATION")),
      matrixtype_(DRT::INPUT::IntegralValue<INPAR::SSI::MatrixType>(
          globaltimeparams.sublist("MONOLITHIC"), "MATRIXTYPE")),
      residual_(Teuchos::null),
      scatrastructuredomain_(Teuchos::null),
      scatrastructureinterface_(Teuchos::null),
      solver_(Teuchos::rcp(
          new LINALG::Solver(DRT::Problem::Instance()->SolverParams(
                                 globaltimeparams.sublist("MONOLITHIC").get<int>("LINEAR_SOLVER")),
              comm, DRT::Problem::Instance()->ErrorFile()->Handle()))),
      strategy_convcheck_(Teuchos::null),
      strategy_assemble_(Teuchos::null),
      meshtying_strategy_s2i_(Teuchos::null),
      structurescatradomain_(Teuchos::null),
      systemmatrix_(Teuchos::null),
      timer_(Teuchos::rcp(new Epetra_Time(comm)))
{
}

/*--------------------------------------------------------------------------*
 | assemble global system of equations                           fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::AssembleMatAndRHS()
{
  // pass scalar transport degrees of freedom to structural discretization
  SetScatraSolution(scatra_->ScaTraField()->Phinp());

  // build system matrix and residual for structure field
  if (iter_ == 1)
    structure_->Evaluate();
  else
  {
    // set up structural increment vector
    const Teuchos::RCP<Epetra_Vector> increment_structure = maps_->ExtractVector(increment_, 1);

    // consider structural meshtying if applicable
    if (SSIInterfaceMeshtying())
    {
      maps_structure_->InsertVector(
          icoup_structure_->MasterToSlave(maps_structure_->ExtractVector(structure_->Dispnp(), 2)),
          1, structure_->WriteAccessDispnp());
      structure_->SetState(structure_->WriteAccessDispnp());
      maps_structure_->InsertVector(
          icoup_structure_->MasterToSlave(maps_structure_->ExtractVector(increment_structure, 2)),
          1, increment_structure);
    }

    // evaluate structural field
    structure_->Evaluate(increment_structure);
  }

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(structure_->Dispnp(), structure_->Velnp());

  // build system matrix and residual for scalar transport field
  scatra_->ScaTraField()->PrepareLinearSolve();

  // evaluate off-diagonal scatra-structure block (domain contributions) of global system matrix
  scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraStructureDomain(
      scatrastructuredomain_);

  // evaluate off-diagonal scatra-structure block (interface contributions) of global system matrix
  if (SSIInterfaceMeshtying())
    scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockScatraStructureInterface(
        scatrastructureinterface_);

  // evaluate off-diagonal structure-scatra block (we only have domain contributions so far) of
  // global system matrix
  scatrastructureOffDiagcoupling_->EvaluateOffDiagBlockStructureScatraDomain(
      structurescatradomain_);

  // assemble scatra block into system matrix
  strategy_assemble_->AssembleScatraDomain(
      systemmatrix_, scatra_->ScaTraField()->SystemMatrixOperator());

  // assemble scatra-strucutre block (domain contributions) into system matrix
  strategy_assemble_->AssembleScatraStructureDomain(systemmatrix_, scatrastructuredomain_);

  // assemble scatra-strucutre block (interface contributions) into system matrix
  if (SSIInterfaceMeshtying())
    strategy_assemble_->AssembleScatraStructureInterface(systemmatrix_, scatrastructureinterface_);

  // assemble structure-scatra block (domain contributions) into system matrix
  strategy_assemble_->AssembleStructureScatraDomain(systemmatrix_, structurescatradomain_);

  // assemble strucutre block into system matrix
  strategy_assemble_->AssembleStructureDomain(systemmatrix_, structure_->SystemMatrix());

  // apply meshtying
  strategy_assemble_->ApplyMeshtyingSystemMatrix(systemmatrix_);

  // finalize global system matrix
  systemmatrix_->Complete();

  // apply scatra Dirichlet
  systemmatrix_->ApplyDirichlet(*scatra_->ScaTraField()->DirichMaps()->CondMap(), true);

  // apply structural Dirichlet conditions
  strategy_assemble_->ApplyStructuralDBCSystemMatrix(systemmatrix_);

  // assemble monolithic RHS
  strategy_assemble_->AssembleRHS(residual_, scatra_->ScaTraField()->Residual(), structure_->RHS());
}

/*-------------------------------------------------------------------------------*
 | build null spaces associated with blocks of global system matrix   fang 01/18 |
 *-------------------------------------------------------------------------------*/
void SSI::SSI_Mono::BuildNullSpaces() const
{
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      // loop over block(s) of global system matrix associated with scalar transport field
      for (int iblock = 0; iblock < scatra_->ScaTraField()->BlockMaps().NumMaps(); ++iblock)
      {
        // store number of current block as string, starting from 1
        std::stringstream iblockstr;
        iblockstr << iblock + 1;

        // equip smoother for current matrix block with empty parameter sublists to trigger null
        // space computation
        Teuchos::ParameterList& blocksmootherparams =
            solver_->Params().sublist("Inverse" + iblockstr.str());
        blocksmootherparams.sublist("Aztec Parameters");
        blocksmootherparams.sublist("MueLu Parameters");

        // equip smoother for current matrix block with null space associated with all degrees
        // of freedom on scalar transport discretization
        scatra_->ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

        // reduce full null space to match degrees of freedom associated with current matrix
        // block
        LINALG::Solver::FixMLNullspace("Block " + iblockstr.str(),
            *scatra_->ScaTraField()->Discretization()->DofRowMap(),
            *scatra_->ScaTraField()->BlockMaps().Map(iblock), blocksmootherparams);
      }

      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      // equip smoother for scatra matrix block with empty parameter sublists to trigger null
      // space computation
      Teuchos::ParameterList& blocksmootherparams = solver_->Params().sublist("Inverse1");
      blocksmootherparams.sublist("Aztec Parameters");
      blocksmootherparams.sublist("MueLu Parameters");

      // equip smoother for scatra matrix block with null space associated with all degrees of
      // freedom on scatra discretization
      scatra_->ScaTraField()->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);

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
  iblockstr << maps_systemmatrix_->NumMaps();

  // equip smoother for structural matrix block with empty parameter sublists to trigger null
  // space computation
  Teuchos::ParameterList& blocksmootherparams =
      solver_->Params().sublist("Inverse" + iblockstr.str());
  blocksmootherparams.sublist("Aztec Parameters");
  blocksmootherparams.sublist("MueLu Parameters");

  // equip smoother for structural matrix block with null space associated with all degrees of
  // freedom on structural discretization
  structure_->Discretization()->ComputeNullSpaceIfNecessary(blocksmootherparams);
}  // SSI::SSI_Mono::BuildNullSpaces

/*----------------------------------------------------------------------------*
 | compute inverse sums of absolute values of matrix row entries   fang 01/18 |
 *----------------------------------------------------------------------------*/
void SSI::SSI_Mono::ComputeInvRowSums(
    const LINALG::SparseMatrix& matrix, const Teuchos::RCP<Epetra_Vector>& invrowsums) const
{
  // compute inverse row sums of matrix
  if (matrix.EpetraMatrix()->InvRowSums(*invrowsums))
    dserror("Inverse row sums of matrix could not be successfully computed!");
}  // SSI::SSI_Mono::ComputeInvRowSums

/*--------------------------------------------------------------------------*
 | return global map of degrees of freedom                       fang 08/17 |
 *--------------------------------------------------------------------------*/
const Teuchos::RCP<const Epetra_Map>& SSI::SSI_Mono::DofRowMap() const { return maps_->FullMap(); }


/*----------------------------------------------------------------------*
 | equilibrate matrix rows                                   fang 01/18 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Mono::EquilibrateMatrixRows(
    LINALG::SparseMatrix& matrix, const Teuchos::RCP<Epetra_Vector>& invrowsums) const
{
  if (matrix.LeftScale(*invrowsums)) dserror("Row equilibration of matrix failed!");
}  // SSI::SSI_Mono::EquilibrateMatrixRows

/*----------------------------------------------------------------------*
 | equilibrate global system of equations if necessary       fang 01/18 |
 *----------------------------------------------------------------------*/
void SSI::SSI_Mono::EquilibrateSystem(const Teuchos::RCP<LINALG::SparseOperator>& systemmatrix,
    const Teuchos::RCP<Epetra_Vector>& residual) const
{
  switch (equilibration_method_)
  {
    case INPAR::SSI::EquilibrationMethod::none:
    {
      // do nothing
      break;
    }

    case INPAR::SSI::EquilibrationMethod::rows_full:
    case INPAR::SSI::EquilibrationMethod::rows_maindiag:
    {
      // initialize vector for inverse sums of absolute values of matrix row entries
      const auto invrowsums = LINALG::CreateVector(*DofRowMap());

      switch (matrixtype_)
      {
        case INPAR::SSI::matrix_block:
        {
          // check matrix
          const auto blocksparsematrix =
              Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix);
          if (blocksparsematrix == Teuchos::null)
            dserror("System matrix is not a block sparse matrix!");

          // perform row equilibration
          for (int i = 0; i < blocksparsematrix->Rows(); ++i)
          {
            // initialize vector for inverse row sums
            const auto invrowsums_block(
                Teuchos::rcp(new Epetra_Vector(blocksparsematrix->Matrix(i, i).RowMap())));

            // compute inverse row sums of current main diagonal matrix block
            if (equilibration_method_ == INPAR::SSI::EquilibrationMethod::rows_maindiag)
              ComputeInvRowSums(blocksparsematrix->Matrix(i, i), invrowsums_block);

            // compute inverse row sums of current row block of global system matrix
            else
            {
              // loop over all column blocks of global system matrix
              for (int j = 0; j < blocksparsematrix->Cols(); ++j)
              {
                // extract current block of global system matrix
                const auto& matrix = blocksparsematrix->Matrix(i, j);

                // loop over all rows of current matrix block
                for (int irow = 0; irow < matrix.RowMap().NumMyElements(); ++irow)
                {
                  // determine length of current matrix row
                  const int length = matrix.EpetraMatrix()->NumMyEntries(irow);

                  if (length > 0)
                  {
                    // extract current matrix row from matrix block
                    int numentries(0);
                    std::vector<double> values(length, 0.0);
                    if (matrix.EpetraMatrix()->ExtractMyRowCopy(
                            irow, length, numentries, &values[0]))
                      dserror(
                          "Cannot extract matrix row with local ID %d from matrix block!", irow);

                    // compute and store current row sum
                    double rowsum(0.0);
                    for (int ientry = 0; ientry < numentries; ++ientry)
                      rowsum += std::abs(values[ientry]);
                    (*invrowsums_block)[irow] += rowsum;
                  }
                }
              }

              // invert row sums of current matrix row block
              if (invrowsums_block->Reciprocal(*invrowsums_block))
                dserror("Vector could not be inverted!");
            }

            // perform row equilibration of matrix blocks in current row block of global system
            // matrix
            for (int j = 0; j < blocksparsematrix->Cols(); ++j)
              EquilibrateMatrixRows(blocksparsematrix->Matrix(i, j), invrowsums_block);

            // insert inverse row sums of current main diagonal matrix block into global vector
            maps_systemmatrix_->InsertVector(invrowsums_block, i, invrowsums);
          }

          break;
        }

        case INPAR::SSI::matrix_sparse:
        {
          // check matrix
          const auto sparsematrix = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix_);
          if (sparsematrix == Teuchos::null) dserror("System matrix is not a sparse matrix!");

          // compute inverse row sums of global system matrix
          ComputeInvRowSums(*sparsematrix, invrowsums);

          // perform row equilibration of global system matrix
          EquilibrateMatrixRows(*sparsematrix, invrowsums);

          break;
        }

        default:
        {
          dserror(
              "Type of global system matrix for scalar-structure interaction not "
              "recognized!");
          break;
        }
      }

      // perform equilibration of global residual vector
      if (residual->Multiply(1.0, *invrowsums, *residual, 0.0))
        dserror("Equilibration of global residual vector failed!");

      break;
    }

    default:
    {
      dserror("Equilibration method not yet implemented!");
      break;
    }
  }
}  // SSI::SSI_Mono::EquilibrateSystem


/*--------------------------------------------------------------------------*
 | finite difference check for global system matrix              fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::FDCheck()
{
  // initial screen output
  if (Comm().MyPID() == 0)
    std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR SSI SYSTEM MATRIX" << std::endl;

  // create global state vector
  Teuchos::RCP<Epetra_Vector> statenp(LINALG::CreateVector(*DofRowMap(), true));
  maps_->InsertVector(scatra_->ScaTraField()->Phinp(), 0, statenp);
  maps_->InsertVector(structure_->Dispnp(), 1, statenp);

  // make a copy of global state vector to undo perturbations later
  Teuchos::RCP<Epetra_Vector> statenp_original = Teuchos::rcp(new Epetra_Vector(*statenp));

  // make a copy of system matrix as Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> sysmat_original = Teuchos::null;
  if (Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(systemmatrix_) != Teuchos::null)
    sysmat_original =
        (new LINALG::SparseMatrix(*(Teuchos::rcp_static_cast<LINALG::SparseMatrix>(systemmatrix_))))
            ->EpetraMatrix();
  else if (Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_) != Teuchos::null)
    sysmat_original =
        (new LINALG::SparseMatrix(
             *(Teuchos::rcp_static_cast<LINALG::BlockSparseMatrixBase>(systemmatrix_)->Merge())))
            ->EpetraMatrix();
  else
    dserror("Type of system matrix unknown!");
  sysmat_original->FillComplete();

  // make a copy of system right-hand side vector
  Teuchos::RCP<Epetra_Vector> rhs_original = Teuchos::rcp(new Epetra_Vector(*residual_));

  // copy and zero out system increment vector if necessary
  Teuchos::RCP<Epetra_Vector> increment_original(Teuchos::null);
  if (iter_ != 1)
  {
    increment_original = Teuchos::rcp(new Epetra_Vector(*increment_));
    increment_->PutScalar(0.);
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
      collid = icoup_structure_->SlaveDofMap()->LID(colgid);
      Comm().MaxAll(&collid, &maxcollid, 1);
      if (maxcollid >= 0) continue;
    }

    // fill global state vector with original state variables
    statenp->Update(1., *statenp_original, 0.);

    // impose perturbation
    if (statenp->Map().MyGID(colgid))
      if (statenp->SumIntoGlobalValue(colgid, 0, scatra_->ScaTraField()->FDCheckEps()))
        dserror("Perturbation could not be imposed on state vector for finite difference check!");
    if (SSIInterfaceMeshtying() and icoup_structure_->PermMasterDofMap()->MyGID(colgid))
      if (statenp->SumIntoGlobalValue(icoup_structure_->SlaveDofMap()->GID(
                                          icoup_structure_->PermMasterDofMap()->LID(colgid)),
              0, scatra_->ScaTraField()->FDCheckEps()))
        dserror("Perturbation could not be imposed on state vector for finite difference check!");
    scatra_->ScaTraField()->Phinp()->Update(1., *maps_->ExtractVector(statenp, 0), 0.);
    structure_->SetState(maps_->ExtractVector(statenp, 1));

    // calculate element right-hand side vector for perturbed state
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
      if (scatra_->ScaTraField()->DirichMaps()->CondMap()->MyGID(rowgid) or
          structure_->GetDBCMapExtractor()->CondMap()->MyGID(rowgid) or
          (SSIInterfaceMeshtying() and icoup_structure_->SlaveDofMap()->MyGID(rowgid)))
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
      const double fdval = -(*residual_)[rowlid] / scatra_->ScaTraField()->FDCheckEps() +
                           (*rhs_original)[rowlid] / scatra_->ScaTraField()->FDCheckEps();

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
      if (abs(relerr1) > scatra_->ScaTraField()->FDCheckTol())
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
        const double left = entry - (*rhs_original)[rowlid] / scatra_->ScaTraField()->FDCheckEps();

        // right-hand side in second comparison
        const double right = -(*residual_)[rowlid] / scatra_->ScaTraField()->FDCheckEps();

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
        if (abs(relerr2) > scatra_->ScaTraField()->FDCheckTol())
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
      printf(
          "--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR "
          "%+12.5e\n\n",
          maxabserrglobal, maxrelerrglobal);
  }

  // undo perturbations of state variables
  scatra_->ScaTraField()->Phinp()->Update(1., *maps_->ExtractVector(statenp_original, 0), 0.);
  structure_->SetState(maps_->ExtractVector(statenp_original, 1));

  // recompute system matrix and right-hand side vector based on original state variables
  AssembleMatAndRHS();

  // restore system increment vector if necessary
  if (increment_original != Teuchos::null) increment_ = increment_original;
}

/*--------------------------------------------------------------------------*
 | initialize monolithic algorithm                               fang 08/17 |
 *--------------------------------------------------------------------------*/
int SSI::SSI_Mono::Init(const Epetra_Comm& comm,     //!< communicator
    const Teuchos::ParameterList& globaltimeparams,  //!< parameter list for time integration
    const Teuchos::ParameterList& scatraparams,      //!< parameter list for scalar transport
    const Teuchos::ParameterList& structparams,      //!< parameter list for structure
    const std::string struct_disname,                //!< name of structural discretization
    const std::string scatra_disname,                //!< name of scalar transport discretization
    bool isAle                                       //!< flag for ALE
)
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
 | output solution to screen and files                           fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Output()
{
  // output scalar transport field
  scatra_->ScaTraField()->Output();

  // output structure field
  structure_->Output();
}

/*--------------------------------------------------------------------------*
 | prepare time step                                             fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::PrepareTimeStep()
{
  // update time and time step
  IncrementTimeAndStep();

  // pass structural degrees of freedom to scalar transport discretization
  SetStructSolution(structure_->Dispnp(), structure_->Velnp());

  // prepare time step for scalar transport field
  scatra_->ScaTraField()->PrepareTimeStep();

  // pass scalar transport degrees of freedom to structural discretization
  // has to be called AFTER scatra_->ScaTraField()->PrepareTimeStep() to ensure
  // consistent scalar transport state vector with valid Dirichlet conditions
  SetScatraSolution(scatra_->ScaTraField()->Phinp());

  // prepare time step for structural field
  structure_->PrepareTimeStep();

  // print time step information to screen
  scatra_->ScaTraField()->PrintTimeStepInfo();
}

/*--------------------------------------------------------------------------*
 | setup monolithic algorithm                                    fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Setup()
{
  // call base class routine
  SSI_Base::Setup();

  // safety checks
  if (scatra_->ScaTraField()->NumScal() != 1)
    dserror(
        "Since the ssi_monolithic framework is only implemented for usage in combination with "
        "volume change laws 'MAT_InelasticDefgradLinScalarIso' or "
        "'MAT_InelasticDefgradLinScalarAniso' so far and these laws are implemented for only "
        "one transported scalar at the moment it is not reasonable to use them with more than one "
        "transported scalar. So you need to cope with it or change implementation! ;-)");

  if (scatra_->ScaTraField()->EquilibrationMethod() != INPAR::SCATRA::EquilibrationMethod::none)
    dserror(
        "You are within the monolithic solid scatra interaction framework but activated a pure "
        "scatra equilibration method. Delete this from 'SCALAR TRANSPORT DYNAMIC' section and set "
        "it in 'SSI CONTROL/MONOLITHIC' instead.");

  if (!scatra_->ScaTraField()->IsIncremental())
    dserror("Must have incremental solution approach for monolithic scalar-structure interaction!");

  // set up scatra-scatra interface meshtying if set in input file
  if (SSIInterfaceMeshtying())
  {
    // extract meshtying strategy for scatra-scatra interface coupling on scatra discretization
    meshtying_strategy_s2i_ = Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(
        scatra_->ScaTraField()->Strategy());

    // safety checks
    if (meshtying_strategy_s2i_ == Teuchos::null)
      dserror("Invalid scatra-scatra interface coupling strategy!");
    if (meshtying_strategy_s2i_->CouplingType() != INPAR::S2I::coupling_matching_nodes)
      dserror(
          "Monolithic scalar-structure interaction only implemented for scatra-scatra "
          "interface "
          "coupling with matching interface nodes!");
  }
}

/*--------------------------------------------------------------------------*
 | setup global system of equations                              fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::SetupSystem()
{
  // merge slave and master side block maps for interface matrix for thermo and scatra
  Teuchos::RCP<Epetra_Map> interface_map_scatra(Teuchos::null);
  Teuchos::RCP<LINALG::MultiMapExtractor> blockmapscatrainterface(Teuchos::null);
  if (SSIInterfaceMeshtying())
  {
    // check whether slave-side degrees of freedom are Dirichlet-free
    std::vector<Teuchos::RCP<const Epetra_Map>> maps(2, Teuchos::null);
    maps[0] = icoup_structure_->SlaveDofMap();
    maps[1] = structure_->GetDBCMapExtractor()->CondMap();
    if (LINALG::MultiMapExtractor::IntersectMaps(maps)->NumGlobalElements() > 0)
      dserror("Must not apply Dirichlet conditions to slave-side structural displacements!");

    interface_map_scatra = LINALG::MultiMapExtractor::MergeMaps(
        {meshtying_strategy_s2i_->CouplingAdapter()->MasterDofMap(),
            meshtying_strategy_s2i_->CouplingAdapter()->SlaveDofMap()});
  }

  // initialize global map extractor
  maps_ = Teuchos::rcp(new LINALG::MapExtractor(
      *LINALG::MergeMap(*scatra_->ScaTraField()->DofRowMap(), *structure_->DofRowMap(), false),
      structure_->DofRowMap(), scatra_->ScaTraField()->DofRowMap()));

  // check global map extractor
  maps_->CheckForValidMapExtractor();

  // initialize global increment vector for Newton-Raphson iteration
  increment_ = LINALG::CreateVector(*DofRowMap(), true);

  // initialize global residual vector
  residual_ = LINALG::CreateVector(*DofRowMap(), true);

  // initialize map extractors associated with blocks of global system matrix
  switch (scatra_->ScaTraField()->MatrixType())
  {
    // one single main-diagonal matrix block associated with scalar transport field
    case INPAR::SCATRA::MatrixType::sparse:
    {
      maps_systemmatrix_ = maps_;
      if (SSIInterfaceMeshtying())
      {
        blockmapscatrainterface = Teuchos::rcp(new LINALG::MultiMapExtractor(*interface_map_scatra,
            std::vector<Teuchos::RCP<const Epetra_Map>>(1, interface_map_scatra)));
        blockmapscatrainterface->CheckForValidMapExtractor();
      }
      break;
    }

    // several main-diagonal matrix blocks associated with scalar transport field
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      // store an RCP to the block maps of the scatra field
      maps_scatra_ = Teuchos::rcpFromRef(scatra_->ScaTraField()->BlockMaps());

      // safety check
      maps_scatra_->CheckForValidMapExtractor();

      // extract maps underlying main-diagonal matrix blocks associated with scalar transport  field
      const int maps_systemmatrix_scatra = maps_scatra_->NumMaps();
      std::vector<Teuchos::RCP<const Epetra_Map>> maps_systemmatrix(maps_systemmatrix_scatra + 1);
      for (int imap = 0; imap < maps_systemmatrix_scatra; ++imap)
        maps_systemmatrix[imap] = maps_scatra_->Map(imap);

      // extract map underlying single main-diagonal matrix block associated with structural
      // field
      maps_systemmatrix[maps_systemmatrix_scatra] = structure_->DofRowMap();

      // initialize map extractor associated with blocks of global system matrix
      maps_systemmatrix_ =
          Teuchos::rcp(new LINALG::MultiMapExtractor(*DofRowMap(), maps_systemmatrix));

      // initialize map extractor associated with all degrees of freedom inside structural
      // field
      map_structure_ =
          Teuchos::rcp(new LINALG::MultiMapExtractor(*structure_->Discretization()->DofRowMap(),
              std::vector<Teuchos::RCP<const Epetra_Map>>(1, structure_->DofRowMap())));

      // safety check
      map_structure_->CheckForValidMapExtractor();
      if (SSIInterfaceMeshtying())
      {
        // build block map for scatra interface by merging slave and master side for each block
        std::vector<Teuchos::RCP<const Epetra_Map>> partial_blockmapscatrainterface(
            maps_systemmatrix_scatra, Teuchos::null);
        for (int iblockmap = 0; iblockmap < maps_systemmatrix_scatra; ++iblockmap)
          partial_blockmapscatrainterface.at(iblockmap) = LINALG::MultiMapExtractor::MergeMaps(
              {meshtying_strategy_s2i_->BlockMapsSlave().Map(iblockmap),
                  meshtying_strategy_s2i_->BlockMapsMaster().Map(iblockmap)});
        blockmapscatrainterface = Teuchos::rcp(
            new LINALG::MultiMapExtractor(*interface_map_scatra, partial_blockmapscatrainterface));
        blockmapscatrainterface->CheckForValidMapExtractor();
      }

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
    case INPAR::SSI::matrix_block:
    {
      // safety check
      if (!solver_->Params().isSublist("AMGnxn Parameters"))
        dserror(
            "Global system matrix with block structure requires AMGnxn block "
            "preconditioner!");

      // initialize global system matrix
      systemmatrix_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *maps_systemmatrix_, *maps_systemmatrix_, 81, false, true));

      // feed AMGnxn block preconditioner with null space information for each block of global
      // block system matrix
      BuildNullSpaces();

      break;
    }

    case INPAR::SSI::matrix_sparse:
    {
      // safety check
      if (scatra_->ScaTraField()->SystemMatrix() == Teuchos::null)
        dserror("Incompatible matrix type associated with scalar transport field!");

      // initialize global system matrix
      systemmatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMap(), 27, false, true));

      break;
    }

    default:
    {
      dserror("Type of global system matrix for scalar-structure interaction not recognized!");
      break;
    }
  }

  // initialize scatra-structure block and structure-scatra block of global system matrix
  switch (scatra_->ScaTraField()->MatrixType())
  {
    case INPAR::SCATRA::MatrixType::block_condition:
    {
      // initialize scatra-structure block
      scatrastructuredomain_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              *map_structure_, scatra_->ScaTraField()->BlockMaps(), 81, false, true));

      // initialize scatra-structure block
      if (SSIInterfaceMeshtying())
        scatrastructureinterface_ =
            Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                *map_structure_, *blockmapscatrainterface, 81, false, true));

      // initialize structure-scatra block
      structurescatradomain_ =
          Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
              scatra_->ScaTraField()->BlockMaps(), *map_structure_, 81, false, true));

      break;
    }

    case INPAR::SCATRA::MatrixType::sparse:
    {
      // initialize scatra-structure block
      scatrastructuredomain_ = Teuchos::rcp(
          new LINALG::SparseMatrix(*scatra_->ScaTraField()->DofRowMap(), 27, false, true));

      // initialize scatra-structure block
      if (SSIInterfaceMeshtying())
        scatrastructureinterface_ =
            Teuchos::rcp(new LINALG::SparseMatrix(*interface_map_scatra, 27, false, true));

      // initialize structure-scatra block
      structurescatradomain_ =
          Teuchos::rcp(new LINALG::SparseMatrix(*structure_->DofRowMap(), 27, false, true));

      break;
    }

    default:
    {
      dserror("Invalid matrix type associated with scalar transport field!");
      break;
    }
  }

  // initialize strategy for assembly
  // create converter
  ADAPTER::CouplingSlaveConverter converter(*icoup_structure_);
  switch (matrixtype_)
  {
    case INPAR::SSI::matrix_block:
    {
      switch (scatra_->ScaTraField()->MatrixType())
      {
        case INPAR::SCATRA::MatrixType::block_condition:
        {
          strategy_assemble_ = Teuchos::rcp(
              new SSI::AssembleStrategyBlockBlock(Teuchos::rcp(this, false), converter));
          break;
        }
        case INPAR::SCATRA::MatrixType::sparse:
        {
          strategy_assemble_ = Teuchos::rcp(
              new SSI::AssembleStrategyBlockSparse(Teuchos::rcp(this, false), converter));
          break;
        }

        default:
        {
          dserror("unknown matrix type");
          break;
        }
      }
      break;
    }
    case INPAR::SSI::matrix_sparse:
    {
      strategy_assemble_ =
          Teuchos::rcp(new SSI::AssembleStrategySparse(Teuchos::rcp(this, false), converter));
      break;
    }
    default:
    {
      dserror("unknown matrix type");
      break;
    }
  }

  // initialize object, that performs evaluations of OD coupling
  scatrastructureOffDiagcoupling_ = Teuchos::rcp(new SSI::ScatraStructureOffDiagCoupling(
      maps_scatra_, maps_structure_, maps_->Map(0), maps_->Map(1), icoup_structure_,
      interface_map_scatra, meshtying_strategy_s2i_, scatra_, structure_));
}

/*---------------------------------------------------------------------------------*
 | set up structural model evaluator for scalar-structure interaction   fang 01/18 |
 *---------------------------------------------------------------------------------*/
void SSI::SSI_Mono::SetupModelEvaluator() const
{
  // construct and register structural model evaluator if necessary
  if (DRT::INPUT::IntegralValue<INPAR::STR::StressType>(
          DRT::Problem::Instance()->IOParams(), "STRUCT_STRESS") != INPAR::STR::stress_none and
      SSIInterfaceMeshtying())
    struct_adapterbase_ptr_->RegisterModelEvaluator("Monolithic Coupling Model",
        Teuchos::rcp(new STR::MODELEVALUATOR::MonolithicSSI(Teuchos::rcp(this, false))));
}

/*--------------------------------------------------------------------------*
 | evaluate time step using Newton-Raphson iteration             fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Solve()
{
  // initialize counter for Newton-Raphson iteration
  iter_ = 0;

  // start Newton-Raphson iteration
  while (true)
  {
    // update iteration counter
    ++iter_;

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
    if (!systemmatrix_->Filled())
      dserror("Complete() has not been called on global system matrix yet!");

    // perform finite difference check on time integrator level
    if ((scatra_->ScaTraField()->FDCheckType() == INPAR::SCATRA::fdcheck_global) and (Step() > 1))
      FDCheck();

    // check termination criterion for Newton-Raphson iteration
    if (strategy_convcheck_->ExitNewtonRaphson(*this)) break;

    // initialize global increment vector
    increment_->PutScalar(0.);

    // store time before solving global system of equations
    time = timer_->WallTime();

    // equilibrate global system of equations if necessary
    EquilibrateSystem(systemmatrix_, residual_);

    // solve global system of equations
    // Dirichlet boundary conditions have already been applied to global system of equations
    solver_->Solve(systemmatrix_->EpetraOperator(), increment_, residual_, true, iter_ == 1);

    // determine time needed for solving global system of equations,
    // and take maximum over all processors via communication
    double mydtsolve = timer_->WallTime() - time;
    Comm().MaxAll(&mydtsolve, &dtsolve_, 1);

    // output performance statistics associated with linear solver into text file if
    // applicable
    if (DRT::INPUT::IntegralValue<bool>(
            *scatra_->ScaTraField()->ScatraParameterList(), "OUTPUTLINSOLVERSTATS"))
      scatra_->ScaTraField()->OutputLinSolverStats(
          *solver_, dtsolve_, Step(), iter_, residual_->Map().NumGlobalElements());

    // update scalar transport field
    scatra_->ScaTraField()->UpdateIter(maps_->ExtractVector(increment_, 0));
    scatra_->ScaTraField()->ComputeIntermediateValues();

    // structure field is updated during the next Newton-Raphson iteration step
  }  // Newton-Raphson iteration
}

/*--------------------------------------------------------------------------*
 | time loop                                                     fang 08/17 |
 *--------------------------------------------------------------------------*/
void SSI::SSI_Mono::Timeloop()
{
  // output initial scalar transport solution to screen and files
  if (Step() == 0)
  {
    SetStructSolution(structure_->Dispnp(), structure_->Velnp());
    scatra_->ScaTraField()->Output();
  }

  // time loop
  while (NotFinished() and scatra_->ScaTraField()->NotFinished())
  {
    // prepare time step
    PrepareTimeStep();

    // store time before calling nonlinear solver
    const double time = timer_->WallTime();

    // evaluate time step
    Solve();

    // determine time spent by nonlinear solver and take maximum over all processors via
    // communication
    double mydtnonlinsolve(timer_->WallTime() - time), dtnonlinsolve(0.);
    Comm().MaxAll(&mydtnonlinsolve, &dtnonlinsolve, 1);

    // output performance statistics associated with nonlinear solver into *.csv file if
    // applicable
    if (DRT::INPUT::IntegralValue<int>(
            *scatra_->ScaTraField()->ScatraParameterList(), "OUTPUTNONLINSOLVERSTATS"))
      scatra_->ScaTraField()->OutputNonlinSolverStats(iter_, dtnonlinsolve, Step(), Comm());

    // prepare structure output
    structure_->PrepareOutput();

    // update scalar transport and structure fields
    Update();

    // output solution to screen and files
    Output();
  }  // while(NotFinished())
}

/*--------------------------------------------------------------------------------------*
 | update scalar transport and structure fields after time step evaluation   fang 08/17 |
 *--------------------------------------------------------------------------------------*/
void SSI::SSI_Mono::Update()
{
  // update scalar transport field
  scatra_->ScaTraField()->Update();

  // update structure field
  structure_->Update();
}
