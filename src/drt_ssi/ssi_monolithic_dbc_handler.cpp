/*----------------------------------------------------------------------*/
/*! \file
\brief Dirichlet boundary condition handler for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/

#include "ssi_monolithic_dbc_handler.H"

#include "ssi_monolithic.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"

#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerBase::DBCHandlerBase(const Teuchos::RCP<const SSI::SSIMono> ssi_mono)
    : ssi_mono_(ssi_mono)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerSparse::DBCHandlerSparse(const Teuchos::RCP<const SSI::SSIMono> ssi_mono)
    : DBCHandlerBase(ssi_mono)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerBlock::DBCHandlerBlock(const Teuchos::RCP<const SSI::SSIMono> ssi_mono)
    : DBCHandlerBase(ssi_mono), position_structure_(-1)
{
  position_structure_ = SSIMono()->GetBlockPositions(SSI::Subproblem::structure)->at(0);
  // safety check
  if (position_structure_ == -1) dserror("Cannot get position of structure block");
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBase::ApplyDBCToRHS(Teuchos::RCP<Epetra_Vector> rhs)
{
  // apply Dirichlet boundary conditions to the structure part of the right hand side
  auto rhs_struct = ssi_mono_->MapsSubProblems()->ExtractVector(
      rhs, SSIMono()->GetProblemPosition(SSI::Subproblem::structure));
  const auto zeros_struct = Teuchos::rcp(new Epetra_Vector(rhs_struct->Map()));
  LINALG::ApplyDirichlettoSystem(
      rhs_struct, zeros_struct, *ssi_mono_->StructureField()->GetDBCMapExtractor()->CondMap());

  // apply Dirichlet boundary conditions to the scatra part of the right hand side
  auto rhs_scatra = ssi_mono_->MapsSubProblems()->ExtractVector(
      rhs, SSIMono()->GetProblemPosition(SSI::Subproblem::scalar_transport));
  const auto zeros_scatra = Teuchos::rcp(new Epetra_Vector(rhs_scatra->Map()));
  LINALG::ApplyDirichlettoSystem(
      rhs_scatra, zeros_scatra, *ssi_mono_->ScaTraField()->DirichMaps()->CondMap());

  // apply Dirichlet boundary conditions to the scatra manifold part of the right hand side
  if (ssi_mono_->IsScaTraManifold())
  {
    auto rhs_scatramanifold = ssi_mono_->MapsSubProblems()->ExtractVector(
        rhs, SSIMono()->GetProblemPosition(SSI::Subproblem::manifold));
    const auto zeros_scatramanifold = Teuchos::rcp(new Epetra_Vector(rhs_scatramanifold->Map()));
    LINALG::ApplyDirichlettoSystem(rhs_scatramanifold, zeros_scatramanifold,
        *ssi_mono_->ScaTraManifold()->DirichMaps()->CondMap());
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBase::ApplyDBCToSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> system_matrix)
{
  // apply the scalar transport Dirichlet boundary conditions to the global system matrix
  system_matrix->ApplyDirichlet(*SSIMono()->ScaTraField()->DirichMaps()->CondMap(), true);

  // apply the scalar transport on manifolds Dirichlet boundary conditions to the global system
  // matrix
  if (ssi_mono_->IsScaTraManifold())
    system_matrix->ApplyDirichlet(*SSIMono()->ScaTraManifold()->DirichMaps()->CondMap(), true);

  // apply the structure Dirichlet boundary conditions to the global system matrix
  ApplyStructuralDBCToSystemMatrix(system_matrix);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerSparse::ApplyStructuralDBCToSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> system_matrix)
{
  // locsys manager of structure
  const auto& locsysmanager_structure = SSIMono()->StructureField()->LocsysManager();

  // map of structural Dirichlet BCs
  const auto& dbcmap_structure = SSIMono()->StructureField()->GetDBCMapExtractor()->CondMap();

  // structural dof row map
  const auto& dofrowmap_structure = SSIMono()->StructureField()->DofRowMap();

  if (locsysmanager_structure == Teuchos::null)
    system_matrix->ApplyDirichlet(*dbcmap_structure);
  else
  {
    auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(system_matrix);

    // extract structural rows of global system matrix
    const Teuchos::RCP<LINALG::SparseMatrix> systemmatrix_structure =
        Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_structure, 27, false, true));
    LINALG::MatrixLogicalSplitAndTransform()(*systemmatrix_sparse, *dofrowmap_structure,
        system_matrix->DomainMap(), 1.0, nullptr, nullptr, *systemmatrix_structure);
    systemmatrix_structure->Complete(system_matrix->DomainMap(), *dofrowmap_structure);

    // apply structural Dirichlet conditions
    locsysmanager_structure->RotateGlobalToLocal(systemmatrix_structure);
    systemmatrix_structure->ApplyDirichletWithTrafo(
        locsysmanager_structure->Trafo(), *dbcmap_structure);
    locsysmanager_structure->RotateLocalToGlobal(systemmatrix_structure);

    // assemble structural rows of global system matrix back into global system matrix
    systemmatrix_sparse->Put(*systemmatrix_structure, 1.0, dofrowmap_structure);
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBlock::ApplyStructuralDBCToSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> system_matrix)
{
  // locsys manager of structure
  const auto& locsysmanager_structure = SSIMono()->StructureField()->LocsysManager();

  // map of structural Dirichlet BCs
  const auto dbcmap_structure = SSIMono()->StructureField()->GetDBCMapExtractor()->CondMap();

  if (locsysmanager_structure == Teuchos::null)
    system_matrix->ApplyDirichlet(*dbcmap_structure);
  else
  {
    auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(system_matrix);

    // apply structural Dirichlet conditions
    for (int iblock = 0; iblock < systemmatrix_block->Cols(); ++iblock)
    {
      locsysmanager_structure->RotateGlobalToLocal(
          Teuchos::rcp(&systemmatrix_block->Matrix(PositionStructure(), iblock), false));
      systemmatrix_block->Matrix(PositionStructure(), iblock)
          .ApplyDirichletWithTrafo(
              locsysmanager_structure->Trafo(), *dbcmap_structure, (iblock == PositionStructure()));
      locsysmanager_structure->RotateLocalToGlobal(
          Teuchos::rcp(&systemmatrix_block->Matrix(PositionStructure(), iblock), false));
    }
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::DBCHandlerBase> SSI::BuildDBCHandler(
    Teuchos::RCP<const SSI::SSIMono> ssi_mono, LINALG::MatrixType matrixtype_ssi)
{
  Teuchos::RCP<SSI::DBCHandlerBase> dbc_handler = Teuchos::null;

  switch (matrixtype_ssi)
  {
    case LINALG::MatrixType::block_field:
    {
      dbc_handler = Teuchos::rcp(new SSI::DBCHandlerBlock(ssi_mono));
      break;
    }
    case LINALG::MatrixType::sparse:
    {
      dbc_handler = Teuchos::rcp(new SSI::DBCHandlerSparse(ssi_mono));
      break;
    }
    default:
    {
      dserror("unknown matrix type of SSI problem");
      break;
    }
  }

  return dbc_handler;
}