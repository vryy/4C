/*----------------------------------------------------------------------*/
/*! \file
\brief Dirichlet boundary condition handler for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/

#include "baci_ssi_monolithic_dbc_handler.hpp"

#include "baci_adapter_str_ssiwrapper.hpp"
#include "baci_lib_locsys.hpp"
#include "baci_linalg_blocksparsematrix.hpp"
#include "baci_linalg_matrixtransform.hpp"
#include "baci_linalg_utils_sparse_algebra_assemble.hpp"
#include "baci_scatra_timint_implicit.hpp"
#include "baci_ssi_monolithic.hpp"
#include "baci_ssi_utils.hpp"

BACI_NAMESPACE_OPEN

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerBase::DBCHandlerBase(const bool is_scatra_manifold,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure)
    : is_scatra_manifold_(is_scatra_manifold),
      scatra_(std::move(scatra)),
      scatra_manifold_(std::move(scatra_manifold)),
      ssi_maps_(std::move(ssi_maps)),
      structure_(std::move(structure))
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerSparse::DBCHandlerSparse(const bool is_scatra_manifold,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure)
    : DBCHandlerBase(is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerBlock::DBCHandlerBlock(const bool is_scatra_manifold,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure)
    : DBCHandlerBase(is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure),
      position_structure_(SSIMaps()->GetBlockPositions(SSI::Subproblem::structure).at(0))
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBase::ApplyDBCToRHS(Teuchos::RCP<Epetra_Vector> rhs)
{
  // apply Dirichlet boundary conditions to the structure part of the right hand side
  const auto& locsysmanager_structure = StructureField()->LocsysManager();
  auto rhs_struct = SSIMaps()->MapsSubProblems()->ExtractVector(
      rhs, UTILS::SSIMaps::GetProblemPosition(SSI::Subproblem::structure));
  const auto zeros_struct =
      Teuchos::rcp(new Epetra_Vector(*StructureField()->GetDBCMapExtractor()->CondMap()));

  if (locsysmanager_structure != Teuchos::null)
    locsysmanager_structure->RotateGlobalToLocal(rhs_struct);
  CORE::LINALG::ApplyDirichletToSystem(
      *rhs_struct, *zeros_struct, *StructureField()->GetDBCMapExtractor()->CondMap());
  if (locsysmanager_structure != Teuchos::null)
    locsysmanager_structure->RotateLocalToGlobal(rhs_struct);

  SSIMaps()->MapsSubProblems()->InsertVector(
      rhs_struct, UTILS::SSIMaps::GetProblemPosition(SSI::Subproblem::structure), rhs);

  // apply Dirichlet boundary conditions to the scatra part of the right hand side
  const auto zeros_scatra =
      Teuchos::rcp(new Epetra_Vector(*ScaTraField()->DirichMaps()->CondMap()));
  CORE::LINALG::ApplyDirichletToSystem(
      *rhs, *zeros_scatra, *ScaTraField()->DirichMaps()->CondMap());

  // apply Dirichlet boundary conditions to the scatra manifold part of the right hand side
  if (IsScaTraManifold())
  {
    const auto zeros_scatramanifold =
        Teuchos::rcp(new Epetra_Vector(*ScaTraManifoldField()->DirichMaps()->CondMap()));
    CORE::LINALG::ApplyDirichletToSystem(
        *rhs, *zeros_scatramanifold, *ScaTraManifoldField()->DirichMaps()->CondMap());
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBase::ApplyDBCToSystemMatrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix)
{
  // apply the scalar transport Dirichlet boundary conditions to the global system matrix
  system_matrix->ApplyDirichlet(*ScaTraField()->DirichMaps()->CondMap(), true);

  // apply the scalar transport on manifolds Dirichlet boundary conditions to the global system
  // matrix
  if (IsScaTraManifold())
    system_matrix->ApplyDirichlet(*ScaTraManifoldField()->DirichMaps()->CondMap(), true);

  // apply the structure Dirichlet boundary conditions to the global system matrix
  ApplyStructureDBCToSystemMatrix(system_matrix);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBase::ApplyStructureDBCToSystemMatrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix)
{
  // locsys manager of structure
  const auto& locsysmanager_structure = StructureField()->LocsysManager();

  // map of structure Dirichlet BCs
  const auto& dbcmap_structure = StructureField()->GetDBCMapExtractor()->CondMap();

  if (locsysmanager_structure == Teuchos::null)
    system_matrix->ApplyDirichlet(*dbcmap_structure);
  else
    ApplyStructureDBCWithLocSysRotationToSystemMatrix(
        system_matrix, dbcmap_structure, locsysmanager_structure);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerSparse::ApplyStructureDBCWithLocSysRotationToSystemMatrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix,
    const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
    Teuchos::RCP<const DRT::UTILS::LocsysManager> locsysmanager_structure)

{
  auto systemmatrix_sparse = CORE::LINALG::CastToSparseMatrixAndCheckSuccess(system_matrix);

  // structure dof row map
  const auto& dofrowmap_structure = StructureField()->DofRowMap();

  // extract structure rows of global system matrix
  const auto systemmatrix_structure =
      Teuchos::rcp(new CORE::LINALG::SparseMatrix(*dofrowmap_structure, 27, false, true));
  CORE::LINALG::MatrixLogicalSplitAndTransform()(*systemmatrix_sparse, *dofrowmap_structure,
      system_matrix->DomainMap(), 1.0, nullptr, nullptr, *systemmatrix_structure);
  systemmatrix_structure->Complete(system_matrix->DomainMap(), *dofrowmap_structure);

  // apply structure Dirichlet conditions
  locsysmanager_structure->RotateGlobalToLocal(systemmatrix_structure);
  systemmatrix_structure->ApplyDirichletWithTrafo(
      *locsysmanager_structure->Trafo(), *dbcmap_structure);
  locsysmanager_structure->RotateLocalToGlobal(systemmatrix_structure);

  // assemble structure rows of global system matrix back into global system matrix
  systemmatrix_sparse->Put(*systemmatrix_structure, 1.0, dofrowmap_structure);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBlock::ApplyStructureDBCWithLocSysRotationToSystemMatrix(
    Teuchos::RCP<CORE::LINALG::SparseOperator> system_matrix,
    const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
    Teuchos::RCP<const DRT::UTILS::LocsysManager> locsysmanager_structure)
{
  auto systemmatrix_block = CORE::LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(system_matrix);

  // apply structure Dirichlet conditions
  for (int iblock = 0; iblock < systemmatrix_block->Cols(); ++iblock)
  {
    locsysmanager_structure->RotateGlobalToLocal(
        Teuchos::rcp(&systemmatrix_block->Matrix(PositionStructure(), iblock), false));
    systemmatrix_block->Matrix(PositionStructure(), iblock)
        .ApplyDirichletWithTrafo(
            *locsysmanager_structure->Trafo(), *dbcmap_structure, (iblock == PositionStructure()));
    locsysmanager_structure->RotateLocalToGlobal(
        Teuchos::rcp(&systemmatrix_block->Matrix(PositionStructure(), iblock), false));
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::DBCHandlerBase> SSI::BuildDBCHandler(const bool is_scatra_manifold,
    CORE::LINALG::MatrixType matrixtype_ssi, Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<ADAPTER::SSIStructureWrapper> structure)
{
  Teuchos::RCP<SSI::DBCHandlerBase> dbc_handler = Teuchos::null;

  switch (matrixtype_ssi)
  {
    case CORE::LINALG::MatrixType::block_field:
    {
      dbc_handler = Teuchos::rcp(new SSI::DBCHandlerBlock(
          is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure));
      break;
    }
    case CORE::LINALG::MatrixType::sparse:
    {
      dbc_handler = Teuchos::rcp(new SSI::DBCHandlerSparse(
          is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure));
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
BACI_NAMESPACE_CLOSE
