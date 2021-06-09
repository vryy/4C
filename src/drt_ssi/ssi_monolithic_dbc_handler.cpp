/*----------------------------------------------------------------------*/
/*! \file
\brief Dirichlet boundary condition handler for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/

#include "ssi_monolithic_dbc_handler.H"

#include "ssi_monolithic.H"
#include "ssi_utils.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"

#include "../drt_lib/drt_locsys.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerBase::DBCHandlerBase(const bool is_scatra_manifold,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure)
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
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure)
    : DBCHandlerBase(is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::DBCHandlerBlock::DBCHandlerBlock(const bool is_scatra_manifold,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure)
    : DBCHandlerBase(is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure),
      position_structure_(-1)
{
  position_structure_ = SSIMaps()->GetBlockPositions(SSI::Subproblem::structure)->at(0);
  // safety check
  if (position_structure_ == -1) dserror("Cannot get position of structure block");
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
  LINALG::ApplyDirichlettoSystem(
      rhs_struct, zeros_struct, *StructureField()->GetDBCMapExtractor()->CondMap());
  if (locsysmanager_structure != Teuchos::null)
    locsysmanager_structure->RotateLocalToGlobal(rhs_struct);

  SSIMaps()->MapsSubProblems()->InsertVector(
      rhs_struct, UTILS::SSIMaps::GetProblemPosition(SSI::Subproblem::structure), rhs);

  // apply Dirichlet boundary conditions to the scatra part of the right hand side
  const auto zeros_scatra =
      Teuchos::rcp(new Epetra_Vector(*ScaTraField()->DirichMaps()->CondMap()));
  LINALG::ApplyDirichlettoSystem(rhs, zeros_scatra, *ScaTraField()->DirichMaps()->CondMap());

  // apply Dirichlet boundary conditions to the scatra manifold part of the right hand side
  if (IsScaTraManifold())
  {
    const auto zeros_scatramanifold =
        Teuchos::rcp(new Epetra_Vector(*ScaTraManifoldField()->DirichMaps()->CondMap()));
    LINALG::ApplyDirichlettoSystem(
        rhs, zeros_scatramanifold, *ScaTraManifoldField()->DirichMaps()->CondMap());
  }
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBase::ApplyDBCToSystemMatrix(Teuchos::RCP<LINALG::SparseOperator> system_matrix)
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
    Teuchos::RCP<LINALG::SparseOperator> system_matrix)
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
    Teuchos::RCP<LINALG::SparseOperator> system_matrix,
    const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
    Teuchos::RCP<const DRT::UTILS::LocsysManager> locsysmanager_structure)

{
  auto systemmatrix_sparse = LINALG::CastToSparseMatrixAndCheckSuccess(system_matrix);

  // structure dof row map
  const auto& dofrowmap_structure = StructureField()->DofRowMap();

  // extract structure rows of global system matrix
  const auto systemmatrix_structure =
      Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap_structure, 27, false, true));
  LINALG::MatrixLogicalSplitAndTransform()(*systemmatrix_sparse, *dofrowmap_structure,
      system_matrix->DomainMap(), 1.0, nullptr, nullptr, *systemmatrix_structure);
  systemmatrix_structure->Complete(system_matrix->DomainMap(), *dofrowmap_structure);

  // apply structure Dirichlet conditions
  locsysmanager_structure->RotateGlobalToLocal(systemmatrix_structure);
  systemmatrix_structure->ApplyDirichletWithTrafo(
      locsysmanager_structure->Trafo(), *dbcmap_structure);
  locsysmanager_structure->RotateLocalToGlobal(systemmatrix_structure);

  // assemble structure rows of global system matrix back into global system matrix
  systemmatrix_sparse->Put(*systemmatrix_structure, 1.0, dofrowmap_structure);
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
void SSI::DBCHandlerBlock::ApplyStructureDBCWithLocSysRotationToSystemMatrix(
    Teuchos::RCP<LINALG::SparseOperator> system_matrix,
    const Teuchos::RCP<const Epetra_Map>& dbcmap_structure,
    Teuchos::RCP<const DRT::UTILS::LocsysManager> locsysmanager_structure)
{
  auto systemmatrix_block = LINALG::CastToBlockSparseMatrixBaseAndCheckSuccess(system_matrix);

  // apply structure Dirichlet conditions
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

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::DBCHandlerBase> SSI::BuildDBCHandler(const bool is_scatra_manifold,
    LINALG::MatrixType matrixtype_ssi, Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,
    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_manifold,
    Teuchos::RCP<SSI::UTILS::SSIMaps> ssi_maps,
    Teuchos::RCP<::ADAPTER::SSIStructureWrapper> structure)
{
  Teuchos::RCP<SSI::DBCHandlerBase> dbc_handler = Teuchos::null;

  switch (matrixtype_ssi)
  {
    case LINALG::MatrixType::block_field:
    {
      dbc_handler = Teuchos::rcp(new SSI::DBCHandlerBlock(
          is_scatra_manifold, scatra, scatra_manifold, ssi_maps, structure));
      break;
    }
    case LINALG::MatrixType::sparse:
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