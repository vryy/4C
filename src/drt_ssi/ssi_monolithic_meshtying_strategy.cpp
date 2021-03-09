/*----------------------------------------------------------------------*/
/*! \file
\brief Mesh tying strategy for monolithic SSI

\level 2

 */
/*----------------------------------------------------------------------*/

#include "ssi_monolithic_meshtying_strategy.H"

#include "Epetra_Map.h"
#include "ssi_monolithic.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"

#include "../drt_lib/drt_locsys.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_matrixtransform.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
SSI::MeshtyingStrategyBase::MeshtyingStrategyBase(const SSI::SSIMono& ssi_mono)
    : meshtying_3_domain_intersection_(
          ssi_mono.SSIInterfaceMeshtying() and ssi_mono.Meshtying3DomainIntersection()),
      ssi_mono_(ssi_mono)
{
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Epetra_Vector SSI::MeshtyingStrategyBase::ApplyMeshtyingToStructureRHS(
    Teuchos::RCP<const Epetra_Vector> structure_rhs)
{
  // make copy of structure right-hand side vector
  Epetra_Vector rhs_structure(*structure_rhs);

  // transform slave-side part of structure right-hand side vector to master side
  const auto rhs_structure_only_slave_dofs =
      SSIMono().MapsCoupStruct()->ExtractVector(rhs_structure, 1);

  const auto rhs_structure_only_master_dofs =
      SSIMono().InterfaceCouplingAdapterStructure()->SlaveToMaster(rhs_structure_only_slave_dofs);

  auto rhs_structure_master =
      SSIMono().MapsCoupStruct()->InsertVector(rhs_structure_only_master_dofs, 2);

  if (Meshtying3DomainIntersection())
  {
    const auto rhs_structure_3_domain_intersection_only_slave_dofs =
        SSIMono().MapsCoupStruct3DomainIntersection()->ExtractVector(rhs_structure, 1);

    const auto rhs_structure_3_domain_intersection_only_master_dofs =
        SSIMono().InterfaceCouplingAdapterStructure3DomainIntersection()->SlaveToMaster(
            rhs_structure_3_domain_intersection_only_slave_dofs);

    const auto rhs_structure_3_domain_intersection_master =
        SSIMono().MapsCoupStruct3DomainIntersection()->InsertVector(
            rhs_structure_3_domain_intersection_only_master_dofs, 2);

    rhs_structure_master->Update(1.0, *rhs_structure_3_domain_intersection_master, 1.0);
  }

  // locsys manager of structure
  const auto& locsysmanager_structure = SSIMono().StructureField()->LocsysManager();

  // apply pseudo Dirichlet conditions to master-side part of structure right-hand side vector
  const auto zeros_structure_master = Teuchos::rcp(new Epetra_Vector(rhs_structure_master->Map()));

  if (locsysmanager_structure != Teuchos::null)
    locsysmanager_structure->RotateGlobalToLocal(rhs_structure_master);

  LINALG::ApplyDirichlettoSystem(rhs_structure_master, zeros_structure_master,
      *SSIMono().StructureField()->GetDBCMapExtractor()->CondMap());

  if (locsysmanager_structure != Teuchos::null)
    locsysmanager_structure->RotateLocalToGlobal(rhs_structure_master);

  // assemble master-side part of structure right-hand side vector
  rhs_structure.Update(1.0, *rhs_structure_master, 1.0);

  // zero out slave-side part of structure right-hand side vector
  SSIMono().MapsCoupStruct()->PutScalar(rhs_structure, 1, 0.0);
  if (Meshtying3DomainIntersection())
    SSIMono().MapsCoupStruct3DomainIntersection()->PutScalar(rhs_structure, 1, 0.0);

  return rhs_structure;
}

/*-------------------------------------------------------------------------*
 *-------------------------------------------------------------------------*/
Teuchos::RCP<SSI::MeshtyingStrategyBase> SSI::BuildMeshtyingStrategy(const SSI::SSIMono& ssi_mono)
{
  return Teuchos::rcp(new SSI::MeshtyingStrategyBase(ssi_mono));
}