/*----------------------------------------------------------------------*/
/*! \file
\file  ad_fbi_constraintbridge_penalty.cpp

\brief Abstract class to be overloaded by different adapter implementations connecting the penalty
constraint enforcement technique with a discretization approach for Fluid-beam interaction.

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/
#include "ad_fbi_constraintbridge_penalty.H"
#include "partitioned_penaltycoupling_assembly_manager_direct.H"
#include <Epetra_FEVector.h>
#include "../linalg/linalg_sparsematrix.H"

void ADAPTER::FBIConstraintBridgePenalty::Setup(
    const Epetra_Map* beam_map, const Epetra_Map* fluid_map)
{
  FBIConstraintBridge::Setup(beam_map, fluid_map);
  fs_ = Teuchos::rcp(new Epetra_FEVector(*beam_map));
  Cmm_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluid_map, 30, true, true,
      LINALG::SparseMatrix::FE_MATRIX));  // todo Is there a better estimator?
  Csm_ = Teuchos::rcp(new LINALG::SparseMatrix(*fluid_map, 30, true, true,
      LINALG::SparseMatrix::FE_MATRIX));  // todo Is there a better estimator?
}
void ADAPTER::FBIConstraintBridgePenalty::Evaluate(
    Teuchos::RCP<const std::vector<Teuchos::RCP<DRT::Discretization>>>
        discretizations)  // todo overload this. Need different assembly
                          // manager for mortar.. dependency injection
                          // or flag?
{
  BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerDirect
      assembly_manager =
          BEAMINTERACTION::SUBMODELEVALUATOR::PartitionedBeamInteractionAssemblyManagerDirect(
              *(GetPairs()));
  assembly_manager.EvaluateForceStiff(discretizations, fm_, fs_, Cmm_, Css_, Csm_, Cms_);
  Cmm_->Complete();
  Csm_->Complete();
}
void ADAPTER::FBIConstraintBridgePenalty::ResetBridge()
{
  ADAPTER::FBIConstraintBridge::ResetBridge();
  fs_->PutScalar(0.0);
  Cmm_->Reset();
  Csm_->Reset();
}
