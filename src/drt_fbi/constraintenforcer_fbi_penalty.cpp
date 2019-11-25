/*----------------------------------------------------------------------*/
/*! \file
\file constraintenforcer_fbi_penalty.cpp

\brief Implements the constraint enforcement technique of a penalty approach (Mortar and GPTS) (for
fluid-beam interaction)

\level 3

\maintainer Nora Hagmeyer
*----------------------------------------------------------------------*/

#include "constraintenforcer_fbi_penalty.H"
#include "constraintenforcer_fbi.H"
#include "ad_fbi_constraintbridge_penalty.H"
#include "beam_to_fluid_meshtying_params.H"
#include "../drt_adapter/ad_str_fbiwrapper.H"
#include "../drt_adapter/ad_fld_moving_boundary.H"
#include <Epetra_Vector.h>

Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::FBIPenaltyConstraintenforcer::AssembleFluidStiffness()
    const
{
  if (Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)
          ->GetCff()
          ->Scale(GetBridge()->GetParams()->GetPenaltyParameter()))
    dserror("Scaling of the penalty stiffness was unsuccessful!\n");

  return Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)
      ->GetCff();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix>
ADAPTER::FBIPenaltyConstraintenforcer::AssembleStructureStiffness() const
{
  return Teuchos::null;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIPenaltyConstraintenforcer::AssembleFluidForce() const
{
  Teuchos::RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(
      (Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)->GetFf())
          ->Map()));
  f->Update(-1.0,
      *(Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)->GetFf()),
      0.0);
  if (f->Scale(GetBridge()->GetParams()->GetPenaltyParameter()))
    dserror("Scaling of the penalty force was unsuccessful!\n");
  std::cout << "The final scaled fluid force looks like " << *f << std::endl;
  return f;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIPenaltyConstraintenforcer::AssembleStructureForce() const
{
  Teuchos::RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(
      (Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)->GetFs())
          ->Map()));
  f->Update(-1.0,
      *(Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)->GetFs()),
      0.0);
  if (f->Scale(GetBridge()->GetParams()->GetPenaltyParameter()))
    dserror("Scaling of the penalty force was unsuccessful!\n");
  return f;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FBIPenaltyConstraintenforcer::PrepareFluidSolve()
{
  GetBridge()->PrepareFluidSolve();
}
