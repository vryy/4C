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

Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::FBIPenaltyConstraintenforcer::AssembleMasterStiffness()
    const
{
  if (Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)
          ->GetCmm()
          ->Scale(GetBridge()->GetParams()->GetPenaltyParameter()))
    dserror("Scaling of the penalty stiffness was unsuccessful!\n");
  Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)
      ->GetCmm()
      ->Scale(1.0);
  /*  std::cout << "Fluid stiffness is "
              << *(Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)
                         ->GetCmm())
              << std::endl;
              */

  return Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)
      ->GetCmm();
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> ADAPTER::FBIPenaltyConstraintenforcer::AssembleSlaveStiffness()
    const
{
  return Teuchos::null;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIPenaltyConstraintenforcer::AssembleMasterForce() const
{
  Teuchos::RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(
      (Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)->GetFm())
          ->Map()));
  f->Update(-1.0,
      *(Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)->GetFm()),
      0.0);
  if (f->Scale(GetBridge()->GetParams()->GetPenaltyParameter()))
    dserror("Scaling of the penalty force was unsuccessful!\n");
  std::cout << "The final scaled fluid force looks like " << *f << std::endl;
  return f;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FBIPenaltyConstraintenforcer::AssembleSlaveForce() const
{
  Teuchos::RCP<Epetra_Vector> f = Teuchos::rcp(new Epetra_Vector(
      (Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)->GetFs())
          ->Map()));
  f->Update(1.0,
      *(Teuchos::rcp_dynamic_cast<ADAPTER::FBIConstraintBridgePenalty>(GetBridge(), true)->GetFs()),
      0.0);
  if (f->Scale(GetBridge()->GetParams()->GetPenaltyParameter()))
    dserror("Scaling of the penalty force was unsuccessful!\n");
  return f;
}
