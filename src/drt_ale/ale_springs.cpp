/*----------------------------------------------------------------------*/
/*!
\file ale_springs.cpp

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#include "ale_springs.H"
#include "ale_utils_mapextractor.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_precond.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_locsys.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ALE::AleSprings::AleSprings(Teuchos::RCP<DRT::Discretization> actdis,
                                Teuchos::RCP<LINALG::Solver> solver,
                                Teuchos::RCP<Teuchos::ParameterList> params,
                                Teuchos::RCP<IO::DiscretizationWriter> output,
                                bool dirichletcond)
  : Ale(actdis,solver,params,output,dirichletcond)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap();
  incr_           = LINALG::CreateVector(*dofrowmap,true);

  xffinterface_ = Teuchos::rcp(new ALE::UTILS::XFluidFluidMapExtractor);
  xffinterface_->Setup(*actdis);

  if (xffinterface_->XFluidFluidCondRelevant())
  {
    // create the toggle vector for fluid-fluid-Coupling
    Teuchos::RCP<Epetra_Vector> dispnp_xff = LINALG::CreateVector(*xffinterface_->XFluidFluidCondMap(),true);
    dispnp_xff->PutScalar(1.0);
    xfftoggle_ = LINALG::CreateVector(*discret_->DofRowMap(),true);
    xffinterface_->InsertXFluidFluidCondVector(dispnp_xff,xfftoggle_);
  }
  // temporay only to ensure that old ale formulation still works with changed
  // ale elements j biehler 10/14
  incremental_=true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleSprings::BuildSystemMatrix(bool full)
{
  if (full)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81,false,true));
  }
  else
  {
    sysmat_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*interface_,*interface_,81,false,true));
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleSprings::Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp, std::string incrementtype)
{
  // We save the current solution here. This will not change the
  // result of our element call, but the next time somebody asks us we
  // know the displacements.
  //
  // Note: What we get here is the sum of all increments in this time
  // step, not just the latest increment. Be careful.


  if (ddisp!=Teuchos::null and incrementtype == "step")
  {
    // Dirichlet -boundaries != 0 are not supported.
    incr_->Update(1.0,*ddisp,1.0,*dispn_,0.0);
    dispnp_->Update(1.0,*incr_,0.0);
  }
  else if (ddisp!=Teuchos::null and incrementtype == "iter")
  {
    if(dispnp_->Update(1.0,*ddisp,1.0))
    {
      dserror("Update of ale displacements not correct");
    }
  }
  else if (ddisp!=Teuchos::null and !(incrementtype=="step") and !(incrementtype=="iter"))
    dserror("Type of increment for ale evaluated not specified correctly. "
            "Choose either 'step' (default) for step increments d^{n+1}_{i+1} = d^n + increment "
            "or choose 'iter' for iteration increments  d^{n+1}_{i+1} = d^{n+1}_{i} + increment "
           );

  if (ddisp!=Teuchos::null)
  {
    // Dirichlet boundaries != 0 are not supported. hahn: Why?!

    incr_->Update(1.0,*ddisp,1.0,*dispn_,0.0);
    dispnp_->Update(1.0,*incr_,0.0);
  }

  EvaluateElements();

  if (LocsysManager() != Teuchos::null) {
    // Transform system matrix and rhs to local coordinate systems
    LocsysManager()->RotateGlobalToLocal(Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_), residual_);

    // When using local systems, a rotated dispnp_ vector needs to be used as dbcval for ApplyDirichlettoSystem
    Teuchos::RCP<Epetra_Vector> dispnp_local = Teuchos::rcp(new Epetra_Vector(*(dispnp_)));
    LocsysManager()->RotateGlobalToLocal(dispnp_local);

    LINALG::ApplyDirichlettoSystem(sysmat_,incr_,residual_,GetLocSysTrafo(),dispnp_local,*(dbcmaps_->CondMap()));
  } else {
    LINALG::ApplyDirichlettoSystem(sysmat_,incr_,residual_,dispnp_,*(dbcmaps_->CondMap()));
  }

  if (xffinterface_->XFluidFluidCondRelevant())
  {
    Teuchos::RCP<Epetra_Vector>  dispnp_ttt = LINALG::CreateVector(*discret_->DofRowMap(),true);
    LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_ttt,xfftoggle_);

    // set dispnp_ of xfem dofs to dispn_
    xffinterface_->InsertXFluidFluidCondVector(xffinterface_->ExtractXFluidFluidCondVector(dispn_), dispnp_);
  }

}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleSprings::Solve()
{
  EvaluateElements();

  // set fixed nodes
  Teuchos::ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);

  ALE::Ale::ApplyDirichletBC(eleparams,dispnp_,Teuchos::null,Teuchos::null,false);

  incr_->Update(1.0,*dispnp_,-1.0,*dispn_,0.0);

  if (LocsysManager() != Teuchos::null) {
    // Transform system matrix and rhs to local coordinate systems
    LocsysManager()->RotateGlobalToLocal(Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_), residual_);

    // When using local systems, a rotated incr_ vector needs to be used as dbcval for ApplyDirichlettoSystem
    Teuchos::RCP<Epetra_Vector> incr_local = Teuchos::rcp(new Epetra_Vector(*(incr_)));
    LocsysManager()->RotateGlobalToLocal(incr_local);

    LINALG::ApplyDirichlettoSystem(sysmat_,incr_,residual_,GetLocSysTrafo(),incr_local,*(dbcmaps_->CondMap()));
  } else {
    LINALG::ApplyDirichlettoSystem(sysmat_,incr_,residual_,incr_,*(dbcmaps_->CondMap()));
  }

  if (xffinterface_->XFluidFluidCondRelevant())
      LINALG::ApplyDirichlettoSystem(sysmat_,incr_,residual_,incr_,xfftoggle_);

  solver_->Solve(sysmat_->EpetraOperator(),incr_,residual_,true);

  incr_->Update(1.0,*dispn_,1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleSprings::SolveBioGr()
{
  EvaluateElements();

  // set fixed nodes
  Teuchos::ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);

  incr_->Update(1.0,*dispnp_,-1.0,*dispn_,0.0);

  if (LocsysManager() != Teuchos::null)
  {
    // Transform system matrix and rhs to local coordinate systems
    LocsysManager()->RotateGlobalToLocal(Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_), residual_);

    // When using local systems, a rotated incr_ vector needs to be used as dbcval for ApplyDirichlettoSystem
    Teuchos::RCP<Epetra_Vector> incr_local = Teuchos::rcp(new Epetra_Vector(*(incr_)));
    LocsysManager()->RotateGlobalToLocal(incr_local);

    LINALG::ApplyDirichlettoSystem(sysmat_,incr_,residual_,GetLocSysTrafo(),incr_local,*(dbcmaps_->CondMap()));
  }
  else
  {
    LINALG::ApplyDirichlettoSystem(sysmat_,incr_,residual_,incr_,*(dbcmaps_->CondMap()));
  }

  solver_->Solve(sysmat_->EpetraOperator(),incr_,residual_,true);

  incr_->Update(1.0,*dispn_,1.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleSprings::Update()
{
  dispn_-> Update(1.0,*incr_,0.0);
  dispnp_->Update(1.0,*incr_,0.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleSprings::EvaluateElements()
{
  // The ALE-springs-algorithm computes its stiffness matrix for the actual
  // displacement! adapter_ale_springs_fixed_ref only uses the initial mesh geometry

  sysmat_->Zero();

  // zero out residual
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // set vector values needed by elements
  discret_->ClearState();

  // action for elements
  eleparams.set("action", "calc_ale_spring");
  eleparams.set("incremental", incremental_);

  // newton-increments in displacements are of no importance for ALE-Evaluate
  discret_->SetState("dispnp", dispn_);
  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  sysmat_->Complete();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ALE::AleSprings::Dispnp() const
{
  return incr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> ALE::AleSprings::Dispn() const
{
  dserror("not implemented!");
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ALE::AleSprings::WriteAccessDispnp() const
{
  return incr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleSprings::SolveAleXFluidFluidFSI()
{
  // At the beginning of the fluid-fluid-fsi step the xfem-dofs are
  // dirichlet values so that they can not change in the next
  // iterations. After the fsi step we put the ALE FSI-dofs to
  // dirichlet and we solve the ALE again to find the real ALE
  // displacement.

  // turn the toggle vector off
  xfftoggle_->PutScalar(0.0);

  // new toggle vector which is on for the fsi-dofs_
  Teuchos::RCP<Epetra_Vector> dispnp_fsicond = LINALG::CreateVector(*interface_->FSICondMap(),true);
  dispnp_fsicond->PutScalar(1.0);
  interface_->InsertFSICondVector(dispnp_fsicond,xfftoggle_);

  BuildSystemMatrix(true);

  Solve();

  // for the next time step set the xfem dofs to dirichlet values
  Teuchos::RCP<Epetra_Vector> dispnp_xff = LINALG::CreateVector(*xffinterface_->XFluidFluidCondMap(),true);
  dispnp_xff->PutScalar(1.0);
  xfftoggle_->PutScalar(0.0);
  xffinterface_->InsertXFluidFluidCondVector(dispnp_xff,xfftoggle_);
}
