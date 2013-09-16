/*----------------------------------------------------------------------*/
/*!
\file adapter_ale_lin.cpp

\brief ALE implementation

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#include "ale_lin.H"
#include "ale_utils_mapextractor.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_precond.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_mapextractor.H"
#include "../drt_io/io.H"

#define scaling_infnorm true

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ALE::AleLinear::AleLinear(RCP<DRT::Discretization> actdis,
                              Teuchos::RCP<LINALG::Solver> solver,
                              Teuchos::RCP<Teuchos::ParameterList> params,
                              Teuchos::RCP<IO::DiscretizationWriter> output,
                              bool incremental,
                              bool dirichletcond)
  : Ale(actdis,solver,params,output,dirichletcond),
    incremental_(incremental)
{
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
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLinear::BuildSystemMatrix(bool full)
{
  // build linear matrix once and for all
  if (full)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81,false,true));
  }
  else
  {
    sysmat_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*interface_,*interface_,81,false,true));
  }

  if (not incremental_)
  {
    EvaluateElements();
    LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,*(dbcmaps_->CondMap()));

   // prepare constant preconditioner on constant matrix

    if (full)
    {
      // partitioned FSI does not use explicit preconditioner objects
    }
    else
    {
      // This is the MFSI case and we need the preconditioner on the inner dofs only
      precond_ = Teuchos::rcp(new LINALG::Preconditioner(LinearSolver()));

      Teuchos::RCP<Epetra_CrsMatrix> A = BlockSystemMatrix()->Matrix(0,0).EpetraMatrix();

      Teuchos::RCP<Epetra_Vector> arowsum;
      Teuchos::RCP<Epetra_Vector> acolsum;

      if (scaling_infnorm)
      {
        arowsum = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
        acolsum = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
        A->InvRowSums(*arowsum);
        A->InvColSums(*acolsum);
        if (A->LeftScale(*arowsum) or
            A->RightScale(*acolsum))
          dserror("ale scaling failed");
      }

      precond_->Setup(A);

      if (scaling_infnorm)
      {
        arowsum->Reciprocal(*arowsum);
        acolsum->Reciprocal(*acolsum);
        if (A->LeftScale(*arowsum) or
            A->RightScale(*acolsum))
          dserror("ale scaling failed");
      }
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLinear::Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp, std::string incrementtype)
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
    dispnp_->Update(1.0,*ddisp,1.0,*dispn_,0.0);
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

  if (incremental_)
  {
    EvaluateElements();
    // dispn_ has zeros at the Dirichlet-entries, so we maintain zeros there
    LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispn_,*(dbcmaps_->CondMap()));

    if (xffinterface_->XFluidFluidCondRelevant())
    {
      Teuchos::RCP<Epetra_Vector>  dispnp_ttt = LINALG::CreateVector(*discret_->DofRowMap(),true);
      LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_ttt,xfftoggle_);

      // set dispnp_ of xfem dofs to dispn_
      xffinterface_->InsertXFluidFluidCondVector(xffinterface_->ExtractXFluidFluidCondVector(dispn_), dispnp_);
    }
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLinear::Solve()
{
  // set fixed nodes
  Teuchos::ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  // the DOFs with Dirchlet BCs are not rebuild, they are assumed to be correct
  if (incremental_)
    EvaluateElements();

  discret_->EvaluateDirichlet(eleparams,dispnp_,Teuchos::null,Teuchos::null,Teuchos::null,Teuchos::null);
  LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,*(dbcmaps_->CondMap()));

  if (xffinterface_->XFluidFluidCondRelevant())
    LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,xfftoggle_);

  solver_->Solve(sysmat_->EpetraOperator(),dispnp_,residual_,true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLinear::SolveBioGr()
{
  solver_->Solve(sysmat_->EpetraOperator(),dispnp_,residual_,true);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLinear::Update()
{
  dispn_->Update(1.0,*dispnp_,0.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLinear::Output()
{
  // We do not need any output -- the fluid writes its
  // displacements itself. But we need restart.

  if (uprestart_ != 0 and step_ % uprestart_ == 0)
  {
    output_->NewStep    (step_,time_);
    output_->WriteVector("dispnp", dispnp_);

    // add restart data
    output_->WriteVector("dispn", dispn_);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLinear::EvaluateElements()
{
  sysmat_->Zero();

  // zero out residual
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // set vector values needed by elements
  discret_->ClearState();

  // action for elements
  eleparams.set("action", "calc_ale_lin_stiff");
  eleparams.set("incremental", incremental_);

  discret_->SetState("dispnp", dispnp_);

  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  sysmat_->Complete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLinear::SolveAleXFluidFluidFSI()
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
