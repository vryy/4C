/*----------------------------------------------------------------------*/
/*!
\file ale_laplace.cpp

\brief ALE implementation

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/


#include "ale_laplace.H"
#include "ale_utils_mapextractor.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_precond.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"

#define scaling_infnorm true

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ALE::AleLaplace::AleLaplace(RCP<DRT::Discretization> actdis,
                              Teuchos::RCP<LINALG::Solver> solver,
                              Teuchos::RCP<ParameterList> params,
                              Teuchos::RCP<IO::DiscretizationWriter> output,
                              bool incremental,
                              bool dirichletcond)
  : Ale(actdis,solver,params,output,dirichletcond),
    incremental_(incremental)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLaplace::BuildSystemMatrix(bool full)
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
void ALE::AleLaplace::Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp)
{
  // We save the current solution here. This will not change the
  // result of our element call, but the next time somebody asks us we
  // know the displacements.
  //
  // Note: What we get here is the sum of all increments in this time
  // step, not just the latest increment. Be careful.

  if (ddisp!=Teuchos::null)
  {
    // Dirichlet boundaries != 0 are not supported.

    dispnp_->Update(1.0,*ddisp,1.0,*dispn_,0.0);
  }

  if (incremental_)
  {
    EvaluateElements();
    // dispn_ has zeros at the Dirichlet-entries, so we maintain zeros there
    LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispn_,*(dbcmaps_->CondMap()));
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLaplace::Solve()
{
  EvaluateElements();

  // set fixed nodes
  ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  // the DOFs with Dirchlet BCs are not rebuild, they are assumed to be correct
  discret_->EvaluateDirichlet(eleparams,dispnp_,null,null,Teuchos::null,Teuchos::null);

  LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,*(dbcmaps_->CondMap()));

  solver_->Solve(sysmat_->EpetraOperator(),dispnp_,residual_,true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLaplace::Update()
{
  dispn_->Update(1.0,*dispnp_,0.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ALE::AleLaplace::Output()
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
void ALE::AleLaplace::EvaluateElements()
{
  sysmat_->Zero();

  // zero out residual
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // set vector values needed by elements
  discret_->ClearState();

  // action for elements
  eleparams.set("action", "calc_ale_laplace");
  eleparams.set("incremental", incremental_);

  discret_->SetState("dispnp", dispnp_);

  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  sysmat_->Complete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ALE::AleLaplace::ExtractDisplacement() const
{
  // We know that the ale dofs are coupled with their original map. So
  // we just return them here.
  return dispnp_;
}

