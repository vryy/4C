/*----------------------------------------------------------------------*/
/*!
\file timint_adjoint.H
\brief Time integration for a hyperleastic quasi static adjoint problem

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

!*/

/*----------------------------------------------------------------------*/
/* headers */
#include <iostream>
#include "../drt_io/io_pstream.H"
#include "Epetra_SerialDenseVector.h"
#include "Teuchos_TimeMonitor.hpp"

#include "timint_adjoint.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_mat/matpar_bundle.H"

#include "../drt_lib/drt_locsys.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"

#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*/
/* constructor */
STR::TimIntAdjoint::TimIntAdjoint(Teuchos::RCP<DRT::Discretization> discret):
discret_(discret),
solver_(Teuchos::null),
dbcmaps_(Teuchos::rcp(new LINALG::MapExtractor())),
dis_(Teuchos::null),
rhs_(Teuchos::null),
stiff_(Teuchos::null),
zeros_(Teuchos::null),
stepn_(0),
isinit_(false)
{
  PrintLogo();

  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  if (not discret_->Filled() || not discret_->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");
  else
    dofrowmap_ = discret_->DofRowMap();

  //setup size of the solution vector
  msteps_ = sdyn.get<int>("NUMSTEP");

  // initialize stiffness matrix
  stiff_ = Teuchos::rcp(new LINALG::SparseMatrix(*(dofrowmap_), 81, true, true));

  //initialize zeros from the dofrowmap
  zeros_ = LINALG::CreateVector(*(dofrowmap_), true);

  // initialize solution at stepn
  disdualn_ = LINALG::CreateVector(*(dofrowmap_), true);
  disn_ = LINALG::CreateVector(*(dofrowmap_), true);
  rhsn_ = LINALG::CreateVector(*(dofrowmap_), true);

  disdual_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_,msteps_,true));

  // get the solver associated with structural dynamic section
  CreateLinearSolver(sdyn);

  // Get Dirichlet Map Extractor
  GetDBCMap();

}

/*----------------------------------------------------------------------*/
/* Setup some quantities */
void STR::TimIntAdjoint::SetupAdjoint(Teuchos::RCP<Epetra_MultiVector> rhs,
                                      Teuchos::RCP<Epetra_MultiVector> dis)
{
  // for now, this is done by the optimizer which has control over the objective function anyway.
  // Take care of the proper oder of the entries in MStep vector!
  rhs_ = rhs;
  dis_ = dis;

  stepn_ = msteps_;
  if (rhs_ != Teuchos::null && dis_ != Teuchos::null)
    isinit_ = true;

}

/*----------------------------------------------------------------------*/
/* integrate */
void STR::TimIntAdjoint::Integrate()
{
  if (!isinit_)
    dserror("rhs is not set up properly for the adjoint equation");

  //integration for the adjoints is reverse! For the quasi static case it doesn't matter
  while ( stepn_ > 0)
  {
    PrepareStep(); // -> stepn_ is stepn_-1 now
    Solve();
    UpdateStep();
  }

  isinit_ = false;

  return;
}

/*----------------------------------------------------------------------*/
/* integrate */
void STR::TimIntAdjoint::EvaluateStiff()
{
  // initialize stiffness matrix to zero
  stiff_->Zero();

  // a dummy internal force vector
  Teuchos::RCP<Epetra_Vector> fintn;
  fintn = LINALG::CreateVector(*(dofrowmap_), true);

  // set the parameters for the discretization
  Teuchos::ParameterList p;
  const std::string action = "calc_struct_nlnstiff";
  p.set("action", action);

  discret_->ClearState();
  discret_->SetState(0,"residual displacement", disn_);
  discret_->SetState(0,"displacement", disn_);

  discret_->Evaluate(p, stiff_, Teuchos::null, fintn, Teuchos::null, Teuchos::null);
  stiff_->Complete();

  discret_->ClearState();

}

/*----------------------------------------------------------------------*/
/* Solve the linear system */
void STR::TimIntAdjoint::Solve()
{
  // Apply Dirichlet
  LINALG::ApplyDirichlettoSystem(stiff_, disdualn_, rhsn_,Teuchos::null,zeros_,*(dbcmaps_->CondMap()));

  // Solve
  solver_->Solve(stiff_->EpetraOperator(),disdualn_,rhsn_,true,true);

  return;
}

/*----------------------------------------------------------------------*/
/* Prepare variables */
void STR::TimIntAdjoint::PrepareStep()
{
  // decrease stepn first
  stepn_ -= 1;

  // extract variables needed for this "time step"
  disn_->Update(1.0,*(*dis_)(stepn_),0.0);
  rhsn_->Update(-1.0,*(*rhs_)(stepn_),0.0);
  // !!! rhs is the optimality condition differentiated w.r.t the displacement; so it has to be
  //multiplied by -1 to be the correct rhsn_ for the dual problem

  // get stiffness matrix according to stepn_
  EvaluateStiff();

  return;
}

/*----------------------------------------------------------------------*/
/* Update step and time */
void STR::TimIntAdjoint::UpdateStep()
{
  (*disdual_)(stepn_)->Update(1.0,*disdualn_,0.0);

  cout << "Finalized step " << stepn_+1 << " / " << msteps_ << endl;

  return;
}

void STR::TimIntAdjoint::CreateLinearSolver(const Teuchos::ParameterList& sdyn)
{
  // get the solver number used for structural problems
  const int linsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL DYNAMIC to a valid number!");

  solver_ = Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                        discret_->Comm(),
                        DRT::Problem::Instance()->ErrorFile()->Handle()));
  discret_->ComputeNullSpaceIfNecessary(solver_->Params());

  return;
}

/*----------------------------------------------------------------------*/
/* Update step and time */
void STR::TimIntAdjoint::GetDBCMap()
{
  // Get Dirichlet Map
  Teuchos::ParameterList p;
  //p.set("total time", timen_);
  discret_->EvaluateDirichlet(p, zeros_, Teuchos::null, Teuchos::null,Teuchos::null, dbcmaps_);
  zeros_->PutScalar(0.0); // paranoia

  return;
}

/*----------------------------------------------------------------------*/
/* UPrint Logo */
void STR::TimIntAdjoint::PrintLogo()
{
  cout << "--------------------------------------------------" << endl;
  cout << "--   Welcome to the adjoint time integration    --" << endl;
  cout << "--------------------------------------------------" << endl;
}
