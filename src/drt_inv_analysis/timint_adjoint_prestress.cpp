/*----------------------------------------------------------------------*/
/*! \file
\brief Time integration for a hyperleastic quasi static adjoint problem including prestressing

\level 3

\maintainer Sebastian Brandstaeter

!*/
#include "timint_adjoint_prestress.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_create.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../drt_lib/drt_locsys.H"

/*----------------------------------------------------------------------*/
/* constructor                                               keh 12/14  */
/*----------------------------------------------------------------------*/
STR::TimIntAdjointPrestress::TimIntAdjointPrestress(Teuchos::RCP<DRT::Discretization> discret)
    : TimIntAdjoint(discret), stiffp_(Teuchos::null)
{
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  stiffp_ = Teuchos::rcp(new LINALG::SparseMatrix(*(dofrowmap_), 81, true, true));

  // initialize solution at stepn
  disdualnp_ = LINALG::CreateVector(*(dofrowmap_), true);
  disdualp_ = Teuchos::rcp(new Epetra_MultiVector(*dofrowmap_, msteps_, true));

  // rhs for the prestress adjoint equation
  rhsnp_ = LINALG::CreateVector(*(dofrowmap_), true);

  // prestress stuff
  pstype_ = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn, "PRESTRESS");
  pstime_ = sdyn.get<double>("PRESTRESSTIME");
}


/*----------------------------------------------------------------------*/
/* bring primal solution and rhs in                          keh 12/14  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjointPrestress::SetupAdjoint(Teuchos::RCP<Epetra_MultiVector> rhs,
    std::vector<double> mtime, Teuchos::RCP<Epetra_MultiVector> dis, std::vector<double> time)
{
  rhs_ = rhs;
  dis_ = dis;
  time_ = time;
  mtime_ = mtime;

  if ((int)time_.size() != msteps_)
    dserror("setup of the timesteps for the adjoint problem messed up");

  if (pstype_ == INPAR::STR::PreStress::mulf)
  {
    // get prestress time
    stepn_ = (int)(pstime_ / dt_);
    timen_ = time_[stepn_ - 1];
    disn_->Scale(0.0);

    // evaluate stiffness matrix
    EvaluateStiff();

    // do locsys already here since the matrix can be reused
    // and must not be re-evaluated before every solve
    if (locsysman_ != Teuchos::null)
      locsysman_->RotateGlobalToLocal(
          Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(stiff_), zeros_);

    zeros_->Scale(0.0);

    // store ready to use stiffness matrix of the prestressed state
    stiffp_->Add(*stiff_, true, 1.0, 0.0);
    stiffp_->Complete();
    stiffp_->ApplyDirichlet(*(dbcmaps_->CondMap()));
  }

  stepn_ = msteps_;
  timen_ = time_[stepn_ - 1];

  if (rhs_ != Teuchos::null && dis_ != Teuchos::null) isinit_ = true;
}


/*----------------------------------------------------------------------*/
/* solve the linear system                                   keh 12/14  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjointPrestress::Solve()
{
  // no need to solve for a linear system with rhs of zeros
  if (steprhsn_ == -1)
  {
    disdualn_->Scale(0.0);
    disdualnp_->Scale(0.0);
    return;
  }
  //-------------------------------------------------------------
  //----- SOLVE for lagrange multiplier of the final state

  // transform rhs to local co-ordinate systems
  if (locsysman_ != Teuchos::null)
    locsysman_->RotateGlobalToLocal(Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(stiff_), rhsn_);

  // Apply Dirichlet
  LINALG::ApplyDirichlettoSystem(
      stiff_, disdualn_, rhsn_, GetLocSysTrafo(), zeros_, *(dbcmaps_->CondMap()));
  solver_->Solve(stiff_->EpetraOperator(), disdualn_, rhsn_, true, true, Teuchos::null);


  //-------------------------------------------------------------
  //----- SOLVE for lagrange multiplier of the prestress state

  // setup of the modified rhs
  stiffn_->SetUseTranspose(true);
  stiffn_->Apply(*disdualn_, *rhsnp_);
  rhsnp_->Update(-1.0, *rhsn_, 1.0);

  LINALG::ApplyDirichlettoSystem(disdualnp_, rhsnp_, zeros_, *(dbcmaps_->CondMap()));
  solver_->Solve(stiffp_->EpetraOperator(), disdualnp_, rhsnp_, true, true, Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*/
/* update state                                              keh 12/14  */
/*----------------------------------------------------------------------*/
void STR::TimIntAdjointPrestress::UpdateStepState()
{
  // Update State
  (*disdual_)(stepn_ - 1)->Update(1.0, *disdualn_, 0.0);
  (*disdualp_)(stepn_ - 1)->Update(1.0, *disdualnp_, 0.0);

  return;
}
