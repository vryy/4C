/*----------------------------------------------------------------------*/
/*!
\file adapter_ale_springs.cpp

\brief

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "adapter_ale_springs.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ADAPTER::AleSprings::AleSprings(RCP<DRT::Discretization> actdis,
                                Teuchos::RCP<LINALG::Solver> solver,
                                Teuchos::RCP<ParameterList> params,
                                Teuchos::RCP<IO::DiscretizationWriter> output,
                                bool dirichletcond)
  : discret_(actdis),
    solver_ (solver),
    params_ (params),
    output_ (output),
    step_(0),
    time_(0.0),
    sysmat_(null),
    restartstep_(0),
    uprestart_(params->get("write restart every", -1))
{
  numstep_ = params_->get<int>("numstep");
  maxtime_ = params_->get<double>("maxtime");
  dt_      = params_->get<double>("dt");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  dispn_          = LINALG::CreateVector(*dofrowmap,true);
  dispnp_         = LINALG::CreateVector(*dofrowmap,true);
  residual_       = LINALG::CreateVector(*dofrowmap,true);
  dirichtoggle_   = LINALG::CreateVector(*dofrowmap,true);
  incr_           = LINALG::CreateVector(*dofrowmap,true);

  DRT::UTILS::SetupNDimExtractor(*actdis,"FSICoupling",interface_);
  DRT::UTILS::SetupNDimExtractor(*actdis,"FREESURFCoupling",freesurface_);

  // set fixed nodes (conditions != 0 are not supported right now)
  ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  discret_->EvaluateDirichlet(eleparams,dispnp_,null,null,dirichtoggle_);

  if (dirichletcond)
  {
    // for partitioned FSI the interface becomes a Dirichlet boundary

    Teuchos::RCP<Epetra_Vector> idisp = LINALG::CreateVector(*interface_.CondMap(),false);
    idisp->PutScalar(1.0);
    interface_.InsertCondVector(idisp,dirichtoggle_);
  }

  if (dirichletcond and freesurface_.Relevant())
  {
    // for partitioned solves the free surface becomes a Dirichlet boundary

    Teuchos::RCP<Epetra_Vector> idisp = LINALG::CreateVector(*freesurface_.CondMap(),false);
    idisp->PutScalar(1.0);
    freesurface_.InsertCondVector(idisp,dirichtoggle_);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::BuildSystemMatrix(bool full)
{
  if (full)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81,false,true));
  }
  else
  {
    sysmat_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(interface_,interface_,81,false,true));
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp)
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

    incr_->Update(1.0,*ddisp,1.0,*dispn_,0.0);
  }

  EvaluateElements();
  LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,dirichtoggle_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::Solve()
{
  EvaluateElements();

  // set fixed nodes
  ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  discret_->EvaluateDirichlet(eleparams,dispnp_,null,null,Teuchos::null);

  incr_->Update(1.0,*dispnp_,-1.0,*dispn_,0.0);

  LINALG::ApplyDirichlettoSystem(sysmat_,incr_,residual_,incr_,dirichtoggle_);

  solver_->Solve(sysmat_->EpetraOperator(),incr_,residual_,true);

  incr_->Update(1.0,*dispn_,1.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::Update()
{
  dispn_-> Update(1.0,*incr_,0.0);
  dispnp_->Update(1.0,*incr_,0.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::Output()
{
  // We do not need any output -- the fluid writes its
  // displacements itself. But we need restart.

  restartstep_ += 1;

  output_->NewStep    (step_,time_);
  output_->WriteVector("dispnp", dispnp_);

  if (restartstep_ == uprestart_)
  {
    restartstep_ = 0;

    // add restart data
    output_->WriteVector("dispn", dispn_);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::Integrate()
{
  while (step_ < numstep_-1 and time_ <= maxtime_)
  {
    PrepareTimeStep();
    Solve();
    Update();
    Output();
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::EvaluateElements()
{
  sysmat_->Zero();

  // zero out residual
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // set vector values needed by elements
  discret_->ClearState();

  // action for elements
  discret_->SetState("dispnp", dispn_);
  eleparams.set("action", "calc_ale_spring");

  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  sysmat_->Complete();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
{
  interface_.InsertCondVector(idisp,dispnp_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::ApplyFreeSurfaceDisplacements(Teuchos::RCP<Epetra_Vector> fsdisp)
{
  freesurface_.InsertCondVector(fsdisp,dispnp_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::AleSprings::ExtractDisplacement() const
{
  return incr_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(dispnp_, "dispnp");
  reader.ReadVector(dispn_,  "dispn");
}


#endif
