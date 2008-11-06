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
#ifdef CCADISCRET

#include "adapter_ale_lin.H"

#define scaling_infnorm true


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ADAPTER::AleLinear::AleLinear(RCP<DRT::Discretization> actdis,
                              Teuchos::RCP<LINALG::Solver> solver,
                              Teuchos::RCP<ParameterList> params,
                              Teuchos::RCP<IO::DiscretizationWriter> output,
                              bool incremental,
                              bool dirichletcond)
  : discret_(actdis),
    solver_ (solver),
    params_ (params),
    output_ (output),
    step_(0),
    time_(0.0),
    incremental_(incremental),
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

  DRT::UTILS::SetupNDimExtractor(*actdis,"FSICoupling",interface_);
  DRT::UTILS::SetupNDimExtractor(*actdis,"FREESURFCoupling",freesurface_);

  // set fixed nodes (conditions != 0 are not supported right now)
  ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  discret_->EvaluateDirichlet(eleparams,dispnp_,null,null,null,dbcmaps_);

  if (dirichletcond)
  {
    // for partitioned FSI the interface becomes a Dirichlet boundary
    std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
    condmaps.push_back(interface_.CondMap());
    condmaps.push_back(dbcmaps_->CondMap());
    Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
    *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged);
  }

  if (dirichletcond and freesurface_.Relevant())
  {
    // for partitioned solves the free surface becomes a Dirichlet boundary
    std::vector<Teuchos::RCP<const Epetra_Map> > condmaps;
    condmaps.push_back(freesurface_.CondMap());
    condmaps.push_back(dbcmaps_->CondMap());
    Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
    *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleLinear::BuildSystemMatrix(bool full)
{
  // build linear matrix once and for all
  if (full)
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap,81,false,true));
  }
  else
  {
    sysmat_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(interface_,interface_,81,false,true));
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
        arowsum = rcp(new Epetra_Vector(A->RowMap(),false));
        acolsum = rcp(new Epetra_Vector(A->RowMap(),false));
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
void ADAPTER::AleLinear::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleLinear::Evaluate(Teuchos::RCP<const Epetra_Vector> ddisp)
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
void ADAPTER::AleLinear::Solve()
{
  // set fixed nodes
  ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);
  discret_->EvaluateDirichlet(eleparams,dispnp_,null,null,Teuchos::null);

  LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,*(dbcmaps_->CondMap()));

  solver_->Solve(sysmat_->EpetraOperator(),dispnp_,residual_,true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleLinear::Update()
{
  dispn_->Update(1.0,*dispnp_,0.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleLinear::Output()
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
void ADAPTER::AleLinear::Integrate()
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
void ADAPTER::AleLinear::EvaluateElements()
{
  sysmat_->Zero();

  // zero out residual
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

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
void ADAPTER::AleLinear::ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
{
  interface_.InsertCondVector(idisp,dispnp_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleLinear::ApplyFreeSurfaceDisplacements(Teuchos::RCP<Epetra_Vector> fsdisp)
{
  freesurface_.InsertCondVector(fsdisp,dispnp_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::AleLinear::ExtractDisplacement() const
{
  // We know that the ale dofs are coupled with their original map. So
  // we just return them here.
  return dispnp_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleLinear::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  reader.ReadVector(dispnp_, "dispnp");
  reader.ReadVector(dispn_,  "dispn");
}


#endif
