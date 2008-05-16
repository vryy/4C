/*----------------------------------------------------------------------*/
/*!
\file

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

#include "adapter_ale.H"

// further includes for AleBaseAlgorithm:
#include "../drt_lib/drt_globalproblem.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

using namespace std;
using namespace Teuchos;


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

  UTILS::SetupNDimExtractor(*actdis,"FSICoupling",interface_);
  UTILS::SetupNDimExtractor(*actdis,"FREESURFCoupling",freesurface_);

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

  EvaluateElements();
  LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,dirichtoggle_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const LINALG::SparseMatrix* ADAPTER::AleLinear::InteriorMatrixBlock() const
{
  LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* bm =
    dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>*>(&*sysmat_);
  if (bm!=NULL)
  {
    return &bm->Matrix(0,0);
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const LINALG::SparseMatrix* ADAPTER::AleLinear::InterfaceMatrixBlock() const
{
  LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* bm =
    dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>*>(&*sysmat_);
  if (bm!=NULL)
  {
    return &bm->Matrix(0,1);
  }
  return NULL;
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

  LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,dirichtoggle_);

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

  UTILS::SetupNDimExtractor(*actdis,"FSICoupling",interface_);
  UTILS::SetupNDimExtractor(*actdis,"FREESURFCoupling",freesurface_);

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
const LINALG::SparseMatrix* ADAPTER::AleSprings::InteriorMatrixBlock() const
{
  LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* bm =
    dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>*>(&*sysmat_);
  if (bm!=NULL)
  {
    return &bm->Matrix(0,0);
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const LINALG::SparseMatrix* ADAPTER::AleSprings::InterfaceMatrixBlock() const
{
  LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>* bm =
    dynamic_cast<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>*>(&*sysmat_);
  if (bm!=NULL)
  {
    return &bm->Matrix(0,1);
  }
  return NULL;
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

    dispnp_->Update(1.0,*ddisp,1.0,*dispn_,0.0);
  }
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
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::AleSprings::Update()
{
  dispn_->Update(1.0,*incr_,1.0);
  dispnp_->Update(1.0,*dispn_,0.0);
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
  // abuse residual vector to sum dispn_ and incr_ (that contains the solution
  // increment.) residual_ will be reset on the next element call.
  residual_->Update(1.0,*dispn_,1.0,*incr_,0.0);
  return residual_;
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::AleBaseAlgorithm::AleBaseAlgorithm()
{
  SetupAle();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::AleBaseAlgorithm::~AleBaseAlgorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::AleBaseAlgorithm::SetupAle()
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("ADAPTER::AleBaseAlgorithm::SetupAle");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numaf,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[genprob.numaf];

  const Teuchos::ParameterList& adyn     = DRT::Problem::Instance()->AleDynamicParams();
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<ParameterList> solveparams = rcp(new ParameterList());
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set<int>("numstep",    fsidyn.get<int>("NUMSTEP"));
  params->set<double>("maxtime", fsidyn.get<double>("MAXTIME"));
  params->set<double>("dt",      fsidyn.get<double>("TIMESTEP"));

  // ----------------------------------------------- restart and output
  // restart
  params->set<int>("write restart every", fsidyn.get<int>("RESTARTEVRY"));

  int aletype = Teuchos::getIntegralValue<int>(adyn,"ALE_TYPE");
  if (aletype==ALE_DYNAMIC::classic_lin)
    ale_ = rcp(new AleLinear(actdis, solver, params, output));
  else if (aletype==ALE_DYNAMIC::springs)
    ale_ = rcp(new AleSprings(actdis, solver, params, output));
  else
    dserror("ale type '%s' unsupported",adyn.get<std::string>("ALE_TYPE").c_str());
}


#endif
