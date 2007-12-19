/*----------------------------------------------------------------------*/
/*!
\file

\brief Solve FSI problems using a Dirichlet-Neumann partitioning approach

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include "fsi_ale.H"

using namespace std;
using namespace Teuchos;


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::AleLinear::AleLinear(RCP<DRT::Discretization> actdis,
                          Teuchos::RCP<LINALG::Solver> solver,
                          Teuchos::RCP<ParameterList> params,
                          Teuchos::RCP<IO::DiscretizationWriter> output)
  : interface_(actdis),
    discret_(actdis),
    solver_ (solver),
    params_ (params),
    output_ (output),
    step_(0),
    time_(0.0),
    maxentriesperrow_(81),
    sysmat_(null),
    haveF_(false),
    haveJacobian_(false),
    restartstep_(0),
    uprestart_(params->get("write restart every", -1))
{
  numstep_   = params_->get<int>("numstep");
  maxtime_ = params_->get<double>("maxtime");
  dt_      = params_->get<double>("dt");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  dispnp_         = LINALG::CreateVector(*dofrowmap,true);

  residual_       = LINALG::CreateVector(*dofrowmap,true);

  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);

  sumdisi_        = LINALG::CreateVector(*dofrowmap,true);

  interface_.SetupCondDofMap("FSICoupling");

  // needed for MFSI
  interface_.SetupOtherDofMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;

  ParameterList eleparams;
  eleparams.set("total time", time_);
  eleparams.set("delta time", dt_);

  sumdisi_->PutScalar(0.);

  discret_->EvaluateDirichlet(eleparams,*dispnp_,*dirichtoggle_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::Evaluate(Teuchos::RCP<const Epetra_Vector> disp)
{
  // We save the current solution here. This will not change the
  // result of our element call, but the next time somebody asks us we
  // know the displacements.
  //
  // Note: What we get here is the sum of all increments in this time
  // step, not just the latest increment. Be careful.
  if (disp!=Teuchos::null)
  {
    dispnp_->Update(1.0,*disp,-1.0,*sumdisi_,1.0);
    sumdisi_->Update(1.0,*disp,0.0);
  }

  EvaluateElements();
  LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,dirichtoggle_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::Solve()
{
  Evaluate(Teuchos::null);
  solver_->Solve(sysmat_,dispnp_,residual_,true,true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::Update()
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::Output()
{
  // We do not need any output -- the fluid writes its
  // displacements itself. But we need restart.

  restartstep_ += 1;

  output_->NewStep    (step_,time_);

  if (restartstep_ == uprestart_)
  {
    restartstep_ = 0;

    // add restart data
    // do we need any variables at all?

    output_->WriteVector("dispnp", dispnp_);
    //output_->WriteVector("dispn",  dispn_);
    //output_->WriteVector("dispnm", dispnm_);
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::Integrate()
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
bool FSI::AleLinear::computeF(const Epetra_Vector &x,
                              Epetra_Vector &F,
                              const FillType fillFlag)
{
  if (!haveF_)
  {
    EvaluateElements();
  }
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool FSI::AleLinear::computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
{
  if (!haveJacobian_)
  {
    EvaluateElements();
  }
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::EvaluateElements()
{
  haveF_ = true;
  haveJacobian_ = true;

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

#if 1
  sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);
#else
  // zero out the stiffness matrix
  if (sysmat_==Teuchos::null)
    sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);
  else
    sysmat_->PutScalar(0.);
#endif

  // zero out residual
  residual_->PutScalar(0.0);

  // create the parameters for the discretization
  ParameterList eleparams;

  // action for elements
  eleparams.set("action", "calc_ale_lin_stiff");

  // other parameters that might be needed by the elements

  // set vector values needed by elements
  discret_->ClearState();
  //discret_->SetState("dispnp", dispnp_);

  discret_->Evaluate(eleparams,sysmat_,residual_);
  discret_->ClearState();

  LINALG::Complete(*sysmat_);
  maxentriesperrow_ = sysmat_->MaxNumEntries();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::SetInterfaceMap(Teuchos::RCP<Epetra_Map> im)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::ApplyInterfaceDisplacements(Teuchos::RCP<Epetra_Vector> idisp)
{
  interface_.InsertCondVector(idisp,dispnp_);

  idisp->PutScalar(1.0);
  interface_.InsertCondVector(idisp,dirichtoggle_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::AleLinear::ExtractDisplacement()
{
  // We know that the ale dofs are coupled with their original map. So
  // we just return them here.
  return dispnp_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::AleLinear::StructCondRHS()
{
  return interface_.ExtractCondVector(dispnp_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::ReadRestart(int step)
{
  IO::DiscretizationReader reader(discret_,step);
  time_ = reader.ReadDouble("time");
  step_ = reader.ReadInt("step");

  // What variables do we need?

  reader.ReadVector(dispnp_, "dispnp");
  //reader.ReadVector(dispn_,  "dispn");
  //reader.ReadVector(dispnm_, "dispnm");
}


#endif
