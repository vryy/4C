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
#ifdef TRILINOS_PACKAGE

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
FSI::AleLinear::AleLinear(RefCountPtr<DRT::Discretization> actdis,
                          Teuchos::RefCountPtr<LINALG::Solver> solver,
                          Teuchos::RefCountPtr<ParameterList> params,
                          Teuchos::RefCountPtr<IO::DiscretizationWriter> output)
  : discret_(actdis),
    solver_ (solver),
    params_ (params),
    output_ (output),
    step_(0),
    time_(0.0),
    maxentriesperrow_(81),
    sysmat_(null),
    haveF_(false),
    haveJacobian_(false)
{
  nstep_   = params_->get<int>("nstep");
  maxtime_ = params_->get<double>("maxtime");
  dt_      = params_->get<double>("dt");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  dispnp_         = LINALG::CreateVector(*dofrowmap,true);
  dispn_          = LINALG::CreateVector(*dofrowmap,true);
  dispnm_         = LINALG::CreateVector(*dofrowmap,true);

  residual_       = LINALG::CreateVector(*dofrowmap,true);

  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);
  //zeros_        = LINALG::CreateVector(*dofrowmap,true);
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

  discret_->EvaluateDirichlet(eleparams,*dispnp_,*dirichtoggle_);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::Solve()
{
  EvaluateElements();

  LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,dirichtoggle_);

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
#if 0
  // currently we do not need any output -- the fluid writes its
  // displacements itself.
  output_->NewStep    (step_,time_);
  output_->WriteVector("dispnp", dispnp_);
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::Integrate()
{
  while (step_ < nstep_-1 and time_ <= maxtime_)
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

  // choose what to assemble
  eleparams.set("assemble matrix 1",true);
  eleparams.set("assemble matrix 2",false);
  eleparams.set("assemble vector 1",true);
  eleparams.set("assemble vector 2",false);
  eleparams.set("assemble vector 3",false);

  // other parameters that might be needed by the elements

  // set vector values needed by elements
  discret_->ClearState();
  //discret_->SetState("dispnp", dispnp_);

  discret_->Evaluate(eleparams,sysmat_,null,residual_,null,null);
  discret_->ClearState();

  LINALG::Complete(*sysmat_);
  maxentriesperrow_ = sysmat_->MaxNumEntries();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::SetInterfaceMap(Teuchos::RefCountPtr<Epetra_Map> im)
{
  imeshmap_ = im;
  extractor_ = rcp(new Epetra_Import(*imeshmap_, *discret_->DofRowMap()));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::AleLinear::ApplyInterfaceDisplacements(Teuchos::RefCountPtr<Epetra_Vector> idisp)
{
  int err = dispnp_->Export(*idisp,*extractor_,Insert);
  if (err)
    dserror("Insert using extractor returned err=%d",err);

  idisp->PutScalar(1.0);
  err = dirichtoggle_->Export(*idisp,*extractor_,Insert);
  if (err)
    dserror("Insert using extractor returned err=%d",err);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector> FSI::AleLinear::ExtractDisplacement()
{
  // We know that the ale dofs are coupled with their original map. So
  // we just return them here.
  return dispnp_;
}


#endif
#endif
