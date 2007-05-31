
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_ale.H"

using namespace std;
using namespace Teuchos;


FSI::AleLinear::AleLinear(RefCountPtr<DRT::Discretization> actdis,
                          LINALG::Solver&       solver,
                          ParameterList&        params,
                          DiscretizationWriter& output)
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
  nstep_   = params_.get<int>("nstep");
  maxtime_ = params_.get<double>("maxtime");
  dt_      = params_.get<double>("dt");

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  dispnp_         = LINALG::CreateVector(*dofrowmap,true);
  dispn_          = LINALG::CreateVector(*dofrowmap,true);
  dispnm_         = LINALG::CreateVector(*dofrowmap,true);

  residual_       = LINALG::CreateVector(*dofrowmap,true);

  dirichtoggle_ = LINALG::CreateVector(*dofrowmap,true);
  zeros_        = LINALG::CreateVector(*dofrowmap,true);
}


void FSI::AleLinear::Integrate()
{
  while (step_ < nstep_-1 && time_ <= maxtime_)
  {
    step_ += 1;
    time_ += dt_;

    ParameterList eleparams;

    eleparams.set("total time", time_);
    eleparams.set("delta time", dt_);

    discret_->EvaluateDirichlet(eleparams,*dispnp_,*dirichtoggle_);

    //double norm;
    //dispnp_->Norm2(&norm);
    //cout << "dispnp norm: " << norm << endl;

    EvaluateElements();

    zeros_->PutScalar(0.0);
    //LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,zeros_,dirichtoggle_);
    LINALG::ApplyDirichlettoSystem(sysmat_,dispnp_,residual_,dispnp_,dirichtoggle_);

    //cout << *dispnp_;
    //cout << *dirichtoggle_;
    //cout << *residual_;

    solver_.Solve(sysmat_,dispnp_,residual_,true,true);

    output_.NewStep    (step_,time_);
    output_.WriteVector("dispnp", dispnp_);
  }
}


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


bool FSI::AleLinear::computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
{
  if (!haveJacobian_)
  {
    EvaluateElements();
  }
  return true;
}


void FSI::AleLinear::EvaluateElements()
{
  haveF_ = true;
  haveJacobian_ = true;

  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // zero out the stiffness matrix
  sysmat_ = LINALG::CreateMatrix(*dofrowmap,maxentriesperrow_);

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

#endif
#endif
