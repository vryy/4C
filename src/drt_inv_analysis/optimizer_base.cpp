/*----------------------------------------------------------------------*/
/*!
\file optimizer_base.cpp

\brief Optimization algorithm base class

<pre>
\level 3
\maintainer Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/
#include "optimizer_base.H"

#include "invana_base.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*/
/* constructor */
INVANA::OptimizerBase::OptimizerBase(const Teuchos::ParameterList& invp):
error_incr_(1.0e6),
restartevry_(invp.get<int>("RESTARTEVRY")),
maxiter_(invp.get<int>("MAXITER")),
stepsize_(invp.get<double>("STEPSIZE")),
runc_(0),
convtol_(invp.get<double>("CONVTOL")),
convcritc_(0),
objval_(0.0),
objval_o_(0.0),
objgrad_(Teuchos::null),
objgrad_o_(Teuchos::null),
sol_(Teuchos::null),
sol_o_(Teuchos::null),
optprob_(Teuchos::null),
solrowmap_(Teuchos::null),
output_(Teuchos::null),
writer_(Teuchos::null),
inputfile_(Teuchos::null),
isinit_(false)
{
  return;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBase::Init(Teuchos::RCP<InvanaBase> optprob)
{
  solrowmap_=optprob->VectorRowLayout();

  sol_ = Teuchos::rcp(new Epetra_MultiVector(*solrowmap_, 1));
  sol_o_ = Teuchos::rcp(new Epetra_MultiVector(*solrowmap_, 1));

  objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*solrowmap_, 1));
  objgrad_o_ = Teuchos::rcp(new Epetra_MultiVector(*solrowmap_, 1));

  optprob_=optprob;

  // output for the optimizer: outputcontrol is "copied"/reproduced to "steal" it from the discretization
  // and give it to the inverse analysis algorithm
  if (DRT::Problem::Instance()->Restart())
    inputfile_ = Teuchos::rcp(new IO::InputControl(DRT::Problem::Instance()->InputControlFile()->FileName(), optprob_->Comm()));

  // setup output
  output_ = Teuchos::rcp(new IO::DiscretizationWriter(optprob_->Discret()));
  output_->SetOutput(DRT::Problem::Instance()->OutputControlFile());

  // wrap output
  writer_ = Teuchos::rcp(new InvanaWriter());
  writer_->Init(output_);

  SetInitialGuess();

  isinit_=true;

}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBase::SetInitialGuess()
{
  sol_->Scale(1.0,optprob_->InitialGuess());
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBase::Evaluate(double* val, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  optprob_->Evaluate(*sol_,val,gradient);
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBase::UpdateSolution(const Epetra_MultiVector& toadd)
{
  sol_o_->Scale(1.0, *sol_);
  sol_->Update(1.0,toadd,1.0);
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBase::UndoUpdateSolution()
{
  sol_->Scale(1.0, *sol_o_);
}

