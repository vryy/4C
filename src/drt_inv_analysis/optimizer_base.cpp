/*----------------------------------------------------------------------*/
/*!
\brief Optimization algorithm base class

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/
#include "optimizer_base.H"

#include "invana_base.H"
#include "initial_guess.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"

#include <Teuchos_ParameterList.hpp>

/*----------------------------------------------------------------------*/
/* constructor */
INVANA::OptimizerBase::OptimizerBase(const Teuchos::ParameterList& invp)
    : error_incr_(1.0e6),
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
      inpar_(invp),
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
  solrowmap_ = optprob->VectorRowLayout();

  sol_ = Teuchos::rcp(new Epetra_MultiVector(*solrowmap_, 1));
  sol_o_ = Teuchos::rcp(new Epetra_MultiVector(*solrowmap_, 1));

  objgrad_ = Teuchos::rcp(new Epetra_MultiVector(*solrowmap_, 1));
  objgrad_o_ = Teuchos::rcp(new Epetra_MultiVector(*solrowmap_, 1));

  optprob_ = optprob;

  // for convenience
  DRT::Problem* problem = DRT::Problem::Instance();

  // restart from this
  if (DRT::Problem::Instance()->Restart())
    inputfile_ = Teuchos::rcp(
        new IO::InputControl(problem->InputControlFile()->FileName(), optprob_->Comm()));

  // setup output to the control file (ignoring the IO section's bin io flag)
  int binio = 1;
  // the restart counter was already increased upon setting the initial control
  // file in the problem! Therefore dont increase counter again here!
  bool adaptname = false;
  Teuchos::RCP<IO::OutputControl> controlfile = Teuchos::rcp(new IO::OutputControl(
      OptProb()->Comm(), problem->ProblemName(), problem->SpatialApproximation(),
      problem->OutputControlFile()->InputFileName(), problem->OutputControlFile()->RestartName(),
      problem->OutputControlFile()->FileName(), problem->NDim(), problem->Restart(),
      problem->OutputControlFile()->FileSteps(), binio, adaptname));


  output_ = Teuchos::rcp(new IO::DiscretizationWriter(optprob_->Discret()));
  output_->SetOutput(controlfile);

  // wrap output
  writer_ = Teuchos::rcp(new InvanaWriter());
  writer_->Init(output_);
  // give a mesh to the output
  writer_->WriteMesh(0, 0.0);

  sol_->Scale(1.0, *optprob_->InitialGuess()->Mean());

  isinit_ = true;
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBase::Evaluate(double* val, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  optprob_->Evaluate(*sol_, val, gradient);
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBase::UpdateSolution(const Epetra_MultiVector& toadd)
{
  sol_o_->Scale(1.0, *sol_);
  sol_->Update(1.0, toadd, 1.0);
}

/*----------------------------------------------------------------------*/
void INVANA::OptimizerBase::UndoUpdateSolution() { sol_->Scale(1.0, *sol_o_); }
