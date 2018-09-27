/*----------------------------------------------------------------------------*/
/*!
\file nln_operator_base.cpp

\brief General nonlinear operator base class

\level 3

\maintainer Matthias Mayr
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_MultiVector.h>

// standard
#include <iostream>
#include <vector>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "nln_operator_base.H"
#include "nln_problem_base.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_dserror.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NLNSOL::NlnOperatorBase::NlnOperatorBase()
    : isinit_(false),
      issetup_(false),
      comm_(Teuchos::null),
      config_(Teuchos::null),
      listname_(""),
      nlnproblem_(Teuchos::null),
      outparams_(Teuchos::null),
      nested_(0)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorBase::Init(const Epetra_Comm& comm,
    Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config, const std::string listname,
    Teuchos::RCP<NLNSOL::NlnProblemBase> nlnproblem, Teuchos::RCP<LINALG::Solver> bacisolver,
    const int nested)
{
  // Enforce to call Setup() after Init()
  SetIsSetup(false);

  // fill member variables with given values
  comm_ = Teuchos::rcp(&comm, false);
  config_ = config;
  listname_ = listname;
  nlnproblem_ = nlnproblem;
  bacisolver_ = bacisolver;
  nested_ = nested;

  // initialize member variables
  outparams_ = Teuchos::rcp(new Teuchos::ParameterList());

  getDefaultOStream()->setOutputToRootOnly(0);
  setVerbLevel(NLNSOL::UTILS::TranslateVerbosityLevelToTeuchos(
      MyGetParameter<std::string>("nonlinear operator: verbosity")));

  // Init() has been called
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const Epetra_Comm& NLNSOL::NlnOperatorBase::Comm() const
{
  if (comm_.is_null()) dserror("Communicator 'comm_' has not been set, yet.");

  return *comm_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> NLNSOL::NlnOperatorBase::Configuration() const
{
  // check if configuration object has already been set
  if (config_.is_null()) dserror("Configuration 'config_' has not been initialized, yet.");

  return config_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<NLNSOL::NlnProblemBase> NLNSOL::NlnOperatorBase::NlnProblem() const
{
  // check if nonlinear problem has already been set
  if (nlnproblem_.is_null())
    dserror(
        "The nonlinear problem 'nlnproblem_' has not been initialized, "
        "yet.");

  return nlnproblem_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorBase::PrintIterSummary(const int iter, const double fnorm2) const
{
  if (getVerbLevel() > Teuchos::VERB_NONE and
      MyGetParameter<bool>("Nonlinear Operator: Print Iterations"))
  {
    *getOStream() << LabelShort() << " iteration " << iter << ":\t|f| = " << std::scientific
                  << std::setprecision(6) << fnorm2 << std::endl;
  }
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnOperatorBase::ContinueIterations(const int iter, const bool converged) const
{
  // initialize return value
  bool retval = false;

  if (iter < GetMaxIter() && not converged) retval = true;

  return retval;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnOperatorBase::CheckSuccessfulConvergence(const int iter, const bool converged) const
{
  // initialize return value
  bool successful = false;

  // make the decision
  if (IsSolver())
  {
    if (converged && iter < GetMaxIter())  // successful convergence
      return true;
    else  // convergence failed
    {
      //      dserror("%s did not converge in %d iterations.", Label(), iter);
      return false;
    }
  }
  else  // no solver, so we don't care for convergence
  {
    return true;
  }

  return successful;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnOperatorBase::IsSolver() const
{
  return MyGetParameter<bool>("Nonlinear Operator: Is Solver");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int NLNSOL::NlnOperatorBase::GetMaxIter() const
{
  return MyGetParameter<int>("Nonlinear Operator: Max Iter");
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NLNSOL::UTILS::OperatorStatus NLNSOL::NlnOperatorBase::ErrorCode(
    const int iter, const bool converged, const bool error, const bool stagnation) const
{
  // initialize error code
  NLNSOL::UTILS::OperatorStatus errorcode = NLNSOL::UTILS::opstatus_undefined;

  // determine error code (order of if-statements is important here!!!)
  if (error)
    errorcode = NLNSOL::UTILS::opstatus_failed;
  else if (stagnation)
    errorcode = NLNSOL::UTILS::opstatus_stagnation;
  else if (CheckSuccessfulConvergence(iter, converged))
    errorcode = NLNSOL::UTILS::opstatus_converged;
  else if (not CheckSuccessfulConvergence(iter, converged))
    errorcode = NLNSOL::UTILS::opstatus_unconverged;
  else
    dserror("Cannot determine error code.");

  // set error code to output parameter list as well
  SetOutParameter("Error Code", errorcode);

  // return error code
  return errorcode;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList> NLNSOL::NlnOperatorBase::GetOutParams() const
{
  if (outparams_.is_null())
    dserror(
        "Output parameter list has not been initialized, yet.\n "
        "It is still Teuchos::null.");

  return outparams_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorBase::SetOutParameterIter(const int iter) const
{
  SetOutParameter("Iterations", iter);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorBase::SetOutParameterConverged(const bool converged) const
{
  SetOutParameter("Converged", converged);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorBase::SetOutParameterResidualNorm(const double norm) const
{
  SetOutParameter("Residual Norm", norm);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorBase::SetOutParameterErrorCode(
    const NLNSOL::UTILS::OperatorStatus errorcode) const
{
  SetOutParameter("Error Code", errorcode);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnOperatorBase::SetOutParameterStagnation(
    Teuchos::RCP<const Teuchos::ParameterList> stagparams) const
{
  if (not stagparams.is_null())
  {
    // copy entries from stagnation sublist to operator output list
    SetOutParameter(
        "Stagnation Detection: status", stagparams->get<bool>("Stagnation Detection: status"));
    SetOutParameter(
        "Stagnation Detection: ratio", stagparams->get<double>("Stagnation Detection: ratio"));
    SetOutParameter("Stagnation Detection: iterations",
        stagparams->get<int>("Stagnation Detection: iterations"));
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> NLNSOL::NlnOperatorBase::BaciLinearSolver() const
{
  if (bacisolver_.is_null()) dserror("No valid Baci linear solver set, yet.");

  return bacisolver_;
}
