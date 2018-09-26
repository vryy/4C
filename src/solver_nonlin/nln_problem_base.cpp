/*----------------------------------------------------------------------------*/
/*!
\file nln_problem_base.cpp

\brief Base class for interface of nonlinear solver to BACI

\maintainer Matthias Mayr

\level 3
*/

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* headers */

// standard

// Epetra
#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Epetra_MultiVector.h>
#include <Epetra_Operator.h>
#include <Epetra_Vector.h>

// NOX
#include <NOX_Abstract_Group.H>
#include <NOX_Epetra_Group.H>
#include <NOX_Epetra_Vector.H>

// Teuchos
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TimeMonitor.hpp>

// baci
#include "nln_problem.H"
#include "nln_utils.H"
#include "nln_utils_debugwriter.H"

#include "../drt_io/io.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NLNSOL::NlnProblemBase::NlnProblemBase()
    : isinit_(false),
      issetup_(false),
      comm_(Teuchos::null),
      config_(Teuchos::null),
      params_(Teuchos::null),
      dofrowmap_(Teuchos::null),
      dbgwriter_(Teuchos::null),
      tolresl2_(0.0),
      lengthscaling_(true)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemBase::Init(const Epetra_Comm& comm,
    Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> config, const std::string listname,
    const Teuchos::ParameterList& params, Teuchos::RCP<const Epetra_Map> dofrowmap,
    Teuchos::RCP<NLNSOL::UTILS::DebugWriterBase> dbgwriter)
{
  // We need to call Setup() after Init()
  issetup_ = false;

  // fill member variables without taking memory ownership
  comm_ = Teuchos::rcp(&comm, false);
  config_ = config;
  listname_ = listname;
  params_ = Teuchos::rcp(&params, false);
  dofrowmap_ = dofrowmap;
  dbgwriter_ = dbgwriter;

  // set model evaluator in derived classes
  SetModelEvaluator();
  SetJacobianOperator();

  // read some parameters from parameter list and store them separately
  tolresl2_ = Configuration()->GetParameter<double>(MyListName(), "Nonlinear Problem: Tol Res L2");
  lengthscaling_ =
      Configuration()->GetParameter<bool>(MyListName(), "Nonlinear Problem: Length Scaled Norms");

  // set the verbosity level
  setVerbLevel(NLNSOL::UTILS::TranslateVerbosityLevelToTeuchos(
      Configuration()->GetParameter<std::string>(MyListName(), "Nonlinear Problem: Verbosity")));

  // Init() has been called
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnProblemBase::ConvergenceCheck(const Epetra_MultiVector& f) const
{
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }
  if (not IsSetup())
  {
    dserror("Setup() has not been called, yet.");
  }

  double fnorm2 = 0.0;

  return ConvergenceCheck(f, fnorm2);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnProblemBase::ConvergenceCheck(const Epetra_MultiVector& f, double& fnorm2) const
{
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }
  if (not IsSetup())
  {
    dserror("Setup() has not been called, yet.");
  }

  if (f.GlobalLength() <= 0) dserror("Cannot compute norm of empty residual vector!");

  // ---------------------------------------------------------------------------
  // compute norm
  // ---------------------------------------------------------------------------
  int err = f.Norm2(&fnorm2);
  if (err != 0)
  {
    dserror("Failed!");
  }

  if (lengthscaling_) fnorm2 /= sqrt(f.GlobalLength());
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Check for convergence
  // ---------------------------------------------------------------------------
  bool converged = false;

  if (fnorm2 < tolresl2_)
    converged = true;
  else
    converged = false;
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Print to screen
  // ---------------------------------------------------------------------------
  if (getVerbLevel() > Teuchos::VERB_MEDIUM)
  {
    *getOStream() << "     *** " << Label() << " residual norm = " << fnorm2;

    if (converged)
      *getOStream() << "  --> Converged!" << std::endl;
    else
      *getOStream() << "  --> Failed to converge!" << std::endl;
  }
  // ---------------------------------------------------------------------------

  return converged;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const Epetra_Comm& NLNSOL::NlnProblemBase::Comm() const
{
  // check if communicator has already been set
  if (comm_.is_null()) dserror("Communicator 'comm_' has not been set, yet.");

  return *comm_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const NLNSOL::UTILS::NlnConfig> NLNSOL::NlnProblemBase::Configuration() const
{
  // check if configuration object has already been set
  if (config_.is_null()) dserror("Configuration 'config_' has not been initialized, yet.");

  return config_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Teuchos::ParameterList> NLNSOL::NlnProblemBase::Params() const
{
  // check if parameter list has already been set
  if (params_.is_null()) dserror("Teuchos::ParameterList 'params_' has not been initialized, yet.");

  return params_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> NLNSOL::NlnProblemBase::DofRowMap() const
{
  // check if Jacobian operator has already been set
  if (dofrowmap_.is_null()) dserror("Epetra_Map 'dofrowmap_' has not been initialized, yet.");

  return dofrowmap_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemBase::WriteVector(Teuchos::RCP<const Epetra_MultiVector> vec,
    const std::string& description, const IO::VectorType vt) const
{
  if (HaveDebugWriter())
  {
    dbgwriter_->WriteVector(vec, description, vt);
  }
  else
  {
    if (getVerbLevel() > Teuchos::VERB_NONE)
    {
      *getOStream() << Label() << ": WARNING: Cant't write debug output of vector '" << description
                    << "', since debug writer 'dbgwriter_' has not been "
                       "set properly, yet."
                    << std::endl;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemBase::WriteVector(
    const Epetra_MultiVector& vec, const std::string& description, const IO::VectorType vt) const
{
  WriteVector(Teuchos::rcp(&vec, false), description, vt);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnProblemBase::HaveDebugWriter() const
{
  return ((not dbgwriter_.is_null()) and dbgwriter_->IsSetup());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<NLNSOL::UTILS::DebugWriterBase> NLNSOL::NlnProblemBase::DebugWriter() const
{
  return dbgwriter_;
}
