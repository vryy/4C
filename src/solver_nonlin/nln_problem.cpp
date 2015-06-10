/*----------------------------------------------------------------------------*/
/*!
\file nln_problem.cpp

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
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

#include "../linalg/linalg_sparseoperator.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NLNSOL::NlnProblem::NlnProblem()
: isinit_(false),
  issetup_(false),
  comm_(Teuchos::null),
  params_(Teuchos::null),
  noxgrp_(Teuchos::null),
  jac_(Teuchos::null),
  dbgwriter_(Teuchos::null),
  tolresl2_(0.0),
  lengthscaling_(true)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblem::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& params, NOX::Abstract::Group& noxgrp,
    Teuchos::RCP<LINALG::SparseOperator> jac,
    Teuchos::RCP<NLNSOL::UTILS::DebugWriterBase> dbgwriter)
{
  // We need to call Setup() after Init()
  issetup_ = false;

  // fill member variables without taking memory ownership
  comm_ = Teuchos::rcp(&comm, false);
  params_ = Teuchos::rcp(&params, false);
  noxgrp_ = Teuchos::rcp(&noxgrp, false);
  jac_ = jac;
  dbgwriter_ = dbgwriter;

  // read some parameters from parameter list and store them separately
  tolresl2_ = Params().get<double>("Nonlinear Problem: Tol Res L2");
  lengthscaling_ = Params().get<bool>("Nonlinear Problem: Length Scaled Norms");

  // set the verbosity level
  setVerbLevel(
      NLNSOL::UTILS::TranslateVerbosityLevel(
          Params().get<std::string>("Nonlinear Problem: Verbosity")));

  // Init() has been called
  SetIsInit();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblem::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblem::ComputeF(const Epetra_MultiVector& x,
    Epetra_MultiVector& f) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnProblem::ComputeF");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() and Setup() have been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // check for correctness of maps
  dsassert(x.Map().PointSameAs(f.Map()), "Maps do not match.");

  // set most recent solution to NOX group
  Teuchos::RCP<Epetra_Vector> xtemp = Teuchos::rcp(new Epetra_Vector(*(x(0))));
  NOX::Epetra::Vector noxvec(xtemp, NOX::Epetra::Vector::CreateCopy); // ToDo (mayr) Change back to CreateView some day?
  NOXGroup().setX(noxvec);

  // ask time integrator to evaluate residual and apply DBCs etc.
  NOX::Abstract::Group::ReturnType ret = NOXGroup().computeF();
  if (ret != NOX::Abstract::Group::Ok)
    dserror("computeF() returned error code %d", ret);

  // extract residual from NOX::Abstract::Group
  Teuchos::RCP<const NOX::Abstract::Vector> noxFRcp = NOXGroup().getFPtr();
  if (noxFRcp.is_null())
    dserror("Could not extract the residual from the NOX Group.");
  Teuchos::RCP<const NOX::Epetra::Vector> noxFRcpEpetra =
      Teuchos::rcp_dynamic_cast<const NOX::Epetra::Vector>(noxFRcp, true);
  Teuchos::RCP<Epetra_MultiVector> frcp =
      Teuchos::rcp(new Epetra_MultiVector(noxFRcpEpetra->getEpetraVector()));
  if (frcp.is_null())
    dserror("Could not extract Epetra_Vector from NOX::Epetra::Vector.");

  // switch sign to obtain a residual in descent direction
  int err = f.Update(-1.0, *frcp, 0.0);
  if (err != 0) { dserror("Update failed."); }

  // check for correctness of maps
  dsassert(x.Map().PointSameAs(f.Map()), "Maps do not match.");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblem::ComputeJacobian() const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time = Teuchos::TimeMonitor::getNewCounter(
      "NLNSOL::NlnProblem::ComputeJacobian");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() and Setup() have been called
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  // Check whether we have a valid residual
  if (not NOXGroup().isF())
    dserror("Cannot compute the Jacobian matrix, since there is no valid "
        "residual.");

  // ask time integrator to evaluate residual and apply DBCs etc.
  NOX::Abstract::Group::ReturnType ret = NOXGroup().computeJacobian();
  if (ret != NOX::Abstract::Group::Ok)
    dserror("computeJacobian() returned error code %d", ret);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnProblem::ConvergenceCheck(const Epetra_MultiVector& f) const
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  double fnorm2 = 0.0;

  return ConvergenceCheck(f, fnorm2);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnProblem::ConvergenceCheck(const Epetra_MultiVector& f,
    double& fnorm2) const
{
  if (not IsInit()) { dserror("Init() has not been called, yet."); }
  if (not IsSetup()) { dserror("Setup() has not been called, yet."); }

  if (f.GlobalLength() <= 0)
    dserror("Cannot compute norm of empty residual vector!");

  // ---------------------------------------------------------------------------
  // compute norm
  // ---------------------------------------------------------------------------
  int err = f.Norm2(&fnorm2);
  if (err != 0) { dserror("Failed!"); }

  if (lengthscaling_)
    fnorm2 /= sqrt(f.GlobalLength());
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
const Epetra_Comm& NLNSOL::NlnProblem::Comm() const
{
  // check if communicator has already been set
  if (comm_.is_null())
    dserror("Communicator 'comm_' has not been set, yet.");

  return *comm_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const Teuchos::ParameterList& NLNSOL::NlnProblem::Params() const
{
  // check if parameter list has already been set
  if (params_.is_null())
    dserror("Parameter list 'params_' has not been initialized, yet.");

  return *params_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NOX::Abstract::Group& NLNSOL::NlnProblem::NOXGroup() const
{
  // check if NOX group has already been set
  if (noxgrp_.is_null())
    dserror("NOX group 'noxgrp_' has not been initialized, yet.");

  return *noxgrp_;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NLNSOL::NlnProblem::GetJacobianOperator() const
{
  // check if Jacobian operator has already been set
  if (jac_.is_null())
    dserror("Jacobian operator 'jac_' has not been initialized, yet.");

//  // check if Jacobian operator is valid // ToDo (mayr) re-introduce safety check
//  if (not NOXGroup().isJacobian())
//    dserror("Jacobian operator is not up-to-date.");

  return jac_->EpetraOperator();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const Epetra_Map& NLNSOL::NlnProblem::DofRowMap() const
{
  // check if Jacobian operator has already been set
  if (jac_.is_null())
    dserror("Jacobian operator 'jac_' has not been initialized, yet.");

  return jac_->OperatorRangeMap();
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblem::WriteVector(Teuchos::RCP<const Epetra_MultiVector> vec,
    const std::string& description,
    const IO::DiscretizationWriter::VectorType vt) const
{
  if (HaveDebugWriter())
  {
    dbgwriter_->WriteVector(vec, description, vt);
  }
  else
  {
    if (getVerbLevel() > Teuchos::VERB_NONE)
    {
      *getOStream() << Label()
          << ": WARNING: Cant't write debug output of vector '"
          << description << "', since debug writer 'dbgwriter_' has not been "
              "set properly, yet."
          << std::endl;
    }
  }

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblem::WriteVector(const Epetra_MultiVector& vec,
    const std::string& description,
    const IO::DiscretizationWriter::VectorType vt) const
{
  WriteVector(Teuchos::rcp(&vec, false), description, vt);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
const bool NLNSOL::NlnProblem::HaveDebugWriter() const
{
  return ((not dbgwriter_.is_null()) and dbgwriter_->IsSetup());
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
Teuchos::RCP<NLNSOL::UTILS::DebugWriterBase> NLNSOL::NlnProblem::DebugWriter() const
{
  return dbgwriter_;
}
