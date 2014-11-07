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

// baci
#include "nln_problem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
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
  tolresl2_(0.0),
  lengthscaling_(true),
  printconvcheck_(false)
{
  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblem::Init(const Epetra_Comm& comm,
    const Teuchos::ParameterList& params, NOX::Abstract::Group& noxgrp,
    Teuchos::RCP<LINALG::SparseOperator> jac)
{
  // We need to call Setup() after Init()
  issetup_ = false;

  // fill member variables without taking memory ownership
  comm_ = Teuchos::rcp(&comm, false);
  params_ = Teuchos::rcp(&params, false);
  noxgrp_ = Teuchos::rcp(&noxgrp, false);
  jac_ = jac;

  // read some parameters from parameter list and store them separately
  tolresl2_ = Params().get<double>("Nonlinear Problem: Tol Res L2");
  lengthscaling_ = Params().get<bool>("Nonlinear Problem: Length Scaled Norms");
  printconvcheck_ =
      Params().get<bool>("Nonlinear Problem: Print Convergence Check");

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

  int err = f.Update(1.0, *frcp, 0.0);
  if (err != 0) { dserror("Update failed."); }

  // check for correctness of maps
  dsassert(x.Map().PointSameAs(f.Map()), "Maps do not match.");

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblem::ComputeJacobian() const
{
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
  if (printconvcheck_ and Comm().MyPID() == 0)
  {
    IO::cout << "     *** " << Label() << " residual norm = " << fnorm2;

    if (converged)
      IO::cout << "  --> Converged!" << IO::endl;
    else
      IO::cout << "  --> Failed to converge!" << IO::endl;
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
Teuchos::RCP<Epetra_Operator> NLNSOL::NlnProblem::GetJacobianOperator()
{
  // check if Jacobian operator has already been set
  if (jac_.is_null())
    dserror("Jacobian operator 'jac_' has not been initialized, yet.");

//  // check if Jacobian operator is valid
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
