/*----------------------------------------------------------------------------*/
/*!
\file nln_problem_nox.cpp

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
#include "nln_problem_nox.H"
#include "nln_utils.H"
#include "nln_utils_debugwriter.H"

#include "../drt_io/io.H"

#include "../drt_lib/drt_dserror.H"

#include "../linalg/linalg_sparseoperator.H"

/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
NLNSOL::NlnProblemNox::NlnProblemNox() : noxgrp_(Teuchos::null), jac_(Teuchos::null) { return; }

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemNox::Setup()
{
  // Make sure that Init() has been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }

  // Setup() has been called
  SetIsSetup();

  return;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnProblemNox::SetModelEvaluator()
{
  if (not Params()->isParameter("NOX Group"))
  {
    dserror("The parameter 'NOX Group' has not been set.");
    return false;
  }

  if (not Params()->isType<Teuchos::RCP<NOX::Abstract::Group>>("NOX Group"))
  {
    dserror(
        "Parameter 'NOX Group' isn't of type "
        "'Teuchos::RCP<NOX::Abstract::Group>'.");
    return false;
  }

  noxgrp_ = Params()->get<Teuchos::RCP<NOX::Abstract::Group>>("NOX Group");

  return true;
}

/*----------------------------------------------------------------------------*/
bool NLNSOL::NlnProblemNox::SetJacobianOperator()
{
  if (not Params()->isParameter("Jacobian Operator"))
  {
    dserror("The parameter 'Jacobian Operator' has not been set.");
    return false;
  }

  if (not Params()->isType<Teuchos::RCP<LINALG::SparseOperator>>("Jacobian Operator"))
  {
    dserror(
        "Parameter 'Jacobian Operator' isn't of type "
        "'Teuchos::RCP<LINALG::SparseOperator>'.");
    return false;
  }

  jac_ = Params()->get<Teuchos::RCP<LINALG::SparseOperator>>("Jacobian Operator");

  return true;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemNox::ComputeF(const Epetra_MultiVector& x, Epetra_MultiVector& f) const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnProblemNox::ComputeF");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() and Setup() have been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }
  if (not IsSetup())
  {
    dserror("Setup() has not been called, yet.");
  }

  // check for correctness of maps
  dsassert(x.Map().PointSameAs(f.Map()), "Maps do not match.");

  // set most recent solution to NOX group
  Teuchos::RCP<Epetra_Vector> xtemp = Teuchos::rcp(new Epetra_Vector(*(x(0))));
  NOX::Epetra::Vector noxvec(xtemp, NOX::Epetra::Vector::CreateView);
  NOXGroup().setX(noxvec);

  // ask time integrator to evaluate residual and apply DBCs etc.
  NOX::Abstract::Group::ReturnType ret = NOXGroup().computeF();
  if (ret != NOX::Abstract::Group::Ok) dserror("computeF() returned error code %d", ret);

  // extract residual from NOX::Abstract::Group
  Teuchos::RCP<const NOX::Abstract::Vector> noxFRcp = NOXGroup().getFPtr();
  if (noxFRcp.is_null()) dserror("Could not extract the residual from the NOX Group.");
  Teuchos::RCP<const NOX::Epetra::Vector> noxFRcpEpetra =
      Teuchos::rcp_dynamic_cast<const NOX::Epetra::Vector>(noxFRcp, true);
  Teuchos::RCP<Epetra_MultiVector> frcp =
      Teuchos::rcp(new Epetra_MultiVector(noxFRcpEpetra->getEpetraVector()));
  if (frcp.is_null()) dserror("Could not extract Epetra_Vector from NOX::Epetra::Vector.");

  // switch sign to obtain a residual in descent direction (=f_ext - f_int)
  int err = f.Update(-1.0, *frcp, 0.0);
  if (err != 0)
  {
    dserror("Update failed.");
  }

  // check for correctness of maps
  if (not x.Map().PointSameAs(f.Map()))
  {
    dserror("Maps do not match.");
  }
  dsassert(x.Map().PointSameAs(f.Map()), "Maps do not match.");

  return;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_MultiVector> NLNSOL::NlnProblemNox::ComputeF(
    const Epetra_MultiVector& x) const
{
  Teuchos::RCP<Epetra_MultiVector> f = Teuchos::rcp(new Epetra_MultiVector(x.Map(), true));
  ComputeF(x, *f);

  return f;
}

/*----------------------------------------------------------------------------*/
void NLNSOL::NlnProblemNox::ComputeJacobian() const
{
  // time measurements
  Teuchos::RCP<Teuchos::Time> time =
      Teuchos::TimeMonitor::getNewCounter("NLNSOL::NlnProblemNox::ComputeJacobian");
  Teuchos::TimeMonitor monitor(*time);

  // Make sure that Init() and Setup() have been called
  if (not IsInit())
  {
    dserror("Init() has not been called, yet.");
  }
  if (not IsSetup())
  {
    dserror("Setup() has not been called, yet.");
  }

  // Check whether we have a valid residual
  if (not NOXGroup().isF())
    dserror(
        "Cannot compute the Jacobian matrix, since there is no valid "
        "residual.");

  // ask time integrator to evaluate residual and apply DBCs etc.
  NOX::Abstract::Group::ReturnType ret = NOXGroup().computeJacobian();
  if (ret != NOX::Abstract::Group::Ok) dserror("computeJacobian() returned error code %d", ret);

  return;
}

/*----------------------------------------------------------------------------*/
NOX::Abstract::Group& NLNSOL::NlnProblemNox::NOXGroup() const
{
  // check if NOX group has already been set
  if (noxgrp_.is_null()) dserror("NOX group 'noxgrp_' has not been initialized, yet.");

  return *noxgrp_;
}

/*----------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NLNSOL::NlnProblemNox::GetJacobianOperator() const
{
  // check if Jacobian operator has already been set
  if (jac_.is_null()) dserror("Jacobian operator 'jac_' has not been initialized, yet.");

  //  // check if Jacobian operator is valid // ToDo (mayr) re-introduce safety check
  //  if (not NOXGroup().isJacobian())
  //    dserror("Jacobian operator is not up-to-date.");

  return jac_->EpetraOperator();
}
