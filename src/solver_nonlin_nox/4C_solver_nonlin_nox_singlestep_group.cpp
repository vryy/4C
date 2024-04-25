/*-----------------------------------------------------------*/
/*! \file

\brief %NOX::NLN implementation of a %::NOX::Epetra::Group
       to use with nonlinear singlestep solver.

\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_solver_nonlin_nox_singlestep_group.hpp"

#include "4C_solver_nonlin_nox_group_prepostoperator.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::SINGLESTEP::Group::Group(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& grpOptionParams,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& i, const ::NOX::Epetra::Vector& x,
    const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys)
    : ::NOX::Epetra::Group(printParams, i, x, linSys),
      NOX::NLN::Group(printParams, grpOptionParams, i, x, linSys)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::NLN::SINGLESTEP::Group::Group(const NOX::NLN::SINGLESTEP::Group& source, ::NOX::CopyType type)
    : ::NOX::Epetra::Group(source, type), NOX::NLN::Group(source, type)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Abstract::Group> NOX::NLN::SINGLESTEP::Group::clone(::NOX::CopyType type) const
{
  Teuchos::RCP<::NOX::Abstract::Group> newgrp =
      Teuchos::rcp(new NOX::NLN::SINGLESTEP::Group(*this, type));
  return newgrp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group& NOX::NLN::SINGLESTEP::Group::operator=(const ::NOX::Epetra::Group& source)
{
  NOX::NLN::Group::operator=(source);
  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::SINGLESTEP::Group::computeX(
    const ::NOX::Abstract::Group& grp, const ::NOX::Abstract::Vector& d, double step)
{
  // Cast to appropriate type, then call the "native" computeX
  const NOX::NLN::SINGLESTEP::Group* nlngrp =
      dynamic_cast<const NOX::NLN::SINGLESTEP::Group*>(&grp);
  if (nlngrp == nullptr) throwError("computeX", "dyn_cast to nox_nln_group failed!");
  const ::NOX::Epetra::Vector& epetrad = dynamic_cast<const ::NOX::Epetra::Vector&>(d);

  computeX(*nlngrp, epetrad, step);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::SINGLESTEP::Group::computeX(
    const NOX::NLN::SINGLESTEP::Group& grp, const ::NOX::Epetra::Vector& d, double step)
{
  prePostOperatorPtr_->runPreComputeX(grp, d.getEpetraVector(), step, *this);

  if (isPreconditioner()) sharedLinearSystem.getObject(this)->destroyPreconditioner();

  resetIsValid();

  step = 1.0;
  xVector.update(-1.0, d, step, grp.xVector);
  prePostOperatorPtr_->runPostComputeX(grp, d.getEpetraVector(), step, *this);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::NLN::SINGLESTEP::Group::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::NLN::SINGLESTEP::Group::" << functionName << " - " << errorMsg << std::endl;

  FOUR_C_THROW(msg.str());
}

FOUR_C_NAMESPACE_CLOSE
