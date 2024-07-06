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
NOX::Nln::SINGLESTEP::Group::Group(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& grpOptionParams,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& i, const ::NOX::Epetra::Vector& x,
    const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys)
    : ::NOX::Epetra::Group(printParams, i, x, linSys),
      NOX::Nln::Group(printParams, grpOptionParams, i, x, linSys)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Nln::SINGLESTEP::Group::Group(const NOX::Nln::SINGLESTEP::Group& source, ::NOX::CopyType type)
    : ::NOX::Epetra::Group(source, type), NOX::Nln::Group(source, type)
{
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Abstract::Group> NOX::Nln::SINGLESTEP::Group::clone(::NOX::CopyType type) const
{
  Teuchos::RCP<::NOX::Abstract::Group> newgrp =
      Teuchos::rcp(new NOX::Nln::SINGLESTEP::Group(*this, type));
  return newgrp;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
::NOX::Abstract::Group& NOX::Nln::SINGLESTEP::Group::operator=(const ::NOX::Epetra::Group& source)
{
  NOX::Nln::Group::operator=(source);
  return *this;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::SINGLESTEP::Group::computeX(
    const ::NOX::Abstract::Group& grp, const ::NOX::Abstract::Vector& d, double step)
{
  // Cast to appropriate type, then call the "native" computeX
  const NOX::Nln::SINGLESTEP::Group* nlngrp =
      dynamic_cast<const NOX::Nln::SINGLESTEP::Group*>(&grp);
  if (nlngrp == nullptr) throw_error("computeX", "dyn_cast to nox_nln_group failed!");
  const ::NOX::Epetra::Vector& epetrad = dynamic_cast<const ::NOX::Epetra::Vector&>(d);

  computeX(*nlngrp, epetrad, step);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::SINGLESTEP::Group::computeX(
    const NOX::Nln::SINGLESTEP::Group& grp, const ::NOX::Epetra::Vector& d, double step)
{
  prePostOperatorPtr_->run_pre_compute_x(grp, d.getEpetraVector(), step, *this);

  if (isPreconditioner()) sharedLinearSystem.getObject(this)->destroyPreconditioner();

  resetIsValid();

  step = 1.0;
  xVector.update(-1.0, d, step, grp.xVector);
  prePostOperatorPtr_->run_post_compute_x(grp, d.getEpetraVector(), step, *this);
  return;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void NOX::Nln::SINGLESTEP::Group::throw_error(
    const std::string& functionName, const std::string& errorMsg) const
{
  std::ostringstream msg;
  msg << "ERROR - NOX::Nln::SINGLESTEP::Group::" << functionName << " - " << errorMsg << std::endl;

  FOUR_C_THROW(msg.str());
}

FOUR_C_NAMESPACE_CLOSE
