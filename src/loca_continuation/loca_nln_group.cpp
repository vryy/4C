/*-----------------------------------------------------------*/
/*!
\file loca_nln_group.cpp

\brief LOCA specific group

\maintainer Michael Hiermeier

\date Nov 16, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "loca_nln_group.H"

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LOCA::NLN::Group::Group(const Teuchos::RCP<LOCA::GlobalData>& loca_gdata,
    Teuchos::ParameterList& grp_options, Teuchos::ParameterList& print_params,
    const Teuchos::RCP<LOCA::Epetra::Interface::Required>& ireq, NOX::Epetra::Vector& initial_guess,
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys, const LOCA::ParameterVector& loca_params)
    : NOX::Epetra::Group(print_params, ireq, initial_guess, linsys),
      NOX::NLN::Group(print_params, grp_options, ireq, initial_guess, linsys),
      LOCA::Abstract::Group(loca_gdata),
      LOCA::Epetra::Group(loca_gdata, print_params, ireq, initial_guess, linsys, loca_params)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LOCA::NLN::Group::Group(const Teuchos::RCP<LOCA::GlobalData>& loca_gdata,
    Teuchos::ParameterList& grp_options, Teuchos::ParameterList& print_params,
    const Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent>& itime,
    NOX::Epetra::Vector& initial_guess, const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys,
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& shifted_linsys,
    const LOCA::ParameterVector& loca_params)
    : NOX::Epetra::Group(print_params, itime, initial_guess, linsys),
      NOX::NLN::Group(print_params, grp_options, itime, initial_guess, linsys),
      LOCA::Abstract::Group(loca_gdata),
      LOCA::Epetra::Group(
          loca_gdata, print_params, itime, initial_guess, linsys, shifted_linsys, loca_params)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
LOCA::NLN::Group::Group(const LOCA::NLN::Group& source, NOX::CopyType type)
    : NOX::Epetra::Group(source, type),
      NOX::NLN::Group(source, type),
      LOCA::Abstract::Group(source, type),
      LOCA::Epetra::Group(source, type)
{
  /* no new variables, pointers or references */
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<NOX::Abstract::Group> LOCA::NLN::Group::clone(NOX::CopyType type) const
{
  return Teuchos::rcp(new LOCA::NLN::Group(*this, type));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group& LOCA::NLN::Group::operator=(const NOX::Abstract::Group& source)
{
  return NOX::NLN::Group::operator=(source);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group& LOCA::NLN::Group::operator=(const NOX::Epetra::Group& source)
{
  return NOX::NLN::Group::operator=(source);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType LOCA::NLN::Group::computeF()
{
  /* Set the parameters prior to computing F
   * We do no direct return, because of the underlying pre/post operator
   * of the NOX::NLN::Group implementation. */
  if (not isF()) userInterface->setParameters(params);

  return NOX::NLN::Group::computeF();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType LOCA::NLN::Group::computeFandJacobian()
{
  /* Set the parameters prior to computing F, if
   * not already done.
   * We do no direct return, because of the underlying pre/post operator
   * of the NOX::NLN::Group implementation. */
  if ((not isF()) and (not isJacobian())) userInterface->setParameters(params);

  return NOX::NLN::Group::computeFandJacobian();
}
