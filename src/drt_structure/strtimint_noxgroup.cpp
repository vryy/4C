/*----------------------------------------------------------------------------*/
/*!
\file strtimint_noxgroup.cpp

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------------*/
/* definitions */

/*----------------------------------------------------------------------------*/
/* headers */
#include "strtimint_noxgroup.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NOX::STR::Group::Group(::STR::TimIntImpl& sti, Teuchos::ParameterList& printParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& i, const NOX::Epetra::Vector& x,
    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys)
    : NOX::Epetra::Group(printParams, i, x, linSys)
{
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::STR::Group::computeF()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Epetra::Group::computeF();

  // Not sure why we call this (historical reasons???)
  sharedLinearSystem.getObject(this);

  // We just evaluated the residual, so its linearization is not valid anymore.
  isValidJacobian = false;

  return ret;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
NOX::Abstract::Group::ReturnType NOX::STR::Group::computeJacobian()
{
  NOX::Abstract::Group::ReturnType ret = NOX::Epetra::Group::computeJacobian();
  if (ret == NOX::Abstract::Group::Ok)
  {
    isValidJacobian = true;

    if (not isValidRHS) isValidRHS = true;
  }

  return ret;
}
